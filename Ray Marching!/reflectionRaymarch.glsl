/*
Van Andrew Nguyen
25/09/21
[Reflection and Refraction]

Followed a tutorial on light refraction for a cube. Uses a 3D texture and redirects the ray in a slightly different direction
so light doesn't pass directly through, giving it a nice look. The index of refraction can be changed with IOR.
*/

// Globals //////////////////////////////////////////////////////////////////////

/*
Max steps is the maximum amount of steps we can make to get a ray length
Max distance is the maximum distance we can shoot a ray
Surf distance is the distance we loop to get to register a hit on the surface
*/
#define MAXSTEPS 100
#define MAXDIS   100.0
#define SURFDIS  0.01

// SDFS //////////////////////////////////////////////////////////////////////

// Capsule Distance
float sdCapsule(vec3 pos, vec3 start, vec3 end, float rad) {
    // We get the projected distance of the origin->start onto start->end
    vec3 ab = end - start;
    vec3 ap = pos - start;
    // Then we normalise that distance and lock it within bounds
    float projection = dot(ab, ap) / dot(ab, ab); // / dot(ab, ab) to normalise 0->1
    projection = clamp(projection, 0.0, 1.0);
    // Then we get the distance from origin to that projected point MINUS radius of capsule
    vec3 len = start + projection * ab;
    float dis = length(pos - len) - rad;
    // Finally, we have the smallest distance to the capsule shape.
    return dis;
}

// Torus Distance
float sdTorus(vec3 pos, vec3 center, vec2 rad) {
    // We subtract the smaller radius from the length of the vector running from
    // origin point to middle of torus
    pos -= center;
    float x = length(pos.xz) - rad.x;
    return length(vec2(x, pos.y)) - rad.y;
}

// Sphere Distance
float sdSphere(vec3 pos, vec3 center, float rad) {
    pos -= center;
    return length(pos) - rad;
}

// Cube Distance
float sdCube(vec3 pos, vec3 size) {
    // We get the length of the position - size (dis) and grab a version > 0
    float val = length(max(abs(pos) - size, 0.0));
    return val;
}

// Shape Operations (From https://www.iquilezles.org/)
float opUnion(float d1, float d2) { return min(d1,d2); }
float opSubtraction(float d1, float d2) { return max(-d1,d2); }
float opIntersection(float d1, float d2) { return max(d1,d2); }

// Rotation Matrix
mat2 rotate(float a) {
    float s = sin(a);
    float c = cos(a);
    return mat2(c, -s, s, c);
}

// Vectors //////////////////////////////////////////////////////////////////////

// Get distance function
float getDist(vec3 pos) {
    
    // Get dist of the various shapes
    vec3 cubePos = pos - vec3(0, 0.8, 5);
    cubePos.xz *= rotate(iTime * 0.5);
    cubePos.yz *= rotate(iTime * 0.5);
    float cubeDis = sdCube(cubePos, vec3(0.6 + 0.05 * sin(iTime))); 
    
    // Final distance to return
    float finalDis = min(cubeDis, cubeDis);
    
    return finalDis;
}

// Return the normal ray
vec3 getNormal(vec3 pos) {
    float dis = getDist(pos);
    vec2 val = vec2(0.01, 0.0);
    
    // To get the slope we give the curve two values super close together
    // Instead of deriving we can do this method via swizzling
    vec3 normal = dis - vec3(getDist(pos-val.xyy), 
                             getDist(pos-val.yxy), 
                             getDist(pos-val.yyx));
    // Same way of doing getDist(p-vec3(0.01, 0.0, 0.0), ... etc
    return normalize(normal);
}

// Ray Marching function
float rayMarch(vec3 rayOrigin, vec3 rayDir) {
    float disFromOrigin = 0.0;
    
    // Now we loop for max steps or until the distance is very small
    for (int i=0;i<MAXSTEPS;i++) {
        // Declare current marching location
        vec3 pos = rayOrigin + rayDir * disFromOrigin;
        float disToScene = getDist(pos);
        disFromOrigin += disToScene;
        
        // Exit condition; we have a hit
        if (disFromOrigin > MAXDIS || disToScene < SURFDIS) { break; }
    }
    
    return disFromOrigin;
}

// Lighting function
float getLight(vec3 pos, vec3 lightOrigin) {
    // Get the light origin
    vec3 lightPos = lightOrigin; // basically a point in 3D space
    //lightPos.xz += vec2(sin(iTime * 8.0), cos(iTime * 8.0)) * 4.0;
    
    // Get the light ray
    vec3 light = normalize(lightPos - pos);
    
    // Get normal
    vec3 normal = getNormal(pos);
    
    // Get diffused lighting 
    // We want 0 if the rays are parallel, 1 if the rays are perpendicular
    // Hence we use dot product
    float dif = clamp(dot(normal, light), 0.0, 1.0);
    float dis = rayMarch(pos + normal * SURFDIS * 2.0, light);
    
    // Getting shadows (if disToSphere < disToLight then we have a shadow)
    if (dis < length(lightPos - pos)) {
        float shadowIntensity = 0.2;
        dif *= shadowIntensity;
    }
    
    return dif;
}

// Main //////////////////////////////////////////////////////////////////////

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab UV Coord
    vec2 uv = (fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;

    // Declare Col
    vec3 col = vec3(0.0);
    vec3 shadowCol = vec3(0.2, 0.1, 0.15);
    vec3 lightCol = vec3(0.8, 0.75, 0.7);

    // Setup Camera
    float camHeight = 1.8;
    float downTilt = -0.2;
    vec3 rayOrigin = vec3(0, camHeight, 0); // this is the camera (origin of vector)
    vec3 rayDir = normalize(vec3(uv.x, uv.y + downTilt, 1));
    
    col = texture(iChannel0, rayDir).rgb;
    
    // Visualise sphere
    float dis = rayMarch(rayOrigin, rayDir); // ray outside object
    if (dis < MAXDIS) {
        vec3 pos = rayOrigin + rayDir * dis;
        float diffuseLight = getLight(pos, vec3(-2, 5, 4));
        vec3 normal = getNormal(pos);
        vec3 mixCol = mix(shadowCol, lightCol, diffuseLight);
        
        // Reflections
        float IOR = 1.1;
        /*
        Index of refraction
        Vacuum  -> 1.00
        Air     -> 1.01
        Water   -> 1.33
        Glass   -> 1.45
        Diamond -> 2.40
        */
        vec3 reflectedRay = reflect(rayDir, normal);
        vec3 refractedRay = refract(rayDir, normal, 1.0 / IOR); // last arg is INDEX of REFRACTION
        float disInterior = rayMarch(pos, refractedRay); // we have to account for the refracted INSIDE the object
        
        // Chromatic Abberation
        float abberation = 0.01;
        vec3 reflectedTex = vec3(0.0); // texture(iChannel0, refractedRay).rgb;
        refractedRay = refract(rayDir, normal, 1.0 / (IOR - abberation));
        reflectedTex.r = texture(iChannel0, refractedRay).r;
        refractedRay = refract(rayDir, normal, 1.0 / IOR);
        reflectedTex.g = texture(iChannel0, refractedRay).g;
        refractedRay = refract(rayDir, normal, 1.0 / (IOR + abberation));
        reflectedTex.b = texture(iChannel0, refractedRay).b;
        
        // Gradual lighting curve
        
        col = vec3(reflectedTex) + mixCol * 0.25;
    }
    
    // Output to screen
    fragColor = vec4(col,1.0);
}