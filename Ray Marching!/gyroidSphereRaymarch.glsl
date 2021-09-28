27/09/21
[Gyroid Sphere]

I experimented with a boolean intersection of a sphere and a gyroid. I added displacement on the sphere surface using 
sine waves and finally got around to properly addressing the issue of lighting an object using it's height as well as the light source.
So now the object evidently has shadows around the low points, and highlights on its raised surface as to be expected in real life.
Paired with the diffuse light; it looks way more realistic. 
The resultant effect is this gross, bubbly raymarched sphere.

Update: I am able to use the light refraction code from another tutorial, so the resultant material looks like a very strange
glass ball, with a lot of bumps.
*/

// Globals //////////////////////////////////////////////////////////////////////

/*
Max steps is the maximum amount of steps we can make to get a ray length
Max distance is the maximum distance we can shoot a ray
Surf distance is the distance we loop to get to register a hit on the surface
*/
#define MAXSTEPS 100
#define MAXDIS 100.0
#define SURFDIS 0.01

#define PI 3.1415

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
float sdCube(vec3 pos, vec3 center, vec3 size) {
    pos -= center;
    
    // We get the length of the position - size (dis) and grab a version > 0
    float val = length(max(abs(pos) - size, 0.0));
    return val;
}

// Gyroid Distance
float sdGyroid(vec3 pos, float scale, float thickness, float offset) {
    pos *= scale;
    float val = abs(dot(sin(pos), cos(pos.zxy)) + offset) / scale - thickness;
    return val;
}

// Shape Operations (From https://www.iquilezles.org/)
float opUnion(float d1, float d2) { return min(d1,d2); }
float opSubtraction(float d1, float d2) { return max(-d1,d2); }
float opIntersection(float d1, float d2) { return max(d1,d2); }

float opDisplacement(vec3 pos, float baseDis) {
    float amp = 0.2;
    float freq = 4.0;
    float t = iTime * 1.0;
    float d1 = baseDis;
    float d2 = amp * sin(pos.x*freq+t) * sin(pos.y*freq+t) * sin(pos.z*freq+t);
    
    return d1 + d2;
}

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
    float sphereDis = opDisplacement(pos, sdSphere(pos, vec3(0, 1, 6), 1.6)); 
    float gyroid1Dis = sdGyroid(pos, 8.0, 0.1, 0.25);
    float gyroid2Dis = sdGyroid(pos, 10.0, 0.2, 0.5);
    float gyroid3Dis = sdGyroid(pos, 8.0, 0.02, 0.8);
    gyroid1Dis += gyroid2Dis * 0.2;
    gyroid1Dis = opSubtraction(gyroid3Dis * 0.8, gyroid1Dis);

    
    // Final distance to return
    float finalDis = opIntersection(sphereDis, gyroid1Dis * 0.8);
    
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
    lightPos.xz += vec2(-abs(sin(iTime)), -abs(cos(iTime))) * 4.0;
    
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
    float downTilt = -0.16;
    vec3 rayOrigin = vec3(0, camHeight, 0); 
    vec3 rayDir = normalize(vec3(uv.x, uv.y + downTilt, 1));
    
    // Assign background
    col = texture(iChannel0, rayDir + vec3(0, 0, 0)).rgb;
    
    // Visualise ray marched items
    float dis = rayMarch(rayOrigin, rayDir);
    if (dis < MAXDIS) {
        vec3 pos = rayOrigin + rayDir * dis;
        vec3 normal = getNormal(pos);
        float diffuseLight = getLight(pos, vec3(0, 5, 6));

        // We multiply the shape by the normal to shade it according to height
        vec3 difCol = mix(shadowCol, lightCol, diffuseLight);
        vec3 mixCol = vec3(0.0);
        
        // Colour banding
        vec3 height = normalize(rayOrigin - pos);
        float dif = clamp(dot(normal, height), 0.0, 1.0);
        mixCol = mix(vec3(0.0), vec3(1.0), dif) * difCol;
        
        // col = vec3(mixCol);
        
        // Reflections
        float IOR = 1.45;
        vec3 reflectedRay = reflect(rayDir, normal);
        vec3 refractedRay = refract(rayDir, normal, 1.0 / IOR); 
        float disInterior = rayMarch(pos, refractedRay);
        
        // Chromatic Abberation
        float abberation = 0.02;
        vec3 reflectedTex = vec3(0.0); // texture(iChannel0, refractedRay).rgb;
        refractedRay = refract(rayDir, normal, 1.0 / (IOR - abberation));
        reflectedTex.r = texture(iChannel0, refractedRay).r;
        refractedRay = refract(rayDir, normal, 1.0 / IOR);
        reflectedTex.g = texture(iChannel0, refractedRay).g;
        refractedRay = refract(rayDir, normal, 1.0 / (IOR + abberation));
        reflectedTex.b = texture(iChannel0, refractedRay).b;
        
        float materialOpaqueness = 0.4;
        col = vec3(reflectedTex) + mixCol * materialOpaqueness;
    }
    
    // Output to screen
    fragColor = vec4(col,1.0);
}
