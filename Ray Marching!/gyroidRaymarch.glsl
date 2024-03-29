/*
Van Andrew Nguyen
24/09/21
[Canopy]

I followed this Art of Code tutorial to create an endless gyroid structure: https://youtu.be/-adHIyjIYgk
Essentially, we create an endlessly repeating shape (a gyroid) and render it using ray marching. To alter the gyroid and shape it into branches,
we can thin out the structure and layer multiple gyroids using boolean operations (subtraction, intersection) to create new shapes.
Then, the final touch is to colour the gyroid and add fog with a 3D camera.
The resultant effect is this incredible scroll through twisting branches.
*/


// Globals //////////////////////////////////////////////////////////////////////

/*
Max steps is the maximum amount of steps we can make to get a ray length
Max distance is the maximum distance we can shoot a ray
Surf distance is the distance we loop to get to register a hit on the surface
*/
#define MAXSTEPS 100
#define MAXDIS 8.0
#define SURFDIS 0.01


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

// Rotation Matrix
mat2 rotate(float a) {
    float s = sin(a);
    float c = cos(a);
    return mat2(c, -s, s, c);
}

// Vectors //////////////////////////////////////////////////////////////////////

// Get distance function
float getDist(vec3 pos) {
    
    // Get dist of shapes
    vec3 cubePos = pos - vec3(0, 1, 4);
    //cubePos.xy *= rotate(iTime);
    float cubeDis = sdCube(cubePos, vec3(0.8)); 
    float gyroid1Dis = sdGyroid(pos, 4.0, 0.03, 1.5); // used as base
    float gyroid2Dis = sdGyroid(pos, 10.0, 0.03, 0.3); // used as bump map
    float gyroid3Dis = sdGyroid(pos, 15.0, 0.02, 0.1); // used as cutout
    
    // Final distance to return
    // Here we multiply gyroid distance by a number < 1.0 because we are decreasing step length in order to not overstep
    gyroid1Dis += gyroid2Dis * 0.15;
    gyroid1Dis = opSubtraction(gyroid3Dis, gyroid1Dis);
    float finalDis = gyroid1Dis * 0.8; //opIntersection(cubeDis, gyroid1Dis * 0.6); // 
    
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
        if (disFromOrigin > MAXDIS || abs(disToScene) < SURFDIS) { break; }
    }
    
    return disFromOrigin;
}

// Lighting function
float getLight(vec3 pos, vec3 lightOrigin) {
    // Get the light origin
    vec3 lightPos = lightOrigin; // basically a point in 3D space
    //lightPos.xz += vec2(sin(iTime), cos(iTime)) * 4.0;
    
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

// Background function; depth perception
vec3 background(vec3 rayDir) {
    vec3 col = vec3(0.0);
    
    vec3 specCol = vec3(0.8, 0.9, 0.1);
    float y = rayDir.y * 0.5 + 0.5; // light is top, dark is bottom. 
    col += y * specCol;
    
    return col;
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
    vec3 specCol = vec3(0.8, 0.9, 0.1);
    float t = iTime * 0.8;
    
    // UV Distortion
    uv += sin(uv * 20.0 + t) * 0.01;
    //uv *= rotate(t / 2.0);
    
    // Setup Camera
    float camHeight = 1.8;
    float downTilt = -0.0;
    vec3 rayOrigin = vec3(0, camHeight, t); // this is the camera (origin of vector)
    vec3 rayDir = normalize(vec3(uv.x, uv.y + downTilt, 1));
    
    // Visualise ray marched items
    float dis = rayMarch(rayOrigin, rayDir);
    if (dis < MAXDIS) {
        vec3 pos = rayOrigin + rayDir * dis;
        vec3 normal = getNormal(pos);
        float diffuseLight = getLight(pos, vec3(1, 4, 2.0 + t));

        // We multiply the shape by the normal to shade it according to height
        vec3 mixCol = mix(shadowCol, lightCol, diffuseLight);
        // Then we darken the shape using the top shadow (normal.y)
        mixCol *= pow(normal.y * 0.5 + 0.5, 2.0);
        
        // Then we can create cheap shadows using the bump map, then highlight it using a smaller smoothstep
        float gyroid2Dis = sdGyroid(pos, 10.0, 0.03, 0.3);
        float gyroid3Dis = sdGyroid(pos, 15.0, 0.02, 0.1);
        mixCol *= smoothstep(-0.1, 0.1, gyroid2Dis); // we darken on either side of 0
        mixCol += smoothstep(-0.01, -0.03, gyroid3Dis) * specCol; // highlight within crack
        
        col = vec3(mixCol);
    }
    
    // Fog
    float fogStart = 0.0;
    float fogEnd = MAXDIS;
    col = mix(col, background(rayDir), smoothstep(fogStart, fogEnd, dis));
    
    // 2D Visualisation
    /*
    col *= 0.0;
    float d = sdGyroid(vec3(uv.x, uv.y, iTime*0.2), 8.0);
    col += abs(d) * 4.0;
    */
    
    // Output to screen
    fragColor = vec4(col,1.0);
}