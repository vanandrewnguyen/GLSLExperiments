/*
Van Andrew Nguyen
21/09/21
[Raymarching Shapes]

I followed "The Art of Code"'s beginner raymarching tutorial and have created a scene with multiple shapes.
Really fun! I also played with colour banding, since the light value is a normalised float, we can band light by 
pitting it wihtin a threshold and assigning it a solid colour, rather than a gradient.
*/

/*
Max steps is the maximum amount of steps we can make to get a ray length
Max distance is the maximum distance we can shoot a ray
Surf distance is the distance we loop to get to register a hit on the surface
*/
#define MAXSTEPS 100
#define MAXDIS 100.0
#define SURFDIS 0.01

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
    return length(pos - center) - rad;
}

// Cube Distance
float sdCube(vec3 pos, vec3 center, vec3 size) {
    pos -= center;
    
    // We get the length of the position - size (dis) and grab a version > 0
    float val = length(max(abs(pos) - size, 0.0));
    return val;
}

// Get distance function
float getDist(vec3 pos) {
    // Here we are using a sphere as collision checker
    /*
    vec3 spherePos = vec3(0, 1, 6); // x , y , z
    int sphereRad = 1;
    vec4 sphere = vec4(spherePos, sphereRad);
    */
    
    // Get dist of the various shapes
    float sphereDis = sdSphere(pos, vec3(0, 1, 6), 0.5); 
    float planeDis = pos.y; // cam height
    float capsuleDis = sdCapsule(pos, vec3(1, 0.25, 4.5), vec3(0, 0.25, 4), 0.25);
    float torusDis = sdTorus(pos, vec3(0, 0.25, 6), vec2(0.8, 0.2));
    float cubeDis = sdCube(pos, vec3(-2, 0.5, 7), vec3(0.5)); 
    
    // Final distance to return
    float finalDis = min(torusDis, planeDis);
    finalDis = min(finalDis, capsuleDis);
    finalDis = min(finalDis, sphereDis);
    finalDis = min(finalDis, cubeDis);
    
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
    lightPos.xz += vec2(sin(iTime), cos(iTime)) * 4.0;
    
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

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab UV Coord
    vec2 uv = (fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;

    // Declare Col
    vec3 col = vec3(0.0);

    // Setup Camera
    float camHeight = 1.8;
    float downTilt = -0.2;
    vec3 rayOrigin = vec3(0, camHeight, 0); // this is the camera (origin of vector)
    vec3 rayDir = normalize(vec3(uv.x, uv.y + downTilt, 1));
    
    // Visualise sphere
    float dis = rayMarch(rayOrigin, rayDir);
    vec3 pos = rayOrigin + rayDir * dis;
    float diffuseLight = getLight(pos, vec3(0, 5, 6));
    
    // Colours ( we can also band colours since diffuseLight is normalised)
    vec3 shadowCol = vec3(0.2, 0.1, 0.15);
    vec3 lightCol = vec3(0.8, 0.75, 0.7);
    // Gradual lighting curve
    vec3 mixCol = mix(shadowCol, lightCol, diffuseLight);
    
    // Sharp two/three tone flat lighting
    /*
    vec3 mixCol = vec3(0.0);
    vec3 middleCol = (shadowCol + lightCol) * 0.5;
    if (diffuseLight < 0.2) { mixCol = shadowCol; 
    } else if (diffuseLight < 0.6) { mixCol = middleCol;  
    } else { mixCol = lightCol; }
    */
    
    col = vec3(mixCol);
    
    // Output to screen
    fragColor = vec4(col,1.0);
}