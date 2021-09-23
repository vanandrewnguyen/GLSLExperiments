/*
Van Andrew Nguyen
23/09/21
[Raymarched Planet]

Today's little experiment is based on transforming the surface of a sphere using a ray marcher. I applied displacement on the distance function (basic sin wave for xyz) and 
tried colouring the planet. I soon realised I messed up; I based the colours off the diffuse lighting rather than a height map, which is completely wrong.
Still, it looks good and I had fun. Will update more soon.
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

#define FBMLAYER 4.0

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

// Shape Operations (From https://www.iquilezles.org/)
float opUnion(float d1, float d2) { 
    return min(d1,d2); 
}

float opSubtraction(float d1, float d2) { 
    return max(-d1,d2); 
}

float opIntersection(float d1, float d2) { 
    return max(d1,d2); 
}

float opSmoothUnion(float a, float b, float c) {
    // a and b are the different values to be smoothed, c is the smoothness level
    float h = clamp(0.5 + 0.5 * (b - a) / c, 0.0, 1.0);
    return mix(b, a, h) - c * h * (1.0 - h);
}

float opDisplacement(float shapeDis, vec3 pos) {
    float t = iTime * 0.1;
    float amp = 0.15 + 0.05 * sin(t);
    float freq = 4.0 + 0.5 * sin(t);
    float d1 = shapeDis;
    float d2 = amp * sin(freq*pos.x+t) * sin(freq*pos.y+t) * sin(freq*pos.z+t);
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
    
    /*
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
    */
    
    float planetDis = sdSphere(pos, vec3(0, 1, 6), 1.6);
    planetDis = opDisplacement(planetDis, pos);
    float finalDis = planetDis;
    
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
    //lightPos.xz += vec2(sin(iTime * 0.1), cos(iTime * 0.1)) * 4.0;
    
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
    vec3 mixCol = vec3(0.0);

    // Setup Camera
    float camHeight = 1.8;
    float downTilt = -0.2;
    vec3 rayOrigin = vec3(0, camHeight, 0); // this is the camera (origin of vector)
    vec3 rayDir = normalize(vec3(uv.x, uv.y + downTilt, 1));
    
    // Visualise sphere
    float dis = rayMarch(rayOrigin, rayDir);
    vec3 pos = rayOrigin + rayDir * dis;
    float diffuseLight = getLight(pos, vec3(1, 5, 3));//vec3(0, 5, 6));
    
    // Base Colours (we can also band colours since diffuseLight is normalised)
    vec3 shadowCol = vec3(0.2, 0.1, 0.15);
    vec3 lightCol = vec3(0.8, 0.75, 0.7);
    vec3 specCol = vec3(0.9, 0.85, 0.8);
    
    // Terrain Colours
    vec3 cloudCol = vec3(0.99, 0.8, 0.6); 
    vec3 mountainCol = vec3(0.3, 0.2, 0.15);
    vec3 snowCol = vec3(0.95, 0.95, 0.99);
    vec3 sandCol = vec3(0.95, 0.91, 0.74);
    vec3 waterCol = vec3(0.45, 0.6, 0.8);
    vec3 grassCol = vec3(0.46, 0.60, 0.47);
    
    mixCol = mix(shadowCol, lightCol, diffuseLight);
    
    // Terrain Colouring (using diffuse light, but what about dis?)
    float ratio = 0.2 * diffuseLight + 0.8;
    
    if (diffuseLight < 0.2) {
        if (diffuseLight > 0.1) {
            mixCol = mix(mixCol, sandCol * (0.75 + diffuseLight), ratio);
        } else {
            if (dis < 6.0) { mixCol = waterCol * 10.0 * (diffuseLight+0.1); }
        }
    } else if (diffuseLight < 0.8) {
        mixCol = mix(mountainCol * diffuseLight, snowCol, diffuseLight * (0.8 - 0.6));
    } else {
        mixCol = mix(mixCol, snowCol * diffuseLight, ratio);
    }
    
    
    // Sharp two/three tone flat lighting
    /*
    vec3 middleCol = (shadowCol + lightCol) * 0.5;
    if (diffuseLight < 0.2) { mixCol = shadowCol; 
    } else if (diffuseLight < 0.6) { mixCol = middleCol;  
    } else if (diffuseLight > 0.95) { mixCol = specCol; 
    } else { mixCol = lightCol; }
    */
    
    col = vec3(mixCol);
    
    // Output to screen
    fragColor = vec4(col,1.0);
}
