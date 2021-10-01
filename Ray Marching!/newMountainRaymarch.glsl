/*
Van Andrew Nguyen
29/09/21
[NEW Raymarched Mountain Terrain]

These mountains are rendered using a single plane. I apply many distortions on this plane to achieve a 'random' 3D terrain effect.

I went back and remade my first ray marched mountain terrain. The previous had many issues with artifacts and edge cases for the noise function.
This new version superimposes a bunch of sine waves at varying frequencies and amplitudes to 'fake' a random terrain.
I use similiar banding techniques to colour the mountain, however this time it's smoothly blended and a LOT faster because of the lack of conditionals.
To top it off I have a moving 3D camera which works wonders.

Resultant effect is not that different however it is faster, cleaner and more accurately lit.
*/

// Globals //////////////////////////////////////////////////////////////////////

/*
Max steps is the maximum amount of steps we can make to get a ray length
Max distance is the maximum distance we can shoot a ray
Surf distance is the distance we loop to get to register a hit on the surface
*/
#define MAXSTEPS 100
#define MAXDIS 12.0
#define SURFDIS 0.01

// Noise //////////////////////////////////////////////////////////////////////

// Return a randomised float using a vec2 (uv coord)
float hash21(in vec2 uv) {
    // Pseudo-random math function to get a random number
    vec2 o = fract(uv * vec2(5859.56, 90.123));
    o += dot(o, o + 554.89);
    return fract(o.x * o.y);	
}

// Return smooth noise 
float smoothNoise(vec2 uv) {
    // Create 1D random value
    vec2 index = uv;
    vec2 localUV = fract(index); 
    vec2 cellID = floor(index); 
    
    localUV = localUV*localUV*(3.0 - 2.0 * localUV); // Hermite Curve
    
    // Get noise values for corners of each cell (bottom/top right + left, then mix it)
    float bl = hash21(cellID);
    float br = hash21(cellID + vec2(1, 0));
    float b = mix(bl, br, localUV.x);
    float tl = hash21(cellID + vec2(0, 1));
    float tr = hash21(cellID + vec2(1, 1));
    float t = mix(tl, tr, localUV.x);
    float noiseCol = mix(b, t, localUV.y);
        
    return noiseCol;
}

// Fractal Brownian Motion 
float fbm(vec2 pos, int iterations) {
    // Declare return value, amplitude of motion, freq of motion
    float val = 0.0;
    float amp = 0.05;
    float freq = 4.0;
    
    // Now loop through layers and return the combined value
    for (int i=0;i<iterations;i++) {
        val += amp * smoothNoise(freq * pos); 
        amp *= 0.5;
        freq * 2.0;
    }
    
    return val;
}


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

// Rotation Matrix
mat2 rotate(float a) {
    float s = sin(a);
    float c = cos(a);
    return mat2(c, -s, s, c);
}

// Vectors //////////////////////////////////////////////////////////////////////

// Get distance function
vec2 getDist(vec3 pos) {
    
    // Get dist of the various shapes
    float planeDis = pos.y;
    float gyroidDis = sdGyroid(pos, 8.0, 0.1, 0.25); // we can add this to top of mountains
    
    float height = 0.0;
    float amp = 0.9;
    float freq = 1.0;
    float iterations = 4.0;
    
    // We superimpose a bunch of 'random' sine waves of varying freq and amp
    for (float i=0.0;i<iterations;i+=1.0) {
        // Create an offset and 2 frequencies
        vec2 offset = pos.xz * 0.4;
        float freq1 = freq * pos.z + 2.0 * pos.x * smoothNoise(offset);
        float freq2 = freq * pos.x;
        // Apply these to the height value
        height += amp * sin(freq1) * sin(freq2); 
        height -= 0.5 * mod(i, 2.0) * amp * sin(freq1 + i * 10.0) * sin(freq2 - i * 10.0); 
        
        amp *= 0.5;
        freq *= 2.0;
    }
    
    planeDis += height;
    
    // Final distance to return (we must decrease the ray march step)
    float finalDis = planeDis * 0.2;
    
    return vec2(finalDis, height);
}

// Return the normal ray
vec3 getNormal(vec3 pos) {
    float dis = getDist(pos).x;
    vec2 val = vec2(0.01, 0.0);
    
    // To get the slope we give the curve two values super close together
    // Instead of deriving we can do this method via swizzling
    vec3 normal = dis - vec3(getDist(pos-val.xyy).x, 
                             getDist(pos-val.yxy).x, 
                             getDist(pos-val.yyx).x);
    // Same way of doing getDist(p-vec3(0.01, 0.0, 0.0), ... etc
    return normalize(normal);
}

// Ray Marching function
vec2 rayMarch(vec3 rayOrigin, vec3 rayDir) {
    float disFromOrigin = 0.0;
    float height = 0.0;
    
    // Now we loop for max steps or until the distance is very small
    for (int i=0;i<MAXSTEPS;i++) {
        // Declare current marching location
        vec3 pos = rayOrigin + rayDir * disFromOrigin;
        float disToScene = getDist(pos).x;
        float height = getDist(pos).y;
        disFromOrigin += disToScene;
        
        // Exit condition; we have a hit
        if (disFromOrigin > MAXDIS || disToScene < SURFDIS) { break; }
    }
    
    // First return arg is distance from origin, 2nd is height of mtn
    return vec2(disFromOrigin, height);
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
    float dis = rayMarch(pos + normal * SURFDIS * 2.0, light).x;
    
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
    
    vec3 sunriseCol = vec3(0.78, 0.90, 0.90); // vec3(0.2, 0.1, 0.15); //
    vec3 specCol = sunriseCol; 
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
    
    vec3 waterCol = vec3(0.13, 0.24, 0.27);
    vec3 lowerFogCol = vec3(0.78, 0.90, 0.90);
    vec3 upperFogCol = vec3(0.93, 0.78, 0.63);
    vec3 snowCol = vec3(0.9, 0.9, 0.95);
    vec3 grassCol = vec3(0.27, 0.33, 0.22);
    vec3 mtnColDark = vec3(0.32, 0.29, 0.26);
    vec3 mtnColLight = vec3(0.76, 0.63, 0.51);
    
    float t = iTime;

    // Setup Camera
    float camHeight = 1.6;
    float downTilt = -0.25;
    vec3 rayOrigin = vec3(0, camHeight, 0. + t); // this is the camera (origin of vector)
    vec3 rayDir = normalize(vec3(uv.x + sin(t * 0.5), uv.y + downTilt, 1));
    
    // Visualise 
    vec2 dis = rayMarch(rayOrigin, rayDir);
    if (dis.x < MAXDIS) {
        vec3 pos = rayOrigin + rayDir * dis.x;
        vec3 normal = getNormal(pos);
        float diffuseLight = getLight(pos, vec3(0, 5, 6. + t));
        float mountainHeight = dis.y + pos.y;
        
        // Gradual lighting curve
        vec3 mixCol = mix(mtnColDark * shadowCol, mtnColLight * lightCol, diffuseLight);
        
        // Grass
        mixCol = mix(grassCol * diffuseLight, mixCol, smoothstep(0.0, 0.4, mountainHeight));
        
        // Water
        mixCol = mix(waterCol, mixCol, smoothstep(0.0, 0.05, mountainHeight));
        
        // Snow
        mixCol = mix(mixCol, snowCol * diffuseLight, smoothstep(0.8, 1.0, mountainHeight));
        
        // Fog
        mixCol = mix(lowerFogCol, mixCol, 0.5 + 0.5 * mountainHeight);
        
        col = vec3(mixCol);
    }
    
    // Fog
    float fogStart = 0.0;
    float fogEnd = MAXDIS;
    col = mix(col, background(rayDir), smoothstep(fogStart, fogEnd, dis.x));
    
    // Output to screen
    fragColor = vec4(col,1.0);
}