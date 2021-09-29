/*
Van Andrew Nguyen
30/09/21
[Meatballs]

I played with displacing a sphere's SDF using superimposed noise and sine waves. I touched on colour banding dependant on camera angle to give
the rocks a low red sheen. Then, I applied a repetition operation to repeatedly draw the balls. The fog came later because I needed to 
keep the fps above 30 frames; rendering everything was way too much strain. Resultant effect is a cool floaty POV through rows and rows of meatballs.
*/

// Globals //////////////////////////////////////////////////////////////////////

/*
Max steps is the maximum amount of steps we can make to get a ray length
Max distance is the maximum distance we can shoot a ray
Surf distance is the distance we loop to get to register a hit on the surface
*/
#define MAXSTEPS 100
#define MAXDIS 42.0
#define SURFDIS 0.01

// Noise //////////////////////////////////////////////////////////////////////

// Return a randomised float using a vec2 (uv coord)
float hash21(in vec2 uv) {
    // Pseudo-random math function to get a random number
    vec2 o = fract(uv * vec2(5859.56, 190.123));
    o += dot(o, o + 54.89);
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

// Rotation Matrix
mat2 rotate(float a) {
    float s = sin(a);
    float c = cos(a);
    return mat2(c, -s, s, c);
}

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
float sdSphere(vec3 pos, float rad) {
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
float opUnion(float d1, float d2) { return min(d1,d2); }
float opSubtraction(float d1, float d2) { return max(-d1,d2); }
float opIntersection(float d1, float d2) { return max(d1,d2); }

float opDisplacement(vec3 pos, float d1, int layers) {
    float d2 = 0.1 * sin(pos.x * 4.0) * sin(pos.y * 4.0) * sin(pos.z * 4.0); 
    d2 += 0.5 * (fbm(pos.xy, layers) + fbm(pos.yz, layers) + fbm(pos.zx, layers));
    return d1 + d2;
}

// Create planet
float sdPlanet(vec3 pos, float rad, float t) {
    pos.xz *= rotate(t);
    pos.xy *= rotate(t);
    float sphereDis = sdSphere(pos, rad); 
    sphereDis = opDisplacement(pos, sphereDis, 4);
    
    return sphereDis;
}

float opRep(vec3 pos, vec3 spacing, float r, float t) {
  vec3 q = mod(pos + 0.5 * spacing, spacing) - 0.5 * spacing;
  return sdPlanet(q, r, t);
}

// Vectors //////////////////////////////////////////////////////////////////////

// Get distance function
float getDist(vec3 pos) {
    
    // Get dist of the various shapes
    vec3 spherePos = pos - vec3(0, -0.5, 12);
    float planetDis = 0.0; //sdPlanet(spherePos, 2.0, iTime * 0.1);
    
    // Parameters
    float rad = 2.0;
    float spacing = 8.0;
    float t = iTime * 0.25;
    
    vec3 finalPos = pos - vec3(4.0, -4.0 + t * 2.0, 0.0);
    
    planetDis = opRep(finalPos, vec3(spacing), rad, (t) + (pos.z / spacing));
    
    // Final distance to return
    float finalDis = planetDis * 0.5;
    
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
    
    vec3 specCol = vec3(0.8, 0.1, 0.2);
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
    vec3 redCol = vec3(0.8, 0.1, 0.2);
    float t = iTime * 4.0 + 2.0 * sin(iTime * 0.5);
    
    // Setup Camera
    float camHeight = 1.8;
    float downTilt = -0.2;
    vec3 rayOrigin = vec3(0, camHeight, t); // this is the camera (origin of vector)
    vec3 rayDir = normalize(vec3(uv.x, uv.y + downTilt, 1));
    
    // Visualise 
    float dis = rayMarch(rayOrigin, rayDir);
    if (dis < MAXDIS) {
        vec3 pos = rayOrigin + rayDir * dis;
        vec3 normal = getNormal(pos);
        float diffuseLight = getLight(pos, vec3(0, 5, 6.0 + t));
        
        // We multiply the shape by the normal to shade it according to height
        vec3 difCol = mix(shadowCol, lightCol, diffuseLight);
        vec3 mixCol = vec3(0.0);
        
        // Height-based lighting (based on normal vector TO camera)
        vec3 height = normalize(rayOrigin - pos);
        float dif = clamp(dot(normal, height), 0.0, 1.0);
        mixCol = mix(shadowCol, lightCol, dif) * difCol;
        
        mixCol = mix(mixCol, redCol, 0.2);
        
        // Colour banding
        if (dif < 0.4) {
            mixCol = mix(mixCol, redCol, clamp(dif - diffuseLight, 0.0, 1.0));
        }
        
        col = vec3(mixCol);
    }
    
    // Fog
    float fogStart = 0.0;
    float fogEnd = MAXDIS;
    col = mix(col, background(rayDir), smoothstep(fogStart, fogEnd, dis));
    
    // Output to screen
    fragColor = vec4(col,1.0);
}