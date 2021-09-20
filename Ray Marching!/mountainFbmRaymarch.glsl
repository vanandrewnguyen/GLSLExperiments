
// Global Constants
#define ZOOM 2.0

#define CLOUDLAYER 4.0
#define CLOUDDENSITY 120.0

#define PI 3.141592653589793

#define MAXSTEPS 100
#define MAXDIS 100.0
#define SURFDIS 0.01

// Get distance function
float getDist(vec3 pos) {
    // Here we are using a sphere as collision checker
    vec3 spherePos = vec3(0, 1, 6); // x , y , z
    int sphereRad = 1;
    vec4 sphere = vec4(spherePos, sphereRad);
    
    // Get dist of the sphere AND plane, then return smaller of the two
    float sphereDis = length(pos - sphere.xyz) - sphere.w; // 4th argument
    float planeDis = pos.y; // plane distance is just camera height; since the plane is even anyways
    float finalDis = min(sphereDis, planeDis);
    
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

// Return a randomised float using a vec2 (uv coord)
float hash21(in vec2 uv) {
    // Pseudo-random math function to get a random number
    vec2 o = fract(uv * vec2(567.567, 123.123));
    o += dot(o, o + 89.89);
    return fract(o.x * o.y);	
}

// Noise function which returns a random float using 4 vertices of a cell
float noise(in vec2 uv) {
    vec2 gridUV = fract(uv);
    vec2 cellID = floor(uv);
    
    // We get the cellID + the corners of a square, then mix it all using a ramp
    vec2 ramp = smoothstep(0.0, 1.0, gridUV);
    
    /*
    Depending on that smoothing function between 0->1 using gridUV, you can change how
    smooth the noise is, how it is interpolated.
    smoothstep(0.0, 1.0, gridUV) -> basic smoothing
    step(0.0, 1.0, gridUV) -> blocky, since there is no smoothing
    */
    
    float bl = hash21(cellID + vec2(0, 0));
    float br = hash21(cellID + vec2(1, 0));
    float tl = hash21(cellID + vec2(0, 1));
    float tr = hash21(cellID + vec2(1, 1));
    float bottom = mix(bl, br, ramp.x);
    float top = mix(tl, tr, ramp.x);
    
    float value = mix(bottom, top, ramp.y); 
    return value;
}

// Fractal Brownian Motion
float fbm(in vec2 pos) {
    // Declare return value, amplitude of motion, freq of motion
    float val = 0.0;
    float amp = 0.65;
    float freq = 2.0;
    
    // Now loop through layers and return the combined value
    for (float i=0.0;i<CLOUDLAYER;i+=1.0) {
        val += amp * noise(freq * pos); 
        // For each layer we always want to half amplitude and increase freq (See Book of Shaders)
        amp *= 0.5;
        freq * 2.0;
    }
    
    
    /*
    Here we can change what kind of shape the clouds take. If we replace val += ... noise()
    we can get different shapes.
    e.g.
    val = 0.5 * sin(iTime + freq * pos.x); gets sine waves
    val = 0.5 * sin(iTime + 4.0 * pos.x) * sin(iTime + 4.0 * pos.y);
    */
    
    return val;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab UV Coord
    vec2 uv = (2.0 * fragCoord - iResolution.xy) / iResolution.y;
    //(fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;
    uv *= ZOOM;
    float t = iTime;
    
    // Declare colours
    vec3 col = vec3(0.0);
    vec3 backCol = vec3(0.9, 0.65, 0.55);
    vec3 cloudCol = vec3(0.6, 0.8, 0.99).bgr; 
    vec3 shadowCol = vec3(0.11, 0.08, 0.05);
    col = backCol;
    
    // Setup Camera
    float camHeight = 1.6;
    float camDownTilt = -0.4;
    vec3 scroll = vec3(0, camHeight, t); // move up, scroll z depth
    vec3 angle = scroll + vec3(0, camDownTilt, 1); // tilt camera down
    
    vec3 camZ = normalize(angle - scroll);
    vec3 camX = normalize(cross(vec3(0, 1, 0), camZ));
    vec3 camY = cross(camZ, camX);
    vec3 rayDir = normalize(uv.x * camX + uv.y * camY + 2.0 * camZ); // we get the normal vector of the camera
    
    // Get sky gradient
    float intensity = 0.5;
    col -= intensity * rayDir.y;
    
    // Sun
    col += smoothstep(0.4, 0.35, length(vec2(uv.x, uv.y - 0.8))) * vec3(1.0, 0.0, 0.0);
    
    // Distance from camera to point
    float dis = 0.0;
    
    // Loop through to draw the clouds
    for (float i=CLOUDDENSITY;i>0.0;i-=1.0) {
        // We grab the position and elevate it using density
        vec3 pos = scroll + 0.05 * i * rayDir;
        vec2 noiseBase = 0.6 * pos.xz;
        float f = pos.y - 1.2 * fbm(noiseBase);
        float density = -f;
 
        // Now we mix the colours together, using the base colour + inverse of cloudCol
        if (density > 0.0) {
            // Here we use 1-(density*col) instead of (1-density)*col
            // because we don't want our colours to ever reach black. 
            float minVal = density * 0.4; // determines low points
            col = mix(col, 1.0 - density * cloudCol, min(1.0, minVal));
            
            /*
            Mixing - we can get a base layer!!! 
            If density is smoothstepped we can get a height value to apply colour banding to :))
            So here I'm applying terrain underneath the cloud/fog layer.
            I can play with the bounds according to height -> less fog on higher areas.
            */
            // Mountain Brown Terrain
            float lowerBound = 0.15;
            float upperBound = 0.4;
            float threshold = smoothstep(lowerBound, upperBound, density);
            col = mix(col, vec3(0.3, 0.2, 0.15), min(1.0, threshold));
            
            // Grass / Snow / Toppings / Bottoms
            // The fbm is using the same uv coords. Cutting it off at a threshold
            // gives you the top half which can be useful!
            float mountainHeight = fbm(noiseBase);
            float snowTheshold = 0.9;
            float waterThreshold = 0.2;
            float grassThreshold = 0.45;
            float stepper = mountainHeight * smoothstep(0.1, 0.5, threshold);
            if (mountainHeight > snowTheshold) {
                col = mix(col, vec3(0.74, 0.82, 0.95), stepper); 
            }
            if (mountainHeight < waterThreshold) {
                col = mix(col, vec3(0.3, 0.2, 0.15), stepper); //vec3(0.27, 0.68, 0.71)
            }
            if (mountainHeight > waterThreshold && mountainHeight < grassThreshold) {
                col = mix(col, vec3(0.46, 0.60, 0.47), stepper); 
            }
            
            // Fog Edge Tint
            float edgeFogRatio = fbm(noiseBase) * 0.08; //fbm(pos.xz * 0.005);
            col = mix(col, vec3(1.0, 0.5, 0.4), edgeFogRatio);
            
            // Darken the edges
            float edgeShadowRatio = fbm(pos.xy * 0.04);
            col = mix(col, shadowCol, edgeShadowRatio);
            
            // Darken the bottom using mountain height - this is fake lighting
            col = mix(col, 0.75 * col + 0.25 * shadowCol, 1.0 - mountainHeight);
            
            // Now we assign the distance of the mountain surface to the ray marcher instead
            dis = rayMarch(pos - f, rayDir);
        }
    }
    
    // Now, we apply the raymarcher to the mountain using the distance variable
    vec3 cPos = scroll + rayDir * dis;
    float diffuseLight = getLight(cPos, scroll + vec3(0, 2.4, 6. + 6.*sin(iTime)));
    
    vec3 mixLightCol = mix(shadowCol, vec3(1.0), diffuseLight);
    col *= vec3(mixLightCol);
    //col = vec3(mixLightCol); <- normalised view of mountains
    
    // Output to screen
    fragColor = vec4(col,1.0);
}