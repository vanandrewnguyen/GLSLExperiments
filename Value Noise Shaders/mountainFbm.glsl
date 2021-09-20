/*
Van Andrew Nguyen
18/09/21
[Snowy Mountains!]

With more experimentation with the forest and cloud fbm shaders, I am able to play with how the colour merges using 
density and height of the noise. I can merge colours and cut them off; for example I've added a brown terrain layer, then
added a red tinge to the white fog. I also used a threshold to limit snow to be mixed in the upper half of the density.
The resultant effect is a scrolling view of snow-tipped mountains covered in fog.
*/

// Global Constants
#define ZOOM 2.0
#define CLOUDLAYER 4.0
#define CLOUDDENSITY 120.0
#define PI 3.141592653589793


// Return a randomised float using a vec2 (uv coord)
float hash21(in vec2 uv) {
    // Pseudo-random math function to get a random number
    vec2 o = fract(uv * vec2(456.456, 123.123));
    o += dot(o, o + 45.45);
    return fract(o.x * o.y);	
}

// Easing Function from: https://github.com/glslify/glsl-easings
float backInOut(float t) {
  float f = t < 0.5
    ? 2.0 * t
    : 1.0 - (2.0 * t - 1.0);

  float g = pow(f, 3.0) - f * sin(f * PI);

  return t < 0.5
    ? 0.5 * g
    : 0.5 * (1.0 - g) + 0.5;
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
    float amp = 0.6;
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
    vec3 col = vec3(1.0);
    vec3 backCol = vec3(0.9, 0.65, 0.55);
    vec3 cloudCol = vec3(0.6, 0.8, 0.99).bgr; 
    col = backCol;
    
    // Setup Camera
    float camHeight = 1.6;
    float camDownTilt = -0.4;
    vec3 scroll = vec3(0, camHeight, t); // move up, scroll z depth
    vec3 angle = scroll + vec3(0, camDownTilt, 1); // tilt camera down
    
    vec3 camZ = normalize(angle - scroll);
    vec3 camX = normalize(cross(vec3(0, 1, 0), camZ));
    vec3 camY = cross(camZ, camX);
    vec3 rd = normalize(uv.x * camX + uv.y * camY + 2.0 * camZ); // we get the normal vector of the camera
    
    // Get sky gradient
    float intensity = 0.5;
    col -= intensity * rd.y;
    
    // Sun
    col += smoothstep(0.4, 0.35, length(vec2(uv.x, uv.y - 0.8))) * vec3(1.0, 0.0, 0.0);

    // Loop through to draw the clouds
    for (float i=CLOUDDENSITY;i>0.0;i-=1.0) {
        // We grab the position and elevate it using density
        vec3 pos = scroll + 0.05 * i * rd;
        float f = pos.y - 1.2 * fbm(0.6 * pos.xz);
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
            col = mix(col, vec3(0.0), edgeShadowRatio);

            // Darken the bottom using mountain height - this is fake lighting
            vec3 shadowCol = vec3(0.11, 0.08, 0.05);
            col = mix(col, 0.6 * col + 0.4 * shadowCol, 1.0 - mountainHeight);
        }
    }
    
    // Output to screen
    fragColor = vec4(col,1.0);
}