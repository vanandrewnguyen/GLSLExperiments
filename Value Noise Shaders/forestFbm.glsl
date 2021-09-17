/*
Van Andrew Nguyen
Van Andrew Nguyen
16/10/21
[Forests!]

I played with the cloud shader using fbm and a 3D camera. Playing with amplitude and frequency you can change the shape of
the noise to more erratic contours - which somewhat resemble trees. Now keep in mind this is super rough; there's probably so much
more you can do. Honestly this attempt is a little scuffed but I had fun playing with values.
*/

#define ZOOM 4.0
#define CLOUDLAYER 4.0
#define CLOUDDENSITY 120.0

// Return a randomised float using a vec2 (uv coord)
float hash21(in vec2 uv) {
    // Pseudo-random math function to get a random number
    vec2 o = fract(uv * vec2(456.456, 123.123));
    o += dot(o, o + 45.45);
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
    float amp = 0.25;
    float freq = 8.0;
    
    // Now loop through layers and return the combined value
    for (float i=0.0;i<CLOUDLAYER;i+=1.0) {
        val += amp * noise(freq * pos);
        val += 0.1 * sin(freq * pos.x) * sin(freq * pos.y); //add onto the 
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
    vec2 uv = (fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;
    uv *= ZOOM;
    
    float t = iTime * 0.5;
    
    // Declare colours
    vec3 col = vec3(1.0);
    vec3 backCol = vec3(0.8, 0.6, 0.8);
    vec3 baseCloudCol = vec3(0.2, 0.3, 0.1).bgr; 
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
    float intensity = 0.7;
    col -= intensity * rd.y;
    
    // Loop through to draw the clouds
    for (float i=CLOUDDENSITY;i>0.0;i-=1.0) {
        vec3 pos = scroll + 0.05 * i * rd;
        float f = pos.y - 1.2 * fbm(0.6 * pos.xz);
        float density = -f;
        // Now we mix the colours together, using the base colour + inverse of cloudCol
        if (density > 0.0) {
            // Here we use 1-(density*col) instead of (1-density)*col
            // because we don't want our colours to ever reach black. 
            float minVal = density * 0.1; // determines low points
            vec3 finalCloudCol = (1.0 - density) * baseCloudCol;
            col = mix(col, finalCloudCol, min(1.0, minVal));
            
            /*
            if (density < 0.01) { col = vec3(0.0); } // <- this gives contour lines
            if (density < 0.01) {
                col = mix(col, vec3(1.0, 0.7, 0.8), i / CLOUDDENSITY * 0.4);
            }
            */
        }
    }
    
    
    // Output to screen
    fragColor = vec4(col,1.0);
}