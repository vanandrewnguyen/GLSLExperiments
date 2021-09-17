/*
Van Andrew Nguyen
16/10/21
[Clouds!]

I looked at Inigo Quilez's cloud shader tried to recreate my own, much simpler version. This one uses fbm much like 
the valueNoise shaders I've written before, only that it creates a 3D camera which gives the smooth noise 'height'.
I am really happy with how this turned out. The resultant effect is a scrolling cloud scene, which is the effect I've 
always wanted to do since day 1.
*/

// Global Constants
#define ZOOM 4.0
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
    float amp = 0.5;
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
    vec2 uv = (fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;
    uv *= ZOOM;
    
    // Declare colours
    vec3 col = vec3(1.0);
    vec3 backCol = vec3(0.45, 0.6, 0.8);
    vec3 cloudCol = vec3(0.6, 0.8, 0.99).bgr; 
    col = backCol;
    
    // Setup Camera
    float camHeight = 1.6;
    float camDownTilt = -0.4;
    vec3 scroll = vec3(0, camHeight, iTime); // move up, scroll z depth
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
        }
    }
    
    // Output to screen
    fragColor = vec4(col,1.0);
}