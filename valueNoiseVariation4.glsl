/*
Van Andrew Nguyen
08/09/21
[Marble]

I added some parameters to the value noise, and changed the colour.
Added smoothstep to colour variable to cutoff a certain threshold. Resultant effect looks much like marble more than oil.
*/

#define PI 3.1415926538

// Return normalised float 0->1
float rand(vec2 o){
    return fract(sin(dot(o, vec2(12.9898, 78.233))) * 43758.5453);
}

// Return float from two inputs to one output
float noise(vec2 p) {
    // These numbers just need to be completely random
    float amp = 1000.0;
    return fract(sin(p.x * 10.0 + p.y * 1234.5) * amp);
}

// Return smooth noise 
float smoothNoise(vec2 uv) {
    // Create 1D random value
    vec2 index = uv;
    vec2 localUV = fract(index); // we split the canvas into 10x10 cells, frac component
    vec2 cellID = floor(index); // <- ID of current cell we are in, integer component
    
    localUV = localUV*localUV*(3.0 - 2.0 * localUV); // Hermite Curve, basically smoothstep.
    //localUV *= vec2(sin(iTime + uv.x));
    
    // Get noise values for corners of each cell (bottom/top right + left, then mix it)
    float bl = noise(cellID);
    float br = noise(cellID + vec2(1, 0));
    float b = mix(bl, br, localUV.x);
    float tl = noise(cellID + vec2(0, 1));
    float tr = noise(cellID + vec2(1, 1));
    float t = mix(tl, tr, localUV.x);
    float noiseCol = mix(b, t, localUV.y);
        
    return noiseCol;
}

// Fractal Brownian Motion
float fbm(in vec2 uv, int iterations) {
    // Generate noise with initial variables
    float value = 0.0;
    float amp = 0.5;
    float freq = 0.0;
    vec2 shift = vec2(10.0);
    
    // Generate layers using a loop instead
    float gain = 0.5;
    
    for (int i=0;i<iterations;i++) {
        float c = smoothNoise(uv);
        value += amp * abs(c);
        uv = uv * 2.0 + shift;
        amp *= gain;
    }
    
    return value;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab uv coord
    vec2 uv = fragCoord.xy / iResolution.xy;
    uv.x *= iResolution.x / iResolution.y;
    vec3 col = vec3(0.0);
    
    // Create a temp
    vec2 xx = vec2(0.0);
    vec2 yy = vec2(0.0);
    xx.x = fbm(uv, 4);
    xx.y = fbm(uv + vec2(1.0), 4); 
    yy.x = fbm(uv + xx + 0.1 * iTime, 4);
    yy.y = fbm(uv + xx + 0.15 * iTime, 4);
    
    vec2 temp = uv + xx + yy;
    float value = fbm(temp, 4); // we want a random iterating value
    
    // Vis
    col = mix(vec3(0.05, 0.02, 0.1), vec3(0.8, 0.4, 1.0), clamp(value*value, 0.0, 1.0));
    col = mix(col, vec3(1.0, 0.3, 0.5), clamp(length(xx - yy), 0.0, 1.0));
    col += fbm(xx * yy * value * 32.0, 4);
    col += smoothstep(0.2 + 0.1*sin(iTime), 0.5, smoothNoise(xx * yy * 4.0));
    
    // Output to screen
    fragColor = vec4(col,1.0);
}
