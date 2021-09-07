/*
Van Andrew Nguyen
06/09/21

This time round I wanted to explore generating and understanding noise functions. Followed this amazing tutorial:
https://youtu.be/zXsWftRdsvU
I'm going to try mess with this to see if I can change the interpolation to get some variation.
*/

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

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab uv coord
    vec2 uv = fragCoord/iResolution.xy;
    
    // Establish smooth noise but then we can add more layers
    // This makes it more complex. One way is to add more layers and double freq, half amp.
    float c = smoothNoise(uv * 8.0);
    c += smoothNoise(uv * 16.0) * 0.5;
    c += smoothNoise(uv * 32.0) * 0.25;
    c += smoothNoise(uv * 64.0) * 0.125;
    // Normalise value
    c /= 1.875;
    
    // Vis
    vec3 col = vec3(c);
    
    // Output to screen
    fragColor = vec4(col,1.0);
}