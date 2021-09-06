/*
Van Andrew Nguyen

Variation one of the basic value noise script. I add and subtract different noise textures of varying thresholds (in smoothstep)
to create a bubbly texture (cutout holes and such, etc). It's also animated.
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

// Generate a layered smooth noise texture
float generateLayeredSmootNoise(vec2 uv) {
    // Establish smooth noise but then we can add more layers
    // This makes it more complex. One way is to add more layers and double freq, half amp.
    float c = smoothNoise(uv * 2.0);
    c += smoothNoise(uv * 4.0) * 0.5;
    c += smoothNoise(uv * 8.0) * 0.25;
    c += smoothNoise(uv * 16.0) * 0.125;
    // Normalise value
    c /= 1.875;
    
    return c;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab uv coord
    vec2 uv = fragCoord/iResolution.xy;
    
    // Use the noise script to get a colour function
    vec2 pos = vec2(uv * 4.0);
    float c = generateLayeredSmootNoise(pos);
    
    // Vis
    float t = abs(1.0 - sin(iTime * 0.1) * 0.5 + 0.25);
    pos += smoothNoise(pos * 1.0)*t;
    
    // Get base colour, then add subtract noise functions to get procedural shapes
    // We play with this to make nice looking textures
    vec3 col = vec3(1.0);
    col *= smoothstep(0.1, 0.11, smoothNoise(pos));
    col += smoothstep(0.3, 0.4, smoothNoise(pos * 4.0));
    col -= smoothstep(0.1, 0.2, smoothNoise(pos * 4.0));
    
    // Output to screen
    fragColor = vec4(col,1.0);
}
