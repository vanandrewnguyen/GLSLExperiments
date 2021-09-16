/*
Van Andrew Nguyen
07/09/21
[Bubble]

Variation two of the noise function. This time I wanted to overlay it on some kind of shape.
The shape is drawn from this resource: https://thebookofshaders.com/11/
The smooth noise is draw onto the shape like a mask.
*/

#define PI 3.1415926538

// Return random normalised float 0->1
float rand(vec2 o){
    return fract(sin(dot(o, vec2(12.9898, 78.233))) * 43758.5453);
}

// Return random vec2
vec2 random2(vec2 o){
    o = vec2(dot(o, vec2(123.4,456.7)), dot(o, vec2(987.6,543.2)));
    return -1.0 + 2.0 * fract(sin(o) * 12345.6);
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

// Draw a shape with ripples 
float drawShape(vec2 uv, float radius, float blurRad, float intensity) {
    /* 
    UV is the coord system
    Radius is shape radius (normalised)
    blurRad is the tightness of the line before it blurs (normalised)
    intensity is the amount of ripple, higher is lower (aim for 4-8)
    */
    // Move the uv to the centre
    uv = vec2(0.5) - uv;
    // Now we grab the distance of current pixel to the origin (centre of canvas)
    float zoom = 2.0;
    float dis = length(uv) * zoom;
    // Grab the angle using x and y
    float angle = atan(uv.y, uv.x);
    // Determine the variation by getting remainder of (angle / 360) - 180 / intensity
    float timeMult = 2.0;
    float variation = (mod(angle + iTime * timeMult, PI * 2.0) - PI) / intensity;
    // Now we add that variation on to the shape edges
    float f = radius;
    float variationWidth = 12.0;
    f += (sin(angle * variationWidth) * pow(variation, 2.0) * 0.1);
    
    return smoothstep(f, f + blurRad, dis);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab uv coord
    vec2 uv = fragCoord/iResolution.xy;
    uv.x *= iResolution.x / iResolution.y;
    uv.x -= 0.4;
    
    // Establish smooth noise 
    float c = smoothNoise(uv * 8.0 + (1.0 * sin(iTime * 0.25)));
    
    // Visualise using animation
    vec2 pos = vec2(uv * 4.0 + (0.5 * sin(iTime * 0.2)));
    float t = abs(1.0 - sin(iTime * 0.2) * 0.5 + 0.25);
    pos += smoothNoise(pos * 1.0)*t;
    
    float shape = (drawShape(uv, 0.3, 0.25, 4.0) - drawShape(uv, 0.5, 0.04, 5.0));
    float bubbles = smoothstep(0.2, 0.5, smoothNoise(pos * 4.0));
    vec3 col = vec3(1.0) * c * shape;
    col *= bubbles;
    
    // Output to screen
    fragColor = vec4(col, 1.0);
}