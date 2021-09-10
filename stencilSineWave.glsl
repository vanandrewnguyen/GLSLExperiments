/*
Van Andrew Nguyen
06/09/21
[Stencil Sine]

This was my first 'shader' using GLSL without tutorial help! Really proud of how this turned out.
We generate a sine wave using y values. Of course, we can use superposition to impose more sine waves on top; they all add up anyways.
Then, we can apply smoothstep and rand(uv) to generate a smooth random stencil effect on the wave.
*/

#define LINETHICKNESS 0.05

// Return a random normalised float
float rand(vec2 co){
    return fract(sin(dot(co, vec2(12.9898, 78.233))) * 43758.5453);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab UV coord
    vec2 uv = fragCoord.xy / iResolution.xy;
    // Push uv.x to bounds [-1, 1]
    uv *= 2.0;
    uv -= 1.0;
    // Make x in ratio with y
    uv.x *= iResolution.x / iResolution.y;
    
    // Generate curve vis
    float amp = 0.25; // amplitude
    float freq = 8.0; // frequency
    float t = -iTime * 2.0 + 0.25*sin(iTime * 2.0); // time dilation
    float y = sin(uv.x * freq + t) * amp; // y value
    
    // Now we have to iterate through a set number
    float i = 0.0;
    float iterations = 4.0;
    for (i=0.0;i<iterations;i+=1.0) {
        // These staggered values are honestly just gained through experimentation
        // The stencil effect is achieved through the rand(uv) function
        // This is because it just multiplies everything by a random float on bounds [0, 1]
        y += sin((uv.x * freq * i * 0.5 * rand(uv)) + t * i) * amp * i * 2.0;
    }
    // We divide the wave by a number to stop it going off-screen
    y *= amp * (2.0 / iterations);
    
    // Visualisation
    vec3 col = vec3(1.0);
    if (uv.y >= y - LINETHICKNESS && uv.y <= y + LINETHICKNESS) {
        // We smooth out the edges
        float dis = abs(uv.y - y)/LINETHICKNESS;
        col = vec3(smoothstep(0.2, 0.8, dis));
    }
    
    // Output to screen
    fragColor = vec4(col,1.0);
}