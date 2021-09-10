/*
van Andrew Nguyen
08/09/21
[Stencil Marble]

I tried playing with creating temp variables for x and y, and using those as parameters. 
Here I was able to get a repeating wave shape, which I morphed into a circle. From there, I wanted to fade out the edges.
Then, I was able to apply intensity which gave it a rough, 'drawn marble' effect.
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

    // Generate xy curve initial values
    float xAmp = 0.35; 
    float xFreq = 4.0;
    float yAmp = 0.45;
    float yFreq = 12.0;
    
    float t = iTime * 2.0 + 0.25*sin(iTime * 2.0); 
    float x = sin(uv.x * xFreq + t) * xAmp;
    float y = sin(uv.y * yFreq + t) * yAmp; 
    
    // Iterate 
    float i = 0.0;
    float iterations = 2.0;
    for (i=0.0;i<iterations;i+=1.0) {
        // These staggered values are honestly just gained through experimentation
        // The stencil effect is achieved through the rand(uv) function
        // This is because it just multiplies everything by a random float on bounds [0, 1]
        x += cos((uv.x * xFreq * i * 0.8) + t * i) * xAmp * i * 2.4;
        y += sin((uv.y * yFreq * i * 0.5) + t * i) * yAmp * i * 2.0;
    }
    // We divide the wave by a number to stop it going off-screen
    x *= xAmp * (2.0 / iterations);
    y *= yAmp * (2.0 / iterations);
    
    // Visualisation
    vec3 col = vec3(1.0);
    float radius = 0.75;
    // Get length(uv) outside of radius and use as gradient curve 
    float patternIntensity = t * 5.0 * (length(uv) / radius); // apply an intensity based on distance from centre
    float gradientOutsideBound = (length(uv) - radius) / radius;
    float stepper = smoothstep(0.0, 0.5, 0.25 * (gradientOutsideBound)); // can be 1.0 - gradient
    // If in circle, draw the pattern, else, we draw the patter * gradient
    if (length(uv) < radius) {
        col += vec3(x * rand(uv * 4.0 * patternIntensity));
        col += vec3(y * rand(uv * 4.0 * patternIntensity));        
    } else {
        col += vec3(x * rand(uv * 4.0 * patternIntensity) * stepper);
        col += vec3(y * rand(uv * 4.0 * patternIntensity) * stepper);
    }
    

    // Output to screen
    fragColor = vec4(col,1.0);
}