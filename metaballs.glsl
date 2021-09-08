/*
Van Andrew Nguyen
08/09/21

Visualisation of metaballs in glsl. We split the canvas up and get neighbouring cells.
Then we can grab distance and create an offset. Then after I had the base effect, 
I applied a basic circle mask on it and radial gradient to make it look like a bacteria petri dish.
*/


#define PI 3.1415926538
#define ZOOM 8.0

// Return normalised float 0->1
float rand(vec2 o){
    return fract(sin(dot(o, vec2(12.9898, 78.233))) * 43758.5453);
}

// Return normalised random vec2
vec2 rand2(vec2 o) {
    return fract(sin(vec2(dot(o, vec2(127.1, 311.7)),
                 dot(o, vec2(269.5, 183.3)))) * 43758.5453);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab uv coord
    vec2 uv = fragCoord.xy / iResolution.xy;
    uv.x *= iResolution.x / iResolution.y;
    uv *= ZOOM;
    
    // Temp coords
    vec2 xx = floor(uv); // we split up the canvas into a grid.
    vec2 yy = fract(uv); // fract and floor return integer and fractional component
    float minDis = 0.8;
    vec3 col = vec3(0.0); // colour of balls
    
    // Now we loop from -1 to 1, for both x and y, hence it is 9 times
    for (int i=-1;i<=1;i++) {
        for (int k=-1;k<=1;k++) {
            // Get neighbouring place in grid
            vec2 neighbor = vec2(float(i), float(k));
            // Get a random offset from neighbour
            vec2 offset = rand2(xx + neighbor);
            offset = 0.5 + 0.5*sin(iTime + (2.0 * PI * offset));
            // Grab position of current and set the distance
            vec2 pos = neighbor + offset - yy;
            minDis = min(minDis, minDis * length(pos));
        }
    }

    // Grab colour to draw
    float lowerBound = 0.05 * rand(uv); // mult by rand(uv) to apply random stencil effect
    float upperBound = 0.075;
    col += smoothstep(lowerBound, upperBound, minDis);
    
    // Get centre position of circle, then calculate radius to draw within it
    vec2 centrePos = vec2(0.5, 0.5);
    float circleRadius = 0.4 * ZOOM + 0.05 * sin(iTime* 2.0);
    float circleWidth = 0.15 * ZOOM + 0.025 * ZOOM * sin(iTime);
    vec2 currentPos = uv - vec2(ZOOM * 0.9, ZOOM / 2.0);
    
    // Band it!
    if (length(currentPos) > circleRadius) {
        col = vec3(1.0);
    }
    col = mix(col, vec3(1.0), length(currentPos) / circleRadius * 0.75);

    // Output to screen
    fragColor = vec4(col, 1.0);
}
