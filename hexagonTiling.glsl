/*
Van Andrew Nguyen
09/09/21

I followed a tutorial from "The Art of Code" in exploring the means to tile not squares, but hexagons.
The resultant effect resembles a glowing beehive. 
*/

#define ZOOM 6.0
#define TIMEMULT 2.0

// Return normalised float 0->1
float rand(vec2 o){
    return fract(sin(dot(o, vec2(12.9898, 78.233))) * 43758.5453);
}

// Calculate distance for hexagon drawing
float hexagonDistance(vec2 uv) {
    // Now, we want to build our hexagon. The vertical lines are easy, just use step to set bounds.
    // For the slanted lines we want to use dot product instead with an angle.
    // dot(vec2 , vec2). e.g. vec(1, 1) gives 45 degree angle. Still need to normalise tho.
    // However a hexagon's slant is not 45 degrees. We use trig to get that ratio.
    
    float c = dot(uv, normalize(vec2(1, sqrt(3.0)))); // the ratio is sqrt(3) by 1
    c = max(c, uv.x); // clamp down the sides
    
    return c;
}

// Return uv coord and id
vec4 hexagonCoord(vec2 uv) {
    // Tile everything
    vec2 gridUVCurrent;
    vec2 rateOfRepetition = vec2(1.0, sqrt(3.0));
    vec2 halfRate = rateOfRepetition * 0.5;
    vec2 gridUVA = mod(uv, rateOfRepetition) - halfRate;
    vec2 gridUVB = mod(uv - halfRate, rateOfRepetition) - halfRate;
    
    // Check which grid current pixel am closest to
    if (length(gridUVA) < length(gridUVB)) {
        gridUVCurrent = gridUVA;
    } else {
        gridUVCurrent = gridUVB;
    }
    
    // Grab an ID for current tile
    float xx = 0.0;
    float yy = hexagonDistance(gridUVCurrent); 
    vec2 id = uv - gridUVCurrent; // difference between the original and new grid gives id
    
    return vec4(xx, yy, id.x, id.y);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab UV Coord and centre
    vec2 uv = (fragCoord - 0.5 * iResolution.xy) / iResolution.y;
    uv *= ZOOM + 1.0 * sin(iTime * TIMEMULT * 0.5);
    
    // Get positive values only
    uv = abs(uv);
    
    // Visualise
    vec3 col = vec3(0.0);
    vec4 hex = hexagonCoord(uv);
    // Change the lower and upper bound of smooth step using sine wave
    float lowerbound = 0.2 + 0.1 * sin(iTime * TIMEMULT);
    float upperbound = 0.8 + 0.1 * cos(iTime * TIMEMULT);
    float stepper = hex.y * sin(hex.z * hex.w + iTime * TIMEMULT);
    float c = smoothstep(lowerbound, upperbound, stepper); // z and w = 3rd, 4th arguments
    
    // Add to the output colour
    col += c;
    col = mix(col, vec3(0.6, 0.6, 0.2), stepper);
    col = mix(col, vec3(0.2, 0.015, 0.01), 1.0 - stepper);
    
    /* Infinite Tiling
    float size=  0.2;
    col += smoothstep(0.1, 0.5, sin(hexagonDistance(uv) * 4.0 + iTime));
    */
    
    // Output to screen
    fragColor = vec4(col,1.0);
}