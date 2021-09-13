/*
Van Andrew Nguyen
13/09/2021
[Cell Membrane]

This was a fast and simple shader which I wrote to simulate the natural movements of microscopic cells. 
I divided the canvas into sections and drew a cirle in each. This circle was morphed in shape and colour using its distance 
from the grid centre and canvas centre. Then, I applied a radial gradient to further the illusion of a lens distortion.
*/

#define ZOOM 5.0

// Return normalised float 0->1
float rand(vec2 o){
    return fract(sin(dot(o, vec2(12.9898, 78.233))) * 43758.5453);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab UV Coord
    vec2 uv = (fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;
    uv *= ZOOM;
    
    vec3 col = vec3(0.0);
    vec3 mergeCol = vec3(0.0);
    
    // Setup grid UV
    vec2 gridUV = fract(uv) - 0.5;
    vec2 cellID = floor(uv);
    
    // Draw shape within cells
    float m = 0.0;
    float t = iTime * 4.0;
    
    // Loop to draw overlap (because we will morph each circle)
    for (float xx=-1.0;xx<=1.0;xx+=1.0) {
        for (float yy=-1.0;yy<=1.0;yy+=1.0) {
            // Grab offset of the nine cells around the centre cell
            vec2 offset = vec2(xx, yy); 
            
            // Then we apply distance to smoothstep to get a circle
            float dis = length(gridUV + offset);
            float ratio = sin(t + length(uv) * 8.0) * 0.5 + 0.5;
            float radius = mix(0.3, 0.4, ratio);
            
            // We merge the outer and inner colours of the cell
            float thickness = 0.06 + 0.02 * sin(t / 4.0 + length(uv)); 
            float colRatio = clamp(dis / length(gridUV) * thickness, 0.0, 1.0);
            mergeCol = mix(vec3(0.4, 0.9, 0.8), vec3(0.85, 0.5, 0.65), colRatio);
            
            // Finally visualise it
            m += smoothstep(radius, radius * 0.75, dis);
            
        }
    }
    
    // Add the colours from the cells, then apply radial gradient
    col += m * mergeCol;
    col = mix(col, vec3(0.0), length(uv) * 0.2);
    col *= 1.0 - vec3(smoothstep(0.6, 0.8, length(uv / ZOOM)));
    
    // Output to screen
    fragColor = vec4(col,1.0);
}