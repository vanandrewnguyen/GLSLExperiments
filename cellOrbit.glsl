/*
Van Andrew Nguyen
13/09/21
[Cell Membrane 2] 

I worked on the original cell.glsl code which only displayed one layer, and 
put all of the drawing code within its own function. I swizzled the col values and mask values through a vec4 and drew layers instead.
I could manipulate a z-depth and so fake a camera moving through a space. Other than that, some tweaking with colours and movement 
made the final result one I am really proud of!
*/

#define ZOOM 5.0

// Return normalised float 0->1
float rand(vec2 o){
    return fract(sin(dot(o, vec2(12.9898, 78.233))) * 43758.5453);
}

// Create a layer of cells
vec4 createCellLayer(vec2 uv, float sizeMult, float t) {
    // Setup grid UV
    vec2 gridUV = fract(uv) - 0.5;
    vec2 cellID = floor(uv);
    float m = 0.0;
    vec3 mergeCol = vec3(0.0);
    
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
            vec3 col1 = vec3(0.4, 0.9, 0.8);
            vec3 col2 = vec3(0.85, 0.5, 0.65);
            // We use the size of the layer (in relation to the camera) to intensify the ratio
            mergeCol = mix(col1, col2, colRatio / sizeMult);
            
            // Finally visualise it
            m += smoothstep(radius, radius * 0.8, dis);
            
        }
    }
    
    return vec4(m, mergeCol);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab UV Coord
    vec2 uv = (fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;
    uv *= ZOOM;
    
    // Declare colour
    vec3 col = vec3(0.0);
    
    // Draw shape within cells
    float t = iTime * 4.0;
    vec4 m;
    float layerNum = 2.0;
    
    // Create rotation matrix
    float s = sin(t * 0.01);
    float c = cos(t * 0.01);
    mat2 rot = mat2(c, -s, s, c);
    uv *= rot;
    
    // Loop through layers
    for (float i=0.0;i<=1.0;i+=1.0/layerNum) {
        float zDepth = fract(i + t * 0.025);
        float maxDis = 1.0;
        float size = mix(maxDis, 0.05, zDepth); 
        
        // Fade Threshold
        float fadeThreshold = 0.3; // when it starts fading
        float fade = smoothstep(0.0, fadeThreshold, zDepth) * 
                     smoothstep(1.0, 0.9, zDepth);
        
        // Final visualisation
        m += createCellLayer(uv * size, size, t) * fade;
        m.yzw = mix(m.yzw, vec3(0.0), 0.6); // apply black fade to all cells
    }
    
    // Add the colours from the cells, then apply radial gradient
    col += m.x * vec3(m.y, m.z, m.w);
    col = mix(col, vec3(0.0), length(uv) * 0.2);
    col *= 1.0 - vec3(smoothstep(0.6, 0.9, length(uv / ZOOM)));
    
    // Output to screen
    fragColor = vec4(col,1.0);
}