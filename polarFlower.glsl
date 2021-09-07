/*
Van Andrew Nguyen
04/09/21

This was made following a tutorial for using polar coordinates in shadertoy. 
Creates a flower which spins and edges which twirl.
*/

#define PI 3.1415926538

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab uv coord
    vec2 uv = ((fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y);
    vec2 colourCoord = vec2(atan(uv.x, uv.y), length(uv));
    
    
    // Shape radius gradient from colourCoord to new uv coords
    float tStagger = 0.45;
    float t = iTime*0.05 + colourCoord.y*(tStagger*sin(iTime));
    uv = vec2(colourCoord.x / (2.0*PI)+0.5+(t), colourCoord.y);
    
    // Create zigzag pattern
    float petalNum = 5.0;
    float xx = uv.x*petalNum;
    float zigzag = min(fract(xx), fract(1.0 - xx));
    float offset = 0.15; // <- stoutness / thickness of flower
    float amplifier = 0.25; // <- length of flower
    float c = smoothstep(0.0, 0.1, zigzag * amplifier + offset - uv.y);
    
    // Final colour
    vec3 col = vec3(0.8 * c + (0.2*sin(iTime*0.5)), 0.3 * c, 0.45)*c; //vec3(colourCoord.x / (2.0*PI) + 0.5);

    // Output to screen
    fragColor = vec4(col, 1.0);
}