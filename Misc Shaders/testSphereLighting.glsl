/*
Van Andrew Nguyen
01/09/21
[Fake Sphere Lighting]

This is my first shader! I played around with creating a sphere by using the distance from the centre of the UV coords to adjust lightness.
Then, I played with shadows in using uv.x in relation to some parameteres to emulate a light being cast.
*/

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Get uv coords
    vec2 uv = fragCoord.xy / iResolution.xy;
    
    // We move it to the centre, then apply ratio to make viewport coords equal.
    uv -= 0.5;
    uv.x *= iResolution.x / iResolution.y;
    
    // Declare variables for length, radius, edges
    float d = length(uv);               // The distance of the pixel from the centre
    float r = 0.4;                      // The radius of the circle
    float c = smoothstep(r, r-0.02, d); // The smoothing function applied to the edge
    
    // We are able to move the centre of light shifting around
    // This is to mimck a real light source being waved around the mock 'sphere'
    if (d > r) { c -= 0.1; }
    
    // Declare variables for the moving shadow, and curved edges on the circle.
    float varyingShadow = 0.2 + 0.1*sin(iTime);
    float centreTop = -smoothstep(0.0, 1.0, c*uv.y - 0.1);
    float centreBottom = -smoothstep(0.0, 1.0, -c*uv.y - 0.1);
    
    // Finally output the colour.
    fragColor = vec4(vec3(c * uv.x + varyingShadow + centreTop + centreBottom), 1.0);
}