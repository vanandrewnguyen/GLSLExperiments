/*
Van Andrew Nguyen
09/16/21
[Heart]

Simple shader using a basic manipulation of x and y coordinates. This time I used the square root of x to get the heart shape.
Then I applied a glow using the y = 1/x but the execution is poor; will fix eventually.
*/


// Draw a heart
float drawHeart(vec2 uv, float radius, float blurRadius) {
    // Stretch coord to form heart
    uv.x *= 0.8;
    uv.y += -sqrt(abs(uv.x)) * 0.4;
    uv.y += 0.1;
    
    // Draw a circle (but with shifted coord it's a heart)
    float dis = length(uv);
    float heart = smoothstep(radius + blurRadius, radius - blurRadius, dis);
    
    return heart;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab UV Coord
    vec2 uv = (fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;

    // Setup col
    vec3 col = vec3(0.0);
    float radius = 0.2 + 0.02 * sin(iTime);
    float blurRadius = 0.01;
    
    float heart = drawHeart(uv, radius, blurRadius);
    vec3 mixCol = mix(vec3(0.9, 0.4, 0.5), vec3(0.3, 0.1, 0.15), length(0.5 * uv / radius));
    
    col = vec3(heart * mixCol);
    
    // Glow
    vec3 glow = clamp(vec3(0.01 / length(uv)), 0.0, 1.0);
    if (col.r < 0.1) {
        col = mix(col, glow * 10.0 * mixCol, radius + length(uv * 0.5));
    }
    
    // Output to screen
    fragColor = vec4(col,1.0);
}