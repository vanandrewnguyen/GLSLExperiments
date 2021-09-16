/*
Van Andrew Nguyen
10/09/21
[Snake Den]

This follows on further experimentataion with truchet tiling. I applied a circular smoothstep by comparing the length of the uv coord
to a constant variable as a radius. This gave a circular motion. I blended the colour to black to smooth out the edges of the tiling.
The resultant effect is one that looks like a snake pit, with constant motion.
*/

#define ZOOM 8.0
#define PI 3.1415926538

// Return normalised float 0->1
float rand(vec2 o) {
    return fract(sin(dot(o, vec2(12.9898, 78.233))) * 43758.5453);
}

// Return random number - noise two to one
float randDot(vec2 o) {
    o = fract(o * vec2(123.45, 456.56));
    o += dot(o, o + 32.21);
    return fract(o.x + o.y); // these are random math operations to make it appear random
}

// Create a truchet pattern (single layer)
float createTruchetPattern(vec2 uv, vec2 gridUV, vec2 cellID) {
    // Grab a random normalised value for each cell
    float n = randDot(cellID); 
    
    // Depending on the random number, we can flip the tile
    if (n < 0.5) { gridUV.x *= -1.0; } 
    
    // Draw diagonal line / circle / shape
    float lineBlur = 0.025;
    float lineWidth = 0.1;
    float radius = 0.5;
    /*
    float dis = abs(gridUV.x + gridUV.y) - lineWidth;                         // straight line from corner to corner
    float dis = abs(abs(gridUV.x + gridUV.y) - 0.5) - lineWidth;              // straight line from middle to middle
    */
    float dis = length(gridUV - 0.5) - radius;                                // gives us a circle at top right corner
    float mask = 1.0 - smoothstep(-lineBlur, lineBlur, abs(dis) - lineWidth);
    dis = length(gridUV + 0.5) - radius;                                      // gives us a circle at bottom left corner
    mask += 1.0 - smoothstep(-lineBlur, lineBlur, abs(dis) - lineWidth);
    
    return mask;
}

// Create truchet layers using pattern script
float createTruchetLayer(vec2 uv, vec2 gridUV, vec2 cellID, float layerNum) {
    float mask;
    
    // Loop through layerNum and add to the mask
    for (float i=0.0;i<layerNum;i+=1.0) {
        mask += createTruchetPattern(uv, gridUV, cellID + i); // add i to offset
    }
    
    return mask;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab UV Coord
    vec2 uv = (fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;
    uv *= ZOOM;
    
    // Declare time and colour
    float t = iTime * 2.0;
    vec3 col = vec3(0.0);
    vec3 snakeColour = vec3(0.75, 0.9, 0.3);
    
    // Divide canvas into tiles
    vec2 gridUV = fract(uv) - 0.5;
    vec2 cellID = floor(uv);
    
    // Get the pattern mask
    float mask = createTruchetLayer(uv, gridUV, cellID, 1.0);
    
    // Visualisation
    // We get the grid corner to draw portions of the snake using trig.
    // Then, we can animate that.
    // 'check' mods x and y and alternates them so it reverses the angle of rotation on alternating tiles
    // For the circular drawing, we can set a radius and ratio using length(uv).
    vec2 gridCorner = gridUV - sign(gridUV.x + gridUV.y + 0.01) * 0.5;
    float angle = atan(gridCorner.x, gridCorner.y); //-PI -> PI from grid corners
    float freq = 0.25;
    float check = mod(cellID.x + cellID.y, 2.0) * 2.0 - 1.0; //-1 to 1
    float flowAnim = sin((check * angle * freq) + t + (length(uv)));
    float radius = 2.0 * ZOOM;
    float circleRatio = clamp(length(uv) / radius, 0.0, 1.0);
    // Finally add the colour and mix it to black (gradient)
    col += snakeColour * mask * flowAnim;
    col = mix(col, vec3(0.0), (0.1 * sin(t)) + mask * flowAnim); 
    
    // Visualise circle
    col -= smoothstep(0.2, 0.6, circleRatio);
    
    // Output to screen
    fragColor = vec4(col,1.0);
}