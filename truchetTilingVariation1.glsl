/*
Van Andrew Nguyen 
10/09/21
[Chessboard]

Following square tiling and hexagonal tiling, I followed a tutorial from "The Art of Code" to truchet tiling.
The idea behind it was to lay a grid full of a simple pattern, then arbitrarily flip some tiles to get a 'maze-like' look. 
In this example I layered multiple truchet tiles to get a more complex pattern, then created a chessboard look with modulo.
Then, I could animate the tiles using an angle variable to draw portions of said tile, and offset it.
The resultant effect is one that looks like a flowing chessboard.
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
    
    vec3 col = vec3(0.0);
    
    // Divide canvas into tiles
    vec2 gridUV = fract(uv) - 0.5;
    vec2 cellID = floor(uv);
    
    float mask = createTruchetLayer(uv, gridUV, cellID, 2.0);
    
    // Visualisation
    float angle = atan((gridUV + 0.5).x, (gridUV + 0.5).y); //-PI -> PI from grid corners
    float freq = 2.0;
    float check = mod(cellID.x + cellID.y, 2.0);
    float flowAnim = angle * freq;
    float offset = cellID.x * 0.8;
    col += mask * sin(flowAnim + iTime + offset);
    col += check;
    
    // Output to screen
    fragColor = vec4(col,1.0);
}