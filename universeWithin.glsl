/*
Van Andrew Nguyen
11/09/21
[The Universe Within Tutorial]

I followed a tutorial from "The Art of Code"; https://youtu.be/3CycKKJiwis
which demonstrates an incredible effect of drawing lines between points on a set grid. These points are put into layers, 
which are moved with parallax towards the camera endlessly. To add my own twist I added a rough grain and blended the edges 
of the composition.
*/


#define ZOOM 1.0

// Return random normalised float from a vec2
float rand(vec2 point) {
    // Mulitply by large random number
    point = fract(point * vec2(234.56, 567.89));
    // Screw around with math to make it seem random (but it's not)
    point += dot(point, point + 12.34);
    return fract(point.x * point.y);
}

// Return random normalised vec2 from vec2
vec2 rand2(vec2 point) {
    float n = rand(point);
    return vec2(n, rand(point+n));
}

// Get the position of a cellID using random position
vec2 grabPosition(vec2 cellID, vec2 offset) {
    // Set random freq and base amp
    vec2 freq = rand2(cellID + offset);
    float amp = 0.4;
    // Declare a x and y points to return
    float x = sin(iTime * freq.x);
    float y = cos(iTime * freq.y);
    
    return offset + vec2(x, y) * amp;
    
}

// Get perpendicular distance to a line from point
float disToLine(vec2 point, vec2 lineStart, vec2 lineEnd) {
    vec2 pointToStart = point - lineStart;
    vec2 pointToEnd = lineEnd - lineStart;
    
    // Basic projection formula
    float t = clamp(dot(pointToStart, pointToEnd) / dot(pointToEnd, pointToEnd), 0.0, 1.0);
    
    return length(pointToStart - pointToEnd * t);
}

// Draw a line 
float drawLine(vec2 point, vec2 lineStart, vec2 lineEnd, float width, float blurWidth) {
    // Grab distance and apply smoothstep
    float dis = disToLine(point, lineStart, lineEnd);
    float m = smoothstep(width, blurWidth, dis);
    
    // Fade line via length (if it's too far away, we fade it away)
    float upperBound = 1.25; // fade line away at 120% length
    float lowerBound = 0.75; // line is at max intensity at 80% length
    m *= smoothstep(upperBound, lowerBound, length(lineStart - lineEnd));
    
    return m;
}

// Create layer of lines
float createLineLayer(vec2 uv, float t) {
    vec2 gridUV = fract(uv) - 0.5;
    vec2 cellID = floor(uv);
    
    // Declare array
    vec2 pos[9];
    int i = 0;
    
    // Loop through neighbouring cells 
    for (float yy=-1.0;yy<=1.0;yy+=1.0) {
        for (float xx=-1.0;xx<=1.0;xx+=1.0) {
            pos[i] = grabPosition(cellID, vec2(xx, yy)); 
            i++;
        }
    }
    
    // Loop through list
    float line = 0.0;
    float lineWidth = 0.02;
    float blurWidth = 0.001;
    for (int i=0;i<9;i++) {
        // Draw a line from the centre (4th tile) to every other tile
        line += drawLine(gridUV, pos[4], pos[i], lineWidth, blurWidth);
        
        // Add circle at point of line
        float ratio = 0.6 + 0.2 * sin(t + cellID.x * cellID.y);
        vec2 falloff = (pos[i] - gridUV) * 64.0 * ratio;
        float circle = 1.0 / (length(falloff) * length(falloff));
        line += circle;
    }
    
    // Draw overlapping lines (1->3, 1->5, 5->7, 7->3)
    line += drawLine(gridUV, pos[1], pos[3], lineWidth, blurWidth);
    line += drawLine(gridUV, pos[1], pos[5], lineWidth, blurWidth);
    line += drawLine(gridUV, pos[5], pos[7], lineWidth, blurWidth);
    line += drawLine(gridUV, pos[7], pos[3], lineWidth, blurWidth);

    return line;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab UV Coord
    vec2 uv = (fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;
    uv *= ZOOM;
    
    // Set colour and new UV
    vec3 col = vec3(0.0);
    float t = iTime;
    
    // Set up multiple layers
    float layerNum = 4.0;
    float line;
    
    // Create rotation matrix
    float s = sin(t * 0.1);
    float c = cos(t * 0.1);
    mat2 rot = mat2(c, -s, s, c);
    uv *= rot;
    // Loop through layers
    for (float i=0.0;i<=1.0;i+=1.0/layerNum) {
        float zDepth = fract(i + t * 0.1);
        float maxDis = 8.0;
        float size = mix(maxDis, 0.1, zDepth); // larger bound = smaller layer
        
        // Fade in and out (get a function that is 0 at start and end)
        float fadeThreshold = 0.4; // when it starts fading
        float fade = smoothstep(0.0, fadeThreshold, zDepth) * 
                     smoothstep(1.0, 1.0 - fadeThreshold, zDepth);
        line += createLineLayer(uv * size + i * 32.0, t * 4.0) * fade;
    }
    
    // We add grain and mix colours
    vec3 shiftColour = vec3(0.9 + 0.25 * sin(iTime), 0.4, 0.1);
    col = vec3(line) * rand(uv) + vec3(line) * 0.6 * shiftColour;
    col = mix(col, vec3(0.0, 0.0, 0.0), length(uv) * 1.2);
    
    
    // Output to screen
    fragColor = vec4(col,1.0);
}
