/*
Van Andrew Nguyen

This is a basic voronoi shader! I followed a tutorial and played with adding new paramters. 
*/

// Function takes vec2 input and outputs random vec2 between 0 and 1.
vec2 giveRandomVec2(vec2 p) {
    // We take output and 
    vec3 o = fract(p.xyx * vec3(123.34, 234.34, 345.56));
    o += dot(o, o+34.45);
    return fract(vec2(o.x*o.y, o.y*o.z));
}

// Main 
void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
    // Declare uv coords
    vec2 uv = fragCoord/iResolution.xy;
    // We times the ratio to make the coords square, the move the origin to the centre.
    uv.x *= iResolution.x / iResolution.y;
    uv.xy -= 0.5;
    
    // Generate random moving points
    float m = 0.0;
    float t = iTime*0.1;
    float minDis = 100.0;
    float cellNumber = 48.0;
    float cellIndex = 0.0;
    
    // Loop through the points
    for (float i=0.0;i<cellNumber;i++) {
        
        // We use the index number to get a random number from 0 - 1.
        vec2 index = giveRandomVec2(vec2(i));
        
        // Make it a moving position
        vec2 point = sin(index*t);
        
        // Calculate distance to draw a circle
        float dis = length(uv - point);
        float radius = 0.05;
        float blurRadius = 0.025;
        m += smoothstep(radius, radius - blurRadius, dis);
        
        // Calculate distances (min, max)
        if (dis < minDis) {
            minDis = dis;
            cellIndex = i;
        }
    }
    
    // Tint colour green and adjust intensity of light
    float intensity = 3.0 + 0.5 * sin(iTime);
    float lightness = 0.25;
  
    //vec3 col = vec3(m); <- dot representation
    //vec3 col = vec3(cellIndex / cellNumber); <- flat cell representation
    //vec3 col = vec3(0.6, 0.3, 0.55)*(intensity*minDis); <- final gradient representation
    // We choose the final representation
    vec3 col = vec3(0.6, 0.3, 0.55)*(intensity*minDis);
    
    // Adjust colour using individual index of the cells.
    col.rb += (cellIndex / cellNumber) * lightness; // Amplifying red and blue
    col.g  -= (1.0 - (cellIndex / cellNumber)) * lightness; // Decreasing green to saturate purple more
    
    // Draw final colour
    fragColor = vec4(col, 1.0);
}
