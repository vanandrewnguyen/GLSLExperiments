/*
Van Andrew Nguyen
09/16/2021
[Starfields]

I followed a tutorial once again from "The Art of Code"; https://www.youtube.com/watch?v=rvDo9LvfoVE
This involves tiling circles on varying points, on multiple layers. A 3D camera is used to move through the layers which loop.
The resultant effect is a really pretty starfield with smooth pans and zooms.
*/

#define PI 3.1415926538

// 2D rotation matrix
mat2 Rot(float angle) {
    float s = sin(angle), c = cos(angle);
    // return matrix
    return mat2(c, -s, s, c);
}

// Random number generator 
float returnRandom(vec2 p) {
    p = fract(p * vec2(123.34, 456.78));
    p += dot(p, p + 45.56);
    return fract(p.x * p.y);
}

// Star function
float Star(vec2 uv, float hasFlare) {
    // Create the centre star
    float dis = length(uv); // <- distance of current pixel to centre of canvas (uv origin)
    float rad = 0.05;
    float m = (rad / dis); // <- light centre
    
    // Create the light rays
    float lightRays = max(0.0, 1.0 - abs(uv.x * uv.y * 1000.0));
    m += lightRays * hasFlare;
    uv *= Rot(PI / 4.0); // <- rotate canvas by 45 degrees
    lightRays = max(0.0, 1.0 - abs(uv.x * uv.y * 1000.0));
    m += lightRays * 0.2 * hasFlare;
    m *= smoothstep(1.0, 0.1, dis); // <- be completely faded at .6, start fading at .1
    
    return m;
}

// Starfield Function
vec3 StarLayer(vec2 uv) {
    vec3 col = vec3(0.0);
    
    // Make a bunch of stars
    vec2 gridUV = fract(uv) - 0.5;
    vec2 cellID = floor(uv);
    
    // Loop nine times
    for (int y=-1;y<=1;y++) {
        for (int x=-1;x<=1;x++) {
            vec2 offset = vec2(x, y);
            float randPos = returnRandom(cellID + offset);
            float size = fract(randPos * 50.0);
            vec2 starPos = (gridUV - offset - vec2(randPos - 0.5, fract(randPos * 5.0) - 0.5));
            float star = Star(starPos, smoothstep(.6, 1.0, size) * 0.5);
            vec3 colour = vec3(0.65 - 0.2 * sin(randPos * PI * iTime), 0.33 - 0.2 * sin(randPos * PI * iTime), 0.80);
            
            // Final output
            star *= sin(iTime * 4.0 + (randPos * 2.0 * PI)) * 0.5 + 1.0;
            col += star * size * colour;
        }
    }
    
    return col;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab uv coord and set the origin to centre
    vec2 uv = (fragCoord - 0.5 * iResolution.xy) / iResolution.y;
    
    // Colour
    vec3 col = vec3(0.0);
    float t = iTime * 0.5 + 0.4*sin(iTime);
    
    // Zoom view out (multiply the uv by a number > 1) and rotate
    uv *= 0.25;
    uv *= 2.0 + 1.0 * sin(iTime);
    uv *= Rot(t * 0.4);
    
    // Create a bunch of star layers and fade them out if they are far away
    float numLayers = 5.0;
    for (float i=0.0;i<1.0;i+=1.0/numLayers) {
        float depth = fract(i+t); // <-
        float scale = mix(20.0, 0.5, depth);
        float fade = depth * smoothstep(1.0, 0.75, depth);
        col += StarLayer(uv * scale + i*100.0 - 4.0*sin(iTime * 0.5)) * fade;
    }
    
    // Final output
    fragColor = vec4(col, 1.0);
}



