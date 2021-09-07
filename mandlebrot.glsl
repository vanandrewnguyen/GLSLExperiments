/*
Van Andrew Nguyen
05/09/21

Have been reading up on the fractals and recursive drawing. Found this amazing resource:
https://darkeclipz.github.io/fractals/
Thanks to it I have followed the tutorial code to render the original mandlebrot set in GLSL.
*/

#define PI 3.1415926538
#define ITERATIONS 64.0
#define BOUND 8.0

// Pick a random value from a vec2 and return it
float random (in vec2 st) {
    return fract(sin(dot(st.xy, vec2(12.989,78.233))) * 43758.543);
}

// Return a random vec2
float rseed = 0.;
vec2 random2() {
    vec2 seed = vec2(rseed++, rseed++);
    return vec2(random(seed + 0.342), random(seed + 0.756));    
}

// Generate a colour palette
vec3 palette( in float t, in vec3 a, in vec3 b, in vec3 c, in vec3 d ) {
    return a + b*cos(2.0 * PI * (c * t + d));
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab uv coord
    vec2 uv = (fragCoord - 0.5 * iResolution.xy) / iResolution.y;
    
    // Zoom out and centre the set
    float zoomMult = 2.0;
    uv *= zoomMult + 0.1 * sin(iTime);
    uv.x -= (1.0 / zoomMult);
    
    // Let's establish some vectors
    vec3 col = vec3(1.0);
    vec2 u = vec2(0.0);
    vec2 c = uv;
    float i;
    
    // Iterate through pixels
    for (i=0.0;i<ITERATIONS;i++) {
        u = mat2(u, -u.y, u.x) * u + c; //exponentiation
        // If |z|^2 > BOUND^2 then sqrt(a^2 + b^2) > B, hence the modulus
        // of the complex number is > bound. 
        if (dot(u, u) > BOUND*BOUND) { break; }
    }
    
    // Make ITERATIONS black
    if (i == ITERATIONS) {
        i = 0.0;
    }
    
    // Grab some arbitrary colours and pass into palette function
    // Creates gradient ramp
    vec3 col1 = vec3(0.0, 0.5, 0.5);
    vec3 col2 = vec3(0.0, 0.5, 0.5);
    vec3 col3 = vec3(0.0, 0.5, 0.33);
    vec3 col4 = vec3(0.0, 0.5, 0.66);
    col = palette((i/ITERATIONS), col1, col2, col3, col4);
    
    /*
    Complex numbers can be expressed as a vec2 variable. E.g. u = [a, b] = a + b*I;
    u + v                  // addition       u+v (vector addition)
    mat2(u, -u.y, u.x) * v // multiplication u*v (matrix multiplication)
    mat2(u, -u.y, u.x) * u // exponentiation u^2 (matrix multiplication)
    
    mat2(-/-) is the same as <<a, -b>|<b, a>>
    */

    // Output to screen
    fragColor = vec4(col,1.0);
}
