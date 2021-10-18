/*
Van Andrew Nguyen
17/10/21
[Himalayas]

This was an extremely satisfying shader which I have to attribute many credits to Inigo Quilez's article on sampling fbm for SDF's.
https://www.iquilezles.org/www/articles/fbmsdf/fbmsdf.htm
This technique allows us to use fbm in ray marching without getting weird artifacts from pushing the marching field too far.
To get a nice terrain I combined this method (of around 2-3 passes) with a rock texture which I used to mold the ridges.
Coloring the mountain; I learnt new techniques! Instead of applying an even banding (e.g. if height exists within a range, colour this), I 
calculated the steepness of the terrain using the normal and in that was able to give nice even coats of snow and grass to the mountain range.
For lighting I used a warm, key light and a cool ambient light to light up the shadows on the left side.
Resultant effect looks, in my opinion, really nice! Very glad to have come from drawing simple circles to this.
*/

// Globals //////////////////////////////////////////////////////////////////////

/*
Max steps is the maximum amount of steps we can make to get a ray length
Max distance is the maximum distance we can shoot a ray
Surf distance is the distance we loop to get to register a hit on the surface
*/
#define MAXSTEPS 100
#define MAXDIS 4.0
#define SURFDIS 0.01

const vec3 shadowCol = vec3(0.2, 0.1, 0.15);
const vec3 lightCol = vec3(0.8, 0.75, 0.7);
const vec3 dirtCol = vec3(0.2, 0.188, 0.180);
const vec3 snowCol = vec3(0.882, 0.921, 0.937);
const vec3 blueFogCol = vec3(0.729, 0.909, 0.992);
const vec3 grassCol = vec3(0.223, 0.270, 0.207);
const vec3 sunRayCol = vec3(0.980, 0.633, 0.480);

// Return a randomised float using a vec2 
float hash21(vec2 uv) {
    // Pseudo-random math function to get a random number
    vec2 o = fract(uv * vec2(5859.56, 90.123));
    o += dot(o, o + 554.89);
    return fract(o.x * o.y);	
}

// Return a randomised float using a vec3 
float hash31(vec3 p) {
    // Again, random math functions, e.g. dot product of vector with random numbers to 
    // generate 'random' result
	p = vec3(dot(p,vec3(127.1,311.7, 74.7)),
			 dot(p,vec3(269.5,183.3,246.1)),
			 dot(p,vec3(113.5,271.9,124.6)));

	return fract(sin(p.x * p.y * p.z)*43758.5453123);
}


// SDFS //////////////////////////////////////////////////////////////////////

// Shape Operations (From https://www.iquilezles.org/)
float opUnion(float d1, float d2) { return min(d1,d2); }
float opSubtraction(float d1, float d2) { return max(-d1,d2); }
float opIntersection(float d1, float d2) { return max(d1,d2); }
float opSmoothUnion( float d1, float d2, float k ) {
    float h = clamp( 0.5 + 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) - k*h*(1.0-h); }
float opSmoothSubtraction( float d1, float d2, float k ) {
    float h = clamp( 0.5 - 0.5*(d2+d1)/k, 0.0, 1.0 );
    return mix( d2, -d1, h ) + k*h*(1.0-h); }
float opSmoothIntersection( float d1, float d2, float k ) {
    float h = clamp( 0.5 - 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) + k*h*(1.0-h); }

// Rotation Matrix
mat2 rotate(float a) {
    float s = sin(a);
    float c = cos(a);
    return mat2(c, -s, s, c);
}

// Sphere Distance
float sdSphere(vec3 pos, vec3 center, float rad) {
    pos -= center;
    return length(pos) - rad;
}

// Sphere Rand Distance
float sdRandSphere(vec3 cellID, vec3 gridUV, vec3 offset) {
    float rad = 0.5 * hash31(cellID + offset);
    return length(gridUV - vec3(offset)) - rad; 
}

// Repetition of random spheres
float sdRandSphereBase(vec3 pos) { 
    /*
    From: https://www.iquilezles.org/www/articles/fbmsdf/fbmsdf.htm
    We are drawing the distance to the 8 corners of a box (surrounding a sphere)
    of varying sizes according to position
    Then, we layer this using fbm 
    */
    
    vec3 cellID = vec3(floor(pos));
    vec3 gridUV = fract(pos);
    
    float sphereDis = min(min(min(sdRandSphere(cellID, gridUV, vec3(0,0,0)),
                                  sdRandSphere(cellID, gridUV, vec3(0,0,1))),
                              min(sdRandSphere(cellID, gridUV, vec3(0,1,0)),
                                  sdRandSphere(cellID, gridUV, vec3(0,1,1)))),
                          min(min(sdRandSphere(cellID, gridUV, vec3(1,0,0)),
                                  sdRandSphere(cellID, gridUV, vec3(1,0,1))),
                              min(sdRandSphere(cellID, gridUV, vec3(1,1,0)),
                                  sdRandSphere(cellID, gridUV, vec3(1,1,1)))));
    return sphereDis;
}

// Use the sdRandSphereBase to add octaves of noise onto a plane and hopefully get smooth natural terrain
float sdTerrainFbm(vec3 pos, int iterations, float dis) {
    float amp = 1.0;
    
    /*
    We generate a bunch of spheres then grab the ones that intersection within a distance
    Then, we smooth union those spheres with the plane surface
    Then repeat with a decreasing space between the spheres and radii
    */
    
    // Loop through layers
    for (int i=0;i<iterations;i++) {
        // Grab new octave
        float oct = amp * sdRandSphereBase(pos);
        
        // Add layer (min to smooth union, max to intersection first)
        oct = opSmoothIntersection(oct, dis - 0.1 * amp, amp * 0.4);
        dis = opSmoothUnion(oct, dis, amp * 0.3);
        
        // Prep next octave by offsetting position
        pos = mat3(0.0, 1.6, 1.2,
                   -1.6, 0.8,-0.9,
                   -1.2,-0.9, 1.3) * pos;
        // Reduce amplitude per octave since we want more detail with less protusion
        amp *= 0.45;
    }
    
    return dis;
}

// Vectors //////////////////////////////////////////////////////////////////////

// Transform terrain position without affecting ray origin
vec3 transformPos(vec3 pos) {
    pos.z += iTime * 0.1;
    return pos;
}

// Get distance function
float getDist(vec3 pos) {
    // Here we move the terrain instead of the camera; because far away from the origin floating points become inaccurate
    // So to scrap rayOrigin.z += iTime etc we move it here instead
    pos = transformPos(pos);
    
    // Get dist of the various shapes
    float planeDis = pos.y;
    float baseDis = sdTerrainFbm(pos, 2, planeDis);
    // Here we are cheating and using a premade texture (grabbing its position vectors)
    float texAmp = 1.0;
    baseDis += texture(iChannel1, pos.xz / 8.0).r * texAmp;
    
    // Final distance to return
    float finalDis = baseDis * 0.4;
    
    return finalDis;
}

// Return the normal ray
vec3 getNormal(vec3 pos) {
    float dis = getDist(pos);
    vec2 val = vec2(0.01, 0.0);
    
    // To get the slope we give the curve two values super close together
    // Instead of deriving we can do this method via swizzling
    vec3 normal = dis - vec3(getDist(pos-val.xyy), 
                             getDist(pos-val.yxy), 
                             getDist(pos-val.yyx));
    // Same way of doing getDist(p-vec3(0.01, 0.0, 0.0), ... etc
    return normalize(normal);
}

// Ray Marching function
float rayMarch(vec3 rayOrigin, vec3 rayDir) {
    float disFromOrigin = 0.0;
    
    // Now we loop for max steps or until the distance is very small
    for (int i=0;i<MAXSTEPS;i++) {
        // Declare current marching location
        vec3 pos = rayOrigin + rayDir * disFromOrigin;
        float disToScene = getDist(pos);
        disFromOrigin += disToScene;
        
        // Exit condition; we have a hit
        if (disFromOrigin > MAXDIS || disToScene < SURFDIS) { break; }
    }
    
    return disFromOrigin;
}

// Lighting function
float getLight(vec3 pos, vec3 lightOrigin) {
    // Get the light origin
    vec3 lightPos = lightOrigin; // basically a point in 3D space
    //lightPos.xz += vec2(sin(iTime), cos(iTime)) * 4.0;
    
    // Get the light ray
    vec3 light = normalize(lightPos - pos);
    
    // Get normal
    vec3 normal = getNormal(pos);
    
    // Get diffused lighting 
    // We want 0 if the rays are parallel, 1 if the rays are perpendicular
    // Hence we use dot product
    float dif = clamp(dot(normal, light), 0.0, 1.0);
    float dis = rayMarch(pos + normal * SURFDIS * 2.0, light);
    
    // Getting shadows (if disToSphere < disToLight then we have a shadow)
    /*
    if (dis < length(lightPos - pos)) {
        float shadowIntensity = 0.3;
        dif *= shadowIntensity;
    }
    */
    float shadowIntensity = 0.5;
    float shadowBlur = 0.4;
    float stepper = clamp(dis / length(lightPos - pos), 0.0, 1.0);
    dif *= smoothstep(shadowIntensity - shadowBlur, shadowIntensity, stepper);
    dif *= 0.6;
    
    return dif;
}

// Background function; depth perception
vec3 background(vec3 rayDir) {
    vec3 col = vec3(0.0);
    
    float y = rayDir.y * 0.5 + 0.5; // light is top, dark is bottom. 
    col += y * blueFogCol;
    float x = rayDir.x * 0.5 + 0.5;
    col += x * sunRayCol * 0.6;
    
    return col;
}

// Main //////////////////////////////////////////////////////////////////////

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab UV Coord
    vec2 uv = (fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;

    // Declare Col
    vec3 col = vec3(0.0);
    
    // Setup Camera
    float camHeight = 0.15;
    float downTilt = -0.25;
    vec3 rayOrigin = vec3(0, camHeight, 0); // this is the camera (origin of vector)
    vec3 rayDir = normalize(vec3(uv.x, uv.y + downTilt, 1));
    vec3 ambLightPos = vec3(0, 10.0, 1.0);
    vec3 keyLightPos = vec3(2, camHeight, 1.5);
    
    // Visualise 
    float dis = rayMarch(rayOrigin, rayDir);
    if (dis < MAXDIS) {
        vec3 pos = rayOrigin + rayDir * dis;
        vec3 normal = getNormal(pos);
        float ambLight = getLight(pos, ambLightPos);
        float diffuseLight = getLight(pos, keyLightPos); 
        
        // Set base dirt colour
        vec3 mixCol = dirtCol * mix(shadowCol, lightCol, diffuseLight);
        
        // Colour based on steepness
        vec3 localUpDir = normalize(pos);
        float steepness = -dot(normal, localUpDir);
        // Grass and Snow
        mixCol = mix(mixCol, grassCol, smoothstep(0.7, 0.4, steepness));
        mixCol = mix(mixCol, snowCol, smoothstep(0.5, 0.3, steepness)); 
        
        // Lighting
        
        // Ambient Light
        mixCol *= blueFogCol * ambLight;
        // Key Light
        mixCol += sunRayCol * diffuseLight * 1.0;
        
        mixCol = clamp(mixCol, 0.0, 1.0);
        col = vec3(mixCol);
        
    }
    
    // Fog
    float fogStart = 0.0;
    float fogEnd = MAXDIS;
    col = mix(col, background(rayDir), smoothstep(fogStart, fogEnd, dis));
    
    // Output to screen
    fragColor = vec4(col,1.0);
}
