/*
Van Andrew Nguyen
01/11/21
[Volumetric Render]

Source blog: https://shaderbits.com/blog/creating-volumetric-ray-marcher
In this shader my end goal was to create realistic clouds. One way outlined was to used volumetric rendering, ray-marching at regular intervals
through a medium and adding absorption at every step. You'd have to ray march for each light within the scene which makes it even slower.
The result is nice and gives a nice opaque look, however once I add noise the ray marcher breaks and freezes. Will have to fix soon. For now,
the scene is just a revolving torus and plane.
*/

// Globals //////////////////////////////////////////////////////////////////////

/*
Max steps is the maximum amount of steps we can make to get a ray length
Max distance is the maximum distance we can shoot a ray
Surf distance is the distance we loop to get to register a hit on the surface
*/
#define MAXSTEPS 100
#define MAXDIS 16.0
#define SURFDIS 0.01
#define PI 3.1415

const vec3 shadowCol = vec3(0.2, 0.1, 0.15);
const vec3 lightCol = vec3(0.8, 0.75, 0.7);
const vec3 fillLightCol = vec3(0.843, 0.976, 0.968);
const vec3 keyLightCol = vec3(0.937, 0.376, 0.101);
const vec3 cloudCol = vec3(0.937, 0.972, 0.984);
const vec3 cloudShadowCol = vec3(0.070, 0.082, 0.109);
const vec3 skyCol = vec3(0.729, 0.917, 0.933);

// Noise //////////////////////////////////////////////////////////////////////

// Flow noise script from https://www.shadertoy.com/view/MtcGRl
vec2 getGradient(vec2 uv, float t) {    
    // Pull random value from a texture rather than generating one
    float rand = texture(iChannel0, uv / 64.0).r;
    
    // Rotate gradient: random starting rotation, random rotation rate
    float angle = 2.0 * PI * rand + 4.0 * t * rand;
    // Rise and run
    return vec2(cos(angle), sin(angle));
}

float tex3DNoise(vec3 pos) {
    vec2 id = floor(pos.xy);
    vec2 f = pos.xy - id;
    vec2 blend = f * f * (3.0 - 2.0 * f); // smoothstep but for a vec2
    float noiseVal = mix(mix(
                            dot(getGradient(id + vec2(0, 0), pos.z), f - vec2(0, 0)),
                            dot(getGradient(id + vec2(1, 0), pos.z), f - vec2(1, 0)),
                            blend.x),
                         mix(
                            dot(getGradient(id + vec2(0, 1), pos.z), f - vec2(0, 1)),
                            dot(getGradient(id + vec2(1, 1), pos.z), f - vec2(1, 1)),
                            blend.x),
                     blend.y);
    return noiseVal / 0.7; 
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

// Torus Distance
float sdTorus(vec3 pos, vec2 rad) {
    float x = length(pos.xz) - rad.x;
    return length(vec2(x, pos.y)) - rad.y;
}

// Sphere Distance
float sdSphere(vec3 pos, float rad) {
    return length(pos) - rad;
}

// Cube Distance
float sdCube(vec3 pos, vec3 size) {
    float val = length(max(abs(pos) - size, 0.0));
    return val;
}

// Vectors //////////////////////////////////////////////////////////////////////

// Get distance function
float getDist(vec3 pos) {
    
    // Get dist of the various shapes
    float planeDis = pos.y + 1.5;

    pos -= vec3(0, 0.5, 6);
    pos.xy *= rotate(iTime);
    pos.zy *= rotate(iTime);
    float torusDis = sdTorus(pos, vec2(1.0, 0.4)); 
    
    // Final distance to return
    float finalDis = opUnion(torusDis, planeDis); //opSmoothUnion(torusDis, sphereDis, 0.5);
    
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
    if (dis < length(lightPos - pos)) {
        float shadowIntensity = 0.2;
        dif *= shadowIntensity;
    }
    
    return dif;
}

// Clouds ////////////////////////////////////////////////////////////////////

float getVolumeOpacity(float volumeDensity, float dis) {
    // Beer Lambert
    return 1.0 / exp(volumeDensity * dis);
}

float lightAttenuation(float len) {
    // Reduces light intensity over distance
    return 1.0 / (1.0 + 0.1 * len + 0.01 * len * len);
}

vec3 getLightExitCol(vec3 colIn, vec3 volumeAlbedo, float opacity) {
    return colIn * (1.0 - opacity) + volumeAlbedo * opacity;
}

float getLightVis(vec3 rayDir, vec3 rayOrigin, float marchSize) {
    float disTravelled = 0.0;
    float lightVis = 1.0;
    float shadowDensity = 0.2;
    
    // March towards the object
    for (int i = 0; i < MAXSTEPS/2; i++) {
        disTravelled += marchSize;
        if (disTravelled > MAXDIS) { break; }
        vec3 pos = rayOrigin + rayDir * disTravelled;
        
        // Check if we are inside the volume
        float disToScene = getDist(pos);
        if (disToScene < SURFDIS) {
            lightVis *= getVolumeOpacity(shadowDensity, marchSize);
        }
    }
    
    return lightVis;
}

// Main //////////////////////////////////////////////////////////////////////

// Background 
vec3 background(vec3 rayDir) {
    vec3 col = vec3(0.0);
    
    float y = rayDir.y * 0.5 + 0.5; 
    col += y * skyCol;
    
    return col;
}


vec3 render(vec3 col, vec3 rayDir, vec3 rayOrigin) {
    vec3 mixCol = col;
    
    // Setup Lights
    vec3 keyLightPos = vec3(2, 4, 5);
    
    // Using modified ray marcher
    float currVis = 1.0;
    float disTravelled = 0.0;
    float marchSize = 0.1;
    float volumeDepth = 0.0;
    float volumeDensity = 0.2;
    float lightVis;
    
    // March towards the object
    for (int i = 0; i < MAXSTEPS; i++) {
        disTravelled += marchSize;
        if (disTravelled > MAXDIS) { break; }
        vec3 pos = rayOrigin + rayDir * disTravelled;
        
        // Check if we are inside the volume
        float disToScene = getDist(pos);
        if (disToScene < SURFDIS) {
            float prevVis = currVis;
            currVis *= getVolumeOpacity(volumeDensity, marchSize);
            float absorption = prevVis - currVis;
            
            // Lighting
            vec3 keyLightDir = pos - keyLightPos;
            lightVis = getLightVis(normalize(keyLightDir), keyLightPos, marchSize);
            vec3 parseKeyLightCol = keyLightCol * lightAttenuation(length(keyLightDir));
            vec3 keyLightExitCol = getLightExitCol(parseKeyLightCol, cloudCol, absorption);
 
            mixCol += keyLightExitCol * absorption * lightVis;
            mixCol += absorption;
        }
    }
    
    mixCol *= lightVis;
    return mixCol;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab UV Coord
    vec2 uv = (fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;

    // Declare Col
    vec3 col = vec3(0.0);

    // Setup Camera
    float camHeight = 1.8;
    float downTilt = -0.2;
    vec3 rayOrigin = vec3(0, camHeight, 0); 
    vec3 rayDir = normalize(vec3(uv.x, uv.y + downTilt, 1));

    vec3 bgCol = background(rayDir);
    col += bgCol; 
    col += render(col, rayDir, rayOrigin);
    
    // Output to screen
    fragColor = vec4(col,1.0);
}