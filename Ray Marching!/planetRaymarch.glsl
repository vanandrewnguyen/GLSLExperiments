/*
Van Andrew Nguyen
22/10/2021
[Planet Shader]

This was a fun shader to build, and it's nice to compare with my previous attempts at a planet sdf. I used two sphere sdf's, one with slight
distortion for the water; the other with rough ridges using a texture. The water sphere is ray marched only twice; once for the nearest distance
and another for depth past the water surface. Was nice introducing specular lighting. The terrain itself needed modification on how I'd usually
shade it -> I needed to get the vector based not on the camera but from the planet origin; hence polar coordinates. 
Resultant effect looks really nice with two diffuse lights.
*/

// Globals //////////////////////////////////////////////////////////////////////

/*
Max steps is the maximum amount of steps we can make to get a ray length
Max distance is the maximum distance we can shoot a ray
Surf distance is the distance we loop to get to register a hit on the surface
*/
#define MAXSTEPS 100
#define MAXDIS 32.0
#define SURFDIS 0.01

#define MATTERRAIN 1
#define MATWATER 2
#define WATERHEIGHT 2.29

const vec3 shadowCol = vec3(0.2, 0.1, 0.15);
const vec3 lightCol = vec3(0.8, 0.75, 0.7);
const vec3 darkSpaceCol = vec3(0.090, 0.129, 0.125);
const vec3 orangeSunray = vec3(0.356, 0.196, 0.125);
const vec3 dirtCol = vec3(0.117, 0.117, 0.101);
const vec3 grassCol = vec3(0.293, 0.350, 0.257);
const vec3 snowCol = vec3(0.882, 0.921, 0.937);
const vec3 waterCol = vec3(0.094, 0.224, 0.279);
const vec3 blueLightCol = vec3(0.729, 0.909, 0.992);
const vec3 orangeLightCol = vec3(0.968, 0.596, 0.431);

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

// Credit to ruojake: https://www.shadertoy.com/view/tlKXDz
float noise(vec3 pos) {
	vec3 id = floor(pos);
    vec3 off = smoothstep(0.0, 1.0, pos - id);
    return mix(
        mix(mix(hash31(id), hash31(id + vec3(1, 0, 0)), off.x),
            mix(hash31(id + vec3(0, 0, 1)), hash31(id + vec3(1, 0, 1)), off.x), off.z),
        mix(mix(hash31(id + vec3(0, 1, 0)), hash31(id + vec3(1, 1, 0)), off.x),
            mix(hash31(id + vec3(0, 1, 1)), hash31(id + 1.0), off.x), off.z), off.y);
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

float sdTerrain(vec3 pos) {
    float dis = 0.0; 
    float amp = 0.4;
    float freq = 0.1;
    float iterations = 4.0;
    for (float i=0.0; i<iterations; i+=1.0) {
        dis += texture(iChannel0, pos.xy * freq).r * amp;
        amp *= 0.5;
        freq *= 2.0;
    }
    dis += texture(iChannel1, pos.xy * 0.1).r * 0.05;
    return sdSphere(pos, vec3(0, 0.8, 6), 2.0 + dis); 
}

float sdWater(vec3 pos) {
    float dis = sdSphere(pos, vec3(0, 0.8, 6), WATERHEIGHT); 
    dis += 0.005 * sin(16.0 * pos.x + iTime + (0.5 + 0.5 * sin(iTime)) );
    dis += 0.005 * sin(24.0 * pos.z + iTime * 0.8);
    dis += 0.01 * noise(iTime + pos * 8.0);
    return dis;
}

// Vectors //////////////////////////////////////////////////////////////////////

// Get distance function
vec2 getDist(vec3 pos, bool includeWater) {
    
    // Get dist of the various shapes
    float terrainDis = sdTerrain(pos);
    float waterDis = (includeWater) ? sdWater(pos) : 2.0;
    
    // Final distance to return
    float finalDis = opUnion(terrainDis, waterDis);
    
    // Material
    int mat = 0;
    mat = (finalDis == terrainDis) ? MATTERRAIN : MATWATER;
    
    return vec2(finalDis, mat);
}

// Return the normal ray
vec3 getNormal(vec3 pos, bool includeWater) {
    float dis = getDist(pos, includeWater).x;
    vec2 val = vec2(0.01, 0.0);
    
    // To get the slope we give the curve two values super close together
    // Instead of deriving we can do this method via swizzling
    vec3 normal = dis - vec3(getDist(pos-val.xyy, includeWater).x, 
                             getDist(pos-val.yxy, includeWater).x, 
                             getDist(pos-val.yyx, includeWater).x);
    // Same way of doing getDist(p-vec3(0.01, 0.0, 0.0), ... etc
    return normalize(normal);
}

// Ray Marching function
vec2 rayMarch(vec3 rayOrigin, vec3 rayDir, bool includeWater) {
    float disFromOrigin = 0.0;
    int mat = 0;
    
    // Now we loop for max steps or until the distance is very small
    for (int i=0;i<MAXSTEPS;i++) {
        // Declare current marching location
        vec3 pos = rayOrigin + rayDir * disFromOrigin;
        vec2 passedDMat = getDist(pos, includeWater);
        float disToScene = passedDMat.x;
        mat = int(passedDMat.y);
        disFromOrigin += disToScene;
        
        // Exit condition; we have a hit
        if (disFromOrigin > MAXDIS || disToScene < SURFDIS) { break; }
    }
    
    return vec2(disFromOrigin, mat);
}

// Lighting function
float getLight(vec3 pos, vec3 lightOrigin, bool includeWater) {
    // Get the light origin
    vec3 lightPos = lightOrigin; // basically a point in 3D space
    //lightPos.xz += vec2(sin(iTime * 0.1), cos(iTime * 0.1)) * 4.0;
    
    // Get the light ray
    vec3 light = normalize(lightPos - pos);
    
    // Get normal
    vec3 normal = getNormal(pos, includeWater);
    
    // Get diffused lighting 
    // We want 0 if the rays are parallel, 1 if the rays are perpendicular
    // Hence we use dot product
    float dif = clamp(dot(normal, light), 0.0, 1.0);
    float dis = rayMarch(pos + normal * SURFDIS * 2.0, light, includeWater).x;
    
    // Getting shadows (if disToSphere < disToLight then we have a shadow)
    if (dis < length(lightPos - pos)) {
        float shadowIntensity = 0.2;
        dif *= shadowIntensity;
    }
    
    return dif;
}

// Main //////////////////////////////////////////////////////////////////////

// Background function; depth perception
vec3 background(vec3 rayDir) {
    vec3 col = vec3(0.0);
    
    float y = rayDir.y * 0.5 + 0.5; // light is top, dark is bottom. 
    col += y * darkSpaceCol;
    float x = rayDir.x * 0.5 + 0.5;
    col += x * orangeSunray * 0.6;
    
    return col;
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
    vec3 rayOrigin = vec3(0, camHeight, 0); // this is the camera (origin of vector)
    vec3 rayDir = normalize(vec3(uv.x, uv.y + downTilt, 1));
    vec3 fillLightPos = vec3(-4, 4, -2);
    vec3 keyLightPos = vec3(4, 4, 2);
    
    vec3 spherePos = vec3(0, 0.8, 6);
    float maxLightDis = 8.0; // used to clamp light sources as ratio
    float IOR = 1.33; // water index of refraction
    
    // Visualise 
    vec2 passedDMat = rayMarch(rayOrigin, rayDir, true);
    float dis = passedDMat.x;
    int mat = int(passedDMat.y);
    if (dis < MAXDIS) {
        vec3 pos = rayOrigin + rayDir * dis;
        vec3 normal = getNormal(pos, true);
        float fillDiffuseLight = getLight(pos, fillLightPos, true);
        float keyDiffuseLight = getLight(pos, keyLightPos, true);
        float finalFill = fillDiffuseLight * clamp(1.0 - (length(fillLightPos - pos) / maxLightDis), 0.0, 1.0);
        float finalKey = keyDiffuseLight * clamp(1.0 - (length(keyLightPos - pos) / maxLightDis), 0.0, 1.0);
        
        // Gradual lighting curve
        vec3 mixCol = vec3(0.0); //mix(shadowCol, lightCol, diffuseLight);
        
        if (mat == MATTERRAIN) {
            mixCol = dirtCol;
            // Steepness
            vec3 localUpDir = normalize(pos);
            float steepness = -dot(normal, localUpDir);
            float radialHeight = abs(length(spherePos - pos)); // origin
            // Snow
            if (radialHeight < WATERHEIGHT + 0.05) {
                mixCol = mix(mixCol, grassCol, smoothstep(0.2, 0.6, steepness)); //grassCol * steepness;
            }
            if (radialHeight > WATERHEIGHT + 0.1) { 
                mixCol = mix(mixCol, snowCol, smoothstep(0.7, 0.85, steepness)); 
            }
        } else {
            mixCol = waterCol;
            // Ray march for depth
            vec3 rayDirIn = refract(rayDir, normal, 1.0 / IOR);
            vec3 posEnter = pos - normal * SURFDIS * 3.0;
            float disIn = rayMarch(posEnter, rayDirIn, false).x;
            float heightStepper = pow(length(disIn), 0.25) * 4.0 + 0.1;
            mixCol = mix(dirtCol, waterCol, heightStepper);
            
            // Specular (take angle diff between rayDir and reflect)
            float specStrength = 0.025;
            vec3 lightReflectRayDir = reflect(keyLightPos - posEnter, normal);
            float spec = pow(max(dot(rayDir, lightReflectRayDir), 0.0), 2.0);
            vec3 specularCol = specStrength * spec * blueLightCol;
            mixCol += specularCol * heightStepper;
        }
        
        mixCol *= fillDiffuseLight;
        mixCol += blueLightCol * finalFill;
        mixCol += orangeLightCol * finalKey;
        
        col = vec3(mixCol);
    }
    
    // Fog
    float fogStart = 0.0;
    float fogEnd = MAXDIS;
    col = mix(col, background(rayDir), smoothstep(fogStart, fogEnd, dis));
    
    // Output to screen
    fragColor = vec4(col,1.0);
}