/*
Van Andrew Nguyen
03/11/2021
[Ocean v1]

This was an attempt to render water from a close-up town down perspective. So of course I had to throw in a sunset for those specular 
reflections. The buoy here is mainly to showcase the reflections in the water. It also bobs and rotates along with the water. 
Water surface SDF uses Perlin Noise FBM.
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

#define TIMEMULT 0.2

#define MATWATER 1
#define MATSAND 2
#define MATBOUY 3

const vec3 shadowCol = vec3(0.2, 0.1, 0.15);
const vec3 lightCol = vec3(0.8, 0.75, 0.7);
const vec3 fillLightCol = vec3(0.843, 0.976, 0.968);
const vec3 keyLightCol = vec3(0.937, 0.376, 0.101);
const vec3 waterDarkCol = vec3(0.188, 0.592, 0.431);
const vec3 waterLightCol = vec3(0.349, 0.909, 0.796);
const vec3 sandCol = vec3(0.6, 0.525, 0.082);
const vec3 bouyCol = vec3(0.274, 0.274, 0.027);
const vec3 bouyTopCol = vec3(0.556, 0.454, 0.043);
const vec3 skyCol = vec3(0.807, 0.517, 0.313);

// Noise /////////////////////////////////////////////////////////////////////

// Smooth out the line using cubic
float smoothFade(float t) {
    return t*t*t*(t*(t*6.0 - 15.0) + 10.0);
}

// Same thing but 2D
vec2 getGradient(vec2 pos) {
	float texW = 256.0;
	vec4 v = texture(iChannel0, vec2(pos.x / texW, pos.y / texW));
    // Map this to [-1, 1]
    return normalize(v.xy * 2.0 - vec2(1.0)); 
}

// Noise function
float getNoise(vec2 pos) {
    // n(p) = (1 - F(p-p0))g(p0)(p-p0) + F(p-p0)g(p1)(p-p1)
    // Source: https://gpfault.net/posts/perlin-noise.txt.html
    // Get four corners as seperate points
    vec2 pos0 = floor(pos);
    vec2 pos1 = pos0 + vec2(1.0, 0.0);
    vec2 pos2 = pos0 + vec2(0.0, 1.0);
    vec2 pos3 = pos0 + vec2(1.0, 1.0);

    // Get gradients for the four corners
    vec2 gradient0 = getGradient(pos0);
    vec2 gradient1 = getGradient(pos1);
    vec2 gradient2 = getGradient(pos2);
    vec2 gradient3 = getGradient(pos3);
    
    // Horizontal Blend
    float hBlend = smoothFade(pos.x - pos0.x); 

    // Vertical Blend
    float vBlend = smoothFade(pos.y - pos0.y); 

    // Get dot product of top two lattice points, then bottom two points
    // Then interpolate both of them
    float p0p1 = (1.0 - hBlend) * dot(gradient0, (pos - pos0)) + 
                 hBlend * dot(gradient1, (pos - pos1));
    float p2p3 = (1.0 - hBlend) * dot(gradient2, (pos - pos2)) + 
                 hBlend * dot(gradient3, (pos - pos3));

    /* Calculate final result */
    return (1.0 - vBlend) * p0p1 + vBlend * p2p3;
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
float sdTorus(vec3 pos, vec3 center, vec2 rad) {
    // We subtract the smaller radius from the length of the vector running from
    // origin point to middle of torus
    pos -= center;
    float x = length(pos.xz) - rad.x;
    return length(vec2(x, pos.y)) - rad.y;
}

// Sphere Distance
float sdSphere(vec3 pos, vec3 center, float rad) {
    pos -= center;
    return length(pos) - rad;
}

// Cube Distance
float sdCube(vec3 pos, vec3 center, vec3 size) {
    pos -= center;
    
    // We get the length of the position - size (dis) and grab a version > 0
    float val = length(max(abs(pos) - size, 0.0));
    return val;
}

// Noise SDF
float sdNoiseFBM(vec2 pos) {
    float dis = 0.0;
    float amp = 0.5;
    float freq = 1.0;
    int iterations = 8;
    // Loop
    for (int i = 0; i < iterations; i++) {
        // Use abs() of the noise to get sharp ridges, for water ripples
        dis += abs(getNoise(pos * freq - iTime * TIMEMULT) * amp);
        amp *= 0.5;
        freq *= 2.0;
    }
    
    return dis;
}

// Vectors //////////////////////////////////////////////////////////////////////

// Get distance function
vec2 getDist(vec3 pos, bool includeWater) {
    
    // Get dist of the various shapes
    float waterPlaneDis = pos.y + sdNoiseFBM(pos.xz * 0.5 + 0.5) * 0.8;
    float oceanFloorPlaneDis = pos.y + 1.5 + getNoise(pos.xz * 0.5 + 0.5 - iTime * TIMEMULT * 0.5) * 1.5;
    
    // Create temporary position variable since we are manipulating it
    vec3 tmpPos = pos;
    tmpPos -= vec3(0, -0.1 + 0.02 * sin(iTime * TIMEMULT), 2);
    tmpPos.xy *= rotate(0.1 * sin(iTime));
    float bouyDis = sdTorus(tmpPos, vec3(0), vec2(0.3, 0.08)) + texture(iChannel2, pos.xz * 0.4).r * 0.004;
    
    // Final distance to return
    float finalDis;
    if (!includeWater) {
        finalDis = oceanFloorPlaneDis;
    } else {
        finalDis = opUnion(waterPlaneDis, bouyDis);
    }
    
    int mat = 0;
    if (finalDis == oceanFloorPlaneDis) {
        mat = MATSAND;
    } else if (finalDis == waterPlaneDis) {
        mat = MATWATER;
    } else if (finalDis == bouyDis) {
        mat = MATBOUY;
    }
    
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
    int mat;
    
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

vec3 background(vec3 rayDir) {
    vec3 col = vec3(0.0);
    
    float y = rayDir.y * 0.2 + 0.8; // light is top, dark is bottom. 
    col += y * skyCol;
    float x = rayDir.x * 0.7 + 0.3;
    col += x * keyLightCol * 0.3;
    
    return col;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab UV Coord
    vec2 uv = (fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;

    // Declare Col
    vec3 col = vec3(0.0);

    // Setup Camera
    float camHeight = 0.8;
    float downTilt = -0.4;
    vec3 rayOrigin = vec3(0, camHeight, 0); // this is the camera (origin of vector)
    vec3 rayDir = normalize(vec3(uv.x, uv.y + downTilt, 1));
    
    // Setup Lights
    vec3 keyDiffuseLightPos = vec3(5, 3, 1);
    vec3 fillDiffuseLightPos = vec3(-2, 9, 4);
    float maxLightDis = 10.0;
    vec3 ambientLight = vec3(0.1);
    
    // Setup Materials
    float IOR = 1.33;
    vec3 refTex = texture(iChannel1, rayDir).rgb;
    
    // Visualise 
    vec2 passedDMat = rayMarch(rayOrigin, rayDir, true);
    float dis = passedDMat.x;
    int mat = int(passedDMat.y);
    if (dis < MAXDIS) {
        // Grab Position of ray
        vec3 pos = rayOrigin + rayDir * dis;
        vec3 normal = getNormal(pos, true);
        vec3 posWaterExit = pos + normal * SURFDIS * 3.0; // slightly above water surface
        vec3 posWaterEnter = pos - normal * SURFDIS * 3.0; // slightly below water surface
        
        // Lighting Variables
        float keyDiffuseLight = getLight(pos, keyDiffuseLightPos, true);
        float fillDiffuseLight = getLight(pos, fillDiffuseLightPos, true);
        float finalFill = fillDiffuseLight * clamp(1.0 - (length(fillDiffuseLightPos - pos) / maxLightDis), 0.0, 1.0);
        float finalKey = keyDiffuseLight * clamp(1.0 - (length(keyDiffuseLightPos - pos) / maxLightDis), 0.0, 1.0);
        float fresnel = dot(normal, (keyDiffuseLightPos - posWaterEnter));
        
        vec3 mixCol; 
        
        if (mat == MATWATER) {
            mixCol = ambientLight;
            // Depth rendering
            vec3 rayDirIn = refract(rayDir, normal, 1.0 / IOR);
            vec2 interiorPassedDMat = rayMarch(posWaterEnter, rayDirIn, false);
            float disIn = interiorPassedDMat.x;
            int matIn = int(interiorPassedDMat.y);
            vec3 posExit = posWaterEnter + rayDirIn * disIn;
            
            float stepper = length(posWaterEnter - posExit) / 1.0;
            if (matIn == MATSAND) {
                mixCol *= mix(vec3(0), waterLightCol, stepper);
            }
            
            // Reflections
            float reflectedVal = 0.2;
            vec3 reflectedRay = reflect(rayDir, normal);
            vec2 exteriorPassedDMat = rayMarch(posWaterExit, reflectedRay, true);
            float disOut = exteriorPassedDMat.x;
            int matOut = int(exteriorPassedDMat.y);
            if (disOut < MAXDIS) {
                // Add colours that are reflected
                vec3 exteriorCol;
                if (matOut == MATBOUY) { 
                    exteriorCol = bouyCol; 
                }
                mixCol += exteriorCol * reflectedVal * fresnel;
            } else {
                mixCol += refTex * reflectedVal * 0.5 * fresnel;
            } 
        } else if (mat == MATSAND) {
            mixCol = sandCol;
        } else if (mat == MATBOUY) {
            mixCol = bouyCol;
            
            // Colour based on steepness
            vec3 localUpDir = vec3(0, 1, 0);
            float steepness = dot(normal, localUpDir);
            mixCol = mix(mixCol, bouyTopCol, smoothstep(0.5, 0.8, steepness));
        }
        
        // Ambient + Diffuse Lighting
        mixCol += ambientLight;
        mixCol += finalKey * keyLightCol;
        mixCol += finalFill * fillLightCol;
        
        // Specular Lighting
        float specStrength = 0.0002;
        vec3 lightReflectRayDir = reflect(keyDiffuseLightPos - posWaterEnter, normal);
        float spec = pow(max(dot(rayDir, lightReflectRayDir), 0.0), 4.0);
        vec3 specularCol = specStrength * spec * keyLightCol * fresnel;
        if (dis < MAXDIS * 0.5) { mixCol += specularCol; }
        
        col = vec3(mixCol);
    }
    
    // Fog / Sky
    float fogStart = 0.0;
    float fogEnd = MAXDIS;
    col = mix(col, background(rayDir), smoothstep(fogStart, fogEnd, dis));
    
    // Output to screen
    fragColor = vec4(col,1.0);
}