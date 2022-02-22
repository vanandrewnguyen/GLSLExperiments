/*
Van Andrew Nguyen
22/02/22
[Wax Cloud]
This is another attempt to render translucent materials. Using the help of this example: https://www.shadertoy.com/view/tsScDG
I used circle tracing and also fixed step marching to render the interior of translucent materials. Was aiming for fluffy clouds but got a wax model instead. I think
it has something to do with the SDF I'm working with. I should be trying to create an inverse value as a density function rather than a distance function. Nonetheless, good 
experimentation with the interactions with light.
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

#define PI 3.1416

#define MATNULL 0
#define MATLIGHT1 1
#define MATLIGHT2 2
#define MATCLOUD 3

const vec3 shadowCol = vec3(0.2, 0.1, 0.15);
const vec3 lightCol = vec3(0.8, 0.75, 0.7);
const vec3 cloudAlbedo = vec3(1.0);
const vec3 fillLightCol = vec3(0.854, 0.364, 0.501);
const vec3 keyLightCol = vec3(0.937, 0.376, 0.101);

const vec3 ambientLight = vec3(0.2);

// Noise /////////////////////////////////////////////////////////////////////

// Return a randomised float using a vec2 (uv coord)
float hash1(float n) {
    return fract( n*17.0*fract( n*0.3183099 ) );
}

float hash21(vec2 uv) {
    return fract(hash1(uv.x * uv.y));
}

vec3 hash33(vec3 p3) {
	p3 = fract(p3 * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yxz+19.19);
    return fract(vec3((p3.x + p3.y)*p3.z, (p3.x+p3.z)*p3.y, (p3.y+p3.z)*p3.x));
}

float noise(in vec3 x) {
    vec3 p = floor(x);
    vec3 f = fract(x);
    
    f = f * f * (3.0 - 2.0 * f);
    
    float n = p.x + p.y * 57.0 + 113.0 * p.z;
    
    float res = mix(mix(mix(hash1(n +   0.0), hash1(n +   1.0), f.x),
                        mix(hash1(n +  57.0), hash1(n +  58.0), f.x), f.y),
                    mix(mix(hash1(n + 113.0), hash1(n + 114.0), f.x),
                        mix(hash1(n + 170.0), hash1(n + 171.0), f.x), f.y), f.z);
    return res;
}

// Return smooth noise 
float smoothNoise(vec2 uv) {
    // Create 1D random value
    vec2 index = uv;
    vec2 localUV = fract(index); 
    vec2 cellID = floor(index); 
    
    localUV = localUV*localUV*(3.0 - 2.0 * localUV); // Hermite Curve
    
    // Get noise values for corners of each cell (bottom/top right + left, then mix it)
    float bl = hash21(cellID);
    float br = hash21(cellID + vec2(1, 0));
    float b = mix(bl, br, localUV.x);
    float tl = hash21(cellID + vec2(0, 1));
    float tr = hash21(cellID + vec2(1, 1));
    float t = mix(tl, tr, localUV.x);
    float noiseCol = mix(b, t, localUV.y);
        
    return noiseCol;
}

mat3 m = mat3( 0.00,  0.80,  0.60,
              -0.80,  0.36, -0.48,
              -0.60, -0.48,  0.64);
float fbm(vec3 p) {
    float f;
    f  = 0.5000 * noise(p); p = m * p * 2.02;
    f += 0.2500 * noise(p); p = m * p * 2.03;
    f += 0.1250 * noise(p);
    return f;
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

// Torus Distance
float sdTorus(vec3 pos, vec3 center, vec2 rad) {
    // We subtract the smaller radius from the length of the vector running from
    // origin point to middle of torus
    pos -= center;
    float x = length(pos.xz) - rad.x;
    return length(vec2(x, pos.y)) - rad.y;
}

// Cube Distance
float sdCube(vec3 pos, vec3 size) { 
    pos = abs(pos) - size;

    // We get the length of the position - size (dis) and grab a version > 0
    float val = length(max(abs(pos) - size, 0.0));
    return val;
}

// Get light pos from getDist
vec3 getKeyLightPos(vec3 pos) {
    // Pos input is a variation of light placement
    float lightRad = 2.0;
    vec3 lightPos = pos + vec3(0.1 * sin(iTime), 2, 0.1 * sin(PI + iTime));
    lightPos.xz += vec2(sin(iTime + PI), cos(iTime + PI)) * (lightRad + 0.5 * sin(iTime));
    return lightPos;
}

vec3 getFillLightPos(vec3 pos) {
    float lightRad = 2.0;
    vec3 lightPos = pos + vec3(0.1 * cos(iTime), 2, 0.1 * cos(PI + iTime));
    lightPos.xz += vec2(sin(iTime), cos(iTime)) * (lightRad + 0.5 * cos(iTime));
    return lightPos;
}

// Cloud SDF (Merge 3 spheres with plane)
float sdCloudVolume(vec3 pos) {
    float disOut = 0.0;
    float cloudRad = 0.8;
    float placeDis = pos.y + 0.1;
    vec3 sphere1Pos = vec3(0, 0.6 + 0.4 * sin(iTime * 1.4), 0);
    vec3 sphere2Pos = vec3(0.1, 1.2 + 0.2 * cos(iTime), -0.9);
    vec3 sphere3Pos = vec3(-0.8, 0.8 + 0.3 * sin(iTime * 0.8), 0.3);
    sphere1Pos.xz += vec2(sin(iTime * 0.9 + PI), cos(iTime * 0.7 + PI)) * (cloudRad + 0.5 * sin(iTime));
    sphere2Pos.xz += vec2(sin(iTime * 0.2), cos(iTime * 1.2)) * (cloudRad + 0.5 * cos(iTime * 1.1));
    sphere3Pos.xz += vec2(sin(iTime * 0.5 - PI), cos(iTime * 1.6 - PI)) * (cloudRad + 0.5 * cos(iTime * 0.8));
    float sphere1Dis = sdSphere(pos - fbm(pos * 1.5) * 2.0, sphere1Pos, 1.0); // - fbm(pos)
    float sphere2Dis = sdSphere(pos - fbm(pos * 2.0) * 1.0, sphere2Pos, 0.7);
    float sphere3Dis = sdSphere(pos - fbm(pos * 1.4) * 1.4, sphere2Pos, 0.8);
    
    disOut = opSmoothUnion(sphere1Dis, sphere2Dis, 0.6);
    disOut = opSmoothUnion(sphere3Dis, disOut, 0.6);
    disOut = opSmoothUnion(disOut, placeDis, 3.4);
    return disOut * 0.8;
}

// Vectors //////////////////////////////////////////////////////////////////////

// Get distance function
vec2 getDist(vec3 pos) {
    
    // Map dist to shapes
    float planeDis = pos.y; 
    // Lights
    float sphereLight1Dis = sdSphere(pos, getKeyLightPos(vec3(0)), 0.15);
    float sphereLight2Dis = sdSphere(pos, getFillLightPos(vec3(0)), 0.15);
    // Cloud 
    float cloudDis = sdCloudVolume(pos);
    cloudDis = opSmoothSubtraction(sphereLight1Dis, cloudDis, 1.4);
    cloudDis = opSmoothSubtraction(sphereLight2Dis, cloudDis, 1.4);
    
    // Final distance to return
    float finalDis = min(sphereLight1Dis, planeDis);
    finalDis = min(sphereLight2Dis, finalDis);
    finalDis = min(cloudDis, finalDis);
    
    // Final material to return
    int mat = MATNULL;
    if (finalDis == sphereLight1Dis) {
        mat = MATLIGHT1;
    } else if (finalDis == sphereLight2Dis) {
        mat = MATLIGHT2;
    } else if (finalDis == cloudDis) {
        mat = MATCLOUD;
    }
    
    return vec2(finalDis, mat);
}

// Return the normal ray
vec3 getNormal(vec3 pos) {
    float dis = getDist(pos).x;
    vec2 val = vec2(0.01, 0.0);
    
    // To get the slope we give the curve two values super close together
    // Instead of deriving we can do this method via swizzling
    vec3 normal = dis - vec3(getDist(pos-val.xyy).x, 
                             getDist(pos-val.yxy).x, 
                             getDist(pos-val.yyx).x);
    // Same way of doing getDist(p-vec3(0.01, 0.0, 0.0), ... etc
    return normalize(normal);
}

// Ray Marching function
vec2 rayMarch(vec3 rayOrigin, vec3 rayDir) {
    float disFromOrigin = 0.0;
    int mat;
    
    // Now we loop for max steps or until the distance is very small
    for (int i=0;i<MAXSTEPS;i++) {
        // Declare current marching location
        vec3 pos = rayOrigin + rayDir * disFromOrigin;
        vec2 passedDMat = getDist(pos);
        float disToScene = passedDMat.x;
        mat = int(passedDMat.y);
        
        disFromOrigin += disToScene;
        
        // Exit condition; we have a hit
        if (disFromOrigin > MAXDIS || disToScene < SURFDIS) { break; }
    }
    
    return vec2(disFromOrigin, mat);
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
    float dis = rayMarch(pos + normal * SURFDIS * 2.0, light).x;
    
    // Getting shadows (if disToSphere < disToLight then we have a shadow)
    if (dis < length(lightPos - pos)) {
        float shadowIntensity = 0.2;
        dif *= shadowIntensity;
    }
    
    return dif;
}

// Main //////////////////////////////////////////////////////////////////////
float beerLambertAttenuation(float absorptionCoefficient, float dis) {
    // Square inverse law
    // In a clear volume like water, you want the absorption to be low since water doesn't absorb light
    return exp(-absorptionCoefficient * dis);
}

float getLightAttenuation(float disToLight) {
    // This is the square inverse law
    float attenuationFactor = 1.4; //1.65;
    return 1.0 / pow(disToLight, attenuationFactor);
}

float getLightVis(vec3 rayOrigin, vec3 rayDir, float maxT, int maxSteps, float marchSize, float absorptionCoefficient) {
    float lightVis = 1.0;
    float dis = 0.0;
    for (int i = 0; i < maxSteps; i++) {
        dis += marchSize;
        if (dis > maxT) { break; }
        vec3 pos = rayOrigin + dis * rayDir;
        lightVis *= beerLambertAttenuation(absorptionCoefficient, marchSize);
    }
    return lightVis * 2.0 + 0.5;
}

float getCloudDensity(vec3 pos) {
    float sdfDis = getDist(pos).x;
    float maxSDFMult = 1.0;
    
    return min(abs(sdfDis), maxSDFMult);
}

vec3 renderCloud(vec3 entryPos, vec3 rayDir) {
    float opaqueVis = 1.0;
    float marchSize = 0.1;
    float shadowMarchSize = 0.1;
    
    float volumeDepth = 0.0;
    int cloudMaxSteps = 32;
    int shadowCloudMaxSteps = 24;
    float absorptionCoefficient = 2.0;
    vec3 volumeCol = vec3(0.0);
    
    for (int i = 0; i < cloudMaxSteps; i++) {
        volumeDepth += marchSize;
        vec3 insidePos = entryPos + rayDir * volumeDepth;
        
        // Get absorption per step
        float prevOpaqueVis = opaqueVis;
        opaqueVis *= beerLambertAttenuation(absorptionCoefficient * getCloudDensity(insidePos), marchSize);
        float currAbsorption = prevOpaqueVis - opaqueVis;
        
        // Compare with diffuse lights (oh boy)
        vec3 keyLightDir = getKeyLightPos(vec3(0)) - insidePos;
        float keyLightDis = length(keyLightDir);
        float keyLightVis = getLightVis(insidePos, keyLightDir, keyLightDis, shadowCloudMaxSteps, shadowMarchSize, absorptionCoefficient);
        vec3 outKeyLightCol = keyLightVis * keyLightCol * getLightAttenuation(keyLightDis);
        volumeCol += currAbsorption * cloudAlbedo * outKeyLightCol;
        
        vec3 fillLightDir = getFillLightPos(vec3(0)) - insidePos;
        float fillLightDis = length(fillLightDir);
        float fillLightVis = getLightVis(insidePos, fillLightDir, fillLightDis, shadowCloudMaxSteps, shadowMarchSize, absorptionCoefficient);
        vec3 outFillLightCol = fillLightVis * fillLightCol * getLightAttenuation(fillLightDis);
        volumeCol += currAbsorption * cloudAlbedo * outFillLightCol;
        
        // Ambient light
        volumeCol += currAbsorption * cloudAlbedo * ambientLight;
    }
    return volumeCol;
}

vec3 render(vec3 rayOrigin, vec3 rayDir) {
    // Assign col
    vec3 col = vec3(0.0);

    // Lighting Setup
    float lightRad = 2.0;
    vec3 keyLightPos = getKeyLightPos(vec3(0, -0.2, 0));
    vec3 fillLightPos = getFillLightPos(vec3(0, -0.2, 0));
    float lightMaxLen = 10.0;
    vec3 mixCol = vec3(ambientLight); 
    
    // Visualise 
    vec2 passedDMat = rayMarch(rayOrigin, rayDir);
    float dis = passedDMat.x;
    int mat = int(passedDMat.y);
    
    if (dis < MAXDIS) {
        vec3 pos = rayOrigin + rayDir * dis;
        vec3 normal = getNormal(pos);
        
        // Lock the lights based on their distance to scene
        float keyDiffuseLight = getLight(pos, keyLightPos);
        float fillDiffuseLight = getLight(pos, fillLightPos);
        float finalFill = fillDiffuseLight * clamp(1.0 - (length(fillLightPos - pos) / lightMaxLen), 0.0, 1.0);
        float finalKey = keyDiffuseLight * clamp(1.0 - (length(keyLightPos - pos) / lightMaxLen), 0.0, 1.0);
        
        // Gradual lighting curve
        //vec3 mixCol = mix(shadowCol, lightCol, diffuseLight);
        
        // Materials
        if (mat == MATNULL) {
            // Checkerboard pattern
            float size = 1.0;
            float total = mod(floor(pos.z * size), 2.0) + 
                          mod(floor(pos.x * size), 2.0);
            bool isEven = mod(total, 2.0) == 0.0;
            mixCol *= (isEven) ? shadowCol : lightCol;
        } else if (mat == MATLIGHT1) {
            mixCol += keyLightCol;
        } else if (mat == MATLIGHT2) {
            mixCol += fillLightCol;
        } 
        
        if (mat == MATCLOUD) {
            mixCol += renderCloud(pos + SURFDIS * 3.0, rayDir);
        } else {
            // Diffuse Lighting
            mixCol += finalKey * keyLightCol;
            mixCol += finalFill * fillLightCol;
        }
        
        // Specular Lighting
            float specStrength = 0.2;
            vec3 lightReflectRayDir = reflect(normalize(keyLightPos - pos), normal);
            float spec = pow(max(dot(rayDir, lightReflectRayDir), 0.0), 32.0);
            vec3 specularCol = specStrength * spec * keyLightCol;
            mixCol += specularCol;
            lightReflectRayDir = reflect(normalize(fillLightPos - pos), normal);
            spec = pow(max(dot(rayDir, lightReflectRayDir), 0.0), 32.0);
            specularCol = specStrength * spec * fillLightCol;
            mixCol += specularCol;
        
        col = vec3(mixCol);
    }
    
    return col;
}

// From BigWings @ Shadertoy
// Update ray direction
vec3 getRayDir(vec2 uv, vec3 p, vec3 l, float z) {
    vec3 f = normalize(l-p),
        r = normalize(cross(vec3(0,1,0), f)),
        u = cross(f,r),
        c = f*z,
        i = c + uv.x*r + uv.y*u,
        d = normalize(i);
    return d;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab UV Coord
    vec2 uv = (fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;
    vec2 m = iMouse.xy/iResolution.xy;
    
    // Declare Col
    vec3 col = vec3(0.0);

    // Setup Camera
    float camX = 0.0;
    float camY = 0.0;
    float camZ = -10.0;
    float camXTilt = 0.0;
    float camYTilt = 0.0;
    float camZoom = 1.2;
    vec3 rayOrigin = vec3(camX, camY, camZ); // this is the camera (origin of vector)
    rayOrigin.yz *= rotate(-m.y * PI + 1.);
    rayOrigin.xz *= rotate(-m.x * 2.0 * PI);
    //vec3 rayDir = normalize(vec3(uv.x + camXTilt, uv.y + camYTilt, camZoom));
    vec3 rayDir = getRayDir(uv, rayOrigin, vec3(camYTilt, 0, camXTilt), camZoom);
    
    col = render(rayOrigin, rayDir);
    
    // Output to screen
    fragColor = vec4(col,1.0);
}
