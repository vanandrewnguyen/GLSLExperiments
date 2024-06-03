// Van-Andrew Nguyen
/*
03/06/24
Gyroid Water Demo

Demo to show SSS in a water material.
*/

////////////// Buffer A /////////////////
// Globals //////////////////////////////////////////////////////////////////////
#define MAXSTEPS 100
#define MAXDIS 24.0
#define SURFDIS 0.01

#define PI 3.1416

#define MATNULL 0
#define MATWATER 1

// Define global settings
const float AMBIENT = 0.2;
const float TIME_MULT = 0.5;

// Define global colours
const vec3 SKY_TOP_COL = vec3(0.941, 0.729, 0.6);
const vec3 SKY_BOTTOM_COL = vec3(0.827, 0.458, 0.211);

const vec3 CAUSTIC_COL = vec3(0.941, 0.729, 0.6);
const vec3 WATER_COL = vec3(0.313, 0.949, 0.847);


// Define lighting
const int NUM_LIGHTS = 2;
const vec3 KEY_LIGHT_POS[] = vec3[](
    vec3(3, 4, 2),
    vec3(-3, 0.5, -3)
);

const vec3 KEY_LIGHT_COL[] = vec3[](
    vec3(1, 0, 0),
    vec3(0, 0, 1)
);

// Noise /////////////////////////////////////////////////////////////////////

// Return a randomised float using a vec2 (uv coord)
float hash1(float n) {
    return fract( n*17.0*fract( n*0.3183099 ) );
}

float hash21(vec2 uv) {
    return fract(hash1(uv.x * uv.y));
}

vec2 hash22(vec2 p) {
    return fract(sin(vec2(dot(p,vec2(127.1,311.7)), dot(p,vec2(269.5,183.3))))*43758.5453);
}

vec3 hash33(vec3 p3) {
	p3 = fract(p3 * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yxz+19.19);
    return fract(vec3((p3.x + p3.y)*p3.z, (p3.x+p3.z)*p3.y, (p3.y+p3.z)*p3.x));
}

// Cellular Noise
float getVoronoi(vec2 uv) {
    float minDis = 1.0;
    vec2 cellID = floor(uv);
    vec2 gridUV = fract(uv);

    for (int y= -1; y <= 1; y++) {
        for (int x= -1; x <= 1; x++) {
            // Neighbor place in the grid
            vec2 neighbor = vec2(x, y);

            // Random position 
            vec2 point = hash22(cellID + neighbor);
            point = 0.5 + 0.5*sin(iTime + 6.2831*point);

			// Vector between the pixel and the point
            vec2 diff = neighbor + point - gridUV;

            // Distance to the point
            float dist = length(diff);
            minDis = min(minDis, dist);
        }
    }

    return minDis;
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

// Cube Distance
float sdCube(vec3 pos, vec3 size) { 
    pos = abs(pos) - size;
    float val = length(max(abs(pos) - size, 0.0));
    return val;
}

float sdTorus(vec3 pos, vec3 center, vec2 rad) {
    pos -= center;
    float x = length(pos.xz) - rad.x;
    return length(vec2(x, pos.y)) - rad.y;
}

// Custom shapes //

// Gerstner Wave Equation
// Source: https://www.shadertoy.com/view/7djGRR
vec3 gerstnerWave(vec2 coord, float wavelength, float steepness, vec2 direction) {
    const float gravitationalConst = 9.81;
    
    vec3 gerstner;
    float k = 2.0 * PI / wavelength;
    float c = sqrt(gravitationalConst / k);
    float a = steepness / k;
    vec2 dir = normalize(direction);
    float f = k * (dot(dir, coord.xy) - c * iTime * (1.0 / TIME_MULT));
    
    gerstner.x += dir.x * (a * cos(f));
    gerstner.y = a * sin(f);
    gerstner.z += dir.y * (a * cos(f));
    
    return gerstner;
}

vec3 gerstnerFBM(vec3 pos) {
    vec3 waves;
    waves += gerstnerWave(pos.xz * 880.0, 60.0, 1.0, vec2(1, 1));
    waves += gerstnerWave(pos.xz * 880.0, 31.0, 1.0, vec2(1, 0.6));
    waves += gerstnerWave(pos.xz * 880.0, 18.0, 1.0, vec2(1, 1.3));
    waves += gerstnerWave(pos.xz * 880.0, 26.0, 1.0, vec2(0.7, 1.0));
    waves += gerstnerWave(pos.xz * 880.0, 22.0, 1.0, vec2(0.8, 0.6));
    
    waves *= 0.025;
    return waves;
}

// Gyroid repeating structure
float sdGyroid(vec3 pos, float scale, float thickness, float bias, float posBias1, float posBias2) {
    // Scale coords
    pos *= scale;
    
    // Gerstner waves
    vec3 waveDis = gerstnerFBM(pos * 0.015);
    pos += waveDis;

    float val = abs(dot(sin(pos * posBias1), cos(pos.zxy * posBias2)) - bias) / (scale * max(posBias1, posBias2)) - thickness;
    
    return val * 0.2;
}

// Vectors //////////////////////////////////////////////////////////////////////

// Get distance function
vec2 getDist(vec3 pos, bool includeWater) {
    
    // Map dist to shapes
    float dummy = sdSphere(pos, vec3(100, 100, 100), 0.1);
    float cubeDis = sdCube(pos + vec3(0, -0.2, 0), vec3(1.4));
    float sphereDis = sdSphere(pos, vec3(0, 1.1, 0), 2.4);
    float sphereSmallDis = sdSphere(pos, vec3(2.5, 1.1, 0), 1.0);
    float torusDis = sdTorus(pos, vec3(0, 1.0, 0), vec2(2.8, 0.5));
    float floorDis = sdCube(pos + vec3(0, 0.5, 0), vec3(1.5, 0.1, 1.5));
    
    // Build gyroid
    float gyroid1Dis = sdGyroid(pos, 2.17, 0.1, 1.2, 1.0, 1.0);
    float gyroid2Dis = sdGyroid(pos, 5.34, 0.02, 0.6, 1.1, 1.5);
    float gyroid3Dis = sdGyroid(pos, 8.65, 0.02, 0.6, 1.1, 1.5);
    float gyroidDis;
    gyroidDis = opIntersection(gyroid1Dis, cubeDis) - gyroid2Dis * 0.1 - gyroid3Dis * 0.05; // bump map + = rocky, - = smooth bubbles
    
    // Final distance
    float finalDis = dummy; // sphereSmallDis;
    if (includeWater) {
        finalDis = gyroidDis; //opUnion(sphereDis, floorDis); // opUnion(sphereDis, sphereSmallDis);
    }
    
    // Final material to return
    int mat = MATNULL;
    if (finalDis == gyroidDis) {
        mat = MATWATER;
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

// Soft shadows
// From https://www.iquilezles.org/www/articles/rmshadows/rmshadows.htm
float getSoftShadow(vec3 pos, vec3 normal, vec3 lightDir, float shadowInt, bool includeWater) {
    float res = 1.0;
    float dis = 0.0;
    float t = SURFDIS;
    vec3 startingPos = pos + normal * SURFDIS * 2.0;
    for(int i = 0; i < MAXSTEPS; i++) {
        dis = getDist(startingPos + lightDir * t, includeWater).x;
        if (dis < SURFDIS) { return 0.0; }
        if (t >= MAXDIS) { return res; }
        res = min(res, shadowInt * dis / t);
        t += dis;
    }
    
    return res;
}

// Hard shadows
void getHardShadow(vec3 pos, vec3 normal, vec3 lightDir, vec3 lightPos, inout float dif, bool includeWater) {
    float dis = rayMarch(pos + normal * SURFDIS * 2.0, lightDir, includeWater).x;
    
    // Getting shadows (if disToSphere < disToLight then we have a shadow)
    if (dis < length(lightPos - pos)) {
        float shadowIntensity = 0.2;
        dif *= shadowIntensity;
    }
}


// Lighting //////////////////////////////////////////////////////////////////////

float getPointLight(vec3 pos, vec3 lightOrigin, bool includeWater) {
    // Get the light origin
    vec3 lightPos = lightOrigin;
    
    // Get the light ray
    vec3 lightDir = normalize(lightPos - pos);
    
    // Get normal
    vec3 normal = getNormal(pos, includeWater);
    
    // Get diffused lighting 
    // We want 0 if the rays are parallel, 1 if the rays are perpendicular
    // Hence we use dot product
    float dif = clamp(dot(normal, lightDir), 0.0, 1.0);
 
    // Hard shadows
    getHardShadow(pos, normal, lightDir, lightPos, dif, includeWater);
    dif *= getSoftShadow(pos, normal, lightDir, 8.0, includeWater);
    
    return dif;
}

float getSpotLight(vec3 pos, vec3 lightOrigin, vec3 endLightDir, float angleMax, float lightDropOff, bool includeWater) {
    // Identical to a point light, however, we restrict the light through a certain degree (use dot prod)
    
    // Get the light origin
    vec3 lightPos = lightOrigin;
    
    // Get angle dif
    vec3 lightDirToIntersect = normalize(lightPos - pos);
    vec3 lightDirToEnd = normalize(lightOrigin - endLightDir);
    float angleDif = dot(lightDirToIntersect, lightDirToEnd);
    
    // Get normal
    vec3 normal = getNormal(pos, includeWater);
    
    float dif = clamp(dot(normal, lightDirToIntersect), 0.0, 1.0);
    dif *= smoothstep(angleMax, angleMax + lightDropOff, angleDif);
    
    // Hard shadows
    getHardShadow(pos, normal, lightDirToIntersect, lightPos, dif, includeWater);
    dif *= getSoftShadow(pos, normal, lightDirToIntersect, 8.0, includeWater);
    
    return dif;
}

float getAmbientOcclusion(vec3 pos, vec3 normal, bool includeWater) {
    float aoSum = 0.0;
    float aoMaxSum = 0.0;
    int aoStepMax = 8;
    float aoStepSize = 0.2;
    
    // Fixed step marching, we march a number of times to see how many objects are nearby
    for (int i = 0; i < aoStepMax; i++) {
        float offset = float(i+1) * aoStepSize;
        vec3 currPos = pos + normal * offset;
        float inc = 1.0 / pow(2.0, float(i));
        
        // Increment our sums
        aoSum += inc * getDist(currPos, includeWater).x;
        aoMaxSum += inc * offset;
    }
    
    // Make sure it's between 0->1
    return aoSum / aoMaxSum;
}

float getLightIntensity(float dif, vec3 lightPos, vec3 pos, float maxLightLen) {
    return dif * clamp(1.0 - (length(lightPos - pos) / maxLightLen), 0.0, 1.0);
}

// Rendering //////////////////////////////////////////////////////////////////////

vec3 getMatCol(int mat) {
    if (mat == MATWATER) {
        return WATER_COL;
    }
    if (mat == MATNULL) {
        return vec3(0);
    }
    return vec3(0);
}

vec3 background(vec3 rayDir) {
    // Texture of cube map
    return texture(iChannel0, rayDir).rgb;
    
    // Gradient background
    vec3 col = vec3(0);
    float y = rayDir.y * 0.4 + 0.6;
    col += (1.0 - y) * vec3(1) * WATER_COL;
    
    // Rays
    float a = atan(rayDir.x, rayDir.z);
    float rays = sin(a * 10.0 + iTime * TIME_MULT) * sin(a * 7.0 - iTime * TIME_MULT) * sin(a * 6.0);
    rays *= smoothstep(0.5, 0.3, y);
    col += rays;
    col = max(col, 0.0);
    col += smoothstep(0.5, 0.0, y);
    
    return col;
}

// https://blog.demofox.org/2017/01/09/raytracing-reflection-refraction-fresnel-total-internal-reflection-and-beers-law/
float fresnelReflectAmount(float n1, float n2, vec3 normal, vec3 incident, float reflectiveness) {
        // Schlick aproximation
        float r0 = (n1-n2) / (n1+n2);
        r0 *= r0;
        float cosX = -dot(normal, incident);
        if (n1 > n2)
        {
            float n = n1/n2;
            float sinT2 = n*n*(1.0-cosX*cosX);
            // Total internal reflection
            if (sinT2 > 1.0)
                return 1.0;
            cosX = sqrt(1.0-sinT2);
        }
        float x = 1.0-cosX;
        float ret = r0+(1.0-r0)*x*x*x*x*x;
 
        // adjust reflect multiplier for object reflectivity
        ret = (reflectiveness + (1.0-reflectiveness) * ret);
        return ret;
}

float beerLambertAttenuation(float absorptionCoefficient, float dis) {
    // Square inverse law
    // In a clear volume like water, you want the absorption to be low since water doesn't absorb light
    return exp(-absorptionCoefficient * dis);
}

float getLightAttenuation(float disToLight) {
    // This is the square inverse law
    float attenuationFactor = 1.5; //1.65;
    return 1.0 / pow(disToLight, attenuationFactor);
}

float getLightVis(vec3 rayOrigin, vec3 rayDir, float maxT, int maxSteps, float marchSize, float absorptionCoefficient) {
    // Get light vis based on distance (exp decay)
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

vec3 renderVolume(vec3 entryPos, vec3 rayDir, float IOR) {    
    // Fixed step ray marcher, calc how much light passes through per step
    float opaqueVis = 1.0;
    float marchSize = 0.1;
    float shadowMarchSize = 0.2;
    
    float volumeDepth = 0.0;
    int volumeMaxSteps = 50;
    int shadowMaxSteps = 20;
    float absorptionCoefficient = 0.5;
    vec3 volumeCol = vec3(0.0);
    
    for (int i = 0; i < volumeMaxSteps; i++) {
        volumeDepth += marchSize;
        vec3 insidePos = entryPos + rayDir * volumeDepth;
        
        // Get absorption per step
        float prevOpaqueVis = opaqueVis;
        opaqueVis *= beerLambertAttenuation(absorptionCoefficient, marchSize);
        float currAbsorption = prevOpaqueVis - opaqueVis;
        
        // Compare with diffuse lights
        for (int i = 0; i < NUM_LIGHTS; i++) {
            vec3 lightDir = KEY_LIGHT_POS[i] - insidePos;
            float lightDis = length(lightDir);
            float lightVis = getLightVis(insidePos, lightDir, lightDis, shadowMaxSteps, shadowMarchSize, absorptionCoefficient);
            vec3 outLightCol = lightVis * KEY_LIGHT_COL[i] * getLightAttenuation(lightDis);
            volumeCol += currAbsorption * WATER_COL * outLightCol;
        }

        // Ambient light
        volumeCol += currAbsorption * WATER_COL * AMBIENT;
        
        // Material underneath (water caustics -> faked with voronoi)
        vec2 passedDMat = getDist(insidePos, false);
        float disToScene = passedDMat.x;
        int matInterior = int(passedDMat.y);
        if (disToScene < SURFDIS) {
            if (matInterior != MATWATER) {
                vec3 intCol = getMatCol(matInterior);
                // First mix the water colour with the interior below
                volumeCol = mix(volumeCol, intCol, 1.0-opaqueVis);
                
                // Then add caustics
                vec3 causticCol = 0.1 * volumeDepth * getVoronoi(rayDir.yz * 20.0) * CAUSTIC_COL;
                volumeCol += causticCol;
                
                return volumeCol;
            }
        }
    }
    
    // If we've reached this stage, this means that we've exited the water 
    // and now need to merge the background col
    // Also add colour dispersion
    float abberation = 0.01;
    
    vec3 posExit = entryPos + rayDir * volumeDepth;
    vec3 exitNormal = -getNormal(posExit, false);
    vec3 rayDirOut;
    
    vec3 backgroundCol = vec3(0.0);
    rayDirOut = refract(rayDir, exitNormal, 1.33 - abberation);
    if (length(rayDirOut) == 0.0) {
        rayDirOut = reflect(rayDir, exitNormal);
    }
    backgroundCol.r = background(rayDirOut).r;
    rayDirOut = refract(rayDir, exitNormal, 1.33);
    if (length(rayDirOut) == 0.0) {
        rayDirOut = reflect(rayDir, exitNormal);
    }
    backgroundCol.g = background(rayDirOut).g;
    rayDirOut = refract(rayDir, exitNormal, 1.33 + abberation);
    if (length(rayDirOut) == 0.0) {
        rayDirOut = reflect(rayDir, exitNormal);
    }
    backgroundCol.b = background(rayDirOut).b;
    
    // Merge colours
    float maxWidth = 0.6;
    float stepper = clamp(volumeDepth / maxWidth, 0.0, 1.0);
    return mix(volumeCol, backgroundCol, maxWidth);
}

vec3 render(vec3 rayOrigin, vec3 rayDir) {
    // Assign col
    vec3 col = background(rayDir);

    // Lighting Setup
    float lightMaxLen = 10.0;
    vec3 ambientLight = vec3(AMBIENT);
    vec3 mixCol = vec3(ambientLight); 
    
    // Visualise 
    vec2 passedDMat = rayMarch(rayOrigin, rayDir, true);
    float dis = passedDMat.x;
    int mat = int(passedDMat.y);
    
    if (dis < MAXDIS) {
        vec3 pos = rayOrigin + rayDir * dis;
        vec3 normal = getNormal(pos, true);
        vec3 posEnter = pos - normal * SURFDIS * 3.0;
        
        // Materials
        float IOR = 1.33;
        if (mat == MATWATER) {
            vec3 rayDirIn = refract(rayDir, normal, 1.0 / IOR);
            mixCol = renderVolume(posEnter, rayDirIn, IOR);
        }
        
        // Ambient Occlusion (check out Alex Evans)
        // Trick is to sample sdf along normal
        if (mat != MATWATER) {
            float ao = getAmbientOcclusion(pos, normal, true);
            mixCol *= ao; // <- video pt 1 has tricks on cheaper AO
        }
        
        // Lighting
        // Diffuse
        if (mat != MATWATER) {
            for (int i = 0; i < NUM_LIGHTS; i++) {
                float keyDiffuse = getPointLight(pos, KEY_LIGHT_POS[i], true);
                mixCol += keyDiffuse * KEY_LIGHT_COL[i] * 0.4;
            }
        }
        
        // Specular
        float specStrength = 0.01;
        for (int i = 0; i < NUM_LIGHTS; i++) {
            vec3 lightReflectRayDir = reflect(KEY_LIGHT_POS[i] - posEnter, normal);
            float spec = pow(max(dot(rayDir, lightReflectRayDir), 0.0), 2.0);
            float fresnel = clamp(dot(normal, (KEY_LIGHT_POS[i] - posEnter)), 0.0, 1.0);
            vec3 specularCol = specStrength * spec * KEY_LIGHT_COL[i] * fresnel;
            mixCol += specularCol;
        }
        
        // Reflections
        if (mat == MATWATER) {
            vec3 reflectedRay = reflect(rayDir, normal);
            float fresnel = fresnelReflectAmount(1.0, IOR, normal, rayDir, 0.1); 
            mixCol += background(reflectedRay) * fresnel;
        }
        
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

// Main //////////////////////////////////////////////////////////////////////

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
    float camZoom = 1.0;
    vec3 rayOrigin = vec3(camX, camY, camZ);
    rayOrigin.yz *= rotate(-m.y * PI + 1.);
    rayOrigin.xz *= rotate(-m.x * 2.0 * PI);
    vec3 rayDir = getRayDir(uv, rayOrigin, vec3(camYTilt, 0, camXTilt), camZoom);
    
    col = render(rayOrigin, rayDir);
    
    // Output to screen
    fragColor = vec4(col,1.0);
}

//////////////   Main   /////////////////
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2 uv = fragCoord.xy / iResolution.xy;
    vec3 col = texture(iChannel0, uv).rgb;
    
    // Vignette
    float vig = uv.x * uv.y * (1.0 - uv.x) * (1.0 - uv.y);
    float vigIntensity = 0.2;
    vig = clamp(pow(16.0 * vig, vigIntensity), 0.0, 1.0);
    col *= 0.5 * vig + 0.5;
    
    // Colour correction
    col = pow(col, vec3(1.0/2.2));
    
    fragColor = vec4(col, 1.0);
}
