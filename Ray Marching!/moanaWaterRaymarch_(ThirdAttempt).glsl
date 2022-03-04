/*
Van Andrew Nguyen
04/03/2022
[Moana Style Water Shader]

Third attempt making a moana-water inspired shader. New additions from v2:
- AO (soft shadows on all SDF's makes it a more natural scene)
- Soft + hard shadows (more realistic lighting)
- Water Caustics (Fake voronoi pattern mapped onto surface below water)
- Gerstner Waves (New mathematical model to represent more realistic wave motion)
- SSS (Water absorbs light realistically with fixed ray marching steps using eqn for light atten.)
- Coloured Coral (colouring is based on gradient of surface)
- Updated camera geometry (can take mouse input)

Overall, super happy with the result. Water is much closer to source material though from a photo-realism standpoint it's still a long ways to go. Will look into 
fluid sim in the future. Next step is to investigate density functions for volumetric lighting and clouds.
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

#define TIMEMULT 0.5

#define AMBIENT 0.2

#define MATNULL 0
#define MATWATER 1
#define MATSAND  2
#define MATCORAL 3

const vec3 shadowCol = vec3(0.2, 0.1, 0.15);
const vec3 lightCol = vec3(0.8, 0.75, 0.7);
const vec3 fillLightCol = vec3(0.854, 0.364, 0.501);
const vec3 keyLightCol = vec3(0.937, 0.376, 0.101);
const vec3 bgCol1 = vec3(0.941, 0.729, 0.6);
const vec3 bgCol2 = vec3(0.827, 0.458, 0.211);

const vec3 waterCol = vec3(0.313, 0.949, 0.847);
const vec3 sandCol = vec3(0.905, 0.854, 0.694);
const vec3 coralCol = vec3(0.243, 0.156, 0.113);
const vec3 algaeCol = vec3(0.427, 0.549, 0.372);
const vec3 skyCol = vec3(0.941, 0.729, 0.6);

const vec3 keyLightPos = vec3(0, 2, 10); //vec3(3, 4, 2);
const vec3 fillLightPos = vec3(-2, 8, 2);

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

// Noise /////////////////////////////////////////////////////////////////////

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

// Gerstner Wave Equation
// Source: https://www.shadertoy.com/view/7djGRR
vec3 gerstnerWave(vec2 coord, float wavelength, float steepness, vec2 direction) {
    const float gravitationalConst = 9.81;
    
    vec3 gerstner;
    float k = 2.0 * PI / wavelength;
    float c = sqrt(gravitationalConst / k);
    float a = steepness / k;
    vec2 dir = normalize(direction);
    float f = k * (dot(dir, coord.xy) - c * iTime * (1.0 / TIMEMULT));
    
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

// Cyclinder Distance
float sdCylinder(vec3 pos, vec3 start, vec3 end, float rad) {
    // We get the projected distance of the origin->start onto start->end
    vec3 ab = end - start;
    vec3 ap = pos - start;
    // Then we normalise that distance and lock it within bounds
    float projection = dot(ab, ap) / dot(ab, ab); // / dot(ab, ab) to normalise 0->1
    // Then we get the distance from origin to that projected point MINUS radius of capsule
    vec3 len = start + projection * ab;
    
    float x = length(pos - len) - rad;
    float y = (abs(projection - 0.5) - 0.5) * length(ab);
    float extLen = length(max(vec2(x, y), 0.0));
    float intLen = min(max(x, y), 0.0);
    
    // Finally, we have the smallest distance to the capsule shape.
    return extLen + intLen;
}

float sdOceanSurf(vec3 pos) {
    // Smooth union a bunch of spheres and a plane (under the rocks surf)
    // This gives a big blobby shape that merges with the ground nicely

    vec3 waveDis = gerstnerFBM(pos * 0.04);
    float waveFloor = pos.y + 1.8;
    float dis;
    
    float rad = 0.8;
    
    // Position + rotation
    vec3 sphere1Pos = vec3(-2.1, -1.5 + 0.4 * sin(iTime * 1.2), 0);
    vec3 sphere2Pos = vec3(1.4, -1.4 + 0.2 * sin(iTime * 0.6), 2.5);
    vec3 sphere3Pos = vec3(0.5, -1.7 + 0.6 * sin(iTime * 1.0), -2.5);
    vec3 sphere4Pos = vec3(0.0, 0.8 + 0.8 * sin(iTime * 1.5), 0.0);
    
    sphere1Pos.xz += vec2(sin(iTime * 0.9 + PI), cos(iTime * 0.7 + PI)) * (rad + 0.5 * sin(iTime));
    sphere2Pos.xz += vec2(sin(iTime * 0.2), cos(iTime * 1.2)) * (rad + 0.5 * cos(iTime * 1.1));
    sphere3Pos.xz += vec2(sin(iTime * 0.5 - PI), cos(iTime * 1.6 - PI)) * (rad + 0.5 * cos(iTime * 0.8));
    sphere4Pos.xz += vec2(sin(iTime * 0.5), cos(iTime * 1.8)) * (rad + 0.5 * cos(iTime * 0.1));
    
    // SDF
    float sphere1Dis = sdSphere(pos + waveDis * 0.18, sphere1Pos, 1.5);
    float sphere2Dis = sdSphere(pos + waveDis * 0.20, sphere2Pos, 2.3);
    float sphere3Dis = sdSphere(pos + waveDis * 0.16, sphere3Pos, 1.3);
    float sphere4Dis = sdSphere(pos + waveDis * 0.12, sphere4Pos, 0.6);

    float k = 2.5;
    dis = opSmoothUnion(sphere1Dis, sphere2Dis, k);
    dis = opSmoothUnion(dis, sphere3Dis, k);
    dis = opSmoothUnion(dis, sphere4Dis, k);
    dis = opSmoothUnion(dis, waveFloor, k);
    
    // Subtract from edge (cut cylindrical platform)
    float cy1Dis = sdCylinder(pos, vec3(0, -20.0, 0), vec3(0, 20.0, 0), 6.5);
    float cy2Dis = sdCylinder(pos, vec3(0, -20.0, 0), vec3(0, 20.0, 0), MAXDIS);
    float finalCyDis = opSubtraction(cy1Dis, cy2Dis);
    
    return opSmoothSubtraction(finalCyDis, dis * 0.8, 0.2);
}

float sdSandFloor(vec3 pos) {
    // Create a sand floor as a cylinder with bump map on y axis
    float bump = 0.0;
    bump += texture(iChannel1, pos.xz * 0.1).r * 0.25;
    bump += texture(iChannel1, 10.0 + pos.xz * 2.0).r * 0.01;

    float dis = sdCylinder(pos, vec3(0, -1.4 - bump, 0), vec3(0, -10.0, 0), 6.5);//pos.y + 1.4;
    return dis;
}

float sdCoral(vec3 pos) {
    // Return a bunch of spheres
    pos.xz += vec2(0.1, 0.8);
    float s1 = sdSphere(pos, vec3(-1.2, -1.2, 3.2), 0.7);
    float s2 = sdSphere(pos, vec3(-1.5, -1.2, 2.1), 0.35);
    float s3 = sdSphere(pos, vec3(1.6, -1.2, 2.3), 0.35);
    float s4 = sdSphere(pos, vec3(-1.1, -1.2, 2.6), 0.5);
    float s5 = sdSphere(pos, vec3(1.3, -3.0, 3.2), 2.0);
    float s6 = sdSphere(pos, vec3(0.0, -2.0, -1.2), 1.1);
    float s7 = sdSphere(pos, vec3(1.3, -1.6, -1.1), 0.5);
    float s8 = sdSphere(pos, vec3(-1.9, -1.5, 0.0), 0.4);
    
    float dis = opUnion(s1, s2);
    dis = opUnion(dis, s3);
    dis = opUnion(dis, s4);
    dis = opUnion(dis, s5);
    dis = opUnion(dis, s6);
    dis = opUnion(dis, s7);
    dis = opUnion(dis, s8);
    
    // Use texture to add bumps
    dis += texture(iChannel1, pos.xz / 2.5).r * 0.10;
    dis += texture(iChannel1, pos.xy / 1.0).r * 0.05;
    dis += texture(iChannel1, pos.yz * 0.5).r * 0.005;
    return dis * 0.9;
}

// Vectors //////////////////////////////////////////////////////////////////////

// Get distance function
vec2 getDist(vec3 pos, bool includeWater) {
    
    // Map dist to shapes
    float oceanFloorDis = sdSandFloor(pos); 
    float oceanSurfDis = sdOceanSurf(pos);
    float oceanCoralDis = sdCoral(pos);
    
    // Final distance to return
    float finalDis = opUnion(oceanFloorDis, oceanCoralDis);
    if (includeWater) {
        finalDis = opUnion(oceanSurfDis, finalDis);
    }
    
    
    // Final material to return
    int mat = MATWATER;
    if (finalDis == oceanFloorDis) {
        mat = MATSAND;
    } else if (finalDis == oceanCoralDis) {
        mat = MATCORAL;
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
    vec3 lightPos = lightOrigin; // basically a point in 3D space
    //lightPos.xz += vec2(sin(iTime), cos(iTime)) * 4.0;
    
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

// Main //////////////////////////////////////////////////////////////////////
vec3 getMatCol(int mat) {
    if (mat == MATSAND) {
        return sandCol;
    } else if (mat == MATCORAL) {
        return coralCol;
    } else if (mat == MATWATER) {
        return waterCol;
    } 
    return vec3(0);
}

vec3 background(vec3 rayDir) {
    vec3 col = vec3(0.0);
    
    float y = rayDir.y * 0.2 + 0.8; // light is top, dark is bottom. 
    col += y * bgCol1;
    float x = rayDir.x * 0.6 + 0.4;
    col += x * bgCol2 * 0.4;
    
    return col;
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

float getLightIntensity(float dif, vec3 lightPos, vec3 pos, float maxLightLen) {
    return dif * clamp(1.0 - (length(lightPos - pos) / maxLightLen), 0.0, 1.0);
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

vec3 renderVolume(vec3 entryPos, vec3 rayDir) {
    // Fixed step ray marcher, calc how much light passes through per step
    float opaqueVis = 1.0;
    float marchSize = 0.1;
    float shadowMarchSize = 0.2;
    
    float volumeDepth = 0.0;
    int volumeMaxSteps = 36;
    int shadowMaxSteps = 20;
    float absorptionCoefficient = 20.0;
    vec3 volumeCol = vec3(0.0);
    
    for (int i = 0; i < volumeMaxSteps; i++) {
        volumeDepth += marchSize;
        vec3 insidePos = entryPos + rayDir * volumeDepth;
        
        // Get absorption per step
        float prevOpaqueVis = opaqueVis;
        opaqueVis *= beerLambertAttenuation(absorptionCoefficient, marchSize);
        float currAbsorption = prevOpaqueVis - opaqueVis;
        
        // Compare with diffuse lights (oh boy)
        vec3 keyLightDir = keyLightPos - insidePos;
        float keyLightDis = length(keyLightDir);
        float keyLightVis = getLightVis(insidePos, keyLightDir, keyLightDis, shadowMaxSteps, shadowMarchSize, absorptionCoefficient);
        vec3 outKeyLightCol = keyLightVis * keyLightCol * getLightAttenuation(keyLightDis);
        volumeCol += currAbsorption * waterCol * outKeyLightCol;
        
        // Compare with diffuse lights (oh boy)
        vec3 fillLightDir = fillLightPos - insidePos;
        float fillLightDis = length(fillLightDir);
        float fillLightVis = getLightVis(insidePos, fillLightDir, fillLightDis, shadowMaxSteps, shadowMarchSize, absorptionCoefficient);
        vec3 outfillLightCol = fillLightVis * fillLightCol * getLightAttenuation(fillLightDis);
        volumeCol += currAbsorption * waterCol * outfillLightCol;

        // Ambient light
        volumeCol += currAbsorption * waterCol * AMBIENT;
        
        // Material underneath (water caustics -> faked with voronoi)
        vec2 passedDMat = getDist(insidePos, false);
        float disToScene = passedDMat.x;
        int matInterior = int(passedDMat.y);
        if (disToScene < SURFDIS) {
            if (matInterior != MATWATER) {
                vec3 intCol = getMatCol(matInterior);
                // First mix the water colour with the interior below
                volumeCol = mix(volumeCol, intCol, opaqueVis + 0.2);
                // Then add caustics
                vec3 causticCol = 0.3 * volumeDepth * getVoronoi(insidePos.xz * 1.0) * skyCol;
                volumeCol += causticCol;
            }
            return volumeCol;
        }
    }
    
    // If we've reached this stage, this means that we've exited the water 
    // and now need to merge the background col
    float maxWidth = 0.6;
    float stepper = clamp(volumeDepth / maxWidth, 0.0, 1.0);
    return mix(volumeCol, background(rayDir), maxWidth);
}

vec3 render(vec3 rayOrigin, vec3 rayDir) {
    // Assign col
    vec3 col = vec3(0.0);

    // Lighting Setup
    float lightMaxLen = 12.0;
    vec3 ambientLight = vec3(AMBIENT);
    vec3 mixCol = vec3(ambientLight); 
    
    // Visualise 
    vec2 passedDMat = rayMarch(rayOrigin, rayDir, true);
    float dis = passedDMat.x;
    int mat = int(passedDMat.y);
    float IOR = 1.33;
    
    if (dis < MAXDIS) {
        vec3 pos = rayOrigin + rayDir * dis;
        vec3 normal = getNormal(pos, true);
        vec3 posEnter = pos - normal * SURFDIS * 3.0;
        
        // Lock the lights based on their distance to scene
        float keyDiffuseLight = getPointLight(pos, keyLightPos, true);
        float fillDiffuseLight = getPointLight(pos, keyLightPos, true); 
        float finalFill = getLightIntensity(fillDiffuseLight, fillLightPos, pos, lightMaxLen);
        float finalKey = getLightIntensity(keyDiffuseLight, keyLightPos, pos, lightMaxLen);

        // Materials
        if (mat == MATNULL) {
            // Checkerboard pattern
            float size = 1.0;
            float total = mod(floor(pos.z * size), 2.0) + 
                          mod(floor(pos.x * size), 2.0);
            bool isEven = mod(total, 2.0) == 0.0;
            mixCol *= (isEven) ? shadowCol : lightCol;
        } else if (mat == MATWATER) {
            vec3 rayDirIn = refract(rayDir, normal, 1.0 / IOR);
            mixCol = renderVolume(posEnter, rayDirIn);
            
        } else if (mat == MATSAND) {
            mixCol *= sandCol;
        } else if (mat == MATCORAL) {
            vec3 finalCoralCol;
            vec3 upDir = vec3(0, 1, 0);
            float angle = dot(upDir, normal);
            finalCoralCol = mix(coralCol, algaeCol, smoothstep(0.7, 0.72, angle));
            
            mixCol *= finalCoralCol;
        }
        
        // Ambient Occlusion (check out Alex Evans)
        // Trick is to sample sdf along normal
        float ao = getAmbientOcclusion(pos, normal, true);
        mixCol *= ao;
        
        // Lighting
        if (mat != MATWATER) {
            mixCol += finalKey * keyLightCol;
            mixCol += finalFill * fillLightCol;
        }
        
        float specStrength = 0.005;
        vec3 lightReflectRayDir = reflect(keyLightPos - posEnter, normal);
        float spec = pow(max(dot(rayDir, lightReflectRayDir), 0.0), 2.0);
        float fresnel = clamp(dot(normal, (keyLightPos - posEnter)), 0.0, 1.0);
        vec3 specularCol = specStrength * spec * keyLightCol * fresnel;
        mixCol += specularCol;
        
        specStrength = 0.002;
        lightReflectRayDir = reflect(fillLightPos - posEnter, normal);
        spec = pow(max(dot(rayDir, lightReflectRayDir), 0.0), 2.0);
        fresnel = clamp(dot(normal, (fillLightPos - posEnter)), 0.0, 1.0);
        specularCol = specStrength * spec * fillLightCol * fresnel;
        mixCol += specularCol;
        
        col = vec3(mixCol);
    }
    
    // Fog / Sky
    float fogStart = MAXDIS / 4.0;
    float fogEnd = MAXDIS;
    col = mix(col, background(rayDir), smoothstep(fogStart, fogEnd, dis));
    
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
    float camZoom = 1.0;
    vec3 rayOrigin = vec3(camX, camY, camZ); // this is the camera (origin of vector)
    rayOrigin.yz *= rotate(-m.y * PI + 1.);
    rayOrigin.xz *= rotate(-m.x * 2.0 * PI);
    vec3 rayDir = getRayDir(uv, rayOrigin, vec3(camYTilt, 0, camXTilt), camZoom);
    
    col = render(rayOrigin, rayDir);
    
    // Output to screen
    fragColor = vec4(col,1.0);
}
