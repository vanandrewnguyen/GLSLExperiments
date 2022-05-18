/*
Van Andrew Nguyen
18/05/22
[Tiny Planet]

http://casual-effects.com/research/McGuire2019ProcGen/McGuire2019ProcGen.pdf
Added some new tricks to an old render -> volumetric atmosphere, stars, proper texturing (triplanar) and coords.
*/

// Globals //////////////////////////////////////////////////////////////////////

/*
Max steps is the maximum amount of steps we can make to get a ray length
Max distance is the maximum distance we can shoot a ray
Surf distance is the distance we loop to get to register a hit on the surface
*/
#define MAXSTEPS 100
#define MAXDIS 24.0
#define SURFDIS 0.01

#define PI 3.1416

#define MATNULL 0
#define MATPLANET 1
#define MATWATER 2

const vec3 shadowCol = vec3(0.2, 0.1, 0.15);
const vec3 lightCol = vec3(0.8, 0.75, 0.7);
const vec3 fillLightCol = vec3(0.854, 0.364, 0.501);
const vec3 keyLightCol = vec3(0.937, 0.376, 0.101);

const vec3 dirtCol = vec3(0.117, 0.117, 0.101);
const vec3 grassCol = vec3(0.293, 0.350, 0.257);
const vec3 snowCol = vec3(0.882, 0.921, 0.937);
const vec3 waterCol = vec3(0.094, 0.224, 0.279);
const vec3 sandCol = vec3(0.890, 0.850, 0.529);
const vec3 atmosphereCol = vec3(0.596, 0.901, 0.870);

const vec3 bgCol1 = vec3(1);
const vec3 bgCol2 = vec3(0.4);

// Noise /////////////////////////////////////////////////////////////////////

// Return a randomised float using a vec2 (uv coord)
float hash11(float n) {
    return fract( n*17.0*fract( n*0.3183099 ) );
}

float hash21(vec2 uv) {
    return fract(hash11(uv.x * uv.y));
}

vec3 hash33(vec3 p3) {
	p3 = fract(p3 * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yxz+19.19);
    return fract(vec3((p3.x + p3.y)*p3.z, (p3.x+p3.z)*p3.y, (p3.y+p3.z)*p3.x));
}

// 3d Noise function from: https://www.shadertoy.com/view/Xsl3Dl
float noise3D( in vec3 p ) {
    vec3 i = floor( p );
    vec3 f = fract( p );
	
	vec3 u = f*f*(3.0-2.0*f);

    return mix( mix( mix( dot( hash33( i + vec3(0.0,0.0,0.0) ), f - vec3(0.0,0.0,0.0) ), 
                          dot( hash33( i + vec3(1.0,0.0,0.0) ), f - vec3(1.0,0.0,0.0) ), u.x),
                     mix( dot( hash33( i + vec3(0.0,1.0,0.0) ), f - vec3(0.0,1.0,0.0) ), 
                          dot( hash33( i + vec3(1.0,1.0,0.0) ), f - vec3(1.0,1.0,0.0) ), u.x), u.y),
                mix( mix( dot( hash33( i + vec3(0.0,0.0,1.0) ), f - vec3(0.0,0.0,1.0) ), 
                          dot( hash33( i + vec3(1.0,0.0,1.0) ), f - vec3(1.0,0.0,1.0) ), u.x),
                     mix( dot( hash33( i + vec3(0.0,1.0,1.0) ), f - vec3(0.0,1.0,1.0) ), 
                          dot( hash33( i + vec3(1.0,1.0,1.0) ), f - vec3(1.0,1.0,1.0) ), u.x), u.y), u.z );
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

// Custom SDF's
float sdPlanet(vec3 pos) {
    float dis;
    dis = sdSphere(pos, vec3(0), 1.95);
    
    // Apply fbm (6 octives, 3rd pow)
    float elevation = 0.0;
    float amp = 1.0;
    float freq = 1.0;
    for (int i = 0; i < 6; i++) {
        elevation += abs(noise3D(pos * freq + float(i) * 10.0) * amp);
        
        amp *= 0.3;
        freq *= 2.0;
    }
    
    elevation *= 1.1;
    
    dis -= elevation;
    return dis * 0.7;
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
    float f = k * (dot(dir, coord.xy) - c * iTime * 1.0);
    
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

// Water SDF
float sdWater(vec3 pos) {
    float dis;
    vec3 waveDis = gerstnerFBM(pos * 0.15) * 0.05;
    dis = sdSphere(pos + waveDis, vec3(0), 2.05);
    //dis = sdSphere(pos, vec3(0), 2.05);
    
    return dis;
}

// Cloud SDF
float sdCloud(vec3 pos) {
    float dis;
    dis = sdSphere(pos, vec3(0), 2.3);
    
    // Apply fbm (6 octives, 3rd pow)
    float elevation = 0.0;
    float amp = 1.0;
    float freq = 1.0;
    for (int i = 0; i < 3; i++) {
        elevation += abs(noise3D(pos * freq + float(i) * 10.0) * amp);
        
        amp *= 0.3;
        freq *= 2.0;
    }
    
    elevation *= 0.5;
    
    dis += elevation;
    
    return dis * 0.9;
}

// Vectors //////////////////////////////////////////////////////////////////////

// Get distance function
vec2 getDist(vec3 pos, bool sphereOnly) {
    
    // Map dist to shapes
    float planetDis = sdPlanet(pos);
    float waterDis = sdWater(pos);
    float sphereDis = sdSphere(pos, vec3(0), 2.2);
    
    // Sphere is for calculating the surface normal of a flat planet, to get the steepness of a mountain
    
    // Final distance to return
    float finalDis = (sphereOnly) ? sphereDis : opUnion(planetDis, waterDis);
    
    // Final material to return
    int mat = MATPLANET;
    if (finalDis == waterDis) {
        mat = MATWATER;
    }
    
    return vec2(finalDis, mat);
}

float getDensity(vec3 pos) {
    // Get sdf, but this time we march inside
    float cloudDis = sdCloud(pos);
    return cloudDis;
}

// Return the normal ray
vec3 getNormalDist(vec3 pos, bool sphereOnly) {
    float dis = getDist(pos, sphereOnly).x;
    vec2 val = vec2(0.01, 0.0);
    
    // To get the slope we give the curve two values super close together
    // Instead of deriving we can do this method via swizzling
    vec3 normal = dis - vec3(getDist(pos-val.xyy, sphereOnly).x, 
                             getDist(pos-val.yxy, sphereOnly).x, 
                             getDist(pos-val.yyx, sphereOnly).x);
    // Same way of doing getDist(p-vec3(0.01, 0.0, 0.0), ... etc
    return normalize(normal);
}

vec3 getNormalDensity(vec3 pos) {
    float dis = getDensity(pos);
    vec2 val = vec2(0.01, 0.0);

    vec3 normal = dis - vec3(getDensity(pos-val.xyy), 
                             getDensity(pos-val.yxy), 
                             getDensity(pos-val.yyx));
    return normalize(normal);
}

// Ray Marching function
vec2 rayMarch(vec3 rayOrigin, vec3 rayDir, bool sphereOnly) {
    float disFromOrigin = 0.0;
    int mat;
    
    // Now we loop for max steps or until the distance is very small
    for (int i=0;i<MAXSTEPS;i++) {
        // Declare current marching location
        vec3 pos = rayOrigin + rayDir * disFromOrigin;
        vec2 passedDMat = getDist(pos, sphereOnly);
        float disToScene = passedDMat.x;
        mat = int(passedDMat.y);
        
        disFromOrigin += disToScene;
        
        // Exit condition; we have a hit
        if (disFromOrigin > MAXDIS || disToScene < SURFDIS) { break; }
    }
    
    return vec2(disFromOrigin, mat);
}

vec4 rayMarchVolume(vec3 rayOrigin, vec3 rayDir) {
    // March to surf dis
    float disFromOrigin = 0.0;
    vec3 surfPos = vec3(0.0);
    for (int i=0;i<MAXSTEPS;i++) {
        vec3 pos = rayOrigin + rayDir * disFromOrigin;
        float disToScene = getDensity(pos);

        disFromOrigin += disToScene;
        
        // Exit condition
        if (disFromOrigin > MAXDIS || disToScene < SURFDIS) { 
            surfPos = pos;
            break; 
        }
    }
    
    return vec4(disFromOrigin, surfPos);
}

// Soft shadows
// From https://www.iquilezles.org/www/articles/rmshadows/rmshadows.htm
float getSoftShadow(vec3 pos, vec3 normal, vec3 lightDir, float shadowInt) {
    float res = 1.0;
    float dis = 0.0;
    float t = SURFDIS;
    vec3 startingPos = pos + normal * SURFDIS * 2.0;
    for(int i = 0; i < MAXSTEPS; i++) {
        dis = getDist(startingPos + lightDir * t, false).x;
        if (dis < SURFDIS) { return 0.0; }
        if (t >= MAXDIS) { return res; }
        res = min(res, shadowInt * dis / t);
        t += dis;
    }
    
    return res;
}

// Hard shadows
void getHardShadow(vec3 pos, vec3 normal, vec3 lightDir, vec3 lightPos, inout float dif) {
    float dis = rayMarch(pos + normal * SURFDIS * 2.0, lightDir, false).x;
    
    // Getting shadows (if disToSphere < disToLight then we have a shadow)
    if (dis < length(lightPos - pos)) {
        float shadowIntensity = 0.2;
        dif *= shadowIntensity;
    }
}

// Lighting //////////////////////////////////////////////////////////////////////

float getPointLight(vec3 pos, vec3 lightOrigin) {
    // Get the light origin
    vec3 lightPos = lightOrigin; // basically a point in 3D space
    //lightPos.xz += vec2(sin(iTime), cos(iTime)) * 4.0;
    
    // Get the light ray
    vec3 lightDir = normalize(lightPos - pos);
    
    // Get normal
    vec3 normal = getNormalDist(pos, false);
    
    // Get diffused lighting 
    // We want 0 if the rays are parallel, 1 if the rays are perpendicular
    // Hence we use dot product
    float dif = clamp(dot(normal, lightDir), 0.0, 1.0);
 
    // Hard shadows
    getHardShadow(pos, normal, lightDir, lightPos, dif);
    dif *= getSoftShadow(pos, normal, lightDir, 8.0);
    
    return dif;
}

float getSpotLight(vec3 pos, vec3 lightOrigin, vec3 endLightDir, float angleMax, float lightDropOff) {
    // Identical to a point light, however, we restrict the light through a certain degree (use dot prod)
    
    // Get the light origin
    vec3 lightPos = lightOrigin;
    
    // Get angle dif
    vec3 lightDirToIntersect = normalize(lightPos - pos);
    vec3 lightDirToEnd = normalize(lightOrigin - endLightDir);
    float angleDif = dot(lightDirToIntersect, lightDirToEnd);
    
    // Get normal
    vec3 normal = getNormalDist(pos, false);
    
    float dif = clamp(dot(normal, lightDirToIntersect), 0.0, 1.0);
    dif *= smoothstep(angleMax, angleMax + lightDropOff, angleDif);
    
    // Hard shadows
    getHardShadow(pos, normal, lightDirToIntersect, lightPos, dif);
    dif *= getSoftShadow(pos, normal, lightDirToIntersect, 8.0);
    
    return dif;
}

float getAmbientOcclusion(vec3 pos, vec3 normal) {
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
        aoSum += inc * getDist(currPos, false).x;
        aoMaxSum += inc * offset;
    }
    
    // Make sure it's between 0->1
    return aoSum / aoMaxSum;
}

float getLightIntensity(float dif, vec3 lightPos, vec3 pos, float maxLightLen) {
    return dif * clamp(1.0 - (length(lightPos - pos) / maxLightLen), 0.0, 1.0);
}

// Main //////////////////////////////////////////////////////////////////////

float beerLambertAttenuation(float absorptionCoefficient, float dis) {
    // Square inverse law
    // In a clear volume like water, you want the absorption to be low since water doesn't absorb light
    return exp(-absorptionCoefficient * dis);
}

vec3 getSkyCol(vec3 rayDir, vec2 uv, vec2 fragCoord) {
    // Shift according to y val of ray dirction
    float k = rayDir.y * 0.5 + 0.5;
    
    // Colour shift based on perspective
    vec3 base1 = bgCol2;
    vec3 base2 = bgCol1;
    vec3 col = mix(base1, base2, k);
    
    // Radial Gradient
    float radius = 0.7;
    float len = length(uv) - radius;
    float atmosphericRadialAttenuation = min(1.0, 0.06 * pow(max(0.0, 1.0 - len / radius), 8.0));
    col = mix(col, atmosphericRadialAttenuation * atmosphereCol, 0.8);
    
    // Radial noise
    float radialNoise = mix(1.0, hash21(normalize(uv) * 40.0 + iTime * 0.00001), 0.05);
    col *= radialNoise;
    
    // Stars
    float scaling = 2.0;
    float clumps = pow(hash21(fragCoord * scaling), 9.0) * 0.5 + pow(hash21(100.0 + fragCoord * scaling * 4.0), 5.0) / 1.5;
    vec3 starCol = vec3(clumps * pow(hash21(fragCoord), 300.0) * 80.0);
    starCol.r *= sqrt(hash21(fragCoord) * 1.2);
    starCol.g *= sqrt(hash21(fragCoord * 4.0));
    starCol *= hash21(iTime * 0.00005 + uv.yx * 10.0);
    
    col += starCol;
    
    return col;
}

vec3 renderVolumetric(vec3 rayOrigin, vec3 rayDir, vec3 inputCol) {
    vec3 col = inputCol;
    
    // March to surf dis
    vec4 density = rayMarchVolume(rayOrigin, rayDir);
    float disFromOrigin = density.x;
    vec3 surfPos = density.yzw;
    
    // Volumetric rendering (kind of scuffed, need to revise)
    if (disFromOrigin < MAXDIS) {
        float opaqueVis = 1.0;
        float prevOpaqueVis = opaqueVis;
        float volumeDis = 0.0;
        float marchSize = 0.025;
        float dis = rayMarch(surfPos, rayDir, false).x;
        if (dis > MAXDIS) {
            vec3 volumeNormal = getNormalDensity(surfPos);
            dis = max(0.0, 0.1 + 100.0 * rayMarchVolume(surfPos - volumeNormal * SURFDIS * 1.0, rayDir).x);
        }
        
        while (volumeDis < dis) {
            opaqueVis *= beerLambertAttenuation(0.2, marchSize);
            float currAbsorption = prevOpaqueVis - opaqueVis;
            volumeDis += marchSize;
            col = mix(col, atmosphereCol, currAbsorption);
        }
        
        
    }
    
    return col;
}

vec3 render(vec3 rayOrigin, vec3 rayDir, vec2 uv, vec2 fragCoord) {
    // Assign col
    vec3 col = vec3(0.0);

    // Lighting Setup
    vec3 keyLightPos = vec3(3, 4, 2);
    vec3 fillLightPos = vec3(-6, 8, 2);
    float lightMaxLen = 10.0;
    vec3 ambientLight = vec3(0.2);
    vec3 mixCol = vec3(ambientLight); 
    
    // Visualise 
    vec2 passedDMat = rayMarch(rayOrigin, rayDir, false);
    float dis = passedDMat.x;
    int mat = int(passedDMat.y);
    
    if (dis < MAXDIS) {
        vec3 pos = rayOrigin + rayDir * dis;
        vec3 normal = getNormalDist(pos, false);
        
        // Lock the lights based on their distance to scene
        float keyDiffuseLight = getPointLight(pos, keyLightPos);
        float fillDiffuseLight = getSpotLight(pos, fillLightPos, vec3(0), PI / 3.5, 0.05); //getPointLight(pos, fillLightPos);
        float finalFill = getLightIntensity(fillDiffuseLight, fillLightPos, pos, lightMaxLen * 1.5);
        float finalKey = getLightIntensity(keyDiffuseLight, keyLightPos, pos, lightMaxLen * 0.7);

        // Materials
        if (mat == MATNULL) {
            // Checkerboard pattern
            float size = 1.0;
            float total = mod(floor(pos.z * size), 2.0) + 
                          mod(floor(pos.x * size), 2.0);
            bool isEven = mod(total, 2.0) == 0.0;
            mixCol *= (isEven) ? shadowCol : lightCol;
        } else if (mat == MATWATER) {
            mixCol = waterCol;
            
            // Ray march for depth
            float IOR = 1.33;
            vec3 rayDirIn = refract(rayDir, normal, 1.0 / IOR);
            vec3 posEnter = pos - normal * SURFDIS * 3.0;
            float disIn = rayMarch(posEnter, rayDirIn, false).x;
            float heightStepper = pow(length(disIn), 0.5) + 2.0;
            mixCol = mix(dirtCol, mixCol, heightStepper);
            
            // Spec
            float specStrength = 0.005;
            vec3 lightReflectRayDir = reflect(keyLightPos - posEnter, normal);
            float spec = pow(max(dot(rayDir, lightReflectRayDir), 0.0), 4.0);
            float fresnel = clamp(dot(normal, (keyLightPos - posEnter)), 0.0, 1.0);
            vec3 specularCol = specStrength * spec * keyLightCol * fresnel;
            mixCol += specularCol;
            
            specStrength = 0.005;
            lightReflectRayDir = reflect(fillLightPos - posEnter, normal);
            spec = pow(max(dot(rayDir, lightReflectRayDir), 0.0), 2.0);
            fresnel = clamp(dot(normal, (fillLightPos - posEnter)), 0.0, 1.0);
            specularCol = specStrength * spec * fillLightCol * fresnel;
            mixCol += specularCol;
            
        } else if (mat == MATPLANET) {
            mixCol = dirtCol;
            
            vec3 colXZ = texture(iChannel0, pos.xz*0.5+0.5).rgb;
            vec3 colXY = texture(iChannel0, pos.xy*0.5+0.5).rgb;
            vec3 colYZ = texture(iChannel0, pos.yz*0.5+0.5).rgb;
            vec3 texCol = colYZ * normal.x + 
                          colXY * normal.z + 
                          colXZ * normal.y;
            
            vec3 planetCentre = vec3(0);
            float disFromCentre = length(pos - planetCentre);
            float landHeight = disFromCentre - 2.05; // water dis
            vec3 sphereNormal = getNormalDist(pos, true);
                float steepness = abs(dot(sphereNormal, normal));
            if (landHeight > 0.13) {
                if (steepness > 0.87) {
                    mixCol = mix(dirtCol, snowCol - (texCol * 0.1), 1.0-smoothstep(1.0, 0.8, steepness));
                }
            } else if (landHeight > 0.02) {
                mixCol = dirtCol * texCol;
                if (steepness > 0.88) {
                    mixCol = mix(mixCol, grassCol * texCol, 1.0-smoothstep(1.0, 0.9, steepness));
                }
            } else {
                mixCol = sandCol - (texCol * 0.2);
            }
            
        }
        
        // Ambient Occlusion (check out Alex Evans)
        // Trick is to sample sdf along normal
        float ao = getAmbientOcclusion(pos, normal);
        mixCol *= ao;
        
        // Lighting
        mixCol += finalKey * keyLightCol;
        mixCol += finalFill * fillLightCol;
        
        col = vec3(mixCol);
    } else {
        // Render bg
        col = getSkyCol(rayDir, uv, fragCoord);
    }
    
    col = renderVolumetric(rayOrigin, rayDir, col);
    
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
    float camZoom = 2.0;
    vec3 rayOrigin = vec3(camX, camY, camZ); // this is the camera (origin of vector)
    rayOrigin.yz *= rotate(-m.y * PI + 1.);
    rayOrigin.xz *= rotate(-m.x * 2.0 * PI);
    vec3 rayDir = getRayDir(uv, rayOrigin, vec3(camYTilt, 0, camXTilt), camZoom);
    
    col = render(rayOrigin, rayDir, uv, fragCoord);
    
    // Output to screen
    fragColor = vec4(col,1.0);
}
