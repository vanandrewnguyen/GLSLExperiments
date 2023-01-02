// Globals //////////////////////////////////////////////////////////////////////
#define MAXSTEPS 200
#define MAXDIS 24.0
#define SURFDIS 0.01

#define PI 3.1416

#define MATNULL 0
#define MATREDCABLE 1
#define MATREDCABLEHIGHLIGHT 2
#define MATCIRCUIT 3
#define MATCIRCUITHIGHLIGHT 4

const vec3 shadowCol = vec3(0.2, 0.1, 0.15);
const vec3 lightCol = vec3(0.8, 0.75, 0.7);
const vec3 fillLightCol = vec3(0.854, 0.364, 0.501);
const vec3 keyLightCol = vec3(0.937, 0.376, 0.101);

// Noise /////////////////////////////////////////////////////////////////////

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

// Return smooth noise 
float smoothNoise(vec2 uv) {
    // Create 1D random value
    vec2 index = uv;
    vec2 localUV = fract(index); 
    vec2 cellID = floor(index); 
    
    localUV = localUV*localUV*(3.0 - 2.0 * localUV); 
    float bl = hash21(cellID);
    float br = hash21(cellID + vec2(1, 0));
    float b = mix(bl, br, localUV.x);
    float tl = hash21(cellID + vec2(0, 1));
    float tr = hash21(cellID + vec2(1, 1));
    float t = mix(tl, tr, localUV.x);
    float noiseCol = mix(b, t, localUV.y);
        
    return noiseCol;
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
    /*
    pos = abs(pos) - size;
    float val = length(max(abs(pos) - size, 0.0));
    return val;
    */
    // tutorial implementation breaks old sdf?
    pos = abs(pos)-size;
	return length(max(pos, 0.))+min(max(pos.x, max(pos.y, pos.z)), 0.);
}

// Cyclinder Distance
float sdCylinder(vec3 pos, vec3 start, vec3 end, float rad) {
    vec3 ab = end - start;
    vec3 ap = pos - start;
    float projection = dot(ab, ap) / dot(ab, ab); 
    vec3 len = start + projection * ab;
    
    float x = length(pos - len) - rad;
    float y = (abs(projection - 0.5) - 0.5) * length(ab);
    float extLen = length(max(vec2(x, y), 0.0));
    float intLen = min(max(x, y), 0.0);
    return extLen + intLen;
}

float sdMinkowskiLen(vec2 pos, float k) {
    // pythag but with an exponent 'k'
    pos = abs(pos);
    return pow(pow(pos.x, k) + pow(pos.y, k), 1.0/k);
}

// Torus Distance
float sdTorus(vec3 pos, float r1, float r2, float coupling) {
    return length(vec2(sdMinkowskiLen(pos.xz, coupling) - r1, pos.y)) - r2;
}


// Vectors //////////////////////////////////////////////////////////////////////

// Get distance function
vec2 getDistTruchetLayer(vec3 pos, float coupling, float wobble, bool cable) {
    
    // Map dist to shapes
    // So trick is using abs(p.y) makes the plane 2d which is viewable from both top and bottom
    // Subtracting a thickness value (0.01) mitigates bugs with an infinitely thin plane intersection
    float gridDis = abs(pos.y) - 0.01; 
    
    // Domain repetition
    vec3 newP = pos;
    vec2 cellID = floor(newP.xz);
    newP.xz = fract(newP.xz) - 0.5;
    float boxDis = sdCube(newP, vec3(0.48));
    
    // Quadrants for repetition of partial torus (1 draw call for 2 quadrants)
    float rand = hash21(cellID);
    if (rand < 0.5) {
        newP.x *= -1.0;
    }
    float s = newP.x > -newP.z ? 1.0 : -1.0;
    newP.xz -= 0.5 * s;
    float angle = atan(newP.x, newP.z); // polar coord from -pi->pi
    // adjust y value based on polar coord, so we get a weave pattern
    float freq = 2.0;
    float amp = 0.05;
    newP.y += cos(angle * freq) * amp * wobble;
    float radius = 0.05;
    radius = cable ? mix(radius, 0.1, smoothstep(0.0, 0.2, cos(angle*freq))) : 
    mix(radius, 0.07, smoothstep(0.0, 0.5, cos(angle*10.0)));
    float torusDis = sdTorus(newP, 0.5, radius, coupling);
    
    // Final distance to return
    // float finalDis = opSubtraction(boxDis, gridDis);
    // finalDis = opUnion(finalDis, torusDis);
    float finalDis = torusDis;
    
    // Final material to return
    int mat = MATNULL;
    
    if (cable) {
        mat = MATREDCABLE;
        if (radius > 0.05) {
            mat = MATREDCABLEHIGHLIGHT;
        }
    } else {
        mat = MATCIRCUIT;
        if (sin(angle * 10.0) > 0.5) {
            mat = MATCIRCUITHIGHLIGHT;
        }
    }
    
    return vec2(finalDis, mat);
}

vec2 getDist(vec3 pos) {
    pos.x += iTime;

    float thickness = 0.05;
    vec2 truchetLayer1 = getDistTruchetLayer(pos, 1.2, 0.0, false);
    truchetLayer1.x -= thickness;
    vec2 truchetLayer2 = getDistTruchetLayer(pos - vec3(0.5, 0, 0.5), 2.8, 2.8, true);
    
    float ground = abs(pos.y + 0.3) - 0.01;
    float finalDis = opUnion(truchetLayer1.x, truchetLayer2.x);
    finalDis = opUnion(ground, finalDis);
    int finalMat = MATNULL;
    if (finalDis == truchetLayer2.x) {
        finalMat = int(truchetLayer2.y);
    } else if (finalDis == truchetLayer1.x) {
        finalMat = int(truchetLayer1.y);
    }
    
    vec2 finalVal = vec2(finalDis * 0.9, finalMat);
    
    return finalVal;
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

// Soft shadows
// From https://www.iquilezles.org/www/articles/rmshadows/rmshadows.htm
float getSoftShadow(vec3 pos, vec3 normal, vec3 lightDir, float shadowInt) {
    float res = 1.0;
    float dis = 0.0;
    float t = SURFDIS;
    vec3 startingPos = pos + normal * SURFDIS * 2.0;
    for(int i = 0; i < MAXSTEPS; i++) {
        dis = getDist(startingPos + lightDir * t).x;
        if (dis < SURFDIS) { return 0.0; }
        if (t >= MAXDIS) { return res; }
        res = min(res, shadowInt * dis / t);
        t += dis;
    }
    
    return res;
}

// Hard shadows
void getHardShadow(vec3 pos, vec3 normal, vec3 lightDir, vec3 lightPos, inout float dif) {
    float dis = rayMarch(pos + normal * SURFDIS * 2.0, lightDir).x;
    
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
    vec3 normal = getNormal(pos);
    
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
    vec3 normal = getNormal(pos);
    
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
        aoSum += inc * getDist(currPos).x;
        aoMaxSum += inc * offset;
    }
    
    // Make sure it's between 0->1
    return aoSum / aoMaxSum;
}

float getLightIntensity(float dif, vec3 lightPos, vec3 pos, float maxLightLen) {
    return dif * clamp(1.0 - (length(lightPos - pos) / maxLightLen), 0.0, 1.0);
}

// Main //////////////////////////////////////////////////////////////////////

vec3 render(vec3 rayOrigin, vec3 rayDir) {
    // Assign col
    vec3 col = vec3(0.0);

    // Lighting Setup
    vec3 keyLightPos = vec3(3, 4, 2);
    vec3 fillLightPos = vec3(-2, 8, 2);
    float lightMaxLen = 10.0;
    vec3 ambientLight = vec3(0.1);
    vec3 mixCol = vec3(ambientLight); 
    
    // Visualise 
    vec2 passedDMat = rayMarch(rayOrigin, rayDir);
    float dis = passedDMat.x;
    int mat = int(passedDMat.y);
    
    if (dis < MAXDIS) {
        vec3 pos = rayOrigin + rayDir * dis;
        vec3 normal = getNormal(pos);
        
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
        } else if (mat == MATREDCABLE) {
            mixCol *= vec3(5, 0, 0);
        } else if (mat == MATREDCABLEHIGHLIGHT) {
            mixCol *= vec3(3, 2, 0);
        } else if (mat == MATCIRCUIT) {
            mixCol *= vec3(0.9);
        } else if (mat == MATCIRCUITHIGHLIGHT) {
            mixCol *= vec3(0);
        }
        
        // Ambient Occlusion (check out Alex Evans)
        // Trick is to sample sdf along normal
        float ao = getAmbientOcclusion(pos, normal);
        mixCol *= ao;
        
        // Lighting
        mixCol += finalKey * keyLightCol;
        mixCol += finalFill * fillLightCol;
        
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
    float camZoom = 2.0;
    vec3 rayOrigin = vec3(camX, camY, camZ); // this is the camera (origin of vector)
    rayOrigin.yz *= rotate(-m.y * PI + 1.);
    rayOrigin.xz *= rotate(-m.x * 2.0 * PI);
    vec3 rayDir = getRayDir(uv, rayOrigin, vec3(camYTilt, 0, camXTilt), camZoom);
    
    col = render(rayOrigin, rayDir);
    
    // Output to screen
    fragColor = vec4(col,1.0);
}
