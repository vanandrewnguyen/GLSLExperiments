/*
Van Andrew Nguyen
9/12/23
[Fur]

Hypertexture using texture lod to extrude sdf's and create a basic fur texture.
*/

// https://www.shadertoy.com/view/3lKyDc

// Globals //////////////////////////////////////////////////////////////////////

/*
Max steps is the maximum amount of steps we can make to get a ray length
Max distance is the maximum distance we can shoot a ray
Surf distance is the distance we loop to get to register a hit on the surface
*/
#define MAXSTEPS 100
#define MAXDIS 48.0
#define SURFDIS 0.001

#define PI 3.1416

#define MATNULL 0
#define MATFUR 1
#define MATGLOBE 2

// Colours
const vec3 shadowCol = vec3(0.2, 0.1, 0.15);
const vec3 lightCol = vec3(0.8, 0.75, 0.7);
const vec3 fillLightCol = vec3(0.854, 0.364, 0.501);
const vec3 keyLightCol = vec3(0.937, 0.376, 0.101);
const vec3 ambientLightCol = vec3(0.50,0.70,1.00);

const vec3 skyCol = vec3(0.6, 0.7, 0.8);
const vec3 sunCol = vec3(1.0, 0.6, 0.3) * 3.0;

// Positions
const vec3 sunPos = vec3(9.0, 4.2, 0.0);

// Noise /////////////////////////////////////////////////////////////////////

// Return a randomised float using a vec2 (uv coord)
float hash1(float n) {
    return fract( n*17.0*fract( n*0.3183099 ) );
}

float hash21(vec2 uv) {
    return fract(hash1(uv.x * uv.y));
}

float hash31(vec3 p) {
	vec3 ip=floor(p);
	p=fract(p);
	p=smoothstep(0.0,1.0,p);
	vec3 st=vec3(7,37,289);
	vec4 pos=dot(ip,st) + vec4(0.0,st.y,st.z,st.y+st.z);
	vec4 val=mix(fract(sin(pos)*7894.552), fract(sin(pos+st.x)*7894.552), p.x);
	vec2 val2=mix(val.xz,val.yw, p.y);
	return mix(val2.x,val2.y, p.z);
}

vec3 hash33(vec3 p) {
	float x=hash31(p);
	float y=hash31(p+sin(hash31(p)));
	float z=hash31(p+cos(hash31(p)));
	return vec3(x,y,z);
}

float dot2(vec3 v) {
    return dot(v,v);
}

float furNoise(in vec3 x)
{
    vec3 i = floor(x);
    vec3 f = fract(x);
	f = f*f*(3.0-2.0*f);
	vec2 uv = (i.xy + vec2(37.0,17.0) * i.z) + f.xy;
	vec2 rg = textureLod( iChannel0, (uv+0.5)/256.0, 1.0).yx;
	return mix(rg.x, rg.y, f.z);
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

// Vectors //////////////////////////////////////////////////////////////////////

// Get distance function
vec2 getDist(vec3 pos) {
    
    // Map dist to shapes
    float planeDis = pos.y;
    float torusDis = sdTorus(pos, vec3(0, 1.0, 3), vec2(1.0, 0.1));
    float sphereDis = sdSphere(pos, vec3(0, 1.5, 0), 1.0);
    float cylinderDis = sdCylinder(pos, vec3(0, 1.0, -3), vec3(0, 2.0, -3), 0.5);
  
    // Final distance to return
    float finalDis = opUnion(planeDis, sphereDis);
    finalDis = opUnion(finalDis, torusDis);
    finalDis = opUnion(finalDis, cylinderDis);
    
    // Final material to return
    int mat = MATNULL;
    if (finalDis == torusDis || finalDis == sphereDis || finalDis == cylinderDis) {
        mat = MATFUR;
    } 
    // if (finalDis == globeDis) {
    //     mat = MATGLOBE;
    // }
    
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
        
        //  Temp fur code... slow
        if (mat == MATFUR || mat == MATGLOBE) {
            vec3 pos = rayOrigin + rayDir * disFromOrigin;
            vec3 normal = getNormal(pos);
            float fre = (mat == MATFUR) ? 1.0 : -dot(normal, rayDir);
            float amp = 0.6;
            disFromOrigin -= fre * amp * furNoise(normal * 60.0);
        }
        
        disFromOrigin += disToScene;
        
        // Exit condition; we have a hit
        if (disFromOrigin > MAXDIS || disToScene < SURFDIS) { 
            break;
        }
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
    
    // Get diff
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

vec3 generateBackgroundCol(in vec3 rayDir) {
    vec3 col = vec3(0.0);
    
    // Fake scattering gradient
    float y = rayDir.y * 0.3 + 0.7; // light is top, dark is bottom. 
    col += y * skyCol;
    
    // Fake sun
    vec3 sunDir = normalize(sunPos);
    col += pow(max(dot(rayDir, normalize(sunDir)), 0.0), 8.0) * sunCol * 0.15;
    col += pow(max(dot(rayDir, normalize(sunDir)), 0.0), 1024.0) * sunCol;
    
    //float x = rayDir.x * 0.5 + 0.5;
    //col += x * sunRayCol * 0.6;
    
    return col;
}

// Main //////////////////////////////////////////////////////////////////////

vec3 render(vec3 rayOrigin, vec3 rayDir) {
    // Lighting Setup
    vec3 keyLightPos = sunPos; // vec3(3, 4, 2);
    vec3 fillLightPos = vec3(-2, 8, 2);
    float lightMaxLen = 10.0;
    vec3 brdf = vec3(0.0); 
    
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
        // float finalKey = getLightIntensity(keyDiffuseLight, keyLightPos, pos, lightMaxLen * 0.7);

        // Lighting
        vec3 diff = keyDiffuseLight * keyLightCol;
        vec3 fill = finalFill * fillLightCol;
        float fre = pow(clamp(1.0 + dot(normal, rayDir), 0.0, 1.0), 2.0);
        float spec = pow(max(dot(rayDir, keyLightPos - pos), 0.0), 1.0);
        float amb = 0.5;
        float ao = getAmbientOcclusion(pos, normal);
          
        // Lighting https://www.shadertoy.com/view/4tSGRw
        brdf += 1.20*diff*vec3(1.00,0.90,0.60);
        brdf += 0.80*spec*vec3(1.00,0.90,0.60)*diff;
        brdf += 0.30*amb*ambientLightCol*ao;
        brdf += 0.40*fre*vec3(1.00,1.00,1.00)*ao;

        // Materials
        if (mat == MATNULL) {
            // Checkerboard pattern
            float size = 1.0;
            float total = mod(floor(pos.z * size), 2.0) + 
                          mod(floor(pos.x * size), 2.0);
            bool isEven = mod(total, 2.0) == 0.0;
            brdf *= (isEven) ? shadowCol : lightCol;
        } else if (mat == MATFUR) {
            // Triplanar mapping
            vec3 colXZ = texture(iChannel1, pos.xz*0.5+0.5).rgb;
            vec3 colYZ = texture(iChannel1, pos.yz*0.5+0.5).rgb;
            vec3 colXY = texture(iChannel1, pos.xy*0.5+0.5).rgb;
            vec3 newNorm = abs(normal);
            newNorm *= pow(newNorm, vec3(2));
            newNorm *= newNorm.x + newNorm.y + newNorm.z; // normalise
            brdf *= 0.1 + (colYZ * newNorm.x + colXZ * newNorm.y + colXY * newNorm.z);
        }

        return brdf;
    }
   
    return generateBackgroundCol(rayDir);
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
