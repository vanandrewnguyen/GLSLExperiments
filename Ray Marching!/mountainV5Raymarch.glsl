/*
Van-Andrew Nguyen
6/1/23
Foggy Mountains V1

Unioned a volumetric fog layer with a basic mountain sdf, both pulling from noise textures. Nothing fancy, although could improve the fog sdf (basically just a plane)
to blend into the mountains better.
*/

// Globals //////////////////////////////////////////////////////////////////////
#define MAXSTEPS 200
#define MAXDIS 24.0
#define SURFDIS 0.01

#define PI 3.1416

#define MATNULL 0
#define MATVOLUMETRIC 1
#define MATMOUNTAIN 2

const vec3 shadowCol = vec3(0.2, 0.1, 0.15);
const vec3 lightCol = vec3(0.8, 0.75, 0.7);
const vec3 skyCol = vec3(0.6, 0.7, 0.8);
const vec3 sunCol = vec3(1.0, 0.6, 0.3) * 3.0;

const vec3 dirtCol = vec3(0.2, 0.188, 0.180);
const vec3 snowCol = vec3(0.882, 0.921, 0.937) * 0.6;
const vec3 blueFogCol = vec3(0.729, 0.909, 0.992);
const vec3 grassCol = vec3(0.223, 0.270, 0.207);
const vec3 sunRayCol = vec3(0.980, 0.633, 0.480);

// Noise /////////////////////////////////////////////////////////////////////
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

// Return smooth noise 
float noise(vec2 uv) {
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

// Credits: stoox https://www.shadertoy.com/view/lttyWl
float noise(vec3 x) {
	vec3 p = floor(x);
	vec3 f = fract(x);

	f = f * f * (3.0 - 2.0 * f);
	float n = p.x + p.y * 55.0 + p.z * 101.0 ;

return mix(
	mix(
		mix(hash11(n), hash11(n + 1.0), f.x),
		mix(hash11(n+55.0), hash11(n + 56.0), f.x),
		f.y),
	mix(
		mix(hash11(n+101.0), hash11(n + 102.0), f.x),
		mix(hash11(n+156.0), hash11(n + 157.0), f.x),
		f.y),
	f.z);
}

float fbm3DNoise(vec3 pos, int octaves) {
    float res = 0.0;
    float freq = 1.0;
    float amp = 1.0;
    for (int i = 0; i < octaves; i++) {
        res += amp * noise(pos * freq);
        
        freq *= 2.0;
        amp *= 0.5;
    }
    
    return res;
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
float sdTorus(vec3 pos, float r1, float r2) {
    return length(vec2(length(pos.xz) - r1, pos.y)) - r2;
}

// Cube Distance
float sdCube(vec3 pos, vec3 size) { 
    pos = abs(pos)-size;
	return length(max(pos, 0.))+min(max(pos.x, max(pos.y, pos.z)), 0.);
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

float sdMountain(vec3 pos) {
    float res = pos.y - 1.8;
    
    float amp = 4.5;
    float freq = 0.3;
    for (int i = 0; i < 3; i++) {
        res += texture(iChannel1, 0.1 * pos.xz * freq).r * amp;
        amp *= 0.3;
        freq *= 1.6;
    }
    
    return res * 0.55;
}

// Vectors //////////////////////////////////////////////////////////////////////

// Get distance function
vec2 getDist(vec3 pos, bool includeFog) {
    
    // Map dist to shapes
    //vec3 movePos = vec3(iTime + pos.x, pos.yz);
    float noiseDis = 0.2 * smoothstep(0.0, 1.5, fbm3DNoise(pos, 3));
    //float cloudDis = pos.y + 0.7; // + noiseDis;
    float cloudDis = sdCube(pos + vec3(0, 1.7, 0), vec3(100.0, 1.0, 100.0));
    float mountainDis = sdMountain(pos); 

    float finalDis = (includeFog) ? opUnion(mountainDis, cloudDis) : mountainDis; 
    
    // Final material to return
    int mat = MATMOUNTAIN; 
    if (finalDis == cloudDis) {
        mat = MATVOLUMETRIC;
    }
    
    return vec2(finalDis, mat);
}

// Return the normal ray
vec3 getNormal(vec3 pos) {
    float dis = getDist(pos, false).x;
    vec2 val = vec2(0.01, 0.0);
    
    // To get the slope we give the curve two values super close together
    // Instead of deriving we can do this method via swizzling
    vec3 normal = dis - vec3(getDist(pos-val.xyy, false).x, 
                             getDist(pos-val.yxy, false).x, 
                             getDist(pos-val.yyx, false).x);
    // Same way of doing getDist(p-vec3(0.01, 0.0, 0.0), ... etc
    return normalize(normal);
}

// Ray Marching function
vec2 rayMarch(vec3 rayOrigin, vec3 rayDir, float side, bool includeFog) {
    float disFromOrigin = 0.0;
    int mat;
    
    // Now we loop for max steps or until the distance is very small
    for (int i=0;i<MAXSTEPS;i++) {
        // Declare current marching location
        vec3 pos = rayOrigin + rayDir * disFromOrigin;
        vec2 passedDMat = getDist(pos, includeFog) * side; // exterior = +1, interior = -1;
        float disToScene = passedDMat.x;
        mat = int(passedDMat.y);
        
        disFromOrigin += disToScene;
        
        // Exit condition; we have a hit
        if (disFromOrigin > MAXDIS || abs(disToScene) < SURFDIS) { break; }
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
    float dis = rayMarch(pos + normal * SURFDIS * 2.0, lightDir, 1.0, false).x;
    
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
        aoSum += inc * getDist(currPos, false).x;
        aoMaxSum += inc * offset;
    }
    
    // Make sure it's between 0->1
    return aoSum / aoMaxSum;
}

float getLightIntensity(float dif, vec3 lightPos, vec3 pos, float maxLightLen) {
    return dif * clamp(1.0 - (length(lightPos - pos) / maxLightLen), 0.0, 1.0);
}

float beerLambertLaw(float absorption, float dis) {
    return exp(-dis * absorption);
}

vec3 background(vec3 rayDir) {
    vec3 col = vec3(0.0);
    
    // Fake scattering gradient
    float y = rayDir.y * 0.3 + 0.7; // light is top, dark is bottom. 
    col += y * skyCol;
    
    // Fake sun
    vec3 sunDir = normalize(vec3(1.5, 0.2, 0.0));
    col += pow(max(dot(rayDir, normalize(sunDir)), 0.0), 8.0) * sunCol * 0.15;
    col += pow(max(dot(rayDir, normalize(sunDir)), 0.0), 1024.0) * sunCol;    
    return col;
}

// Main //////////////////////////////////////////////////////////////////////

vec3 render(vec3 inCol, vec3 rayOrigin, vec3 rayDir, bool includeFog) {
    // Assign col
    vec3 col = inCol;
  
    // Lighting Setup    
    vec3 ambLightPos = vec3(0, 10.0, 1.0);
    vec3 keyLightPos = vec3(3, 2.5, 5.0);
    
    float lightMaxLen = 10.0;
    vec3 ambientLight = vec3(0.05);
    vec3 mixCol = vec3(ambientLight); 
    
    // Visualise 
    vec2 passedDMat = rayMarch(rayOrigin, rayDir, 1.0, includeFog);
    float dis = passedDMat.x;
    int mat = int(passedDMat.y);
    
    // Base background
    //col = background(rayDir);
    
    // If hit? render
    if (dis < MAXDIS) {
        vec3 pos = rayOrigin + rayDir * dis;
        vec3 normal = getNormal(pos);
        
        // Lock the lights based on their distance to scene
        float ambLight = getPointLight(pos, ambLightPos);
        float diffuseLight = getPointLight(pos, keyLightPos); 
        
        // Materials
        if (mat == MATNULL) {
            // Checkerboard pattern
            const float size = 1.0;
            float total = mod(floor(pos.z * size), 2.0) + 
                          mod(floor(pos.x * size), 2.0);
            bool isEven = mod(total, 2.0) == 0.0;
            mixCol *= (isEven) ? shadowCol : lightCol;
        } else if (mat == MATVOLUMETRIC) {            
            // Raymarch into the medium
            float volumeDepth = 0.0f;
            
            vec3 fogCol = vec3(0.88,0.96,0.98);
            float stepSize = 0.01f;
            vec3 posEnter = pos - normal * SURFDIS * 3.0;

            for(int i = 0; i < MAXSTEPS; i++) {
                vec3 currPos = posEnter + volumeDepth * rayDir;
                // 3 octaves, stepped to raise colour ceiling
                float currDensity = smoothstep(0.4, 1.5, fbm3DNoise(currPos, 3));
                
                float disToScene = getDist(currPos, includeFog).x * -1.0;
                volumeDepth += stepSize * currDensity; // if put before currpos changes edges
                
                if (volumeDepth > MAXDIS || abs(disToScene) < SURFDIS) { break; }
            }
            
            float absorption = 0.2;
            float volumeOpacity = 1.0-beerLambertLaw(absorption, volumeDepth);

            col = mix(col, fogCol, volumeOpacity);

        } else if (mat == MATMOUNTAIN) {
            col = dirtCol;
            
            // Colour based on steepness
            vec3 localUpDir = vec3(0, 1, 0);
            float steepness = dot(normal, localUpDir);
            // Grass and Snow
            vec3 grassMix = mix(col, grassCol, smoothstep(-1.4, -2.0, pos.y));
            col = mix(col, grassMix, smoothstep(0.5, 1.0, steepness));
            vec3 snowMix = mix(col, snowCol, smoothstep(-1.5, -0.8, pos.y));
            col = mix(col, snowMix, smoothstep(0.5, 1.0, steepness));
        }
        
        if (mat != MATVOLUMETRIC) {
            // Ambient Occlusion (check out Alex Evans)
            // Trick is to sample sdf along normal
            float ao = getAmbientOcclusion(pos, normal);
            mixCol *= ao;

            // Lighting
            mixCol += blueFogCol * ambLight * 0.1;
            mixCol += sunRayCol * diffuseLight * 0.5;

            col += vec3(mixCol);
        }
    }
    
    // Fog
    float fogStart = 0.0;
    float fogEnd = MAXDIS * 1.1;
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
    float camZoom = 0.8;
    vec3 rayOrigin = vec3(camX, camY, camZ); // this is the camera (origin of vector)
    rayOrigin.yz *= rotate(-m.y * PI + 1.);
    rayOrigin.xz *= rotate(-m.x * 2.0 * PI);
    vec3 rayDir = getRayDir(uv, rayOrigin, vec3(camYTilt, 0, camXTilt), camZoom);
    
    col = render(vec3(0), rayOrigin, rayDir, false);
    col = render(col, rayOrigin, rayDir, true);
    
    // Output to screen
    fragColor = vec4(col,1.0);
}
