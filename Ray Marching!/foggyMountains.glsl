// Van Andrew Nguyen
// [12/05/23]
// Foggy Mountains Scene
// Attempt to follow along with Inigo Quilez: https://www.youtube.com/watch?v=BFld4EBO2RE
// Learn a lot on how to use fog to block out unneeded features. Sky features a gradient with clouds and a fake sun.
// Terrain uses fbm perlin noise with textures to surface, and trees! using domain repetition.


// Main
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2 uv = fragCoord.xy / iResolution.xy;
    vec3 col = texture(iChannel0, uv).rgb;
    
    // Colour correction
    col.r = pow(col.r, 1.0);
    col.g = pow(col.g, 0.9);
    col.b = pow(col.b, 1.0);
        
    // Vignette
    float vig = uv.x * uv.y * (1.0 - uv.x) * (1.0 - uv.y);
    float vigIntensity = 0.2;
    vig = clamp(pow(16.0 * vig, vigIntensity), 0.0, 1.0);
    col *= 0.5 * vig + 0.5;
    
    fragColor = vec4(col, 1.0);
}

// Buffer A
// Max steps is the maximum amount of steps we can make to get a ray length
// Max distance is the maximum distance we can shoot a ray
// Surf distance is the distance we loop to get to register a hit on the surface
#define MAXSTEPS 200
#define MAXDIS 50.0
#define SURFDIS 0.01

#define PI 3.1416

// Define global materials
const int MATNULL = 0;
const int MATGROUND = 1;
const int MATTREE = 2;

// Define global colours
const vec3 SKY_TOP_COL = vec3(0.77,0.56,0.93);
const vec3 SKY_BOTTOM_COL = vec3(0.96,0.62,0.26);

const vec3 FOG_TOP_COL = SKY_TOP_COL; // vec3(0.78,0.94,0.98);
const vec3 FOG_BOTTOM_COL = vec3(0.4,0.37,0.41); // vec3(0.28,0.36,0.42);
const vec3 DISTANCE_FOG_COL = vec3(0.32,0.29,0.38); // vec3(0.3);

const vec3 CLOUD_TOP_COL = vec3(0.94,0.89,0.71);
const vec3 CLOUD_BOTTOM_COL = vec3(0.58,0.51,0.57);

const vec3 SUN_RAY_COL = vec3(0.980, 0.633, 0.480);
const vec3 SUN_COL = vec3(1.0, 0.6, 0.3) * 3.0;

const vec3 SAND_COL = vec3(0.98,0.95,0.81);
const vec3 DIRT_COL = vec3(0.32,0.29,0.25);
const vec3 GRASS_COL = vec3(0.36,0.46,0.11);
const vec3 TREE_COL = vec3(0.08,0.68,0.08);


// Define global positions
const vec3 SUN_POS = vec3(2.0, -0.4, -8.0);
const vec3 FILL_POS = vec3(-5.0, 0.5, 10.0);

// Pseudo-random hash functions
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

// Perlin Noise
float noise(vec2 uv) {
    // Create 1D random value
    vec2 index = uv;
    vec2 localUV = fract(index); 
    vec2 cellID = floor(index); 
    
    localUV = localUV*localUV*(3.0 - 2.0 * localUV); 
    
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

float opUnion(float d1, float d2) { 
    return min(d1, d2); 
}

mat2 rotate(float a) {
    float s = sin(a);
    float c = cos(a);
    return mat2(c, -s, s, c);
}

float sdSphere(vec3 pos, float rad) {
    return length(pos) - rad;
}

float sdTerrain(in vec3 pos) {
    float dis = pos.y;
    
    // Bigger landscape curves
    // 'pos' is needs to be offset because it is centered for rayOrigin
    pos.x += 50.0;
    pos.z -= 10.0;
    dis += noise(pos.xz * 0.25) * 3.0;
    
    // Smaller details using smaller curves
    float amp = 1.0;
    float freq = 0.5;
    int octaves = 19;
    
    for (int i = 0; i < octaves; i++) {
        if (i <= 6 || i >= 14) {
            dis += noise(pos.xz * freq) * amp;
        
            // Offset
            pos.x += 5.0;
            pos.z -= 3.3;
            freq *= 2.0;
            amp *= 0.5;
        }
    }
    
    // Cheat and use textures for surfacing
    dis += texture(iChannel1, pos.xz * 0.05).r * 0.05;
    dis += texture(iChannel0, pos.xz * 0.05).r * 0.05;
    
    return dis * 0.6;
}

float sdTrees(vec3 pos, float dis) {
    float finalDis = pos.y + 100.0;

    // Trees!
    vec2 uv = (pos.xz - vec2(12.42, 5.05)) * 15.0;
    vec2 cellID = floor(uv);
    vec2 localUV = fract(uv) - 0.5;
    
    float rad = 0.2 + (hash21(cellID) * 2.0 - 1.0) * 0.07;
    float yOffset = 0.2;
    vec3 localCoord = vec3(localUV.x, 
                           dis + yOffset, 
                           localUV.y);
    
    // Segment trees based on noise pattern
    float n = noise(pos.xz * 2.0 - vec2(-15.29, 4.59));
    if (n > 0.4) {
        float trees = sdSphere(localCoord, rad);
        trees -= texture(iChannel1, cellID + localUV * 10.0).r * 0.05;
        finalDis = trees;
    }

    return finalDis * 0.5;
}

vec2 getDist(vec3 pos) {
    // Map dist to virtual scene
    float terrainDis = sdTerrain(pos);
    float treeDis = sdTrees(pos, terrainDis);
    
    // Group final distance and material
    float finalDis = opUnion(terrainDis, treeDis);
    int mat = MATNULL;
    if (finalDis == terrainDis) {
        mat = MATGROUND;
    } else if (finalDis == treeDis) {
        mat = MATTREE;
    }
    
    return vec2(finalDis, mat);
}

vec3 getNormal(vec3 pos) {
    float dis = getDist(pos).x;
    vec2 val = vec2(0.01, 0.0);
    
    // Swizzle
    vec3 normal = dis - vec3(getDist(pos-val.xyy).x, 
                             getDist(pos-val.yxy).x, 
                             getDist(pos-val.yyx).x);
    // Same way of doing getDist(p-vec3(0.01, 0.0, 0.0), ... etc
    return normalize(normal);
}

vec2 rayMarch(vec3 rayOrigin, vec3 rayDir) {
    float disFromOrigin = 0.0;
    int mat;

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

void getHardShadow(vec3 pos, vec3 normal, vec3 lightDir, vec3 lightPos, inout float dif) {
    float dis = rayMarch(pos + normal * SURFDIS * 2.0, lightDir).x;
    
    // Getting shadows (if disToSphere < disToLight then we have a shadow)
    if (dis < length(lightPos - pos)) {
        float shadowIntensity = 0.2;
        dif *= shadowIntensity;
    }
}

float getPointLight(vec3 pos, vec3 lightOrigin) {
    vec3 lightPos = lightOrigin;
    vec3 lightDir = normalize(lightPos - pos);
    vec3 normal = getNormal(pos);
    
    // Get diff (0->1)
    float dif = clamp(dot(normal, lightDir), 0.0, 1.0);
 
    // Shadows
    getHardShadow(pos, normal, lightDir, lightPos, dif);
    dif *= getSoftShadow(pos, normal, lightDir, 8.0);
    
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
    
    // Normalize 0->1
    return aoSum / aoMaxSum;
}

float getLightIntensity(float dif, vec3 lightPos, vec3 pos, float maxLightLen) {
    return dif * clamp(1.0 - (length(lightPos - pos) / maxLightLen), 0.0, 1.0);
}

// Render 2D Clouds
vec3 renderClouds(vec3 rayDir, vec3 originCol) {
    vec3 col = vec3(0.0);
    
    rayDir.y -= 0.15;
    
    // Create varying uv coords to print from
    vec2 uv1 = rayDir.xy;
    vec2 uv2 = rayDir.xy;
    vec2 uv3 = rayDir.xy;
    vec2 uv4 = rayDir.xy;
    uv1.x += iTime / 10.0;
    uv2.x += iTime / 15.0;
    uv3.x -= iTime / 20.0;
    uv4.x += iTime / 25.0;
    
    // FBM but manual
    col += vec3(noise(uv1 * 4.0)) * 1.0;
    col += vec3(noise(uv1 * 8.0)) * 0.5;
    col += vec3(noise(uv2 * 16.0)) * 0.25;
    col += vec3(noise(uv3 * 32.0)) * 0.125;
    col += vec3(noise(uv4 * 64.0)) * 0.0625;
    col /= 1.6; 
    
    // Get rid of the darker areas
    col *= smoothstep(-1.0, 1.0, col);
    vec3 newSkyCol = mix(mix(CLOUD_BOTTOM_COL, CLOUD_TOP_COL, col * 2.0), originCol, 1.0 - col);
    
    return newSkyCol;
}

vec3 renderSun(vec3 rayDir, vec3 originCol, float intensity) {
    // Fake sun
    vec3 col = originCol;
    vec3 sunDir = normalize(SUN_POS);
    col += pow(max(dot(rayDir, normalize(sunDir)), 0.0), 2.0) * SUN_RAY_COL * 0.04 * intensity;
    col += pow(max(dot(rayDir, normalize(sunDir)), 0.0), 8.0) * SUN_COL * 0.01 * intensity;
    col += pow(dot(rayDir, normalize(sunDir)), 1024.0) * SUN_COL * 0.9 * intensity;
    return col;
}

vec3 background(vec3 rayDir) {
    float grad = 16.0 * (rayDir.y * 0.5 + 0.5) - 6.95;
    vec3 col = mix(SKY_BOTTOM_COL, SKY_TOP_COL, grad);
    
    col = renderClouds(rayDir, col);
    col = renderSun(rayDir, col, 1.0);
    return col;
}

vec3 render(vec3 rayOrigin, vec3 rayDir) {
    // Lighting Setup
    vec3 keyLightPos = SUN_POS;
    vec3 fillLightPos = FILL_POS;
    float lightMaxLen = 10.0;
    vec3 ambientLight = vec3(0.6);
    vec3 mixCol = vec3(ambientLight); 

    vec2 passedDMat = rayMarch(rayOrigin, rayDir);
    float dis = passedDMat.x;
    int mat = int(passedDMat.y);

    if (dis < MAXDIS) {
        vec3 pos = rayOrigin + rayDir * dis;
        vec3 normal = getNormal(pos);
        
        // Lock the lights based on their distance to scene
        float keyDiffuseLight = getPointLight(pos, keyLightPos);
        float fillDiffuseLight = getPointLight(pos, keyLightPos);
        
        // Ambient Occlusion (check out Alex Evans)
        // Trick is to sample sdf along normal
        float ao = getAmbientOcclusion(pos, normal);
        mixCol *= ao;
        
        // Lighting
        mixCol += keyDiffuseLight * SUN_COL * 0.4;
        mixCol += fillDiffuseLight * SKY_TOP_COL * 0.3;
        
        // Materials
        if (mat == MATGROUND) {
            float upN = dot(vec3(0,1,0), normal);
            vec3 materialCol = mix(DIRT_COL, GRASS_COL, pow(upN, 6.0));
            mixCol *= materialCol;
            
            // Triplanar mapping
            vec3 colXZ = texture(iChannel1, pos.xz*0.5+0.5).rgb;
            vec3 colYZ = texture(iChannel1, pos.yz*0.5+0.5).rgb;
            vec3 colXY = texture(iChannel1, pos.xy*0.5+0.5).rgb;
            vec3 newNorm = abs(normal);
            newNorm *= pow(newNorm, vec3(2));
            newNorm *= newNorm.x + newNorm.y + newNorm.z; // normalise
            mixCol += 0.2 * (colYZ * newNorm.x + colXZ * newNorm.y + colXY * newNorm.z);
            
        } else if (mat == MATTREE) {
            float n = noise(pos.xz * 4.0 + vec2(19.21, 43.12));
            vec3 treeCol = mix(TREE_COL, GRASS_COL, n);
            mixCol *= treeCol;
        }

        // Vertically based fog
        float groundYCoord = -2.2;
        float diff = pos.y - groundYCoord;
        float n1 = noise(pos.xz * 0.75 + vec2(-15.29, 4.59)); // foggy area
        
        // Use FBM on fog to distort
        vec2 nOffset = vec2(-25.29, 42.59);
        float nAmp = 0.8;
        float nFreq = 1.25;
        float n2 = 0.0; // texture
        for (int i = 0; i < 3; i++) {
            n2 += noise(pos.xz * nFreq + nOffset) * nAmp;
            nFreq *= 2.0;
            nAmp *= 0.3;
            nOffset += vec2(4.96, 17.32);
        }
        vec3 fogCol = mix(FOG_BOTTOM_COL, FOG_TOP_COL, max(abs(n2), 0.1) * (diff / -2.0));
        mixCol = mix(mixCol, fogCol, n1 * (1.0 - diff));
        
        // View distance Fog
        float fogStart = 0.0;
        float fogEnd = MAXDIS;
        float disStep = smoothstep(fogStart, fogEnd, dis);
        float deltaR = exp(-1.0 * disStep * 1.);
        float deltaG = exp(-1.0 * disStep * 2.);
        float deltaB = exp(-1.0 * disStep * 4.);
        mixCol.r = (vec3(deltaR) * mixCol + vec3((1.0 - deltaR)) * DISTANCE_FOG_COL).r;
        mixCol.g = (vec3(deltaG) * mixCol + vec3((1.0 - deltaG)) * DISTANCE_FOG_COL).g;
        mixCol.b = (vec3(deltaB) * mixCol + vec3((1.0 - deltaB)) * DISTANCE_FOG_COL).b;

        mixCol = renderSun(rayDir, mixCol, 0.7);
    } else {
        mixCol = background(rayDir);
    }

    return mixCol;
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

void mainImage (out vec4 fragColor, in vec2 fragCoord) {
    // Grab UV Coord
    vec2 uv = (fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;
    vec2 m = iMouse.xy/iResolution.xy;

    vec3 col = vec3(0.0);

    // Setup Camera
    float camX = 0.0;
    float camY = 0.0;
    float camZ = -6.0;
    float camXTilt = 0.0;
    float camYTilt = 0.0;
    float camZoom = 2.0;
    vec3 rayOrigin = vec3(camX, camY, camZ);
    rayOrigin.yz *= rotate(-m.y * PI + 1.);
    rayOrigin.xz *= rotate(-m.x * 2.0 * PI);
    vec3 rayDir = getRayDir(uv, rayOrigin, vec3(camYTilt, 0, camXTilt), camZoom);
    
    // normalize(vec3(uv.x + 2.0, uv.y - 0.5, 2));
    col = render(rayOrigin, rayDir);
    
    // Output to screen
    fragColor = vec4(col,1.0);
}
