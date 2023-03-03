/*
Van-Andrew Nguyen
03/03/23
[Sahara Desert]

Wanted to play with gradient noise functions to create ripples in the sand. Was originally influenced by 
Journey's sand. 
*/


// Main
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2 uv = fragCoord.xy / iResolution.xy;
    vec3 col; //texture(iChannel0, uv).rgb;
    
    // Chromatic Abberation
    vec2 offset = (uv - 0.5) * 0.004;
    col.x = texture(iChannel0, uv - offset).x;
    col.y = texture(iChannel0, uv).y;
    col.z = texture(iChannel0, uv + offset).z;
    
    // Vignette
    float vig = uv.x * uv.y * (1.0 - uv.x) * (1.0 - uv.y);
    float vigIntensity = 0.2;
    vig = clamp(pow(16.0 * vig, vigIntensity), 0.0, 1.0);
    col *= 0.5 * vig + 0.5;
    
    fragColor = vec4(col, 1.0);
}


// Buffer A
#define MAXSTEPS 100
#define MAXDIS 4.0
#define SURFDIS 0.01
#define PI 3.14159

// Define global colours
const vec3 BLUE_FOG_COL = vec3(0.729, 0.909, 0.992);
const vec3 SUN_RAY_COL = vec3(0.980, 0.633, 0.480);
const vec3 SUN_COL = vec3(1.0, 0.6, 0.3) * 3.0;
const vec3 SAND_COL = vec3(0.98,0.95,0.81);
const vec3 BASE_CLOUD_COL = vec3(0.937, 0.376, 0.101); 
const vec3 HIGHLIGHT_CLOUD_COL = vec3(0.968, 0.894, 0.929);

// Define global positions
const vec3 sunPos = vec3(6.0, 0.2, 10.0);

// Hash noise functions
float hash11(float n) {
    return fract( n*17.0*fract( n*0.3183099 ) );
}

float hash21(vec2 uv) {
    return fract(hash11(uv.x * uv.y));
}

vec2 hash22(vec2 p) {
    float n = sin(dot(p, vec2(113, 1)));
    p = fract(vec2(2097152, 262144)*n)*2. - 1.;
    return p;
}

float hash31(vec3 p) {
	p = vec3(dot(p,vec3(127.1,311.7, 74.7)),
			 dot(p,vec3(269.5,183.3,246.1)),
			 dot(p,vec3(113.5,271.9,124.6)));

	return fract(sin(p.x * p.y * p.z)*43758.5453123);
}

float smoothNoise(vec2 uv) {
    // Create 1D random value
    vec2 index = uv;
    vec2 localUV = fract(index); 
    vec2 cellID = floor(index); 
    
    localUV = localUV*localUV*(3.0 - 2.0 * localUV); 
    
    // Get noise values for corners of each cell
    float bl = hash21(cellID);
    float br = hash21(cellID + vec2(1, 0));
    float b = mix(bl, br, localUV.x);
    float tl = hash21(cellID + vec2(0, 1));
    float tr = hash21(cellID + vec2(1, 1));
    float t = mix(tl, tr, localUV.x);
    float noiseCol = mix(b, t, localUV.y);
        
    return noiseCol;
}

// Smooth out the line using cubic
float smoothFade(float t) {
    return t*t*t*(t*(t*6.0 - 15.0) + 10.0);
}

// Use texture to generate gradient vector for noise
vec2 getGradientTexture(vec2 pos) {
	float texW = 256.0;
	vec4 v = texture(iChannel0, vec2(pos.x / texW, pos.y / texW));
    return normalize(v.xy * 2.0 - vec2(1.0)); 
}

// Noise generation from texture gradient
float getNoise(vec2 pos) {
    // n(p) = (1 - F(p-p0))g(p0)(p-p0) + F(p-p0)g(p1)(p-p1)
    // Source: https://gpfault.net/posts/perlin-noise.txt.html
    // Get four corners as seperate points
    vec2 pos0 = floor(pos);
    vec2 pos1 = pos0 + vec2(1.0, 0.0);
    vec2 pos2 = pos0 + vec2(0.0, 1.0);
    vec2 pos3 = pos0 + vec2(1.0, 1.0);

    // Get gradients for the four corners
    vec2 gradient0 = getGradientTexture(pos0);
    vec2 gradient1 = getGradientTexture(pos1);
    vec2 gradient2 = getGradientTexture(pos2);
    vec2 gradient3 = getGradientTexture(pos3);
    
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

    // Calculate final result
    return (1.0 - vBlend) * p0p1 + vBlend * p2p3;
}

// Frankensteined from https://www.shadertoy.com/view/ld3BzM
// Gradient function without texture
float getGradient21(in vec2 f){
   const vec2 e = vec2(0, 1);
    vec2 p = floor(f);
    f -= p; 

    // Hermite curve (cubic smoothing)
    vec2 w = f*f*(3.0 - 2.0 * f); 

    // Four corners, mixed (see above)
    float c = mix(mix(dot(hash22(p + e.xx), f - e.xx), dot(hash22(p + e.yx), f - e.yx), w.x),
                  mix(dot(hash22(p + e.xy), f - e.xy), dot(hash22(p + e.yy), f - e.yy), w.x), w.y);
    
    // Change bounds to 0->1 from -1->1
    c = c * 0.5 + 0.5;
    return c;
}

float getGradient11(float x, float offs){
    x = abs(fract(x / 6.283 + offs - 0.25) - 0.5) * 2.0;
    
    // Triangle wave from changed bounds (square then multiply by higher bound)
    float x2 = clamp(x * x * (-1.0 + 2.0 * x), 0.0, 1.0); 
    x = smoothstep(0.0, 1.0, x);
    return mix(x, x2, 0.15);
}

// SDFS //////////////////////////////////////////////////////////////////////

// Shape Operations (From https://www.iquilezles.org/)
float opUnion(float d1, float d2) { return min(d1,d2); }
float opSubtraction(float d1, float d2) { return max(-d1,d2); }
float opIntersection(float d1, float d2) { return max(d1,d2); }

// Rotation Matrix
mat2 rotate(float a) {
    float s = sin(a);
    float c = cos(a);
    return mat2(c, -s, s, c);
}

float sdSandSurf(vec2 pos) {
    // Rotate position layer
    vec2 newPos = rotate(PI / 18.0) * pos;
    // Add noise to new layer to create wavy effect
    newPos.y += (getGradient21(newPos * 18.0) - 0.5) * 0.05; 
    // Apply gradient lines to repeat them
    float grad1 = getGradient11(newPos.y * 80.0, 0.0); 
    
    // Do the same as above for second gradient
    newPos = rotate(-PI / 20.0) * pos; 
    newPos.y += (getGradient21(newPos * 12.0) - 0.5) * 0.05;
    float grad2 = getGradient11(newPos.y * 80.0, 0.5); 

    // Rotate again
    newPos = rotate(PI / 4.0) * pos;

    // Repeating the period by switching xy and yx
    float a2 = dot(sin(newPos * 12.0 - cos(newPos.yx * 12.0)), vec2(0.25)) + 0.5;
    float a1 = 1.0 - a2;

    // Blend gradients
    float c = 1.0 - (1.0 - grad1 * a1) * (1.0 - grad2 * a2);
   
    return c;
}

float sdSandBase(vec2 pos, float base) { 
    // Shape the general landscape
    // Some domain warping here
    float layer1 = getNoise(getNoise(pos) + pos) * 0.1; 
    pos = rotate(PI / 8.0) * pos;
    float layer2 = abs(getNoise(pos * 0.2));
    float dune = layer1 + layer2 * 0.7;
    
    // Ripples (change position from above to avoid repetitive landscapes)
    vec2 newPos = vec2(pos.y - pos.x, pos.x + pos.y) * 0.5; 
    float c1 = sdSandSurf(newPos);
    vec2 rotatedPos = rotate(PI / 12.0) * newPos;
    float c2 = sdSandSurf(rotatedPos * 1.25);
    c1 = mix(c1, c2, smoothstep(0.1, 0.9, getGradient21(p * vec2(4))));
    
    // Noise (using texture to add grainy noise heightmap on top)
    c1 += texture(iChannel0, pos / 8.0).r * 0.1;
    c1 += texture(iChannel0, pos * 2.0 + vec2(-0.1, -0.2) * iTime * 0.5 * dune).r * 2.0;
    
    float final = base + dune + c1 * 0.004;
    return final;
}

// Vectors //////////////////////////////////////////////////////////////////////

// Move camera forwards (don't move camera position because we'll get floating point error after a while)
vec3 transformPos(vec3 pos) {
    pos.z += iTime * 0.04;
    return pos;
}

float getDist(vec3 pos) {
    pos = transformPos(pos);
    float baseDis = pos.y;
    float duneDis = sdSandBase(pos.xz * 3.0, pos.y) * 0.8;
    float finalDis = baseDis + duneDis;
    
    return finalDis;
}

vec3 getNormal(vec3 pos) {
    float dis = getDist(pos);
    vec2 val = vec2(0.01, 0.0);
    
    // Swizzle to get slope
    vec3 normal = dis - vec3(getDist(pos-val.xyy), 
                             getDist(pos-val.yxy), 
                             getDist(pos-val.yyx));
    return normalize(normal);
}

float rayMarch(vec3 rayOrigin, vec3 rayDir) {
    float disFromOrigin = 0.0;
    for (int i=0;i<MAXSTEPS;i++) {
        vec3 pos = rayOrigin + rayDir * disFromOrigin;
        float disToScene = getDist(pos);
        disFromOrigin += disToScene;
        
        // Exit condition; we have a hit
        if (disFromOrigin > MAXDIS || disToScene < SURFDIS) { break; }
    }
    
    return disFromOrigin;
}

float getLight(vec3 pos, vec3 lightOrigin) {
    // Get the light origin, light ray, normal
    vec3 lightPos = lightOrigin;
    vec3 light = normalize(lightPos - pos);
    vec3 normal = getNormal(pos);
    
    // Get diffused lighting 
    float dif = clamp(dot(normal, light), 0.0, 1.0);
    float dis = rayMarch(pos + normal * SURFDIS * 2.0, light);
    
    // Fake lighting, based on light distance
    float shadowIntensity = 0.5;
    float shadowBlur = 0.4;
    float stepper = clamp(dis / length(lightPos - pos), 0.0, 1.0);
    dif *= smoothstep(shadowIntensity - shadowBlur, shadowIntensity, stepper);
    dif *= 0.6;
    
    /*
    Old lighting equation; don't need this since harsh shadows make the dreamy effect of the dunes redundant
    if (dis < length(lightPos - pos)) {
        float shadowIntensity = 0.3;
        dif *= shadowIntensity;
    }
    */
    return dif;
}

vec3 renderClouds(vec3 rayDir, vec3 originCol) {
    vec3 col = vec3(0.0);
    
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
    col += vec3(smoothNoise(uv1 * 4.0)) * 1.0;
    col += vec3(smoothNoise(uv1 * 8.0)) * 0.5;
    col += vec3(smoothNoise(uv2 * 16.0)) * 0.25;
    col += vec3(smoothNoise(uv3 * 32.0)) * 0.125;
    col += vec3(smoothNoise(uv4 * 64.0)) * 0.0625;
    col /= 1.8;
    
    // Get rid of the darker areas by forcing values into threshold
    col *= smoothstep(-0.2, 1.4, col);
    vec3 output = mix(mix(BASE_CLOUD_COL, HIGHLIGHT_CLOUD_COL, col * 1.35), originCol, 1.0 - col);
    return output;
}

vec3 background(vec3 rayDir, bool includeClouds) {
    vec3 col = vec3(0.0);
    
    // Get gradient: light is top, dark is bottom. 
    float y = rayDir.y * 0.5 + 0.5; 
    col += y * BLUE_FOG_COL;
    
    // Used for different passes (before or after rendering landscape)
    if (includeClouds) {
        col = renderClouds(4.0 * (rayDir * 0.5 + 0.5), col);
    }
    
    // Fake sun
    vec3 sunDir = normalize(sunPos);
    col += pow(max(dot(rayDir, normalize(sunDir)), 0.0), 2.0) * SUN_RAY_COL * 0.3;
    col += pow(max(dot(rayDir, normalize(sunDir)), 0.0), 8.0) * SUN_COL * 0.1;
    col += pow(max(dot(rayDir, normalize(sunDir)), 0.0), 1024.0) * SUN_COL * 0.9;
    
    return col;
}

// Main //////////////////////////////////////////////////////////////////////

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
    vec2 uv = (fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;
    vec3 col = vec3(0.0);
    
    // Setup Camera
    float camHeight = 0.15;
    float downTilt = -0.25;
    vec3 rayOrigin = vec3(0, camHeight, 0);
    vec3 rayDir = normalize(vec3(uv.x, uv.y + downTilt, 1));
    vec3 ambLightPos = vec3(0, 10.0, 1.0);
    vec3 keyLightPos = vec3(2, 0.15, 1.5);
    
    // Fog
    float fogStart = 0.0;
    float fogEnd = MAXDIS;
    col = background(rayDir, true);
    
    // March 
    float dis = rayMarch(rayOrigin, rayDir);
    if (dis < MAXDIS) {
        vec3 pos = rayOrigin + rayDir * dis;
        vec3 normal = getNormal(pos);
        float ambLight = getLight(pos, ambLightPos);
        float diffuseLight = getLight(pos, keyLightPos); 
        
        vec3 mixCol = SAND_COL;
        
        // Specular Lighting
        float fresnel = dot(normal, (sunPos - pos));
        float specStrength = 0.05;
        vec3 lightReflectRayDir = reflect(sunPos - pos, normal);
        float spec = pow(max(dot(rayDir, rayDir), 0.0), 16.0);
        vec3 specularCol = specStrength * spec * SUN_COL * fresnel;
        mixCol += specularCol * (1.0 - (dis / MAXDIS));
        
        // Lighting
        mixCol *= BLUE_FOG_COL * ambLight;
        mixCol += SUN_RAY_COL * diffuseLight * 1.0;
        
        // Output
        mixCol = clamp(mixCol, 0.0, 1.0);
        col = vec3(mixCol);
        col = mix(col, background(rayDir, false), smoothstep(fogStart, fogEnd, dis));
    } else {
        // Render clouds behind landscape
        col = mix(col, background(rayDir, true), smoothstep(fogStart, fogEnd, dis));
    }
    
    fragColor = vec4(col,1.0);
}
