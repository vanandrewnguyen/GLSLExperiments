/*
Van Andrew Nguyen
26/10/21
[Moana Water Shader]

Went back and redid my first attempt at ray marching water. This time I followed the blog;
https://wallisc.github.io/rendering/2020/12/08/Making-Of-Moana-the-shadertoy.html
more closely and ended up with a much closer result, though I am still missing caustics and shadows. I've yet to learn how to create
volumetric clouds for the background as well.
Still, am very happy with the result.
*/

// Globals //////////////////////////////////////////////////////////////////////

/*
Max steps is the maximum amount of steps we can make to get a ray length
Max distance is the maximum distance we can shoot a ray
Surf distance is the distance we loop to get to register a hit on the surface
*/
#define MAXSTEPS 100
#define MAXDIS 8.0
#define SURFDIS 0.01

#define MATWATER 1
#define MATSAND  2
#define MATCORAL 3

const vec3 shadowCol = vec3(0.2, 0.1, 0.15);
const vec3 lightCol = vec3(0.8, 0.75, 0.7);
const vec3 fillLightCol = vec3(0.843, 0.976, 0.968);
const vec3 keyLightCol = vec3(0.937, 0.376, 0.101);
const vec3 waterCol = vec3(0.313, 0.949, 0.847);
const vec3 sandCol = vec3(0.905, 0.854, 0.694);
const vec3 coralCol = vec3(0.243, 0.156, 0.113);
const vec3 skyCol = vec3(0.729, 0.917, 0.933);

// Noise //////////////////////////////////////////////////////////////////////

// Return a randomised float using a vec2 (uv coord)
float hash1(float n) {
    return fract( n*17.0*fract( n*0.3183099 ) );
}

float hash21(vec2 uv) {
    return fract(hash1(uv.x * uv.y));
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

// Fractal Brownian Motion 
float fbm(vec2 pos, int iterations) {
    // Declare return value, amplitude of motion, freq of motion
    float val = 0.0;
    float amp = 0.05;
    float freq = 4.0;
    
    // Now loop through layers and return the combined value
    for (int i=0;i<iterations;i++) {
        val += amp * smoothNoise(freq * pos); 
        amp *= 0.5;
        freq * 2.0;
    }
    
    return val;
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

// Water SDF
float sdWater(vec3 pos) {
    float dis = pos.y + 0.1;
    float amp = 0.5;
    float freq = 1.0;
    float iterations = 4.0;
    float t = iTime;
    
    // Use FBM for main wave shapes
    for (float i = 0.0; i < iterations; i += 1.0) {
        dis += amp * sin(pos.x * freq + i * amp + t * i) * 
                     sin(pos.z * freq + -i * amp + t * i);
        
        amp *= 0.25;
        freq *= 2.0;
    }
    
    // Use noise fbm to add irregular waves
    dis += fbm(pos.xz * 1.25, 2);
    
    // Use texture for finer wave details
    float texAmp = 0.2 * sin(pos.x * pos.z);
    float texFreq = 4.0 + 0.2 * pos.z;
    float texOffset = iTime * 0.05;
    dis += texture(iChannel0, pos.xz / texFreq + texOffset).r * texAmp;
    
    return dis * 0.8;
}

// Coral SDF
float sdCoral(vec3 pos) {
    // Return a bunch of spheres
    float s1 = sdSphere(pos, vec3(-1.2, -1.2, 3.2), 0.7);
    float s2 = sdSphere(pos, vec3(-1.5, -1.2, 2.1), 0.35);
    float s3 = sdSphere(pos, vec3(1.6, -1.2, 2.3), 0.35);
    float s4 = sdSphere(pos, vec3(-1.1, -1.2, 2.6), 0.5);
    float s5 = sdSphere(pos, vec3(1.3, -3.0, 3.2), 2.0);
    
    float dis = opUnion(s1, s2);
    dis = opUnion(dis, s3);
    dis = opUnion(dis, s4);
    dis = opUnion(dis, s5);
    
    // Use texture to add bumps
    dis += texture(iChannel0, pos.xy / 2.5).r * 0.10;
    dis += texture(iChannel0, pos.xy / 1.0).r * 0.05;
    dis += texture(iChannel1, pos.xz * 0.5).r * 0.005;
    return dis * 0.9;
}

// Sand SDF
float sdSand(vec3 pos) {
    float dis = pos.y + 1.2;
    dis += texture(iChannel1, pos.xz * 0.1).r * 0.04;
    dis += texture(iChannel1, 10.0 + pos.xz * 2.0).r * 0.005;
    
    return dis;
}

// Vectors //////////////////////////////////////////////////////////////////////

// Get distance function
vec2 getDist(vec3 pos, bool includeWater) {
    
    // Get dist of the various shapes
    float viewerSphereDis = sdSphere(pos, vec3(0, 0.5, 1.5), 2.0 + 0.05 * sin(iTime)); 
    float waterSurfaceDis = sdWater(pos); 
    float oceanFloorDis = sdSand(pos); 
    float coralDis = sdCoral(pos);
    
    float finalDis;
    if (includeWater) {
        finalDis = opSmoothSubtraction(viewerSphereDis, waterSurfaceDis, 0.8);
        finalDis = opUnion(finalDis, oceanFloorDis);
        finalDis = opUnion(finalDis, coralDis);
    } else {
        finalDis = opSmoothUnion(oceanFloorDis, coralDis, 0.02);
    }
    
    int mat = 0;
    if (finalDis == oceanFloorDis) {
        mat = MATSAND;
    } else if (finalDis == coralDis) {
        mat = MATCORAL;
    } else {
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
    int mat = 0;
    
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
    //lightPos.xz += vec2(sin(iTime), cos(iTime)) * 4.0;
    
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
        float shadowIntensity = 0.5;
        dif *= shadowIntensity;
    }
    
    return dif;
}

// Main //////////////////////////////////////////////////////////////////////

// Background function; depth perception
vec3 background(vec3 rayDir) {
    vec3 col = vec3(0.0);
    
    float y = rayDir.y * 0.2 + 0.8; // light is top, dark is bottom. 
    col += y * skyCol;
    float x = rayDir.x * 0.6 + 0.4;
    col += x * keyLightCol * 0.4;
    
    return col;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab UV Coord
    vec2 uv = (fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;

    // Declare Col
    vec3 col = vec3(0.0);

    // Setup Camera
    float camHeight = -0.6;
    float downTilt = 0.05;
    vec3 rayOrigin = vec3(0, camHeight, 0); // this is the camera (origin of vector)
    vec3 rayDir = normalize(vec3(uv.x, uv.y + downTilt, 1));
    
    // Lights
    vec3 keyDiffuseLightPos = vec3(3, 2, 6);
    vec3 fillDiffuseLightPos = vec3(-1, 4, 1);
    float maxLightDis = 10.0;
    vec3 ambientLight = vec3(0.2);
    
    // Materials
    float IOR = 1.33;
    vec3 mixCol = ambientLight; 
    
    // Visualise 
    vec2 passedDMat = rayMarch(rayOrigin, rayDir, true);
    float dis = passedDMat.x;
    int mat = int(passedDMat.y);
    
    if (dis < MAXDIS) {
        vec3 pos = rayOrigin + rayDir * dis;
        vec3 normal = getNormal(pos, true);
        vec3 posEnter = pos - normal * SURFDIS * 3.0;
        float keyDiffuseLight = getLight(pos, keyDiffuseLightPos, true);
        float fillDiffuseLight = getLight(pos, fillDiffuseLightPos, true);
        float finalFill = fillDiffuseLight * clamp(1.0 - (length(fillDiffuseLightPos - pos) / maxLightDis), 0.0, 1.0);
        float finalKey = keyDiffuseLight * clamp(1.0 - (length(keyDiffuseLightPos - pos) / maxLightDis), 0.0, 1.0);
        
        // Materials
        if (mat == MATWATER) {            
            // Depth rendering
            vec3 rayDirIn = refract(rayDir, normal, 1.0 / IOR);
            vec2 interiorPassedDMat = rayMarch(posEnter, rayDirIn, false);
            float disIn = interiorPassedDMat.x;
            int matIn = int(interiorPassedDMat.y);
            vec3 posExit = posEnter + rayDirIn * disIn;
            
            // Mix the colour we are hitting
            vec3 interiorCol = sandCol;
            if (matIn == MATSAND) { 
                interiorCol = sandCol; 
            } else if (matIn == MATCORAL) {
                interiorCol = coralCol;
            }
            float stepper = abs(length(posEnter - posExit));
            mixCol *= mix(interiorCol, waterCol, stepper);
            if (stepper > MAXDIS * 0.5) { 
                // We haven't hit the sea bed, hence display the sky
                mixCol = mix(mixCol, skyCol, stepper / (MAXDIS * 10.0));
            }
            
            // Reflections
            float reflectedVal = 0.1;
            vec3 reflectedRay = reflect(rayDir, normal);
            vec2 exteriorPassedDMat = rayMarch(posEnter, rayDirIn, false);
            float disOut = exteriorPassedDMat.x;
            int matOut = int(exteriorPassedDMat.y);
            if (disOut < MAXDIS) {
                // Add colours that are reflected
                vec3 exteriorCol = sandCol;
                if (matOut == MATSAND) { 
                    exteriorCol = sandCol; 
                } else if (matOut == MATWATER) {
                    exteriorCol = waterCol;
                }
                mixCol += exteriorCol * reflectedVal;
            } else {
                mixCol += keyLightCol * reflectedVal * 2.0;
            }
            
        } else if (mat == MATSAND) {
            mixCol *= sandCol;
        } else if (mat == MATCORAL) {
            mixCol *= coralCol;
        }
        
        // Lighting
        mixCol += finalKey * keyLightCol;
        mixCol += finalFill * fillLightCol;
        
        // Specular
        float specStrength = 0.000015;
        float fresnel = dot(normal, (keyDiffuseLightPos - posEnter));
        vec3 lightReflectRayDir = reflect(keyDiffuseLightPos - posEnter, normal);
        float spec = pow(max(dot(rayDir, lightReflectRayDir), 0.0), 4.0);
        vec3 specularCol = specStrength * spec * keyLightCol * fresnel;
        mixCol += specularCol;
    }
    
    col = vec3(mixCol);
    
    // Fog / Sky
    float fogStart = 0.0;
    float fogEnd = MAXDIS;
    col = mix(col, background(rayDir), smoothstep(fogStart, fogEnd, dis));
    
    // Output to screen
    fragColor = vec4(col,1.0);
}
