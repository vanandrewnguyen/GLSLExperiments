/*
Van Andrew Nguyen
19/10/21
[Mountain Ranges v4]

This is my fourth attempt at ray marching mountains, and it was a joy to combine my previous works to make this new shader.
The resultant effect are mountain ridges with pools of water which reflect and refract light. 
I did not use any fancy fbm techniques to make the mountains, rather I relied on textures (layers and layers and layers of them)
to shape the mountain side. I divided the mountain and water into two materials which are coloured seperately.
The mountains are coloured based on steepness rather than colour-banding. I ray march twice for the water; once to see how deep the water is
and mix it with reference to the mountain terrain below it, once to get the reflections of the normal.
The lighting is identical to how the Himalayas shader is setup, one fill and one key light.
*/

// Globals //////////////////////////////////////////////////////////////////////

/*
Max steps is the maximum amount of steps we can make to get a ray length
Max distance is the maximum distance we can shoot a ray
Surf distance is the distance we loop to get to register a hit on the surface
*/
#define MAXSTEPS 100
#define MAXDIS 16.0
#define SURFDIS 0.01

#define MATMOUNTAIN 1
#define MATWATER 2

const vec3 shadowCol = vec3(0.2, 0.1, 0.15);
const vec3 lightCol = vec3(0.8, 0.75, 0.7);
const vec3 dirtCol = vec3(0.117, 0.117, 0.101);
const vec3 grassCol = vec3(0.223, 0.270, 0.207);
const vec3 snowCol = vec3(0.882, 0.921, 0.937);
const vec3 waterCol = vec3(0.094, 0.274, 0.329);
const vec3 blueLightCol = vec3(0.729, 0.909, 0.992);
const vec3 orangeLightCol = vec3(0.968, 0.596, 0.431);

// Return a randomised float using a vec2 
float hash21(vec2 uv) {
    // Pseudo-random math function to get a random number
    vec2 o = fract(uv * vec2(5859.56, 90.123));
    o += dot(o, o + 554.89);
    return fract(o.x * o.y);	
}

// Return a randomised float using a vec3 
float hash31(vec3 p) {
    // Again, random math functions, e.g. dot product of vector with random numbers to 
    // generate 'random' result
	p = vec3(dot(p,vec3(127.1,311.7, 74.7)),
			 dot(p,vec3(269.5,183.3,246.1)),
			 dot(p,vec3(113.5,271.9,124.6)));

	return fract(sin(p.x * p.y * p.z)*43758.5453123);
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

// Displacement
float sdDisplace(vec3 pos, float iterations, float dis) {
    float height = 0.0;
    
    // Loop through and add height based on a texture - play with freq
    for (float i=0.0;i<iterations;i+=1.0) {
        vec3 newPos = pos += i * 10.0;
        float texAmp = 1.2 - 0.2 * i;
        float texFreq = 16.0 * i;
        height += texture(iChannel0, newPos.xz / texFreq).r * texAmp;
        height += texture(iChannel0, newPos.xz / texFreq).r * texAmp;
        height += texture(iChannel0, newPos.xz / texFreq).r * texAmp;
    }
        
    return height + dis;
}

float sdMountains(vec3 pos) {
    // Generate plane and texture height map
    float iterations = 4.0;
    float mtnDis = pos.y - iterations;
    mtnDis = sdDisplace(pos, iterations, mtnDis);
    mtnDis += texture(iChannel1, pos.xz / 10.0).r * 0.1;
    mtnDis *= 0.25;
    
    return mtnDis;
}

float sdWater(vec3 pos) {
    // Generate plane and waves
    float waterDis = pos.y + 2.1;
    waterDis += 0.005 * sin(10.0 * pos.x + iTime);
    waterDis += 0.005 * sin(15.0 * pos.z + iTime);
    
    return waterDis;
}

// Vectors //////////////////////////////////////////////////////////////////////

// Transform terrain position without affecting ray origin
vec3 transformPos(vec3 pos) {
    pos.z += iTime * 2.0;
    return pos;
}

// Get distance function
vec2 getDist(vec3 pos, bool includeWater) {
    pos = transformPos(pos);
    
    // Get dist
    float mtnDis = sdMountains(pos);
    float waterDis = (!includeWater) ? MAXDIS : sdWater(pos);
    
    // Final distance to return
    float finalDis = opUnion(waterDis, mtnDis);
    
    // Final material to return
    int mat = 0;
    if (finalDis == mtnDis * 0.25) {
        mat = MATMOUNTAIN;
    } else if (finalDis == waterDis) {
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
    vec2 passedDMat = vec2(0.0);
    
    // Now we loop for max steps or until the distance is very small
    for (int i=0;i<MAXSTEPS;i++) {
        // Declare current marching location
        vec3 pos = rayOrigin + rayDir * disFromOrigin;
        passedDMat = getDist(pos, includeWater);
        float disToScene = passedDMat.x;
        
        disFromOrigin += disToScene;
        
        // Exit condition; we have a hit
        if (disFromOrigin > MAXDIS || disToScene < SURFDIS) { break; }
    }
    
    return vec2(disFromOrigin, passedDMat.y);
}

// Lighting function
float getLight(vec3 pos, vec3 lightOrigin, bool includeWater) {
    // Get the light origin
    vec3 lightPos = lightOrigin; // basically a point in 3D space
    
    // Get the light ray
    vec3 light = normalize(lightPos - pos);
    
    // Get normal
    vec3 normal = getNormal(pos, includeWater);
    
    // Get diffused lighting 
    // We want 0 if the rays are parallel, 1 if the rays are perpendicular
    // Hence we use dot product
    float dif = clamp(dot(normal, light), 0.0, 1.0);
    vec2 passedDMat = rayMarch(pos + normal * SURFDIS * 2.0, light, includeWater);
    float dis = passedDMat.x;
    int mat = int(passedDMat.y);
    
    // Getting shadows (if disToSphere < disToLight then we have a shadow)
    if (dis < length(lightPos - pos)) {
        float shadowIntensity = 0.2;
        dif *= shadowIntensity;
    }
    
    return dif;
}

// Main //////////////////////////////////////////////////////////////////////

// Background function; depth perception
vec3 background(vec3 rayDir) {
    vec3 col = vec3(0.0);
    
    float y = rayDir.y * 0.5 + 0.5; // light is top, dark is bottom. 
    col += y * blueLightCol;
    float x = rayDir.x * 0.5 + 0.5;
    col += x * orangeLightCol * 0.6;
    
    return col;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab UV Coord
    vec2 uv = (fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;

    // Declare Col
    vec3 col = vec3(0.0);

    // Setup Camera
    float camHeight = 0.0;
    float downTilt = -0.25;
    vec3 rayOrigin = vec3(0, camHeight, 0); // this is the camera (origin of vector)
    vec3 rayDir = normalize(vec3(uv.x, uv.y + downTilt, 1));
    vec3 fillLightPos = vec3(-1, 6, 1);
    vec3 keyLightPos = vec3(4, 4, 4);
    
    float maxLightDis = 8.0; // used to clamp light sources as ratio
    float IOR = 1.33; // water index of refraction
    
    // Visualise 
    vec2 passedDMat = rayMarch(rayOrigin, rayDir, true);
    float dis = passedDMat.x;
    int mat = int(passedDMat.y);
    
    if (dis < MAXDIS) {
        vec3 pos = rayOrigin + rayDir * dis;
        vec3 normal = getNormal(pos, true);
        float fillDiffuseLight = getLight(pos, fillLightPos, true);
        float keyDiffuseLight = getLight(pos, keyLightPos, true);
        float finalFill = fillDiffuseLight * clamp(1.0 - (length(fillLightPos - pos) / maxLightDis), 0.0, 1.0);
        float finalKey = keyDiffuseLight * clamp(1.0 - (length(keyLightPos - pos) / maxLightDis), 0.0, 1.0);
        
        // Gradual lighting curve
        vec3 mixCol = mix(shadowCol, lightCol, fillDiffuseLight); //shadowCol; //
        
        // Colour based on steepness
        vec3 localUpDir = normalize(pos);
        float steepness = -dot(normal, localUpDir);
        
        // Materials
        if (mat == MATMOUNTAIN-1) {
            mixCol = dirtCol;
            // Grass (we clamp it to the bottom half of the mountains)
            vec3 grassMix = mix(mixCol, grassCol, smoothstep(0.6, 0.4, steepness));
            mixCol = mix(mixCol, grassMix, 1.0 - clamp(pos.y + 2.4, 0.0, 1.0)); 
            // Snow
            vec3 snowMix = mix(mixCol, snowCol, smoothstep(0.7, 0.6, steepness));
            mixCol = mix(mixCol, snowMix, clamp(pos.y + 1.6, 0.0, 1.0)); 
        } else if (mat == MATWATER) {
            // Ray march for depth
            vec3 rayDirIn = refract(rayDir, normal, 1.0 / IOR);
            vec3 posEnter = pos - normal * SURFDIS * 3.0;
            float disIn = rayMarch(posEnter, rayDirIn, false).x;
            float stepper = length(disIn);
            mixCol = mix(dirtCol, waterCol, stepper);
            
            // Reflections
            vec3 reflectRayDir = reflect(rayDir, normal);
            float reflectDis = rayMarch(posEnter, reflectRayDir, false).x;
            if (reflectDis < MAXDIS) {
                mixCol = mix(mixCol, mix(dirtCol, snowCol, length(reflectDis) * 0.5), 0.02 * length(reflectDis));
            }
        }
        
        // Lighting
        mixCol *= fillDiffuseLight;
        mixCol += blueLightCol * finalFill;
        mixCol += orangeLightCol * finalKey;
        
        col = vec3(mixCol);
    }
    
    // Fog
    float fogStart = 0.0;
    float fogEnd = MAXDIS;
    col = mix(col, background(rayDir), smoothstep(fogStart, fogEnd, dis));
    
    // Output to screen
    fragColor = vec4(col,1.0);
}