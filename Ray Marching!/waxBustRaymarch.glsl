/*
Van Andrew Nguyen
24/10/2021
[Wax Bust]

Source for this shader: https://www.alanzucconi.com/2017/08/30/fast-subsurface-scattering-1/
I followed a tutorial on implementing subsurface scattering to create a wax-like material seen above. I used this to render a scene with
a wax bust, a glass plate and a wooden table, all with different materials. Shader features phong lighting model + subsurface scattering + realistic material textures.
I played with refraction and reflection with the wax and glass materials to some effect. Final result isn't too believeable yet, but is close.
*/

// Globals //////////////////////////////////////////////////////////////////////

/*
Max steps is the maximum amount of steps we can make to get a ray length
Max distance is the maximum distance we can shoot a ray
Surf distance is the distance we loop to get to register a hit on the surface
*/
#define MAXSTEPS 100
#define MAXDIS 100.0
#define SURFDIS 0.01

#define MATWAX 1
#define MATGLASS 2
#define MATWOOD 3

const vec3 shadowCol = vec3(0.2, 0.1, 0.15);
const vec3 lightCol = vec3(0.8, 0.75, 0.7);
const vec3 fillLightCol = vec3(0.843, 0.976, 0.968);
const vec3 keyLightCol = vec3(0.937, 0.376, 0.101);

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

float sdCappedCylinder(vec3 pos, float rad, float height) {
  vec2 dis = abs(vec2(length(pos.xz), pos.y)) - vec2(rad, height);
  return min(max(dis.x, dis.y),0.0) + length(max(dis, 0.0));
}

// Capsule Distance
float sdCapsule(vec3 pos, vec3 start, vec3 end, float rad) {
    // We get the projected distance of the origin->start onto start->end
    vec3 ab = end - start;
    vec3 ap = pos - start;
    // Then we normalise that distance and lock it within bounds
    float projection = dot(ab, ap) / dot(ab, ab); // / dot(ab, ab) to normalise 0->1
    projection = clamp(projection, 0.0, 1.0);
    // Then we get the distance from origin to that projected point MINUS radius of capsule
    vec3 len = start + projection * ab;
    float dis = length(pos - len) - rad;
    // Finally, we have the smallest distance to the capsule shape.
    return dis;
}

// Torus Distance
float sdTorus(vec3 pos, vec3 center, vec2 rad) {
    // We subtract the smaller radius from the length of the vector running from
    // origin point to middle of torus
    pos -= center;
    float x = length(pos.xz) - rad.x;
    return length(vec2(x, pos.y)) - rad.y;
}

// Sphere Distance
float sdSphere(vec3 pos, float rad) {
    return length(pos) - rad;
}

// Cube Distance
float sdCube(vec3 pos, vec3 center, vec3 size, float bevel) {
    pos -= center;
    
    // We get the length of the position - size (dis) and grab a version > 0
    float val = length(max(abs(pos) - size, 0.0)) - bevel;
    return val;
}

// Bust
float sdBust(vec3 pos) {
    // Base: - vec3(0, 1, 6)
    // Head Shape
    float sphereTopHeadDis = sdSphere(pos, 0.8); 
    float sphereBotHeadDis = sdSphere(pos - vec3(0, -0.4, 0), 0.55); 
    float sphereBotHeadCutoutDis = sdSphere(pos - vec3(0, -0.6, 0.5), 0.5);
    float sphereLeftCheekCutoutDis = sdSphere(pos - vec3(0.8, -0.6, -0.1), 0.15);
    float sphereRightCheekCutoutDis = sdSphere(pos - vec3(-0.8, -0.6, -0.1), 0.15);
    
    // Eyes + Nose
    float sphereLeftEyeCutoutDis = sdSphere(pos - vec3(0.45, -0.2, -0.9), 0.25); 
    float sphereRightEyeCutoutDis = sdSphere(pos - vec3(-0.45, -0.2, -0.9), 0.25); 
    float lineSegNoseDis = sdCapsule(pos, vec3(0, -0.4, -0.8), vec3(0, 0, -0.7), 0.1);
    float sphereLeftEyeDis = sdSphere(pos - vec3(0.25, -0.2, -0.5), 0.25);
    float sphereRightEyeDis = sdSphere(pos - vec3(-0.25, -0.2, -0.5), 0.25);
    
    // Neck
    float lineSegNeckDis = sdCapsule(pos, vec3(0, -0.2, 0.35), vec3(0, -1.5, 0.3), 0.15);
    float cubeBaseDis = sdCube(pos, vec3(0, -1.5, 0.2), vec3(0.65, 0.1, 0.65), 0.04);
    float ringBaseDis = sdTorus(pos, vec3(0, -1.5, 0.2), vec2(1.0, 0.3));
    
    // Final distance to return
    float finalDis = opSmoothUnion(sphereTopHeadDis, sphereBotHeadDis, 0.4);
    finalDis = opSmoothUnion(finalDis, lineSegNoseDis, 0.1);
    finalDis = opSmoothSubtraction(sphereBotHeadCutoutDis, finalDis, 0.4);
    finalDis = opSmoothSubtraction(sphereLeftCheekCutoutDis, finalDis, 0.4);
    finalDis = opSmoothSubtraction(sphereRightCheekCutoutDis, finalDis, 0.4);
    finalDis = opSmoothSubtraction(sphereLeftEyeCutoutDis, finalDis, 0.2);
    finalDis = opSmoothSubtraction(sphereRightEyeCutoutDis, finalDis, 0.2);
    finalDis = opSmoothUnion(finalDis, sphereLeftEyeDis, 0.025);
    finalDis = opSmoothUnion(finalDis, sphereRightEyeDis, 0.025);
    finalDis = opSmoothUnion(finalDis, lineSegNeckDis, 0.2);
    finalDis = opSmoothUnion(finalDis, cubeBaseDis, 0.5);
    finalDis = opSubtraction(ringBaseDis * 0.8, finalDis * 0.8);
    return finalDis;
}

// Vectors //////////////////////////////////////////////////////////////////////

// Get distance function
vec2 getDist(vec3 pos) {
    
    // Transformations
    vec3 originPos = pos - vec3(0, 1, 6);
    originPos.xz *= rotate(iTime * 0.5);
    
    // Bust -> Glass Plate -> Wooden Table
    float bustDis = sdBust(originPos);
    float plateDis = sdCappedCylinder(originPos - vec3(0, -1.7, 0), 1.5, 0.05);
    plateDis += texture(iChannel1, pos.xz).r * 0.001;
    float tableDis = sdCube(pos, vec3(0, -1.45, 6), vec3(3.0, 0.4, 3.0), 0.1);
    tableDis += texture(iChannel2, pos.xz).r * 0.08;
    
    // Final distance to return
    float finalDis = opUnion(bustDis, plateDis);
    finalDis = opUnion(finalDis, tableDis);
    
    // Final material to return
    int mat = 0;
    if (finalDis == bustDis) {
        mat = MATWAX;
    } else if (finalDis == plateDis) {
        mat = MATGLASS;
    } else if (finalDis == tableDis) {
        mat = MATWOOD;
    }
    
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
    int mat = 0;
    
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

// Lighting function
vec4 getLight(vec3 pos, vec3 lightOrigin) {
    // Get the light origin
    vec3 lightPos = lightOrigin; // basically a point in 3D space
    lightPos.xz += vec2(sin(iTime), cos(iTime)) * 8.0;
    //lightPos.x += 3.0 * sin(iTime + lightOrigin.z);
    
    // Get the light ray
    vec3 lightRay = normalize(lightPos - pos);
    
    // Get normal
    vec3 normal = getNormal(pos);
    
    // Get diffused lighting 
    // We want 0 if the rays are parallel, 1 if the rays are perpendicular
    // Hence we use dot product
    float dif = clamp(dot(normal, lightRay), 0.0, 1.0);
    float dis = rayMarch(pos + normal * SURFDIS * 2.0, lightRay).x;
    
    // Getting shadows (if disToSphere < disToLight then we have a shadow)
    if (dis < length(lightPos - pos)) {
        float shadowIntensity = 0.2;
        dif *= shadowIntensity;
    }
    
    return vec4(dif, lightRay);
}

// Main //////////////////////////////////////////////////////////////////////

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab UV Coord
    vec2 uv = (fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;

    // Declare Col
    vec3 col = vec3(0.0);

    // Setup Camera
    float camHeight = 1.8;
    float downTilt = -0.2;
    vec3 rayOrigin = vec3(0, camHeight, 0); // this is the camera (origin of vector)
    vec3 rayDir = normalize(vec3(uv.x, uv.y + downTilt, 1));
    
    // Lighting Setup
    float IOR = 1.45;
    vec3 fillDiffuseLightPos = vec3(0, 0.4, 12); //vec3(-1.5, 0.4, 12);
    vec3 keyDiffuseLightPos = vec3(0, 0.5, 2); //vec3(1.5, 0.5, 2);
    float maxLightDis = 8.0;
    
    // Texture Setup
    vec3 texCol = texture(iChannel0, rayDir).rgb;
    col = texCol;
    vec3 woodTex = texture(iChannel2, rayDir.xy).rgb;
    
    // Visualise 
    vec2 passedDMat = rayMarch(rayOrigin, rayDir);
    float dis = passedDMat.x;
    int mat = int(passedDMat.y);
    
    if (dis < MAXDIS) {
        vec3 pos = rayOrigin + rayDir * dis;
        vec3 normal = getNormal(pos);
        vec4 keyDif = getLight(pos, keyDiffuseLightPos);
        float keyDiffuseLightVal = keyDif.x;
        vec4 fillDif = getLight(pos, fillDiffuseLightPos);
        float fillDiffuseLightVal = fillDif.x;
        float finalFill = fillDiffuseLightVal * clamp(1.0 - (length(fillDiffuseLightPos - pos) / maxLightDis), 0.0, 1.0);
        float finalKey = keyDiffuseLightVal * clamp(1.0 - (length(keyDiffuseLightPos - pos) / maxLightDis), 0.0, 1.0);
        
        // Gradual lighting curve
        float lightStep = clamp(keyDiffuseLightVal + fillDiffuseLightVal, 0.0, 1.0);
        vec3 mixCol = vec3(0.0); //mix(shadowCol, lightCol, lightStep);
        
        // Textures
        if (mat == MATWOOD) {
            mixCol += woodTex * lightStep;
        }
        
        // Local Thickness (raymarch again for thickness of material)
        vec3 rayDirIn = refract(rayDir, normal, 1.0 / IOR);
        vec3 posEnter = pos - normal * SURFDIS * 3.0;
        float disIn = rayMarch(posEnter, rayDirIn).x;
        
        // Refraction
        if (mat == MATGLASS) {
            vec3 posExit = posEnter + rayDirIn * disIn;
            vec3 normalExit = -getNormal(posExit);
            vec3 rayDirOut = refract(rayDirIn, normalExit, IOR / 1.0);
            // Total Internal Reflection
            if (length(rayDirOut) == 0.0) {
                rayDirOut = rayDirIn;
            }
            vec3 refractedTex = texture(iChannel2, rayDirOut.xy).rgb;
            mixCol += refractedTex * lightStep;
        }
        
        // Back Translucency + Surface Scatter 
        /*
        Source: https://www.alanzucconi.com/2017/08/30/fast-subsurface-scattering-1/
        Here we want to vary the translucency of the sdf using the direction of the light and viewport.
        lightRay -> light vector from light origin to surface
        delta -> higher values cause more intense light scattering, lower values reuslt in a 'thinner' looking material
        power -> intensity of the curve dropoff of the back translucency
        */
        if (mat == MATWAX) {
            float delta, backPower; 
            
            // Key Light
            backPower = 8.0;
            delta = 0.7; 
            vec3 keyLightRay = -keyDif.yzw; 
            float keyBackIntensity = dot(rayDir, -normalize(keyLightRay + normal * delta));
            keyBackIntensity = pow(clamp(keyBackIntensity, 0.0, 1.0), backPower); 
            mixCol = mix(mixCol, keyLightCol, keyBackIntensity + length(disIn));
            
            // Fill Light
            backPower = 2.0;
            delta = 0.2;
            vec3 fillLightRay = -fillDif.yzw; 
            float fillBackIntensity = dot(rayDir, -normalize(fillLightRay + normal * delta));
            fillBackIntensity = pow(clamp(fillBackIntensity, 0.0, 1.0), backPower);
            mixCol = mix(mixCol, fillLightCol, 0.25 * fillBackIntensity + length(disIn));
        }
        
        // Reflections
        float reflectiveVal = 0.1;
        vec3 reflectedRay = reflect(rayDir, normal);
        vec3 reflectiveTex = texture(iChannel0, reflectedRay).rgb;
        mixCol += reflectiveTex * reflectiveVal;
        
        // Specular Lighting
        float specStrength = 0.25;
        float specPower = 8.0;
        vec3 reflectedLightRay = reflect(keyDif.yzw, normal);
        float spec = pow(max(dot(rayDir, reflectedLightRay), 0.0), specPower);
        vec3 specCol = specStrength * spec * keyLightCol;
        mixCol += specCol;
        
        // Lighting
        mixCol += finalKey * keyLightCol;
        mixCol += finalFill * fillLightCol;
        
        // Final colour
        col = vec3(mixCol);
    }
    
    // Output to screen
    fragColor = vec4(col,1.0);
}
