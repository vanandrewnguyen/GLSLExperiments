/*
Van Andrew Nguyen
10/10/2021
[Glass Marbles]

This was an adventurous exploration of refraction and reflection. I created five spheres using sdf's and played with how much light
is absorbed, the material roughness, and index of refraction. The resultant effect looks really nice! I had to tweak my ray marching
template to account for different 'materials' within the shader and it was great to finally understand my code on my own without tutorial help.
*/

// Globals //////////////////////////////////////////////////////////////////////

/*
Max steps is the maximum amount of steps we can make to get a ray length
Max distance is the maximum distance we can shoot a ray
Surf distance is the distance we loop to get to register a hit on the surface
*/
#define MAXSTEPS 100
#define MAXDIS 50.0
#define SURFDIS 0.01
#define MAXBOUNCES 2

#define MATBALL1 1
#define MATBALL2 2
#define MATBALL3 3
#define MATBALL4 4
#define MATBALL5 5

#define MATSMUDGE 10
#define MATSCRATCH 11

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

// Minimum function for vec2 (check x component)
vec2 vec2Min(vec2 a, vec2 b) {
    return a.x < b.x ? a : b;
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
float sdSphere(vec3 pos, vec3 center, float rad) {
    pos -= center;
    return length(pos) - rad;
}

// Cube Distance
float sdCube(vec3 pos, vec3 center, vec3 size) {
    pos -= center;
    
    // We get the length of the position - size (dis) and grab a version > 0
    float val = length(max(abs(pos) - size, 0.0));
    return val;
}

// Vectors //////////////////////////////////////////////////////////////////////

// Get distance function
vec2 getDist(vec3 pos) {
    
    // Get dist of the various shapes
    float rad = 0.7;
    float offs = 1.8;
    float t = iTime * 2.0;
    float amp = 0.4;
    float sphere1Dis = sdSphere(pos, vec3(-offs*2.0, 1.0+amp*sin(t-2.0), 6), rad); 
    float sphere2Dis = sdSphere(pos, vec3(-offs, 1.0+amp*sin(t-1.0), 6), rad); 
    float sphere3Dis = sdSphere(pos, vec3(0, 1.0+amp*sin(t), 6), rad); 
    float sphere4Dis = sdSphere(pos, vec3(offs, 1.0+amp*sin(t+1.0), 6), rad); 
    float sphere5Dis = sdSphere(pos, vec3(offs*2.0, 1.0+amp*sin(t+2.0), 6), rad); 
    
    // Final distance to return
    float finalDis = sphere1Dis;
    finalDis = opUnion(finalDis, sphere2Dis);
    finalDis = opUnion(finalDis, sphere3Dis);
    finalDis = opUnion(finalDis, sphere4Dis);
    finalDis = opUnion(finalDis, sphere5Dis);
    
    // Materials (index to return)
    int mat = 0;
    if (finalDis == sphere1Dis) {
        mat = MATBALL1;
    } else if (finalDis == sphere2Dis) {
        mat = MATBALL2;
    } else if (finalDis == sphere3Dis) {
        mat = MATBALL3;
    } else if (finalDis == sphere4Dis) {
        mat = MATBALL4;
    } else if (finalDis == sphere5Dis) {
        mat = MATBALL5;
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
vec2 rayMarch(vec3 rayOrigin, vec3 rayDir, float side) {
    float disFromOrigin = 0.0;
    int mat = 0;
    
    // Now we loop for max steps or until the distance is very small
    for (int i=0;i<MAXSTEPS;i++) {
        // Declare current marching location
        vec3 pos = rayOrigin + rayDir * disFromOrigin;
        vec2 passedDist = getDist(pos);
        float disToScene = passedDist.x * side;
        mat = int(passedDist.y);
        
        disFromOrigin += disToScene;
        
        // Exit condition; we have a hit
        if (disFromOrigin > MAXDIS || disToScene < SURFDIS) { break; }
    }
    
    return vec2(disFromOrigin, mat);
}

// Lighting function
float getLight(vec3 pos, vec3 lightOrigin) {
    // Get the light origin
    vec3 lightPos = lightOrigin; // basically a point in 3D space
    lightPos.xz += vec2(sin(iTime), cos(iTime)) * 4.0;
    
    // Get the light ray
    vec3 light = normalize(lightPos - pos);
    
    // Get normal
    vec3 normal = getNormal(pos);
    
    // Get diffused lighting 
    // We want 0 if the rays are parallel, 1 if the rays are perpendicular
    // Hence we use dot product
    float dif = clamp(dot(normal, light), 0.0, 1.0);
    float dis = rayMarch(pos + normal * SURFDIS * 2.0, light, 1.0).x;
    
    // Getting shadows (if disToSphere < disToLight then we have a shadow)
    if (dis < length(lightPos - pos)) {
        float shadowIntensity = 0.2;
        dif *= shadowIntensity;
    }
    
    return dif;
}

// Main //////////////////////////////////////////////////////////////////////

vec3 renderObject(inout vec3 rayOrigin, inout vec3 rayDir, inout vec3 refVal, inout float roughness, bool last) {
    // Init Colours
    vec3 col = texture(iChannel0, rayDir).rgb;
    vec3 shadowCol = vec3(0.2, 0.1, 0.15);
    vec3 lightCol = vec3(0.8, 0.75, 0.7);
    vec3 marbleCol = vec3(0.2, 0.043, 0);
    
    // Visualise 
    vec2 passedDist = rayMarch(rayOrigin, rayDir, 1.0);
    float dis = passedDist.x;
    int mat = int(passedDist.y);
    
    // Reset reflective value
    refVal = vec3(0.0);
    
    if (dis < MAXDIS) {
        // Get position and normal
        vec3 pos = rayOrigin + rayDir * dis;
        vec3 normal = getNormal(pos);
        float diffuseLight = getLight(pos, vec3(0, 5, 6));
        
        // Gradual lighting curve
        vec3 mixCol = vec3(0.0); //mix(shadowCol, lightCol, diffuseLight);
        
        // Reflection (setup reflected ray, texture to draw, and fresnel value)
        vec3 reflectedRay = reflect(rayDir, normal);
        vec3 refTex = texture(iChannel0, reflectedRay).rgb;
        float fresnel = 1.0 - dot(normal, -rayDir); // calculate how many rays are hitting (depending on angle of camera)
        fresnel = clamp(fresnel, 0.0, 1.0);
        fresnel = pow(fresnel, 4.0);
        
        // Refraction
        float IOR = 1.0;
        float baseIOR = 1.0;
        float matDensity = 1.0;
        /*
        Water is 1.33 / Air is 1.01 / Diamond is 2.4
        Each ball material is slightly different, from being clear and transparent to dense and opaque.
        */
        if (mat == MATBALL1) {
            IOR = 1.01;
            matDensity = 1.0;
            roughness = 0.4;
        } else if (mat == MATBALL2) {
            IOR = 1.33;
            matDensity = 0.85;
            roughness = 0.2;
        } else if (mat == MATBALL3) {
            IOR = 1.50;
            matDensity = 0.65;
            roughness = 0.15;
        } else if (mat == MATBALL4) {
            IOR = 2.13;
            matDensity = 0.45;
            roughness = 0.05;
        } else if (mat == MATBALL5) {
            IOR = 2.40;
            matDensity = 0.25;
            roughness = 0.01;
        }
        
        // March inside the object (in then out, similiar to above march)
        vec3 refractedRayIn = refract(rayDir, normal, baseIOR / IOR);
        vec3 posEnter = pos - normal * SURFDIS * 3.0;
        float disIn = rayMarch(posEnter, refractedRayIn, -1.0).x;
        vec3 posExit = posEnter + refractedRayIn * disIn;
        vec3 normalExit = -getNormal(posExit);
        vec3 refractedRayOut = refract(refractedRayIn, normalExit, IOR / baseIOR);
        
        // Total Internal Reflection
        if (length(refractedRayOut) == 0.0) {
            refractedRayOut = reflectedRay;
        }
        
        // Change reftex from reflected to refracted, then add back reflection
        refTex = texture(iChannel0, refractedRayOut).rgb;
        refTex += texture(iChannel0, reflectedRay).rgb * roughness;
        
        // Materials
        float altLight = 0.4 + 0.6 * diffuseLight;
        float lightAbsorption = 1.2 * clamp((1.0 - fresnel), 0.0, 1.0);
        
        if (mat == MATBALL1 || mat == MATBALL2 || mat == MATBALL3 || mat == MATBALL4 || mat == MATBALL5) {
            mixCol = mix(refTex, vec3(0.0), clamp(lightAbsorption * (1.0 - matDensity), 0.0, 1.0));
            refVal = vec3(0.8);
            if (last) {
                col += mixCol;//refVal * refTex;
            }
        }
        
        // Ray march again to get reflections of other sdfs 
        rayOrigin = pos + normal * SURFDIS * 3.0;
        rayDir = reflectedRay;
        
        col = vec3(mixCol);
    }
    
    return col;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab UV Coord
    vec2 uv = (fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;

    // Declare Col
    vec3 col = vec3(0.0);
    vec3 bounce = vec3(0.0);
    
    // Setup Camera
    float camHeight = 1.8;
    float downTilt = -0.2;
    float camZoom = 1.0;
    vec3 rayOrigin = vec3(0, camHeight, 0); // this is the camera (origin of vector)
    vec3 rayDir = normalize(vec3(uv.x, uv.y + downTilt, camZoom));
    
    // Render image
    vec3 refVal = vec3(0.0); //reflectivity
    vec3 refFilter = vec3(1.0); 
    float roughness = 1.0;
    col = renderObject(rayOrigin, rayDir, refVal, roughness, false);
    // Loop for as many reflection bounces as we need 
    for (int i=0;i<MAXBOUNCES;i++) {
        refFilter *= refVal;
        bounce = refFilter * renderObject(rayOrigin, rayDir, refVal, roughness, i==MAXBOUNCES-1);
        col += bounce * refVal * roughness;
    }
    
    // Output to screen
    fragColor = vec4(col,1.0);
}