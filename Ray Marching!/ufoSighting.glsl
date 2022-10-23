/*
Van Andrew Nguyen
23/10/22
[UFO-Sighting]

Wanted to give voronoi-textures a shot with using them as a heightmap for landscapes. Played around with a basic scene showing a UFO sdf shining a light onto the terrain.
*/


vec3 getShipPos(float time) {
    return vec3(0, 0.2 * sin(time) + 1.0, 0);
}

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
#define MATSHIP 1
#define MATGROUND 2

const vec3 shadowCol = vec3(0.2, 0.1, 0.15);
const vec3 lightCol = vec3(0.8, 0.75, 0.7);
const vec3 fillLightCol = vec3(0.854, 0.364, 0.501);
const vec3 keyLightCol = vec3(0.937, 0.376, 0.101);
const vec3 dirtCol = vec3(0.17,0.14,0.12);
const vec3 snowCol = vec3(0.44,0.5,0.55);
const vec3 grassCol = vec3(0.8, 0.1, 0.2);
const vec3 beamCol = vec3(0.24,0.83,0.95);

// Noise /////////////////////////////////////////////////////////////////////

// Return a randomised float using a vec2 (uv coord)
float hash1(float n) {
    return fract( n*17.0*fract( n*0.3183099 ) );
}

float hash21(vec2 uv) {
    return fract(hash1(uv.x * uv.y));
}

vec2 hash22(vec2 p) {
    vec3 o = fract(p.xyx * vec3(123.34, 234.34, 345.56));
    o += dot(o, o+34.45);
    return fract(vec2(o.x*o.y, o.y*o.z));
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

float voronoi(vec2 uv) {
    // Voronoi
    vec2 cellID = floor(uv);
    vec2 cellUV = fract(uv);
    float minDist = 1.0;
    
    for (int y= -1; y <= 1; y++) {
        for (int x= -1; x <= 1; x++) {
            // Get neighbours in 3x3 kernel
            vec2 neighbour = vec2(float(x),float(y));
            // Get animated random point and use difference to set intensity
            vec2 point = hash22(cellID + neighbour);
            point = 0.5 + 0.5 * sin(iTime * 0.1 + PI * 2.0 * point);
            vec2 diff = neighbour + point - cellUV;
            float dist = length(diff);
            minDist = min(minDist, dist);
        }
    }
    
    return minDist;
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

float sdFloor(vec3 pos) {
    vec2 uv = vec2(pos.x, pos.z);
    float roundness = 1.25;
    float amp = 1.0 + texture(iChannel0, uv * 0.22).r;
    float heightMap = pow(voronoi(uv * 0.5), roundness) * amp;
    return (pos.y + heightMap) * 0.5; 
}

float sdShip(vec3 pos) {
    vec3 shipPos = getShipPos(iTime);
    float torus1Dis = sdTorus(pos, vec3(shipPos.x, shipPos.y + 0.3, shipPos.z), vec2(0.8, 0.02));
    float torus2Dis = sdTorus(pos, vec3(shipPos.x, shipPos.y + 0.4, shipPos.z), vec2(0.6, 0.03));
    float torus3Dis = sdTorus(pos, vec3(shipPos.x, shipPos.y + 0.5, shipPos.z), vec2(0.4, 0.04));
    float sphereDis = sdSphere(pos, vec3(shipPos.x, shipPos.y + 0.6, shipPos.z), 0.1);
    
    float dis = opSmoothUnion(torus3Dis, sphereDis, 0.5);
    dis = opSmoothUnion(dis, torus2Dis, 0.3);
    dis = opSmoothUnion(dis, torus1Dis, 0.4);

    return dis;
}

// Get distance function
vec2 getDist(vec3 pos) {
    
    // Map dist to shapes
    float planeDis = sdFloor(pos);
    float torusDis = sdShip(pos);
    
    // Final distance to return
    float finalDis = opUnion(torusDis, planeDis);
    
    // Final material to return
    int mat = MATNULL;
    if (finalDis == planeDis) {
        mat = MATGROUND;
    } else if (finalDis == torusDis) {
        mat = MATSHIP;
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

vec3 background(vec3 rayDir) {
    vec3 col = vec3(0.0);
    
    float y = rayDir.y * 0.2 + 0.8; // light is top, dark is bottom. 
    col += y * fillLightCol;
    float x = rayDir.x * 0.6 + 0.4;
    col += x * keyLightCol * 0.4;
    
    return col;
}

vec3 render(vec3 rayOrigin, vec3 rayDir) {
    // Assign col
    vec3 col = vec3(0.0);

    // Lighting Setup
    vec3 keyLightPos = vec3(3, 4, 2);
    vec3 fillLightPos = vec3(-2, 8, 2);
    float lightMaxLen = 10.0;
    vec3 ambientLight = vec3(0.2);
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
        vec3 sP = getShipPos(iTime);
        vec3 shipRimLightPos = vec3(sP.x, sP.y + 0.5, sP.z);
        vec3 shipBeamLightPos = vec3(sP.x, sP.y - 0.7, sP.z);
        float shipSpotLight = getSpotLight(pos, shipBeamLightPos, vec3(0, -1, 0), PI / 6.0, 0.02);
        
        
        // Materials
        if (mat == MATNULL) {
            // Checkerboard pattern
            float size = 1.0;
            float total = mod(floor(pos.z * size), 2.0) + 
                          mod(floor(pos.x * size), 2.0);
            bool isEven = mod(total, 2.0) == 0.0;
            mixCol *= (isEven) ? shadowCol : lightCol;
        } else if (mat == MATGROUND) {
            vec3 localUpDir = vec3(0,1,0);
            float steepness = dot(normal, localUpDir);
            // Grass and Snow
            mixCol = mix(dirtCol, grassCol, smoothstep(0.7, 0.4, steepness));
            mixCol = mix(mixCol, snowCol, smoothstep(0.6, 1.0, steepness));
        } else if (mat == MATSHIP) {
            float fresnel = clamp(1.0 - dot(normal, -rayDir), 0.0, 2.0);
            mixCol = vec3(fresnel);
            vec3 reflectedRay = reflect(rayDir, normal);
            vec3 reflectiveTex = texture(iChannel1, reflectedRay).rgb;
            mixCol = 0.8 * (reflectiveTex * 0.5 + (beamCol * 0.5));
        }
        
        // Ambient Occlusion (check out Alex Evans)
        // Trick is to sample sdf along normal
        float ao = getAmbientOcclusion(pos, normal);
        mixCol *= ao;
        
        // Lighting
        mixCol += finalKey * keyLightCol;
        mixCol += finalFill * fillLightCol;
        mixCol += shipSpotLight * beamCol;
        
        col = vec3(mixCol);
    }
    
    // Fog / Sky
    float fogStart = MAXDIS / 4.0;
    float fogEnd = MAXDIS;
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
    float camZoom = 2.0;
    vec3 rayOrigin = vec3(camX, camY, camZ); 
    rayOrigin.yz *= rotate(-m.y * PI + 1.);
    rayOrigin.xz *= rotate(-m.x * 2.0 * PI);
    vec3 rayDir = getRayDir(uv, rayOrigin, vec3(camYTilt, 0, camXTilt), camZoom);
    
    col = render(rayOrigin, rayDir);
    
    // Output to screen
    fragColor = vec4(col,1.0);
}
