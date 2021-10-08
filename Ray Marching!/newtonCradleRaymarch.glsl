/*
Van Andrew Nguyen
07/10/21
[Newton's Cradle]

I followed the Art of Code's tutorial to create a newton's cradle model. Currently it just uses reflections from the cube map.
This is going to change; I am attempting to add reflections between the objects themselves through ray marching.
Furthermore, I can experiment with surface imperfects by dulling and etching on the surfaces of the sdfs.

Update: I have added reflections.
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
#define MAXBOUNCES 1

#define MATBASE   1
#define MATBAR    2
#define MATBALL   3
#define MATSTRING 4

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
float sdTorus(vec3 pos, vec2 rad, float thickness) {
    // We subtract the smaller radius from the length of the vector running from
    // origin point to middle of torus
    float x = length(pos.xy) - rad.x; //xz
    return length(vec2(x, pos.z)) - rad.y + thickness; //y
}

// Sphere Distance
float sdSphere(vec3 pos, float rad) {
    return length(pos) - rad;
}

// Cube Distance
float sdCube(vec3 pos, vec3 size) { 
    pos = abs(pos) - size;

    // We get the length of the position - size (dis) and grab a version > 0
    float val = length(max(abs(pos) - size, 0.0));
    return val;
}

// 2D Square Distance
float sdSquare(vec2 pos, vec2 size) {
    pos = abs(pos) - size;

    // We get the length of the position - size (dis) and grab a version > 0
    float val = length(max(pos, 0.0)) + min(max(pos.x, pos.y), 0.0);
    return val;
}

// Newton Ball (Combine a sphere with a torus ring)
vec2 sdNewtonBall(vec3 pos, vec3 origin, float rad, float angle) {
    // Local variables
    float ropeHeight = 1.3 + (rad/2.0);
    float ropeWidth = 0.65;
    
    // Set new drawing origin
    pos.y -= ropeHeight + rad; // move the origin up for rotation
    vec3 newPos = pos - origin;
    newPos.xy *= rotate(angle); // rotation
    newPos.y += ropeHeight + rad; // then reset origin
    vec3 newPosElevated = newPos - vec3(0, rad, 0); 
    
    
    // Draw the ball, hook and rope
    float sphereDis = sdSphere(newPos, rad);
    float torusDis = sdTorus(newPosElevated, vec2(0.05), 0.04);
    newPosElevated.z = abs(newPosElevated.z);
    float lineDis = sdCapsule(newPosElevated, vec3(0), vec3(0, ropeHeight, ropeWidth), 0.002);
    
    // Return value
    float val = opSmoothUnion(sphereDis, torusDis, 0.025);
    val = opUnion(val, lineDis);
        
    return vec2(val, val == lineDis ? MATSTRING : MATBALL);
}

// Vectors //////////////////////////////////////////////////////////////////////

// Get distance function
vec2 getDist(vec3 pos) {
    
    // Setup origin coords so we can apply transformations easily
    float t = sin(iTime*4.0);
    vec3 origin = vec3(0, 0, 5);
    pos -= origin;
    
    pos.xz *= rotate(iTime * 0.4);
    
    // Base (3D Box)
    float cubeBevel = 0.08;
    float cubeDis = sdCube(pos, vec3(0.8, 0.1, 0.4)) - cubeBevel; 
    cubeDis = opSubtraction(-cubeDis, -pos.y); // cut the box in half
    
    // Bars (2D Box)
    float barBevel = 0.2;
    float barThickness = 0.05;
    float barWidthDis = 0.65;
    float barDis = length(vec2(sdSquare(pos.xy, vec2(1.3, 2.4)) - barBevel, 
                               abs(pos.z) - barWidthDis)) - barThickness;
    barDis = opSubtraction(-barDis, -pos.y); // cut the bar in half
    
    // Balls (Sphere), we create five copies
    pos += origin;
    float rad = 0.2;
    float timingLeft = min(0.0, t);
    float timingRight = max(0.0, t);
    
   vec2 ballOneDis = sdNewtonBall(pos - origin, vec3(rad*4.0, 1, 0), rad, timingLeft),
         ballTwoDis = sdNewtonBall(pos - origin, vec3(rad*2.0, 1, 0), rad, (t+timingLeft)*0.05),
         ballThreeDis = sdNewtonBall(pos - origin, vec3(0, 1, 0), rad, t*0.05),
         ballFourDis = sdNewtonBall(pos - origin, vec3(-rad*2.0, 1, 0), rad, (t+timingRight)*0.05),
         ballFiveDis = sdNewtonBall(pos - origin, vec3(-rad*4.0, 1, 0), rad, timingRight);
    vec2 ballCollectionDis = vec2Min(ballOneDis, 
                              vec2Min(ballTwoDis, 
                              vec2Min(ballThreeDis, 
                              vec2Min(ballFourDis, ballFiveDis))));
    
    // Surface Imperfections
    cubeDis += sin(pos.x * 8.0) * 0.002; 
    
    // Final distance to return
    float finalDis = cubeDis;
    finalDis = opUnion(finalDis, barDis);
    finalDis = opUnion(finalDis, ballCollectionDis.x);
    
    // Materials (index to return)
    int mat = 0;
    if (finalDis == cubeDis) {
        mat = MATBASE;
    } else if (finalDis == barDis) {
        mat = MATBAR;
    } else if (finalDis == ballCollectionDis.x) {
        mat = int(ballCollectionDis.y); 
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
    
    vec2 grabbedDist = vec2(0.0);
    
    // Now we loop for max steps or until the distance is very small
    for (int i=0;i<MAXSTEPS;i++) {
        // Declare current marching location
        vec3 pos = rayOrigin + rayDir * disFromOrigin;
        grabbedDist = getDist(pos);
        
        float disToScene = grabbedDist.x;
        disFromOrigin += disToScene;
        
        // Exit condition; we have a hit
        if (disFromOrigin > MAXDIS || disToScene < SURFDIS) { break; }
    }
    
    return vec2(disFromOrigin, grabbedDist.y);
}

// Lighting function
float getLight(vec3 pos, vec3 lightOrigin) {
    // Get the light origin
    vec3 lightPos = lightOrigin; // basically a point in 3D space
    //lightPos.xz += vec2(sin(iTime), cos(iTime)) * 4.0;
    
    // Get the light ray
    vec3 light = normalize(lightPos - pos);
    
    // Get normal
    vec3 normal = getNormal(pos);
    
    // Get diffused lighting 
    // We want 0 if the rays are parallel, 1 if the rays are perpendicular
    // Hence we use dot product
    float dif = clamp(dot(normal, light), 0.0, 1.0);
    float dis = rayMarch(pos + normal * SURFDIS * 2.0, light).x;
    
    // Getting shadows (if disToSphere < disToLight then we have a shadow)
    if (dis < length(lightPos - pos)) {
        float shadowIntensity = 0.2;
        dif *= shadowIntensity;
    }
    
    return dif;
}

// Main //////////////////////////////////////////////////////////////////////

/*
Little note; you can set variables as 'inout' so if you change it within the function,
it changes outside of the function as well like a normal pointer.
*/
vec3 renderObject(inout vec3 rayOrigin, inout vec3 rayDir, inout vec3 refVal, bool last) {
    // Init Colours
    vec3 col = texture(iChannel0, rayDir).rgb;
    vec3 shadowCol = vec3(0.2, 0.1, 0.15);
    vec3 lightCol = vec3(0.8, 0.75, 0.7);
    
    // Visualise the object
    vec2 grabbedDist = rayMarch(rayOrigin, rayDir);
    float dis = grabbedDist.x;
    
    refVal = vec3(0.0); // reset reflection before render
    
    if (dis < MAXDIS) {
        vec3 pos = rayOrigin + rayDir * dis;
        vec3 normal = getNormal(pos);
        float diffuseLight = getLight(pos, vec3(0, 5, 6));
        
        // Reflection
        vec3 reflectedRay = reflect(rayDir, normal);
        vec3 refTex = texture(iChannel0, reflectedRay).rgb;
        float fresnel = 1.0 - dot(normal, -rayDir); // calculate how many rays are hitting (depending on angle of camera)
        fresnel = clamp(fresnel, 0.0, 1.0);
        fresnel = pow(fresnel, 4.0);
        
        // Gradual lighting curve
        vec3 mixCol = mix(shadowCol, lightCol, diffuseLight);
        
        // Materials
        // Here we change the colour and reflectivity of each material index
        int mat = int(grabbedDist.y);
        float brighterLight = (0.4 + 0.6 * diffuseLight);
        if (mat == MATBASE) {
            mixCol = vec3(0.313, 0.105, 0.145) * diffuseLight; 
            refVal = vec3(mix(0.01, 0.5, fresnel)); //fresnel;
        } else if (mat == MATBAR) {
            mixCol = (0.5 + 0.5 * refTex) * vec3(0.992, 0.890, 0.866) * brighterLight; 
            refVal = vec3(0.8);
        } else if (mat == MATBALL) {
            mixCol = refTex * vec3(0.992, 0.890, 0.866) * brighterLight;
            refVal = vec3(0.8);
            if (last) {
                col += refVal * refTex;
            }
        } else if (mat == MATSTRING) {
            mixCol = vec3(0.968, 0.945, 0.929) * diffuseLight;
            refVal = vec3(0.0);
        }
        
        // Ray march again to get reflections of other sdfs (called in 'bounce')
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
    float camHeight = 2.4;
    float t = 0.5 + 0.5 * sin(iTime);
    float downTilt = mix(-0.5, -0.3, t); //-0.5; //-0.3;
    float camZoom = mix(2.4, 1.2, t);    //2.4; //1.2;  //2.0 + sin(iTime);
    vec3 rayOrigin = vec3(0, camHeight, 0); 
    vec3 rayDir = normalize(vec3(uv.x, uv.y + downTilt, camZoom));
    
    // Render image
    vec3 refVal = vec3(0.0); //reflectivity
    vec3 refFilter = vec3(1.0); 
    col = renderObject(rayOrigin, rayDir, refVal, false);
    // Loop for as many reflection bounces as we need 
    for (int i=0;i<MAXBOUNCES;i++) {
        refFilter *= refVal;
        bounce = refFilter * renderObject(rayOrigin, rayDir, refVal, i==MAXBOUNCES-1);
        col += bounce * refVal;
    }
    
    // Output to screen
    fragColor = vec4(col,1.0);
}