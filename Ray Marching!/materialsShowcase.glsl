/*
Van Andrew Nguyen
29/11/2021
[Materials Showcase: Chessboard]

My most technical shader to date. I combined everything I've learnt on crafting material properties; namely reflection, 
refraction and sub-surface scattering to add multiple bounces for each. The shader is a chessboard with eight different 
pawns of different materials; from left to right: brushed metal, clay, acrylic / clear plastic, glass, hard plastic, 
metal, wax, rock.
I'm happy with how it looks however, the shader itself takes quite a while to render (multiple ray marches per pawn, 
three passes, times eight pawns, so my laptop gets toasty). The board itself is reflective which plays off on the pawn's reflections.

The only missing piece I haven't discovered how to write is how light sources interact with objects, i.e. how light passes
through a transparent or bounces of a metallic surface. So far the pawns interact amongst themselves but never with the light.
*/

// Globals //////////////////////////////////////////////////////////////////////

/*
Max steps is the maximum amount of steps we can make to get a ray length
Max distance is the maximum distance we can shoot a ray
Surf distance is the distance we loop to get to register a hit on the surface
*/
#define MAXSTEPS 100
#define MAXDIS 32.0
#define SURFDIS 0.01

#define MATNULL 0
#define MATGLASS 1
#define MATCLEARPLASTIC 2
#define MATSOLIDPLASTIC 3
#define MATMETAL 4
#define MATMATTEMETAL 5
#define MATWAX 6
#define MATROCK 7
#define MATCLAY 9

const vec3 shadowCol = vec3(0.2, 0.1, 0.15);
const vec3 lightCol = vec3(0.8, 0.75, 0.7);
const vec3 fogCol = vec3(0.615, 0.533, 0.666);
const vec3 fillLightCol = vec3(0.854, 0.364, 0.501);
const vec3 keyLightCol = vec3(0.937, 0.376, 0.101);

const vec3 planeCol1 = vec3(0.921, 0.847, 0.698);
const vec3 planeCol2 = vec3(0.254, 0.235, 0.227);
const vec3 plasticInteriorCol = vec3(0.145, 0.894, 0.894);
const vec3 baseWaxCol = vec3(0.639, 0.525, 0.458);
const vec3 baseClayCol = vec3(0.807, 0.796, 0.713);
const vec3 baseRockCol = vec3(0.356, 0.349, 0.431);

// Noise /////////////////////////////////////////////////////////////////////

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

// SDF for Chess Pawn
float sdChessPiece(vec3 pos, vec3 start) {
    float foot = sdCylinder(pos, start, start + vec3(0, 0.15, 0), 0.3);
    float base = sdTorus(pos, start + vec3(0, 0.2, 0), vec2(0.2, 0.05));
    float upperBase = sdTorus(pos, start + vec3(0, 0.3, 0), vec2(0.12, 0.04));
    float stem = sdCylinder(pos, start, start + vec3(0, 0.7, 0), 0.1);
    float neck = sdTorus(pos, start + vec3(0, 0.65, 0), vec2(0.12, 0.04));
    float head = sdSphere(pos, start + vec3(0, 0.8, 0), 0.18);
    
    float dis = opSmoothUnion(base, foot, 0.05);
    dis = opSmoothUnion(dis, upperBase, 0.15);
    dis = opSmoothUnion(dis, stem, 0.1);
    dis = opSmoothUnion(dis, neck, 0.05);
    dis = opUnion(dis, head);
    
    return dis;
}

// Noise using textures
float sdRockTex(vec3 pos, float amp1, float amp2) {
    float dis = 0.0;
    dis += texture(iChannel1, pos.xy / 2.5).r * amp1;
    dis += texture(iChannel1, pos.xy / 1.0).r * amp2;
    return dis;
}

float sdClayTex(vec3 pos, float amp1, float amp2) {
    float dis = 0.0;
    dis += texture(iChannel2, pos.xy / 10.0).r * amp1;
    dis += texture(iChannel2, pos.xy / 4.0).r * amp2;
    return dis * 0.8;
}

// Vectors //////////////////////////////////////////////////////////////////////

// Get distance function
vec2 getDist(vec3 pos, bool includeWater) {
    
    // Get SDF dis
    float planeHeight = -0.5;
    float planeDis = sdCube(pos - vec3(0, planeHeight, 6), vec3(4, 0.02, 4));

    float pawnClay = sdChessPiece(pos, vec3(-1, planeHeight, 4.5)) + sdClayTex(pos, 0.02, 0.01);
    float pawnGlass = sdChessPiece(pos, vec3(0, planeHeight, 4.5));
    float pawnMetal = sdChessPiece(pos, vec3(1, planeHeight, 4.5));
    float pawnRock = sdChessPiece(pos, vec3(2, planeHeight, 4.5)) + sdRockTex(pos, 0.03, 0.015);
    
    float pawnMatteMetal = sdChessPiece(pos, vec3(-1.5, planeHeight, 3.5));
    float pawnClearPlastic = sdChessPiece(pos, vec3(-0.5, planeHeight, 3.5));
    float pawnSolidPlastic = sdChessPiece(pos, vec3(0.5, planeHeight, 3.5)) + sdRockTex(pos, 0.01, 0.005);
    float pawnWax = sdChessPiece(pos, vec3(1.5, planeHeight, 3.5)); 
    
    // Final distance to return
    float finalDis = planeDis;
    finalDis = opUnion(finalDis, pawnClay);
    finalDis = opUnion(finalDis, pawnGlass);
    finalDis = opUnion(finalDis, pawnMetal);
    finalDis = opUnion(finalDis, pawnMatteMetal);
    finalDis = opUnion(finalDis, pawnClearPlastic);
    finalDis = opUnion(finalDis, pawnSolidPlastic);
    finalDis = opUnion(finalDis, pawnRock);
    finalDis = opUnion(finalDis, pawnWax);
    
    // Final material to return
    int mat = -1;
    if (finalDis == planeDis) {
        mat = MATNULL;
    } else if (finalDis == pawnClay) {
        mat = MATCLAY;
    } else if (finalDis == pawnGlass) {
        mat = MATGLASS;
    } else if (finalDis == pawnMetal) {
        mat = MATMETAL;
    } else if (finalDis == pawnMatteMetal) {
        mat = MATMATTEMETAL;
    } else if (finalDis == pawnClearPlastic) {
        mat = MATCLEARPLASTIC;
    } else if (finalDis == pawnSolidPlastic) {
        mat = MATSOLIDPLASTIC;
    } else if (finalDis == pawnWax) {
        mat = MATWAX;
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
vec2 rayMarch(vec3 rayOrigin, vec3 rayDir, bool includeWater, float side) {
    float disFromOrigin = 0.0;
    int mat = 0;
    
    // Now we loop for max steps or until the distance is very small
    for (int i=0;i<MAXSTEPS;i++) {
        // Declare current marching location
        vec3 pos = rayOrigin + rayDir * disFromOrigin;
        vec2 passedDMat = getDist(pos, includeWater);
        float disToScene = passedDMat.x * side;
        mat = int(passedDMat.y);
        
        disFromOrigin += disToScene;
        
        // Exit condition; we have a hit
        if (disFromOrigin > MAXDIS || disToScene < SURFDIS) { break; }
    }
    
    return vec2(disFromOrigin, mat);
}

// Lighting function
vec4 getLight(vec3 pos, vec3 lightOrigin, bool includeWater) {
    // Get the light origin
    vec3 lightPos = lightOrigin; // basically a point in 3D space
    
    // Get the light ray
    vec3 lightRay = normalize(lightPos - pos);
    
    // Get normal
    vec3 normal = getNormal(pos, includeWater);
    
    // Get diffused lighting 
    // We want 0 if the rays are parallel, 1 if the rays are perpendicular
    // Hence we use dot product
    float dif = clamp(dot(normal, lightRay), 0.0, 1.0);
    float dis = rayMarch(pos + normal * SURFDIS * 2.0, lightRay, includeWater, 1.0).x;
    
    // Getting shadows (if disToSphere < disToLight then we have a shadow)
    if (dis < length(lightPos - pos)) {
        float shadowIntensity = 0.2;
        dif *= shadowIntensity;
    }
    
    return vec4(dif, lightRay);
}

// Main //////////////////////////////////////////////////////////////////////

vec3 getReflectivity(int mat, vec3 normal, vec3 rayDir) {
    vec3 ref = vec3(0.0);
    // Don't forget fresnel (non-metals need it)
    float fresnel = clamp(1.0 - dot(normal, -rayDir), 0.0, 1.0);
    
    // Switch up ref based on material
    // We can add colour by changing ref to become a col, e.g. 1,0,0
    // Of course this only changes the colour of a reflective material (e.g. clay must be changed outside)
    if (mat == MATNULL) {
        ref = vec3(0.2 * fresnel);
    } else if (mat == MATGLASS) {
        ref = vec3(0.0);
    } else if (mat == MATCLEARPLASTIC) {
        ref = vec3(0.0);
    } else if (mat == MATSOLIDPLASTIC) {
        ref = vec3(0.4 * fresnel);
    } else if (mat == MATMETAL) {
        ref = vec3(0.9);
    } else if (mat == MATMATTEMETAL) {
        ref = vec3(0.5);
    } else if (mat == MATWAX) {
        ref = vec3(0.1 * fresnel);
    } else if (mat == MATROCK) {
        ref = vec3(0.0);
    } else if (mat == MATCLAY) {
        ref = vec3(0.0);
    }
    
    return ref;
}

float getIOR(int mat) {
    float IOR = 0.0;
    
    // Switch up ref based on material
    // We can add colour by changing ref to become a col, e.g. 1,0,0
    // Of course this only changes the colour of a reflective material (e.g. clay must be changed outside)
    if (mat == MATGLASS) {
        IOR = 1.52;
    } else if (mat == MATCLEARPLASTIC) {
        IOR = 1.46;
    }
    
    return IOR;
}

vec3 renderImage(inout vec3 rayOrigin, inout vec3 rayDir, 
                 inout vec3 reflectiveness, inout float IOR,
                 inout bool doBounce, inout bool doRefract) {
    // Assign cubemap
    vec3 col = texture(iChannel0, rayDir).rgb;

    // Lighting Setup
    vec3 keyLightPos = vec3(6, 4, 6);
    vec3 fillLightPos = vec3(-2, 8, 2);
    float lightMaxLen = 10.0;
    vec3 ambientLight = vec3(0.15);
    vec3 mixCol = vec3(ambientLight); 
    
    // Visualise 
    vec2 passedDMat = rayMarch(rayOrigin, rayDir, false, 1.0);
    float dis = passedDMat.x;
    int mat = int(passedDMat.y);
    
    if (dis < MAXDIS) {
        vec3 pos = rayOrigin + rayDir * dis;
        vec3 normal = getNormal(pos, false);
        vec3 posExtEnter = pos + normal * SURFDIS * 3.0;
        vec3 posInEnter = pos - normal * SURFDIS * 3.0;
        
        // Lock the lights based on their distance to scene
        vec4 keyDiff = getLight(pos, keyLightPos, true);
        vec4 fillDiff = getLight(pos, fillLightPos, true);
        float keyDiffuseLight = keyDiff.x;
        float fillDiffuseLight = fillDiff.x;
        float finalFill = fillDiffuseLight * clamp(1.0 - (length(fillLightPos - pos) / lightMaxLen), 0.0, 1.0);
        float finalKey = keyDiffuseLight * clamp(1.0 - (length(keyLightPos - pos) / lightMaxLen), 0.0, 1.0);

        // Materials!
        reflectiveness = getReflectivity(mat, normal, rayDir);
        IOR = getIOR(mat);
        if (mat == MATNULL) {
            // Checkerboard pattern
            float size = 1.0;
            float total = mod(floor(pos.z * size), 2.0) + 
                          mod(floor(pos.x * size), 2.0);
            bool isEven = mod(total, 2.0) == 0.0;
            mixCol *= (isEven) ? planeCol1 : planeCol2;
        } else {
            mixCol *= vec3(finalKey + finalFill);
        }
        
        // Reflections!
        if (reflectiveness != vec3(0.0)) {
            doBounce = true;
            vec3 reflectRay = reflect(rayDir, normal);
            if (mat == MATMATTEMETAL) {
                // If the material is matte metal we want to randomise the reflected rays for an inaccurate reflection)
                float spacing = 1.6; // lower is more refined mix
                float amp = 0.25;
                float rand1 = amp * (hash1(normal.x * spacing) - 0.5);
                float rand2 = amp * (hash1(normal.y * spacing) - 0.5);
                float rand3 = amp * (hash1(normal.z * spacing) - 0.5);
                vec3 rayDirRand = rayDir + vec3(rand1, rand2, rand3);
                reflectRay = reflect(rayDirRand, normal);
            }
            
            // Change the ray variables (since we are marching from surface);
            rayOrigin = posExtEnter;
            rayDir = reflectRay;
        }
        
        // Refraction!
        if (IOR != 0.0) {
            vec3 refractRayIn = refract(rayDir, normal, 1.0 / IOR);
            float disIn = rayMarch(posInEnter, refractRayIn, false, -1.0).x;
            vec3 posExit = rayOrigin + disIn * refractRayIn;
            vec3 normalExit = -getNormal(posExit, false);
            vec3 refractRayOut = refract(refractRayIn, normalExit, IOR);
            
            // Total internal reflection
            if (length(refractRayOut) == 0.0) {
                refractRayOut = reflect(refractRayIn, normalExit);
            }
            
            // Clear plastic (interior colour based on thickness of material)
            if (mat == MATCLEARPLASTIC) {
                float maxThickness = 0.8;
                vec3 interiorCol = mix(mixCol, plasticInteriorCol, clamp(abs(disIn / maxThickness), 0.0, 1.0));
                mixCol = interiorCol;
            }
            
            // Change variables
            doRefract = true;
            rayOrigin = posExit - normalExit * SURFDIS * 3.0;
            rayDir = refractRayOut;
        }
        
        // Subsurface scattering!
        if (mat == MATWAX) {
            mixCol = baseWaxCol;
            float delta, backPower; 
            backPower = 16.0;
            delta = 0.1; 
            
            // Key Light
            vec3 refractRayIn = refract(rayDir, normal, 1.0 / 1.45);
            float disIn = rayMarch(posInEnter, refractRayIn, false, -1.0).x;
            vec3 keyLightRay = -keyDiff.yzw; 
            float keyBackIntensity = dot(rayDir, -normalize(keyLightRay + normal * delta));
            keyBackIntensity = pow(clamp(keyBackIntensity, 0.0, 1.0), backPower); 
            mixCol = mix(mixCol, keyLightCol, length(disIn));
        }
        
        // Colours
        if (mat == MATCLAY) { mixCol *= baseClayCol; }
        if (mat == MATROCK) { mixCol *= baseRockCol; }

        // Lighting
        mixCol += finalKey * keyLightCol;
        mixCol += finalFill * fillLightCol;
        
        col = vec3(mixCol);
    } 
    
    return col;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab UV Coord
    vec2 uv = (fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;
    
    // Setup Camera
    float camXOffset = 0.5;
    float camYOffset = 1.8;
    float downTilt = -0.6;
    float sideTilt = -0.1;
    float camZoom = 1.1;
    vec3 rayOrigin = vec3(camXOffset, camYOffset, 0); // this is the camera (origin of vector)
    vec3 rayDir = normalize(vec3(uv.x + sideTilt, uv.y + downTilt, camZoom));

    // Render
    bool doBounce = false;
    bool doRefract = false;
    vec3 reflectiveness = vec3(0.0);
    float IOR = 0.0;
    vec3 col = renderImage(rayOrigin, rayDir, reflectiveness, IOR, doBounce, doRefract);
    if (doRefract) {
        vec3 refractCol = IOR * renderImage(rayOrigin, rayDir, reflectiveness, IOR, doBounce, doRefract);
        col += refractCol;
    }
    if (doBounce) {
        vec3 bounce = reflectiveness * renderImage(rayOrigin, rayDir, reflectiveness, IOR, doBounce, doRefract);
        col += bounce;
    }
    
    // Gamma Correction
    col = pow(col, vec3(0.4545));
    
    // Output to screen
    fragColor = vec4(col,1.0);
}