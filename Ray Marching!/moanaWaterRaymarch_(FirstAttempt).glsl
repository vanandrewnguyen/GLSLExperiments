/*
Van Andrew Nguyen
10/01/21
[Raymarched Water]

I honestly struggled with this one a lot and deviated from the 'correct' way of volumetric lighting.
Saw this blog post: https://wallisc.github.io/rendering/2020/12/08/Making-Of-Moana-the-shadertoy.html
I tried to follow it; here's what I got up to.
1. Build a ray marcher
2. Render two planes, one for the ocean bed, one for the water surface
3. Apply fbm + superimposed sine waves on water surface to get 'natural look'
4. Cut out a sphere using smooth subtraction and place camera within it
5. Apply light absorption onto water - the code is sampled from the blog
6. Ray march for a second time after we render the water; to get the distances of the ocean bed (this deviates from marching inside the volume)
7. Render the ocean bed and light it in accordance to the water
8. Refract light rays through the water to get that warping effect of the ocean bed

That's all I managed to get up; in the blog there are a dozen more steps, which I feel is way over my depth.
Still, it was a fun and fustrating project. Resultant effect is, well, a less detailed version of Chris's shader.
*/

// Globals //////////////////////////////////////////////////////////////////////

/*
Max steps is the maximum amount of steps we can make to get a ray length
Max distance is the maximum distance we can shoot a ray
Surf distance is the distance we loop to get to register a hit on the surface
*/
#define MAXSTEPS 100
#define MAXDIS 40.0
#define SURFDIS 0.01

// Noise //////////////////////////////////////////////////////////////////////

// Return a randomised float using a vec2 (uv coord)
float hash1(float n) {
    return fract( n*17.0*fract( n*0.3183099 ) );
}

float hash21(vec2 uv) {
    /*
    // Pseudo-random math function to get a random number
    vec2 o = fract(uv * vec2(145.56, 7498.12));
    o += dot(o, o + 159.89);
    return fract(o.x * o.y);
    */
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
float fbm(vec2 pos, int iterations, float t) {
    // Declare return value, amplitude of motion, freq of motion
    float val = 0.0;
    float amp = 0.05;
    float freq = 4.0;
    
    // Now loop through layers and return the combined value
    for (int i=0;i<iterations;i++) {
        val += amp * smoothNoise(freq * pos + t); 
        amp *= 0.5;
        freq * 2.0;
    }
    
    return val;
}

// SDFS //////////////////////////////////////////////////////////////////////

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

// Water SDF
float sdWater(vec3 pos) {
    float height = pos.y;
    
    // Superimpose sine waves on water surface
    float amp = 0.2;
    float freq = 1.0;
    float iterations = 2.0;
    float t = iTime;
    
    for (float i=0.0;i<iterations;i+=1.0) {
        vec2 offset = pos.xz * (0.5 - 0.5 * (i / iterations));
        float displacement = amp * 
                             sin(pos.x * smoothNoise(offset) + (pos.x + i + offset.x) * freq + t) * 
                             sin(pos.z + (pos.z + i + offset.y) * freq - t);
        height += displacement;
        
        amp *= 0.5;
        freq *= 2.0;
    }
    
    // Finally add noise as small touch on water surface
    height += fbm(pos.xz, 4, t);
    
    // Return
    return height;
}

// Sand SDF
float sdSand(vec3 pos) {
    // Return something lower than pos.y + a bunch of sine waves
    return pos.y + 1.0 + 0.005 * sin(pos.x * 16.0) * sin(pos.z * 16.0);
}

// Coral SDF
float sdCoral(vec3 pos) {
    // Return a bunch of spheres
    float s1 = sdSphere(pos, vec3(-1.2, -1.0, 3.6), 0.6);
    float s2 = sdSphere(pos, vec3(-1.5, -1.0, 2.4), 0.3);
    float s3 = sdSphere(pos, vec3(1.6, -1.0, 2.8), 0.25);

    float dis = opUnion(s1, s2);
    dis = opUnion(dis, s3);
    return dis;
}

// Vectors //////////////////////////////////////////////////////////////////////

// Get distance function
vec4 getDist(vec3 pos, bool includeWater) {
    // Colours
    vec3 finalCol = vec3(0.0);
    vec3 waterCol = vec3(0.431, 0.968, 0.968);
    vec3 sandCol = vec3(1.0, 0.96, 0.68);
    vec3 coralCol = vec3(0.28, 0.25, 0.);
    
    // Create two layers, one for the water surface, other for ocean bottom
    float sphereDis = sdSphere(pos, vec3(0, 1.3, 1.2), 2.5); 
    float coralDis = sdCoral(pos);
    float planeBottomDis = sdSand(pos); 
    float planeTopDis = sdWater(pos);
    
    // Final distance to return
    float finalDis = 0.0;
    if (includeWater == true) {
        finalDis = opSmoothSubtraction(sphereDis, planeTopDis * 0.8, 0.9); //planeBottomDis; 
        finalDis = opUnion(finalDis, planeBottomDis);
        finalDis = opUnion(finalDis, coralDis);
    } else {
        finalDis = opUnion(planeBottomDis, coralDis);
    }
    
    // Final colour to return
    if (finalDis == planeBottomDis) { 
        finalCol = sandCol; 
    } else if (finalDis == coralDis) {
        finalCol = coralCol;
    } else { finalCol = waterCol; }
    
    //float finalDis = planeTopDis * 0.9;
    
    return vec4(finalDis, finalCol);
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
vec4 rayMarch(vec3 rayOrigin, vec3 rayDir, bool includeWater) {
    float disFromOrigin = 0.0;
    vec3 passedColour = vec3(0.0);
    
    // Now we loop for max steps or until the distance is very small
    for (int i=0;i<MAXSTEPS;i++) {
        // Declare current marching location
        vec3 pos = rayOrigin + rayDir * disFromOrigin;
        vec4 passedDist = getDist(pos, includeWater);
        float disToScene = passedDist.x;
        passedColour = passedDist.yzw;
        disFromOrigin += disToScene;
        
        // Exit condition; we have a hit
        if (disFromOrigin > MAXDIS || disToScene < SURFDIS) { break; }
    }
    
    return vec4(disFromOrigin, passedColour);
}

// See absorption rate of light through a medium
float beerLambert(float absorptionCoefficient, float distanceTraveled) {
    // High is less dense medium like water, Low is dense medium like milk
    return exp(-absorptionCoefficient * distanceTraveled);
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
    vec4 passedDist = rayMarch(pos + normal * SURFDIS * 2.0, light, includeWater);
    float dis = passedDist.x;
    
    // Getting shadows (if disToSphere < disToLight then we have a shadow)
    if (dis < length(lightPos - pos)) {
        float shadowIntensity = 0.2;
        dif *= shadowIntensity;
    }
    
    return dif;
}

// Background function; depth perception
vec3 background(vec3 rayDir) {
    vec3 skyCol = vec3(0.807, 0.850, 0.992);
    
    vec3 col = vec3(0.0);

    vec3 specCol = skyCol; //vec3(0.2, 0.1, 0.15);
    float ratio = 0.7;
    float y = rayDir.y * (1.0 - ratio) + ratio; // light is top, dark is bottom. 
    col += y * specCol;
    
    return col;
}

// Main //////////////////////////////////////////////////////////////////////

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Grab UV Coord
    vec2 uv = (fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;

    // Declare Col for Lighting + Material
    vec3 col = vec3(0.0);
    vec3 mixCol = vec3(0.0);
    vec3 shadowCol = vec3(0.2, 0.1, 0.15);
    vec3 lightCol = vec3(0.8, 0.75, 0.7);
    
    vec3 waterCol = vec3(0.431, 0.968, 0.968); //vec3(0.43, 0.89, 0.97);
    vec3 sandCol = vec3(1.0, 0.96, 0.68);
    vec3 coralCol = vec3(0.28, 0.25, 0.);

    // Setup Camera
    float camHeight = -0.3;
    float downTilt = 0.02;
    vec3 rayOrigin = vec3(0, camHeight, 0); // this is the camera (origin of vector)
    vec3 rayDir = normalize(vec3(uv.x, uv.y + downTilt, 1));
    vec3 lightPos = vec3(0, 5, 6);
    
    // Visualise 
    vec4 passedDist = rayMarch(rayOrigin, rayDir, true);
    float dis = passedDist.x;
    vec3 passedCol = passedDist.yzw;
    
    // Water 
    if (dis < MAXDIS) {
        // Get position, normal and lighting
        vec3 pos = rayOrigin + rayDir * dis;
        vec3 normal = getNormal(pos, true);
        float diffuseLight = getLight(pos, lightPos, true);
        
        // Refract the dang light ray
        float IORAir = 1.00;
        float IORWater = 1.04; //0.96;
        rayDir = refract(rayDir, normal, IORAir / IORWater);
        
        // Gradual lighting curve
        mixCol = mix(shadowCol, lightCol, diffuseLight);
        
        /*
        1. Ray march inside a volume
        2. Calculate light absorption at each step
        */
        
        // We are calculating light absorption here - we haven't checked what's behind the water
        float opaqueVisiblity = 1.0;
        float marchSize = 0.6;
        float volumeDepth = 0.0;
        float opaqueDepth = 100.0;
        float volumeAlbedo = 0.1;
        vec3 posEnter = pos - normal * SURFDIS * 4.0;
        
        for(int i=0;i<MAXSTEPS;i++) {
            volumeDepth += marchSize;
            if(volumeDepth > opaqueDepth) break;

            vec3 posIn = posEnter + volumeDepth * rayDir;
            
            // This is a bad check but it works for now
            if (passedCol == waterCol) {
                float previousOpaqueVisiblity = opaqueVisiblity;
                opaqueVisiblity *= beerLambert(0.05, marchSize);
                float absorptionFromMarch = previousOpaqueVisiblity - opaqueVisiblity;
                float lightDis = length((lightPos - posIn));
                mixCol += absorptionFromMarch * volumeAlbedo * (diffuseLight * (lightDis / dis));
            }
        }
        
        // Add it to the output colour
        col += vec3(mixCol) * passedCol;
    }
    
    vec4 passedDistBelow = rayMarch(rayOrigin, rayDir, false);
    float disBelow = passedDistBelow.x;
    vec3 passedColBelow = passedDistBelow.yzw;
    
    // Ocean Bed
    if (disBelow < MAXDIS) {
        vec3 posBelow = rayOrigin + rayDir * disBelow;
        vec3 normalBelow = getNormal(posBelow, false);
        float diffuseLightBelow = getLight(posBelow, lightPos, false);
        
        vec3 pos = rayOrigin + rayDir * dis;
        float maxOpaqueDis = MAXDIS;
        float waterDensity = 0.4;
        
        // MODULATE is the number we need to grab from above position in relevance to the camera pos
        // 0.0 is solid water, 1.0 is completely transparent
        // We need to strike a balance.
        // So we kinda need the distance to the dis of the sea bed - dis to nearest ocean pixel NORMALISED
        // So here I get the length, absolute value, clamped
        
        float modulate = 1.0 - clamp(abs(length(posBelow - pos) / maxOpaqueDis), 0.0, 1.0); //0.1;
        mixCol = mix(shadowCol, lightCol, diffuseLightBelow) * modulate * waterDensity;
        
        // Output colour
        col += vec3(mixCol) * passedColBelow * (diffuseLightBelow - modulate * 0.5);
    }
    
    // Fog
    float fogStart = 0.0;
    float fogEnd = MAXDIS;
    col = mix(col, background(rayDir), smoothstep(fogStart, fogEnd, dis));
    
    // Output to screen
    fragColor = vec4(col, 1.0);
}
