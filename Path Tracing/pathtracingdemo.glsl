/*
Van Andrew Nguyen
22/4/22
[Path Tracing]

Dabbled into path-tracing, I wanted to try learn sample theory and other methods of rendering scenes.
*/

// Utility
// Credits: https://blog.demofox.org/2020/06/06/casual-shadertoy-path-tracing-2-image-improvement-and-glossy-reflections/
vec3 LessThan(vec3 f, float value) {
    return vec3(
        (f.x < value) ? 1.0f : 0.0f,
        (f.y < value) ? 1.0f : 0.0f,
        (f.z < value) ? 1.0f : 0.0f);
}
 
vec3 LinearToSRGB(vec3 rgb) {
    rgb = clamp(rgb, 0.0f, 1.0f);
     
    return mix(
        pow(rgb, vec3(1.0f / 2.4f)) * 1.055f - 0.055f,
        rgb * 12.92f,
        LessThan(rgb, 0.0031308f)
    );
}
 
vec3 SRGBToLinear(vec3 rgb) {
    rgb = clamp(rgb, 0.0f, 1.0f);
     
    return mix(
        pow(((rgb + 0.055f) / 1.055f), vec3(2.4f)),
        rgb / 12.92f,
        LessThan(rgb, 0.04045f)
    );
}

// ACES tone mapping curve fit to go from HDR to LDR
//https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
vec3 ACESFilm(vec3 x) {
    float a = 2.51f;
    float b = 0.03f;
    float c = 2.43f;
    float d = 0.59f;
    float e = 0.14f;
    return clamp((x*(a*x + b)) / (x*(c*x + d) + e), 0.0f, 1.0f);
}

// Shlick approx
float FresnelReflectAmount(float n1, float n2, vec3 normal, vec3 incident, float f0, float f90) {
        // Schlick aproximation
        float r0 = (n1-n2) / (n1+n2);
        r0 *= r0;
        float cosX = -dot(normal, incident);
        if (n1 > n2)
        {
            float n = n1/n2;
            float sinT2 = n*n*(1.0-cosX*cosX);
            // Total internal reflection
            if (sinT2 > 1.0)
                return f90;
            cosX = sqrt(1.0-sinT2);
        }
        float x = 1.0-cosX;
        float ret = r0+(1.0-r0)*x*x*x*x*x;
 
        // adjust reflect multiplier for object reflectivity
        return mix(f0, f90, ret);
}

// Buffer A
#define SURFDIS 0.01
#define MAXDIS 100.0
#define MAXBOUNCES 8

#define RENDERPERFRAME 8

#define PI 3.1416
#define FOV 90.0

// Credits to https://blog.demofox.org/2020/05/25/casual-shadertoy-path-tracing-1-basic-camera-diffuse-emissive/

struct materialData {
    // emissive = how much light is cast from this object
    // albedo = base colour of the object
    // chance = percentage of light reflected / refracted
    // roughness = blurriness of reflection / refraction
    // col = colour of the specular reflection, e.g. coloured metals, coloured glass
    // IOR = index of refraction
    vec3 emissive;
    vec3 albedo;
    float specChance;
    float specRoughness;
    vec3 specCol;
    float IOR;
    float refractionChance;
    float refractionRoughness;
    vec3 refractionCol;
};

// Keeps data of the ray in a struct
struct rayHitData {
    float dist;
    vec3 normal;
    bool fromInside;
    materialData matData;
};

// Default material (like python fixture)
materialData getDefaultMat() {
    materialData mat;
    mat.albedo = vec3(0.0);
    mat.emissive = vec3(0.0);
    mat.specChance = 0.0;
    mat.specRoughness = 0.0;
    mat.specCol = vec3(0.0);
    mat.IOR = 1.0;
    mat.refractionChance = 0.0;
    mat.refractionRoughness = 0.0;
    mat.refractionCol = vec3(0.0);
    return mat;
}

// Randomises the seed used to bounce rays off surfaces
uint wangHash(inout uint seed) {
    // Hash func by Thomas Wang
    seed = uint(seed ^ uint(61)) ^ uint(seed >> uint(16));
    seed *= uint(9);
    seed = seed ^ (seed >> 4);
    seed *= uint(0x27d4eb2d);
    seed = seed ^ (seed >> 15);
    return seed;
}

float hash11(inout uint seed) {
    return float(wangHash(seed)) / 4294967296.0;
}

// Generates a random unit vector along the positive hemisphere
vec3 generateRandomUnitVector(inout uint seed) {
    // Remember we are influence the seed for every ray bounce, hence it's an inout var
    float z = hash11(seed) * 2.0 - 1.0;
    float a = hash11(seed) * PI * 2.0;
    float r = sqrt(1.0 - z * z);
    float x = r * cos(a);
    float y = r * sin(a);
    return vec3(x, y, z);
}

float scalarTriple(vec3 u, vec3 v, vec3 w) {
    // Find normal of u, v vectors and get mag on vector w
    return dot(cross(u, v), w);
}

// For indepth analysis on the ray-sphere intersections check out: https://www.youtube.com/watch?v=8fWZM8hCX5E
// Intersection funct for quad
bool testQuadTrace(in vec3 rayPos, in vec3 rayDir, inout rayHitData data, in vec3 a, in vec3 b, in vec3 c, in vec3 d) {
    // calculate normal and flip vertices order if needed
    vec3 normal = normalize(cross(c-a, c-b));
    if (dot(normal, rayDir) > 0.0f) {
        normal *= -1.0f;
        
		vec3 temp = d;
        d = a;
        a = temp;
        
        temp = b;
        b = c;
        c = temp;
    }
    
    vec3 p = rayPos;
    vec3 q = rayPos + rayDir;
    vec3 pq = q - p;
    vec3 pa = a - p;
    vec3 pb = b - p;
    vec3 pc = c - p;
    
    // determine which triangle to test against by testing against diagonal first
    vec3 m = cross(pc, pq);
    float v = dot(pa, m);
    vec3 intersectPos;
    if (v >= 0.0f) {
        // test against triangle a,b,c
        float u = -dot(pb, m);
        if (u < 0.0f) return false;
        float w = scalarTriple(pq, pb, pa);
        if (w < 0.0f) return false;
        float denom = 1.0f / (u+v+w);
        u*=denom;
        v*=denom;
        w*=denom;
        intersectPos = u*a+v*b+w*c;
    } else {
        vec3 pd = d - p;
        float u = dot(pd, m);
        if (u < 0.0f) return false;
        float w = scalarTriple(pq, pa, pd);
        if (w < 0.0f) return false;
        v = -v;
        float denom = 1.0f / (u+v+w);
        u*=denom;
        v*=denom;
        w*=denom;
        intersectPos = u*a+v*d+w*c;
    }
    
    float dist;
    if (abs(rayDir.x) > 0.1f) {
        dist = (intersectPos.x - rayPos.x) / rayDir.x;
    } else if (abs(rayDir.y) > 0.1f) {
        dist = (intersectPos.y - rayPos.y) / rayDir.y;
    } else {
        dist = (intersectPos.z - rayPos.z) / rayDir.z;
    }
    
	if (dist > SURFDIS && dist < data.dist) {
        data.dist = dist;        
        data.normal = normal; 
        data.fromInside = false;
        return true;
    }    
    
    return false;
}

// Intersection funct for box
bool testBoxTrace(in vec3 rayPos, in vec3 rayDir, inout rayHitData data, in vec3 verts[8]) {
    // Technically could call testQuadTrace() for 6 faces, but that's slow
    /*
    0 1 2 3 = bottom vert in ccw
    4 5 6 7 = top vert in ccw
            7--------6
	       /|       /|
	      4--------5 |
	      | |      | |
	      | 3------|-2
	      |/       |/
	      0--------1    
    */
    if (testQuadTrace(rayPos, rayDir, data, verts[0], verts[1], verts[5], verts[4])) { data.normal *= -1.0f; return true; }
    if (testQuadTrace(rayPos, rayDir, data, verts[1], verts[2], verts[6], verts[5])) { data.normal *= -1.0f; return true; }
    if (testQuadTrace(rayPos, rayDir, data, verts[0], verts[3], verts[7], verts[4])) { data.normal *= -1.0f; return true; }
    if (testQuadTrace(rayPos, rayDir, data, verts[3], verts[2], verts[6], verts[7])) { data.normal *= -1.0f; return true; }
    if (testQuadTrace(rayPos, rayDir, data, verts[4], verts[5], verts[6], verts[7])) { data.normal *= -1.0f; return true; }
    if (testQuadTrace(rayPos, rayDir, data, verts[0], verts[1], verts[2], verts[3])) { data.normal *= -1.0f; return true; }

    return false;
}

// Intersection funct for sphere
bool testSphereTrace(in vec3 rayPos, in vec3 rayDir, inout rayHitData data, in vec4 sphere) {
	//get the vector from the center of this sphere to where the ray begins.
	vec3 m = rayPos - sphere.xyz;

    //get the dot product of the above vector and the ray's vector
	float b = dot(m, rayDir);

	float c = dot(m, m) - sphere.w * sphere.w;

	//exit if r's origin outside s (c > 0) and r pointing away from s (b > 0)
	if(c > 0.0 && b > 0.0) {
		return false;
    }

	//calculate discriminant
	float discr = b * b - c;

	//a negative discriminant corresponds to ray missing sphere
	if(discr < 0.0) {
		return false;
    }
    
	//ray now found to intersect sphere, compute smallest t value of intersection
    bool fromInside = false;
	float dist = -b - sqrt(discr);
    if (dist < 0.0f) {
        fromInside = true;
        dist = -b + sqrt(discr);
    }
    
	if (dist > SURFDIS && dist < data.dist) {
        data.fromInside = fromInside;
        data.dist = dist;        
        data.normal = normalize((rayPos+rayDir*dist) - sphere.xyz) * (fromInside ? -1.0f : 1.0f);
        return true;
    }
    
    return false;
}

// Get scene 
void scene(in vec3 rayPos, in vec3 rayDir, inout rayHitData hitInfo) {
    float wallLen = 10.0;
    float farPoint = 24.0;
    float nearPoint = 18.0;
    
    // Right wall
    {
        vec3 A = vec3(wallLen, -wallLen, farPoint);
        vec3 B = vec3(wallLen, -wallLen, nearPoint);
        vec3 C = vec3(wallLen,  wallLen, nearPoint);
        vec3 D = vec3(wallLen,  wallLen, farPoint);
        if (testQuadTrace(rayPos, rayDir, hitInfo, A, B, C, D)) {
            hitInfo.matData = getDefaultMat();
            hitInfo.matData.albedo = vec3(0.1f, 0.7f, 0.1f);
        }
    } 
    
    // Left wall
    {
        vec3 A = vec3(-wallLen, -wallLen, farPoint);
        vec3 B = vec3(-wallLen, -wallLen, nearPoint);
        vec3 C = vec3(-wallLen,  wallLen, nearPoint);
        vec3 D = vec3(-wallLen,  wallLen, farPoint);
        if (testQuadTrace(rayPos, rayDir, hitInfo, A, B, C, D)) {
            hitInfo.matData = getDefaultMat();
            hitInfo.matData.albedo = vec3(0.7f, 0.1f, 0.1f);
        }
    } 
    
    // Top wall
    {
        vec3 A = vec3(-wallLen, wallLen, farPoint);
        vec3 B = vec3( wallLen, wallLen, farPoint);
        vec3 C = vec3( wallLen, wallLen, nearPoint);
        vec3 D = vec3(-wallLen, wallLen, nearPoint);
        if (testQuadTrace(rayPos, rayDir, hitInfo, A, B, C, D)) {
            hitInfo.matData = getDefaultMat();
            hitInfo.matData.albedo = vec3(0.7f, 0.7f, 0.7f);
        }
    } 
    
    // Bottom wall
    {
        vec3 A = vec3(-wallLen, -wallLen, farPoint);
        vec3 B = vec3( wallLen, -wallLen, farPoint);
        vec3 C = vec3( wallLen, -wallLen, nearPoint);
        vec3 D = vec3(-wallLen, -wallLen, nearPoint);
        if (testQuadTrace(rayPos, rayDir, hitInfo, A, B, C, D)) {
            hitInfo.matData = getDefaultMat();
            hitInfo.matData.albedo = vec3(0.7f, 0.7f, 0.7f);
        }
    } 
    
    // Back wall
    {
        vec3 A = vec3(-wallLen, -wallLen, farPoint);
        vec3 B = vec3( wallLen, -wallLen, farPoint);
        vec3 C = vec3( wallLen,  wallLen, farPoint);
        vec3 D = vec3(-wallLen,  wallLen, farPoint);
        if (testQuadTrace(rayPos, rayDir, hitInfo, A, B, C, D)) {
            hitInfo.matData = getDefaultMat();
            hitInfo.matData.albedo = vec3(0.7f, 0.7f, 0.7f);
        }
    } 
    
    // Lighting
    float pad = 2.0;
    {
        vec3 A = vec3(-wallLen + pad*2.0, wallLen - 0.1, farPoint - pad);
        vec3 B = vec3( wallLen - pad*2.0, wallLen - 0.1, farPoint - pad);
        vec3 C = vec3( wallLen - pad*2.0, wallLen - 0.1, nearPoint + pad);
        vec3 D = vec3(-wallLen + pad*2.0, wallLen - 0.1, nearPoint + pad);
        if (testQuadTrace(rayPos, rayDir, hitInfo, A, B, C, D)) {
            hitInfo.matData = getDefaultMat();
            hitInfo.matData.albedo = vec3(0.7f, 0.7f, 0.7f);
            hitInfo.matData.emissive = vec3(1.0f, 0.9f, 0.7f) * 20.0;
        }
    } 
    
    /*
    // Balls
    for (int i = 0; i < 5; i++) {
        if (testSphereTrace(rayPos, rayDir, hitInfo, vec4(-6.0 + float(i*3), -7.5, 18.0, 1.5))) {
            hitInfo.matData = getDefaultMat();
            hitInfo.matData.albedo = vec3(1.0f, 0.8f, 0.7f);
            hitInfo.matData.specChance = float(i) / 5.0;
            hitInfo.matData.specRoughness = 0.0;
            hitInfo.matData.specCol = hitInfo.matData.albedo;
            hitInfo.matData.IOR = 1.33;
        } 
    }
    
    for (int i = 0; i < 5; i++) {
        if (testSphereTrace(rayPos, rayDir, hitInfo, vec4(-6.0 + float(i*3), -7.5 + 4.0, 20.0, 1.5))) {
            hitInfo.matData = getDefaultMat();
            hitInfo.matData.albedo = vec3(0.8f, 1.0f, 0.6f);
            hitInfo.matData.specChance = 0.9;
            hitInfo.matData.specRoughness = 1.0 - float(i) / 5.0;
            hitInfo.matData.specCol = hitInfo.matData.albedo;
            hitInfo.matData.IOR = 1.33;
        } 
    }
    
    for (int i = 0; i < 5; i++) {
        if (testSphereTrace(rayPos, rayDir, hitInfo, vec4(-6.0 + float(i*3), -7.5 + 4.0 + 4.0, 22.0, 1.5))) {
            hitInfo.matData = getDefaultMat();
            hitInfo.matData.albedo = vec3(0.8f, 0.7f, 1.0f);
            hitInfo.matData.specChance = 0.1;
            hitInfo.matData.specRoughness = 0.5;
            hitInfo.matData.specCol = hitInfo.matData.albedo;
            hitInfo.matData.IOR = 1.33;
            hitInfo.matData.refractionChance = float(i) / 5.0;
            hitInfo.matData.refractionRoughness = 1.0 - float(i) / 5.0;
            hitInfo.matData.refractionCol = hitInfo.matData.albedo;
            
        } 
    }
    */
    
    // Box
    {
        /*
            7--------6
	       /|       /|
	      4--------5 |
	      | |      | |
	      | 3------|-2
	      |/       |/
	      0--------1
        */
        float boxT = -2.0;
        float boxB = -10.0;
        float boxL = -6.0;
        float boxR = -2.0;
        float boxFar = 23.5;
        float boxClose = 19.0;
        
        vec3 verts[8];
        verts[0] = vec3(boxL, boxB, boxClose);
        verts[1] = vec3(boxR, boxB, boxClose);
        verts[2] = vec3(boxR, boxB, boxFar);
        verts[3] = vec3(boxL, boxB, boxFar);
        verts[4] = vec3(boxL, boxT, boxClose);
        verts[5] = vec3(boxR, boxT, boxClose);
        verts[6] = vec3(boxR, boxT, boxFar);
        verts[7] = vec3(boxL, boxT, boxFar);
        if (testBoxTrace(rayPos, rayDir, hitInfo, verts)) {
            hitInfo.matData = getDefaultMat();
            hitInfo.matData.albedo = vec3(0.7);
        }
    }
    
    if (testSphereTrace(rayPos, rayDir, hitInfo, vec4(-8.0f, -8.5f, 19.0f, 1.5f))) {
        hitInfo.matData = getDefaultMat();
        hitInfo.matData.albedo = vec3(1.0f, 0.8f, 0.8f);
        hitInfo.matData.specChance = 0.9;
        hitInfo.matData.specRoughness = 0.0;
        hitInfo.matData.specCol = hitInfo.matData.albedo;
    } 
    
    if (testSphereTrace(rayPos, rayDir, hitInfo, vec4(0.4f, -8.2f, 21.0f, 1.8f))) {
        hitInfo.matData = getDefaultMat();
        hitInfo.matData.albedo = vec3(0.9, 0.6, 0.7);
        hitInfo.matData.specChance = 0.5;
        hitInfo.matData.specRoughness = 0.5;
        hitInfo.matData.specCol = hitInfo.matData.albedo;
    }
    
    if (testSphereTrace(rayPos, rayDir, hitInfo, vec4(5.0f, -8.0f, 20.0f, 2.0f))) {
        hitInfo.matData = getDefaultMat();
        hitInfo.matData.albedo = vec3(0.8f, 0.8f, 1.0f);
        hitInfo.matData.specChance = 0.2;
        hitInfo.matData.specRoughness = 0.3;
        hitInfo.matData.specCol = hitInfo.matData.albedo;
        hitInfo.matData.refractionChance = 1.0;
        hitInfo.matData.refractionRoughness = 0.0;
        hitInfo.matData.IOR = 1.33;
        hitInfo.matData.refractionCol = vec3(0.1,0.2,0);
    }
    
    /*
    if (testSphereTrace(rayPos, rayDir, hitInfo, vec4(-5.0f, -4.0f, 20.0f, 2.0f))) {
        hitInfo.matData = getDefaultMat();
        hitInfo.matData.albedo = vec3(1.0f, 0.8f, 0.8f);
        hitInfo.matData.specChance = 0.9;
        hitInfo.matData.specRoughness = 0.0;
        hitInfo.matData.specCol = hitInfo.matData.albedo;
    } 
     
    if (testSphereTrace(rayPos, rayDir, hitInfo, vec4(0.0f, -6.0f, 20.0f, 2.0f))) {
        hitInfo.matData = getDefaultMat();
        hitInfo.matData.albedo = vec3(0.8f, 1.0f, 0.8f);
        hitInfo.matData.specChance = 0.5;
        hitInfo.matData.specRoughness = 0.5;
        hitInfo.matData.specCol = hitInfo.matData.albedo;
    }    
     
    if (testSphereTrace(rayPos, rayDir, hitInfo, vec4(5.0f, -8.0f, 20.0f, 2.0f))) {
        hitInfo.matData = getDefaultMat();
        hitInfo.matData.albedo = vec3(0.8f, 0.8f, 1.0f);
        hitInfo.matData.specChance = 0.5;
        hitInfo.matData.specRoughness = 0.3;
        hitInfo.matData.specCol = hitInfo.matData.albedo;
        hitInfo.matData.refractionChance = 1.0;
        hitInfo.matData.refractionRoughness = 0.0;
        hitInfo.matData.IOR = 1.33;
        hitInfo.matData.refractionCol = vec3(0.1,0.2,0);
    }
    */
      
}

// Equivalent to scene() in template
vec3 getRayColour(in vec3 startRayPos, in vec3 startRayDir, inout uint randomSeed) {
    // Init
    // throughput is the colour which is bounced within the scene
    vec3 col = vec3(0.0f, 0.0f, 0.0f);
    vec3 throughput = vec3(1.0f, 1.0f, 1.0f);
    vec3 rayPos = startRayPos;
    vec3 rayDir = startRayDir;
     
    for (int bounceIndex = 0; bounceIndex <= MAXBOUNCES; ++bounceIndex) {
        // Shoot ray from the current pixel
        rayHitData hitInfo;
        hitInfo.matData = getDefaultMat();
        hitInfo.dist = MAXDIS;
        hitInfo.fromInside = false;
        
        // Now grab the scene info
        scene(rayPos, rayDir, hitInfo);
         
        // If the ray hits nothing, exit
        if (hitInfo.dist == MAXDIS) {
            // We can also grab data from the cubemap
            //col += SRGBToLinear(texture(iChannel1, rayDir).rgb * throughput);
            break;
        }
        
        // If inside the object, add absorption (beer lambert, it's inverse exp)
        if (hitInfo.fromInside) {
            throughput *= exp(-hitInfo.matData.refractionCol * hitInfo.dist);
        }
        
        // Get chances for bounced ray
        float specChance = hitInfo.matData.specChance;
        float refractionChance = hitInfo.matData.refractionChance;
        
        // Apply fresnel
        float rayProb = 1.0f;
        if (specChance > 0.0) {
            // Material actually has spec, get n1, n2 as index of air (1.0) and mat.IOR
            specChance = FresnelReflectAmount(
            hitInfo.fromInside ? hitInfo.matData.IOR : 1.0,
            !hitInfo.fromInside ? 1.0 : hitInfo.matData.IOR,
            rayDir, hitInfo.normal, hitInfo.matData.specChance, 1.0f);
            
            // Chance mult keeps diff / refraction ratio the same
            float chanceMult = (1.0 - specChance) / (1.0 - hitInfo.matData.specChance);
            refractionChance *= chanceMult;
        }
        
        // Now we can choose whether we do a diff, spec, refractive ray
        float doSpec = 0.0;
        float doRefraction = 0.0;
        float selectRay = hash11(randomSeed);
        if (specChance > 0.0 && selectRay < specChance) {
            doSpec = 1.0;
            rayProb = specChance;
        } else if (refractionChance > 0.0 && selectRay < specChance + refractionChance) {
            doRefraction = 1.0;
            rayProb = refractionChance;
        } else {
            rayProb = 1.0 - (specChance + refractionChance);
        }
        
        // Don't let rayProb be 0.0 since we'll divide by it later
        rayProb = max(rayProb, 0.001);
        
        // Update ray pos (notice how we push away from surface on EITHER side)
        if (doRefraction == 1.0) {
            rayPos = (rayPos + rayDir * hitInfo.dist) - hitInfo.normal * SURFDIS;
        } else {
            rayPos = (rayPos + rayDir * hitInfo.dist) + hitInfo.normal * SURFDIS;
        }
        
        // Now we need a new ray direction
        // Diff uses normal pointing hemisphere to shoot new ray
        // Spec just reflects the incident ray
        // Rough surfaces lerps diff and spec ray (so is more random)
        vec3 diffRayDir = normalize(hitInfo.normal + generateRandomUnitVector(randomSeed));
         
        // Notice how we mix the rays depending on the roughness of the material
        vec3 specRayDir = reflect(rayDir, hitInfo.normal);
        float specRoughness = hitInfo.matData.specRoughness * hitInfo.matData.specRoughness;
        specRayDir = normalize(mix(specRayDir, diffRayDir, specRoughness));

        vec3 refractionRayDir = refract(rayDir, hitInfo.normal, hitInfo.fromInside ? hitInfo.matData.IOR : 1.0f / hitInfo.matData.IOR);
        float refractionRoughness = hitInfo.matData.refractionRoughness * hitInfo.matData.refractionRoughness;
        refractionRayDir = normalize(mix(refractionRayDir, normalize(-hitInfo.normal + generateRandomUnitVector(randomSeed)), refractionRoughness));

        rayDir = mix(diffRayDir, specRayDir, doSpec);
        rayDir = mix(rayDir, refractionRayDir, doRefraction);
        
        // Add in lighting from sources
        col += hitInfo.matData.emissive * throughput;
        
        // Update col mult
        if (doRefraction == 0.0) {
            throughput *= mix(hitInfo.matData.albedo, hitInfo.matData.specCol, doSpec);
        }
        // Divide throughput because we only choose 1/3 rays to shoot (think of FBM)
        throughput /= rayProb;
        
        // Apparently it's called Russian Roulette
        // As throughput approaches 0, the ray is more likely to be terminated early
        // Hence straggling rays are amplified
        float p = max(throughput.r, max(throughput.g, throughput.b)); // max of any channel
        if (hash11(randomSeed) > p) {
            break;
        }
        // Add the energy we 'lose' by randomly terminating paths
        throughput *= 1.0f / p;            
    }

    return col;
}


// Main
void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
    // Init random seed (because we sample random dirs and grab the avg in the 2nd buffer)
    uint randomSeed = uint(uint(fragCoord.x) * uint(1973) + uint(fragCoord.y) * uint(9277) + uint(iFrame) * uint(26699)) | uint(1);
    vec2 uv = fragCoord/iResolution.xy;

    // Setting up camera (target is on -1->1 map)
    float camDis = 1.0 / tan(FOV * 0.5 * PI / 180.0);
    vec3 rayPos = vec3(0.0f, 0.0f, 0.0f);
    // Anti aliasing -> adding a random vector on rayTarget on the imaginary plane
    vec2 targetOffset = vec2(hash11(randomSeed), hash11(randomSeed)) - 0.5;
    vec3 rayTarget = vec3(((fragCoord + targetOffset) / iResolution.xy) * 2.0f - 1.0f, camDis);
    
    // Changing aspect ratio
    float ratio = iResolution.x / iResolution.y;
    rayTarget.y /= ratio;
    vec3 rayDir = normalize(rayTarget - rayPos);
    
    // We need to average the pixel values because we're sampling a bunch of rays that bounce randomly 
    vec3 col = vec3(0);
    if (iFrame > 2) {
        col = getRayColour(rayPos, rayDir, randomSeed);
        vec3 prevCol = texture(iChannel0, uv).rgb;
        col = mix(prevCol, col, 1.0 / float(iFrame + 1));
    }
    
    fragColor = vec4(col, 1.0f);
}

// Image Out
void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
    
    vec2 uv = fragCoord / iResolution.xy;
    vec3 col = texture(iChannel0, uv).rgb;
    
    // Colour conversions, because rgb is not true to human eyes
    // First we convert high dynamic range to a normalised value, since lights are unbounded 1->?
    // Then we convert our linear colours to SRGB
    float exposure = 0.5;
    col *= exposure;
    
    col = ACESFilm(col);
    col = LinearToSRGB(col);
    
    fragColor = vec4(col, 1.0);
}
