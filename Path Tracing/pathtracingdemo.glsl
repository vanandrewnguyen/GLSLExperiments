/*
Van Andrew Nguyen
22/4/22
[Path Tracing]

Dabbled into path-tracing, I wanted to try learn sample theory and other methods of rendering scenes.
*/

// Buffer A
#define SURFDIS 0.01
#define MAXDIS 100.0
#define MAXBOUNCES 8

#define PI 3.1416
#define FOV 90.0

// Credits to https://blog.demofox.org/2020/05/25/casual-shadertoy-path-tracing-1-basic-camera-diffuse-emissive/
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

vec3 generateRandomUnitVector(inout uint seed) {
    // Remember we are influence the seed for every ray bounce, hence it's an inout var
    float z = hash11(seed) * 2.0 - 1.0;
    float a = hash11(seed) * PI * 2.0;
    float r = sqrt(1.0 - z * z);
    float x = r * cos(a);
    float y = r * sin(a);
    return vec3(x, y, z);
}

// Keeps data of the ray in a struct
struct rayHitData {
    float dist;
    vec3 normal;
    vec3 emissive;
    vec3 albedo;
};

float scalarTriple(vec3 u, vec3 v, vec3 w) {
    // Find normal of u, v vectors and get mag on vector w
    return dot(cross(u, v), w);
}

// Intersection funct for sphere
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
        return true;
    }    
    
    return false;
}

// Intersection funct for quad
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
        data.dist = dist;        
        data.normal = normalize((rayPos+rayDir*dist) - sphere.xyz) * (fromInside ? -1.0f : 1.0f);
        return true;
    }
    
    return false;
}

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
        if (testQuadTrace(rayPos, rayDir, hitInfo, A, B, C, D))
        {
            hitInfo.albedo = vec3(0.1f, 0.7f, 0.1f);
            hitInfo.emissive = vec3(0.0f, 0.0f, 0.0f);
        }
    } 
    
    // Left wall
    {
        vec3 A = vec3(-wallLen, -wallLen, farPoint);
        vec3 B = vec3(-wallLen, -wallLen, nearPoint);
        vec3 C = vec3(-wallLen,  wallLen, nearPoint);
        vec3 D = vec3(-wallLen,  wallLen, farPoint);
        if (testQuadTrace(rayPos, rayDir, hitInfo, A, B, C, D))
        {
            hitInfo.albedo = vec3(0.7f, 0.1f, 0.1f);
            hitInfo.emissive = vec3(0.0f, 0.0f, 0.0f);
        }
    } 
    
    // Top wall
    {
        vec3 A = vec3(-wallLen, wallLen, farPoint);
        vec3 B = vec3( wallLen, wallLen, farPoint);
        vec3 C = vec3( wallLen, wallLen, nearPoint);
        vec3 D = vec3(-wallLen, wallLen, nearPoint);
        if (testQuadTrace(rayPos, rayDir, hitInfo, A, B, C, D))
        {
            hitInfo.albedo = vec3(0.7f, 0.7f, 0.7f);
            hitInfo.emissive = vec3(0.0f, 0.0f, 0.0f);
        }
    } 
    
    // Bottom wall
    {
        vec3 A = vec3(-wallLen, -wallLen, farPoint);
        vec3 B = vec3( wallLen, -wallLen, farPoint);
        vec3 C = vec3( wallLen, -wallLen, nearPoint);
        vec3 D = vec3(-wallLen, -wallLen, nearPoint);
        if (testQuadTrace(rayPos, rayDir, hitInfo, A, B, C, D))
        {
            hitInfo.albedo = vec3(0.7f, 0.7f, 0.7f);
            hitInfo.emissive = vec3(0.0f, 0.0f, 0.0f);
        }
    } 
    
    // Back wall
    {
        vec3 A = vec3(-wallLen, -wallLen, farPoint);
        vec3 B = vec3( wallLen, -wallLen, farPoint);
        vec3 C = vec3( wallLen,  wallLen, farPoint);
        vec3 D = vec3(-wallLen,  wallLen, farPoint);
        if (testQuadTrace(rayPos, rayDir, hitInfo, A, B, C, D))
        {
            hitInfo.albedo = vec3(0.7f, 0.7f, 0.7f);
            hitInfo.emissive = vec3(0.0f, 0.0f, 0.0f);
        }
    } 
    
    // Lighting
    float pad = 2.0;
    {
        vec3 A = vec3(-wallLen + pad, wallLen - 0.1, farPoint - pad);
        vec3 B = vec3( wallLen - pad, wallLen - 0.1, farPoint - pad);
        vec3 C = vec3( wallLen - pad, wallLen - 0.1, nearPoint + pad);
        vec3 D = vec3(-wallLen + pad, wallLen - 0.1, nearPoint + pad);
        if (testQuadTrace(rayPos, rayDir, hitInfo, A, B, C, D))
        {
            hitInfo.albedo = vec3(0.7f, 0.7f, 0.7f);
            hitInfo.emissive = vec3(1.0f, 0.9f, 0.7f) * 20.0;
        }
    } 

    // Balls
    if (testSphereTrace(rayPos, rayDir, hitInfo, vec4(-5.0f, -5.0f, 20.0f, 2.0f))) {
        hitInfo.albedo = vec3(1.0f, 0.8f, 0.8f);
        hitInfo.emissive = vec3(0.0f, 0.0f, 0.0f);        
    } 
     
    if (testSphereTrace(rayPos, rayDir, hitInfo, vec4(0.0f, -5.0f, 20.0f, 2.0f))) {
        hitInfo.albedo = vec3(0.8f, 1.0f, 0.8f);
        hitInfo.emissive = vec3(0.0f, 0.0f, 0.0f);        
    }    
     
    if (testSphereTrace(rayPos, rayDir, hitInfo, vec4(5.0f, -5.0f, 20.0f, 2.0f))) {
        hitInfo.albedo = vec3(0.8f, 0.8f, 1.0f);
        hitInfo.emissive = vec3(0.0f, 0.0f, 0.0f);
    }
     
    // Light
    {
        vec3 A = vec3(-wallLen,  wallLen, farPoint);
        vec3 B = vec3( wallLen,  wallLen, farPoint);
        vec3 C = vec3( wallLen,  wallLen, nearPoint);
        vec3 D = vec3(-wallLen,  wallLen, nearPoint);
        if (testQuadTrace(rayPos, rayDir, hitInfo, A, B, C, D))
        {
            hitInfo.albedo = vec3(0.7f, 0.7f, 0.7f);
            hitInfo.emissive = vec3(1.0f, 0.9f, 0.7f) * 20.0;
        }
    }    
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
        hitInfo.dist = MAXDIS;
        scene(rayPos, rayDir, hitInfo);
         
        // If the ray hits nothing, exit
        if (hitInfo.dist == MAXDIS) {
            break;
        }
         
        // Update ray position and update its direction along the normal
        rayPos = (rayPos + rayDir * hitInfo.dist) + hitInfo.normal * SURFDIS;
        rayDir = normalize(hitInfo.normal + generateRandomUnitVector(randomSeed));        
         
        // Add in emissive colour of the hit, then update the bounced colour
        col += hitInfo.emissive * throughput;
        throughput *= hitInfo.albedo;      
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
    vec3 rayTarget = vec3(uv * 2.0f - 1.0f, camDis);
    
    // Changing aspect ratio
    float ratio = iResolution.x / iResolution.y;
    rayTarget.y /= ratio;
    vec3 rayDir = normalize(rayTarget - rayPos);
    
    // We need to average the pixel values because we're sampling a bunch of rays that bounce randomly 
    vec3 col = getRayColour(rayPos, rayDir, randomSeed);
    vec3 prevCol = texture(iChannel0, uv).rgb;
    col = mix(prevCol, col, 1.0 / float(iFrame + 1));
    
    fragColor = vec4(col, 1.0f);
}

// Image Out
void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
    vec2 uv = fragCoord / iResolution.xy;
    vec3 col = texture(iChannel0, uv).rgb;
    fragColor = vec4(col, 1.0);
}
