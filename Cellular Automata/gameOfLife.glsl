

// MAIN
/*
We use buffers to check the previous state, since in shaders you only calculate each pixel
per frame instead of frame by frame. Has no record of previous colour.
*/

#define ZOOM 2.0
#define CYCLE 100.0

const vec3 bgCol = vec3(0.058, 0.082, 0.066);
const vec3 col1 = vec3(0.152, 0.647, 0.427);
const vec3 col2 = vec3(0.215, 0.345, 0.627);
const vec3 col3 = vec3(0.588, 0.317, 0.211);
const vec3 col4 = vec3(0.772, 0.203, 0.254);

// We need to make a colour gradient which cycles between 0.0->1.0
vec3 getGradient(float value) {
    // Value should be the cycle of the state / CYCLE (between 0 -> 1)
    // Seperate into four bands for four colours
    if (value < 0.25) {
        return mix(bgCol, col1, value * 4.0);
    } else if (value < 0.5) {
        return mix(col1, col2, (value - 0.25) * 4.0);
    } else if (value < 0.75) {
        return mix(col2, col3, (value - 0.50) * 4.0);
    } else {
        return mix(col3, col4, (value - 0.75) * 4.0);
    }
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord/iResolution.xy / ZOOM;
    
    // Read from buffer (previous state)
    vec4 state = texture(iChannel0, uv, 0.0);
    
    vec3 col = getGradient(state.y / 100.0);

    if (state.x > 0.0) {
	    fragColor = vec4(col, 1.0);
    } else {
        // Reduce colour intensity if it's not alive
        const float blackMix = 0.5;
        fragColor = vec4(mix(col, bgCol, blackMix), 1.0);
    }
}

// BUFFER A
#define THRES 0.8
#define ZOOM 2.0
#define BRUSHSIZE 1.0
#define CYCLE 100.0

int getNeighbours(ivec2 pos) {
    int num = 0;
    
    for (int y = -1; y <= 1; y++) {
        for (int x = -1; x <= 1; x++) {
            if (x == 0 && y == 0) continue;
            num += texelFetch(iChannel1, pos + ivec2(x, y), 0).r > THRES ? 1 : 0;
        }
    }
    
    return num;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord / iResolution.xy;
    vec4 col = vec4(0); 
    
    // Init using frame data
    if (iFrame < 1) {
        col = texture(iChannel0, uv);
    } else {
        // Game of life rules //
        
        // Here we're pulling from the previous frame, since iChannel writes to buffer A and reads from it
        vec4 previousState = texelFetch(iChannel1, ivec2(fragCoord), 0);
        float isAlive = previousState.x;//texelFetch(iChannel1, ivec2(fragCoord), 0).r > THRES;
        int numNeighbours = getNeighbours(ivec2(fragCoord));
        
        // Apply the three rules
        if (isAlive > 0.0 && (numNeighbours == 2 || numNeighbours == 3)) {
            isAlive = 1.0;
        } else if (isAlive <= 0.0 && numNeighbours == 3) {
            isAlive = 1.0;
        } else {
            isAlive = 0.0;
        }
        
        // Mouse control
        if (distance(fragCoord, iMouse.xy / ZOOM) < BRUSHSIZE) {
            isAlive = 1.0;
        }
        
        // Create a new state
        vec4 newState = previousState;
        newState.x = isAlive;
        
        if (isAlive > 0.0) {
            // Cycle is how long it takes to cycle between colours
            newState.y = mod(newState.y + isAlive, CYCLE);
        } else {
            // You'd want the board to decay if there are no alive cells (no slime trail)
            // So newstate.y will eventually reach 0.0 if there are no populating cells
            const float decayAmount = 0.05;
            newState.y = clamp(newState.y - decayAmount, 0.0, CYCLE);
        }
        
        col = newState;
    }
    
    fragColor = col;
}
