/*
Van Andrew Nguyen
13/12/21
[Falling Sand Simulation]

I have been playing a lot of Noita and decided to write my own falling sand simulation. Because instead of using a big array I am writing in a pixel shader, I can't keep 
track of previous states or position, so I need to use buffers to get a pixel's previous state, which is a bit finicky.
*/

// MAIN
/*
We use buffers to check the previous state, since in shaders you only calculate each pixel
per frame instead of frame by frame. Has no record of previous colour.
*/

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord/iResolution.xy;
    
    // Read from buffer (previous state)
    vec3 col = texture(iChannel0, uv).rgb;

    fragColor = vec4(col,1.0);
}

// BUFFER A
#define THRES 0.6
#define BRUSHSIZE 16.0

const vec3 sandCol = vec3(0.917, 0.776, 0.380);
const vec3 bgCol = vec3(0.067, 0.051, 0.104);

vec3 background(vec2 uv) {
    vec3 col = bgCol;
    float y = uv.y * 0.8 + 0.2;
    col += y * vec3(0.2, 0, 0);
    
    return col;
}

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

bool isEmptyPosition(ivec2 pos) {
    return texelFetch(iChannel1, pos, 0).r > THRES ? false : true;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord/iResolution.xy;
    vec4 col = vec4(0); 
    
    // Init using frame data
    if (iFrame < 1) {
        col = texture(iChannel0, uv);
    } else {
        // Falling Sand Rules // 
        
        // Here we're pulling from the previous frame, since iChannel writes to buffer A and reads from it
        bool isAlive = texelFetch(iChannel1, ivec2(fragCoord), 0).r > THRES;
        int numNeighbours = getNeighbours(ivec2(fragCoord));
        
        /* 
        Apply rules
        If there is a space below, fall down
        If there is a space below and left, fall left
        If there is a space below and right, fall right
        */
        vec3 currCol = background(uv);
        if (isAlive) {
            // Check if pixel below is not empty, then we stay alive
            ivec2 pos = ivec2(fragCoord) + ivec2(0, -1);
            if (!isEmptyPosition(pos)) {
                currCol = sandCol;
            } 
            // Now check if there are empty pixels to the bottom left/right
            if (isEmptyPosition(pos + ivec2(-1, 0)) || isEmptyPosition(pos + ivec2(1, 0))) {
                currCol = background(uv);
            }
            // Now finally check if we are at the bottom, if so stay alive
            // Same with walls
            if (pos.y < 1 || 
                (pos + ivec2(-1, 1)).x < 0 || 
                (pos + ivec2(1, 1)).x > int(iResolution.x) - 1) {
                currCol = sandCol;
            }
        } else {
            // If there is a pixel above us, we are alive in the next iteration
            if (!isEmptyPosition(ivec2(fragCoord) + ivec2(0, 1))) {
                currCol = sandCol;
            } 
            // If there is sand below us (we shouldn't be falling)
            // and there are gaps on either side of us, we need to move off and fall to the side.
            if (!isEmptyPosition(ivec2(fragCoord) + ivec2(0, -1))) {
                if (!isEmptyPosition(ivec2(fragCoord) + ivec2(-1, 1)) ||
                    !isEmptyPosition(ivec2(fragCoord) + ivec2(1, 1))) {
                    currCol = sandCol;
                }
            }
        }
        
        // Mouse control
        if (iMouse.z > 0.0 && length(iMouse.xy - fragCoord.xy) < BRUSHSIZE) {
            currCol = sandCol;
        }
        
        // Edit sand colour
        if (currCol == sandCol) {
            float height = 0.0;
            while (texelFetch(iChannel1, ivec2(fragCoord) + ivec2(0, height), 0).r > THRES) {
                height += 1.0;
            }
            currCol += height * 0.005;
        }
        
        col = vec4(currCol, 1);
    }
    
    fragColor = col;
}
