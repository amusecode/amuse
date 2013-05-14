//Copyright 2010 JogAmp Community. All rights reserved.

#version 140
#extension GL_OES_standard_derivatives : enable

uniform vec3    ColorStatic;
uniform float   Alpha;
uniform float   Weight;
uniform sampler2D RBOTexture;
uniform vec2    TextureSize;

in vec2    tCoordsFS;

out vec4 fragColor;

//
// 2-pass shader w/o weight
//

const vec4 tex_weights = vec4(0.075, 0.06, 0.045, 0.025);

void main (void)
{    
    vec2 rtex = vec2(abs(tCoordsFS.x),abs(tCoordsFS.y));
    vec3 c = ColorStatic.rgb;

    float alpha = 0.0;
    float enable = 1.0;
    
    if((tCoordsFS.x == 0.0) && (tCoordsFS.y == 0.0)) {
         alpha = Alpha;
    }
    else if((tCoordsFS.x >= 5.0)) {
        vec2 dfx = dFdx(tCoordsFS);
        vec2 dfy = dFdy(tCoordsFS);
        
        vec2 size = 1.0/TextureSize;

        rtex -= 5.0;
        vec4 t = texture(RBOTexture, rtex)* 0.18;

        t += texture(RBOTexture, rtex + size*(vec2(1, 0)))*tex_weights.x;
        t += texture(RBOTexture, rtex - size*(vec2(1, 0)))*tex_weights.x;
        t += texture(RBOTexture, rtex + size*(vec2(0, 1)))*tex_weights.x;
        t += texture(RBOTexture, rtex - size*(vec2(0, 1)))*tex_weights.x;
        
        t += texture(RBOTexture, rtex + 2.0*size*(vec2(1, 0)))*tex_weights.y;
        t += texture(RBOTexture, rtex - 2.0*size*(vec2(1, 0)))*tex_weights.y;
        t += texture(RBOTexture, rtex + 2.0*size*(vec2(0, 1)))*tex_weights.y; 
        t += texture(RBOTexture, rtex - 2.0*size*(vec2(0, 1)))*tex_weights.y;
        
        t += texture(RBOTexture, rtex + 3.0*size*(vec2(1, 0)))*tex_weights.z;
        t += texture(RBOTexture, rtex - 3.0*size*(vec2(1, 0)))*tex_weights.z;
        t += texture(RBOTexture, rtex + 3.0*size*(vec2(0, 1)))*tex_weights.z;
        t += texture(RBOTexture, rtex - 3.0*size*(vec2(0, 1)))*tex_weights.z;
        
        t += texture(RBOTexture, rtex + 4.0*size*(vec2(1, 0)))*tex_weights.w;
        t += texture(RBOTexture, rtex - 4.0*size*(vec2(1, 0)))*tex_weights.w;
        t += texture(RBOTexture, rtex + 4.0*size*(vec2(0, 1)))*tex_weights.w;
        t += texture(RBOTexture, rtex - 4.0*size*(vec2(0, 1)))*tex_weights.w;
        
        /** discard freezes NV tegra2 compiler
        if(t.w == 0.0){
            discard;
        } */
        
        c = t.xyz;
        alpha = Alpha * t.w;
    }
    else if ((tCoordsFS.x > 0.0) && (rtex.y > 0.0 || rtex.x == 1.0)) {
        rtex.y -= 0.1;
          
        if(rtex.y < 0.0 && tCoordsFS.y < 0.0) {
            // discard; // freezes NV tegra2 compiler
            alpha = 0.0;
        } else {
            rtex.y = max(rtex.y, 0.0);

            vec2 dtx = dFdx(rtex);
            vec2 dty = dFdy(rtex);
              
            vec2 f = vec2((dtx.y - dtx.x + 2.0*rtex.x*dtx.x), (dty.y - dty.x + 2.0*rtex.x*dty.x));
            float position = rtex.y - (rtex.x * (1.0 - rtex.x));

            // FIXME: will we ever set gcu_Alpha != 1.0 ? If not, a==alpha!
            float a = clamp(0.5 - ( position/length(f) ) * sign(tCoordsFS.y), 0.0, 1.0);
            alpha = Alpha * a;
        }
    }
    fragColor = vec4(c, alpha);
    //fragColor = vec4(1,1,1,1);
}