//Copyright 2010 JogAmp Community. All rights reserved.

#version 140
#extension GL_OES_standard_derivatives : enable

uniform vec3    ColorStatic;
uniform float   Alpha;
uniform float   Weight;
uniform sampler2D RBOTexture;
uniform vec2    TextureSize;
uniform int     secondPass;

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
    } else if(secondPass != 0) {
        vec4 t = texture(RBOTexture, rtex);
        
        c = vec3(1,1,1);//t.xyz;
        alpha = 1;//Alpha * t.w;    
    }

    fragColor = vec4(c, alpha);
}