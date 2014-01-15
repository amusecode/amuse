#version 140

uniform sampler2D Texture;

out vec4 fragColor;

uniform int blurDirection;

uniform int scrWidth;
uniform int scrHeight;

uniform float offset[3] = float[]( 0.0, 1.3846153846, 3.2307692308 );
uniform float weight[3] = float[]( 0.2270270270, 0.3162162162, 0.0702702703 );

void main(void)
{
	fragColor = texture2D( Texture, vec2(gl_FragCoord.x/scrWidth, gl_FragCoord.y/scrHeight)) * weight[0];

	vec2 direction;
	if (blurDirection == 0) {		
	    for (int i=1; i<3; i++) {
	        fragColor += texture2D( Texture, vec2((gl_FragCoord.x + offset[i])/float(scrWidth), gl_FragCoord.y/float(scrHeight))) * weight[i];	        
	        fragColor += texture2D( Texture, vec2((gl_FragCoord.x + offset[i])/float(scrWidth), gl_FragCoord.y/float(scrHeight))) * weight[i];
	    }	
    } else {
	    for (int i=1; i<3; i++) {
	        fragColor += texture2D( Texture, vec2(gl_FragCoord.x/float(scrWidth), (gl_FragCoord.y + offset[i])/float(scrHeight))) * weight[i];
	        fragColor += texture2D( Texture, vec2(gl_FragCoord.x/float(scrWidth), (gl_FragCoord.y - offset[i])/float(scrHeight))) * weight[i];
	    }		
	}
}
