#version 140

uniform sampler2D texture_map;

in vec2 tCoord;

void main() { 		
	vec4 color = vec4(texture(texture_map, tCoord).rgb, 1.0);
    gl_FragColor = vec4(color.rgb, 1.0);
} 
