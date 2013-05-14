#version 140

uniform vec4 Color;

in vec3 vertex_normal;
in vec3 eye_direction;

out vec4 fragColor;

void main() {
	vec3 matColor = vec3(Color.r, Color.g, Color.b);	
	
	float dotP = dot(vertex_normal, vec3(0,0,1));	

	//fragColor = vec4(matColor,0.5-(dotP/2.0));
	
	fragColor = vec4(matColor, dotP*dotP*dotP);
}
