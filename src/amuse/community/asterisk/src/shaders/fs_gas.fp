#version 140

//uniform float gas_opacity_factor;
//uniform vec4 node_color;

in vec4 pColor;

out vec4 fragColor;

void main() {
	//fragColor = node_color * gas_opacity_factor;
	
	fragColor = pColor; 
	
	//vec4(1.0, 1.0, 1.0, 1.0);
}
