#version 140

uniform vec4 Color;

out vec4 fragColor;

void main() {	
	fragColor = Color;
	
	if (Color.a > 1.0) {
		fragColor.rgb = Color.rgb * Color.a;
	}
}
