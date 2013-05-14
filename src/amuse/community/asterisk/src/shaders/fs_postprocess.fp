#version 140

uniform sampler2D axesTexture;
uniform sampler2D pointGasTexture;
uniform sampler2D octreeGasTexture;
uniform sampler2D sphereTexture;
uniform sampler2D starTexture;
uniform sampler2D starHaloTexture;
uniform sampler2D hudTexture;

uniform float axesBrightness;
uniform float pointGasBrightness;
uniform float octreeGasBrightness;
uniform float sphereBrightness;
uniform float starBrightness;
uniform float starHaloBrightness;
uniform float hudBrightness;
uniform float overallBrightness;

uniform int scrWidth;
uniform int scrHeight;

out vec4 fragColor;

void main() {
	vec2 tCoord   = vec2(gl_FragCoord.x/float(scrWidth), gl_FragCoord.y/float(scrHeight));
	
	vec4 starColor = vec4(texture(starTexture, tCoord).rgb, 1.0) * starBrightness * 0.1;
	vec4 starHaloColor = vec4(texture(starHaloTexture, tCoord).rgb, 1.0) * starHaloBrightness * 0.1;	
	
	vec4 combinedStarColor = mix(starColor, starHaloColor, 0.5); 
		
	vec4 axesColor = vec4(texture(axesTexture, tCoord).rgb, 1.0) * axesBrightness * 0.1;
  	vec4 pointGasColor  = vec4(texture(pointGasTexture, tCoord).rgb, 1.0) * pointGasBrightness * 0.1;	
  	vec4 octreeGasColor  = vec4(texture(octreeGasTexture, tCoord).rgb, 1.0) * octreeGasBrightness * 0.1;	
	vec4 sphereColor = vec4(texture(sphereTexture, tCoord).rgb, 1.0) * sphereBrightness * 0.1;
	vec4 hudColor = vec4(texture(hudTexture, tCoord).rgb, 1.0) * hudBrightness * 0.1;
    
    vec4 color = axesColor + pointGasColor + octreeGasColor + sphereColor + combinedStarColor + hudColor;
    
    //vec4 color = mix(starColor * starBrightness, starHaloColor * starHaloBrightness, 0.5); 
    //	 color = mix(color, sphereColor * sphereBrightness * 2.5, 0.1);
    //	 color = mix(color, pointGasColor * pointGasBrightness * 2.5, 0.1);
    //	 color = mix(color, octreeGasColor * octreeGasBrightness * 2.5, 0.1);
    //	 color = mix(color, axesColor * axesBrightness * 2.5, 0.1);
    //	 color = mix(color, hudColor * hudBrightness * 2.5, 0.1);
    
    fragColor = vec4(color.rgb * overallBrightness, 1.0);
}
