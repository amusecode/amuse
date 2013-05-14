#version 140

uniform sampler2D Texture;

uniform int scrWidth;
uniform int scrHeight;

out vec4 fragColor;

void main(void) {
  vec2 tCoord   = vec2(gl_FragCoord.x/float(scrWidth), gl_FragCoord.y/float(scrHeight));

  vec4 fragmentColor = vec4(0.0, 0.0, 0.0, 0.0);
  vec2 blurSize = vec2(0.002, 0.002);

  fragmentColor += texture2D(Texture, tCoord - 4.0 * blurSize) * 0.05;
  fragmentColor += texture2D(Texture, tCoord - 3.0 * blurSize) * 0.09;
  fragmentColor += texture2D(Texture, tCoord - 2.0 * blurSize) * 0.12;
  fragmentColor += texture2D(Texture, tCoord - 1.0 * blurSize) * 0.15;
  fragmentColor += texture2D(Texture, tCoord                 ) * 0.16;
  fragmentColor += texture2D(Texture, tCoord + 1.0 * blurSize) * 0.15;
  fragmentColor += texture2D(Texture, tCoord + 2.0 * blurSize) * 0.12;
  fragmentColor += texture2D(Texture, tCoord + 3.0 * blurSize) * 0.09;
  fragmentColor += texture2D(Texture, tCoord + 4.0 * blurSize) * 0.05;
  
  fragColor = vec4(fragmentColor.rgb, 1.0);
}
