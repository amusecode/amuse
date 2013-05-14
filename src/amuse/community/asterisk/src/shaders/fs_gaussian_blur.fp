#version 140

uniform sampler2D Texture;

uniform int scrWidth;
uniform int scrHeight;

uniform int blurDirection;
uniform float blurSize;
uniform int blurType;

uniform float Sigma;
uniform float NumPixelsPerSide;
uniform float Alpha;

const float pi = 3.14159265;

const vec2 vertical = vec2(0.0, 1.0);
const vec2 horizontal = vec2(1.0, 0.0);

out vec4 fragColor;

// The sigma value for the gaussian function: higher value means more blur
// A good value for 9x9 is around 3 to 5
// A good value for 7x7 is around 2.5 to 4
// A good value for 5x5 is around 2 to 3.5
// ... play around with this based on what you need :)

// blurSize should usually be equal to
// 1.0f / texture_pixel_width for a horizontal blur, and
// 1.0f / texture_pixel_height for a vertical blur.

vec4 gaussianBlur(sampler2D tex, vec2 tCoord, vec2 multiplyVec, int maxTexSize, float blurSize, float numPixelsPerSide, float sigma) {
	// Incremental Gaussian Coefficent Calculation (See GPU Gems 3 pp. 877 - 889)
	vec3 incrementalGaussian;
  	incrementalGaussian.x = 1.0 / (sqrt(2.0 * pi) * sigma);
  	incrementalGaussian.y = exp(-0.5 / (sigma * sigma));
  	incrementalGaussian.z = incrementalGaussian.y * incrementalGaussian.y;

  	vec4 avgValue = vec4(0.0, 0.0, 0.0, 0.0);
  	float coefficientSum = 0.0;

  	// Take the central sample first...
  	avgValue += texture(tex, tCoord) * incrementalGaussian.x;
  	coefficientSum += incrementalGaussian.x;
  	incrementalGaussian.xy *= incrementalGaussian.yz;

  	// Go through the remaining 8 vertical samples (4 on each side of the center)
  	
  	for (float i = 1.0; i <= numPixelsPerSide; i++) {
  		vec2 offset = vec2((i * blurSize * multiplyVec/float(maxTexSize)));
  		if (
  		tCoord.x - offset.x < 0.0 || tCoord.x + offset.x > 1.0 ||
  			tCoord.y - offset.y < 0.0 || tCoord.y + offset.y > 1.0) {
  			avgValue += 2 * texture(tex, tCoord) * incrementalGaussian.x;
  		} else {
  			avgValue += texture(tex, tCoord - offset) * incrementalGaussian.x;  		         
	    	avgValue += texture(tex, tCoord + offset) * incrementalGaussian.x;
	    }
	             
	    coefficientSum += 2.0 * incrementalGaussian.x;
	    incrementalGaussian.xy *= incrementalGaussian.yz;
  	}
  	
  	return vec4(avgValue / coefficientSum);
}

void main() {
	vec2 tCoord   = vec2(gl_FragCoord.x/float(scrWidth), gl_FragCoord.y/float(scrHeight));
	
	float sigma;
	float numPixelsPerSide;
	
	if (blurType == 0) {
		sigma = Sigma;
		numPixelsPerSide = NumPixelsPerSide;
	} else if (blurType == 1) {
		sigma = 1.0;
		numPixelsPerSide = 1.0; // 3x3
	} else if (blurType == 2) {
		sigma = 2.0;
		numPixelsPerSide = 2.0; // 5x5
	} else if(blurType == 3) {
		sigma = 2.5;
		numPixelsPerSide = 3.0; // 7x7
	} else if(blurType == 4) {
		sigma = 3.0;
		numPixelsPerSide = 4.0; // 9x9
	} else if(blurType == 6) {
		sigma = 4.0;
		numPixelsPerSide = 6.0; // 13x13
	} else if(blurType == 8) {
		sigma = 6.0;
		numPixelsPerSide = 8.0; // 17x17
	} else {
		sigma = 7.0;
		numPixelsPerSide = 10.0; // 21x21
	}
	
	vec2 direction;
	if (blurDirection == 0) {
		direction = horizontal;
		fragColor = vec4(gaussianBlur(Texture, tCoord, direction, scrWidth, blurSize, numPixelsPerSide, sigma).rgb, Alpha);
	} else {
		direction = vertical;
		fragColor = vec4(gaussianBlur(Texture, tCoord, direction, scrHeight, blurSize, numPixelsPerSide, sigma).rgb, Alpha);
	}
	
  	
}
