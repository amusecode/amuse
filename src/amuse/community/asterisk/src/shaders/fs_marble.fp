#version 140

in vec3  MCposition;

uniform sampler3D Noise;
uniform vec3 Color2;
uniform float Offset;

const vec3 Color1 = vec3(0.8, 0.7, 0.0);
const float NoiseScale = 1.2;
const float LightIntensity = 1.0; 

float undulate(float x) {
    if(x < -0.4) {
        return 0.15 + 2.857 * (x + 0.75)*(x + 0.75);
    } else if(x < 0.4) {
        return 0.95 - 2.8125 * (x*x);
    } else {
        return 0.26 + 2.666 * (x - 0.7)*(x - 0.7);
    }
}
 
float marble(vec3 p) {
    float PI = 3.1415;
    float s = 0.25;
    float a = 2.0;
    
    // this computes 3 octaves of turbulence
    vec4 t1 = 1.0 *texture(Noise, 0.25*p);
    vec4 t2 = 0.5 *texture(Noise, 0.50*p);
    vec4 t3 = 0.25*texture(Noise, 1.00*p);
 
    float m = undulate(sin(p.y*2.0*PI + a*(t1.x + t2.x + t3.x)));
   
    return m;
}
 
void main() {   
    float intensity = marble(MCposition);
    gl_FragColor = vec4(mix(Color1, Color2, intensity) * LightIntensity, 1.0);
}