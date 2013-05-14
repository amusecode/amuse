#version 150
#extension GL_EXT_geometry_shader4 : enable

layout(points) in;
//layout(triangle_strip, max_vertices = 256) out;
layout(points) out;

//uniform float radius;
//const float radius = 1.0;
//uniform int ndiv;

//uniform vec3 vdata[12];
//uniform int tindices[20][3];

//uniform mat4 PMatrix;
//uniform mat4 global_MVMatrix;

/*
struct Parameters {
	vec3 a;
	vec3 b;
	vec3 c;
	
	int div;	
};


void makeVertices(vec3 a, vec3 b, vec3 c, int div) {	
	struct Parameters stack[256];
	int index = 0;
	
	stack[0] = Parameters(a, b, c, div);
	index = 1;
	
	while (index > 0) {
		Parameters current = stack[index-1];
		index--;
		
		if (current.div <= 0) {
	        vec3 ra = current.a * radius;
	        vec3 rb = current.b * radius;
	        vec3 rc = current.c * radius;
	
			gl_Position = gl_PositionIn[0] + vec4(ra, 1);
			EmitVertex();
			gl_Position = gl_PositionIn[0] + vec4(rb, 1);
			EmitVertex();
			gl_Position = gl_PositionIn[0] + vec4(rc, 1);
			EmitVertex();
	
			EndPrimitive();
	    } else {
	        vec3 ab = normalize(vec3(current.a.x+current.b.x, current.a.y+current.b.y, current.a.z+current.b.z));
	        vec3 ac = normalize(vec3(current.a.x+current.c.x, current.a.y+current.c.y, current.a.z+current.c.z));
	        vec3 bc = normalize(vec3(current.b.x+current.c.x, current.b.y+current.c.y, current.b.z+current.c.z));
	        
	        stack[index] = Parameters(current.a,  ab, ac, current.div - 1);
	        index++;
	        
	        stack[index] = Parameters(current.b,  bc, ab, current.div - 1);
	        index++;
	        
	        stack[index] = Parameters(current.c,  ac, bc, current.div - 1);
	        index++;
	        
	        stack[index] = Parameters(ab, bc, ac, current.div - 1);
	        index++;
	    }
	}
}
*/

void main()
{
	//vec3 vdata = vec3[12](	vec3(-X, 0.0, Z), 	vec3(X, 0.0, Z), 	vec3(-X, 0.0, -Z), 	vec3(X, 0.0, -Z),
	//						vec3(0.0, Z, X), 	vec3(0.0, Z, -X), 	vec3(0.0, -Z, X), 	vec3(0.0, -Z, -X),
	//						vec3(Z, X, 0.0), 	vec3(-Z, X, 0.0), 	vec3(Z, -X, 0.0), 	vec3(-Z, -X, 0.0) );
								
	//int tindices[20][3] = 	int[20](int[3]( 1, 4, 0 ), 	int[3]( 4, 9, 0 ), 	int[3]( 4, 5, 9 ), 	int[3]( 8, 5, 4 ), 
	//								int[3]( 1, 8, 4 ), 	int[3]( 1, 10, 8 ),	int[3]( 10, 3, 8 ), int[3]( 8, 3, 5 ), 
	//								int[3]( 3, 2, 5 ), 	int[3]( 3, 7, 2 ), 	int[3]( 3, 10, 7 ), int[3]( 10, 6, 7 ), 
	//								int[3]( 6, 11, 7 ),	int[3]( 6, 0, 11 ), int[3]( 6, 1, 0 ), 	int[3]( 10, 1, 6 ), 
	//								int[3]( 11, 0, 9 ),	int[3]( 2, 11, 9 ), int[3]( 5, 2, 9 ), 	int[3]( 11, 2, 7 ) );
				
	//for (int i = 0; i < 20; i++) {
	//	makeVertices(vdata[tindices[i][0]], vdata[tindices[i][1]], vdata[tindices[i][2]], ndiv);
	//}
	
	gl_Position = gl_PositionIn[0]; //vec4(gl_PositionIn[0].x, 			gl_PositionIn[0].y, 		gl_PositionIn[0].z, 	1.0);	
	EmitVertex();

	EndPrimitive();

	//gl_Position = vec4(gl_PositionIn[0].x+radius, 	gl_PositionIn[0].y, 		gl_PositionIn[0].z, 		1.0);
	//EmitVertex();
	//gl_Position = vec4(gl_PositionIn[0].x, 			gl_PositionIn[0].y+radius, 	gl_PositionIn[0].z, 		1.0);
	//EmitVertex();
	//gl_Position = vec4(gl_PositionIn[0].x, 			gl_PositionIn[0].y, 		gl_PositionIn[0].z+radius, 	1.0);
	//EmitVertex();

	//EndPrimitive();
	
}
