#version 140
#extension GL_EXT_geometry_shader4 : enable

layout(triangles) in;
layout(triangle_strip, max_vertices = 3) out;

in vec3 EyespaceNormal[3];
out vec3 EyespaceNormalOut;
 
// a passthrough geometry shader for color and position
void main()
{
	vec4 normal;
	
	normal = vec4(normalize(cross(gl_PositionIn[2].xyz - gl_PositionIn[0].xyz, gl_PositionIn[1].xyz - gl_PositionIn[0].xyz)), 1.0);
    gl_Position = gl_PositionIn[0] + normal;
    EyespaceNormalOut = EyespaceNormal[0];
    EmitVertex();   
    
    normal = vec4(normalize(cross(gl_PositionIn[0].xyz - gl_PositionIn[1].xyz, gl_PositionIn[2].xyz - gl_PositionIn[1].xyz)), 1.0);
    gl_Position = gl_PositionIn[1] + normal;
    EyespaceNormalOut = EyespaceNormal[1];
    EmitVertex();   
    
    normal = vec4(normalize(cross(gl_PositionIn[1].xyz - gl_PositionIn[2].xyz, gl_PositionIn[0].xyz - gl_PositionIn[2].xyz)), 1.0);
    gl_Position = gl_PositionIn[2] + normal;
    EyespaceNormalOut = EyespaceNormal[2];
    EmitVertex();
    
    EndPrimitive();
    
    //vec3 normal0 = gl_PositionIn[0].xyz + 2 * -cross(gl_PositionIn[1].xyz - gl_PositionIn[0].xyz, gl_PositionIn[2].xyz - gl_PositionIn[0].xyz);
    //vec3 normal1 = gl_PositionIn[1].xyz + 2 * -cross(gl_PositionIn[2].xyz - gl_PositionIn[1].xyz, gl_PositionIn[0].xyz - gl_PositionIn[1].xyz);
    //vec3 normal2 = gl_PositionIn[2].xyz + 2 * -cross(gl_PositionIn[0].xyz - gl_PositionIn[2].xyz, gl_PositionIn[1].xyz - gl_PositionIn[2].xyz);
    
    //float tipX = (normal0.x + normal1.x + normal2.x) / 3.0;  
    //float tipY = (normal0.y + normal1.y + normal2.y) / 3.0;
    //float tipZ = (normal0.z + normal1.z + normal2.z) / 3.0;
    
    //vec4 tip = vec4(tipX, tipY, tipZ, 1.0);
    
    //gl_Position = gl_PositionIn[0];
    //EmitVertex(); 
    //gl_Position = gl_PositionIn[1];
    //EmitVertex(); 
    //gl_Position = tip;
    //EmitVertex(); 
    
    //EndPrimitive();
    
    //gl_Position = gl_PositionIn[1];
    //EmitVertex(); 
    //gl_Position = gl_PositionIn[2];
    //EmitVertex(); 
    //gl_Position = tip;
    //EmitVertex(); 
    
    //EndPrimitive();
    
    //gl_Position = gl_PositionIn[2];
    //EmitVertex(); 
    //gl_Position = gl_PositionIn[0];
    //EmitVertex(); 
    //gl_Position = tip;
    //EmitVertex(); 
    
    //EndPrimitive();
  //for(int i = 0; i < gl_VerticesIn; ++i)
  //{
    // copy color
    //gl_FrontColor = gl_FrontColorIn[i];
 
    // copy position
    //gl_Position = gl_PositionIn[i];
    //EyespaceNormalOut = EyespaceNormal[i];
 
    // done with the vertex
    //EmitVertex();
  //}
  
  //EndPrimitive();
}
