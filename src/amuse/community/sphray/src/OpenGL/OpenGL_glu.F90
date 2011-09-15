MODULE OpenGL_GLU
USE, INTRINSIC :: ISO_C_BINDING
USE OpenGL_kinds
IMPLICIT NONE
PRIVATE
!  Boolean values
INTEGER(GLboolean), PARAMETER, PUBLIC :: GLU_FALSE                = 0 ! 0
INTEGER(GLboolean), PARAMETER, PUBLIC :: GLU_TRUE                 = 1 ! 1
!  Version values
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_VERSION_1_1          = 1 ! 1
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_VERSION_1_2          = 1 ! 1
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_VERSION_1_3          = 1 ! 1
!  StringName values
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_VERSION              = 100800 ! 100800
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_EXTENSIONS           = 100801 ! 100801
!  ErrorCode values
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_INVALID_ENUM         = 100900 ! 100900
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_INVALID_VALUE        = 100901 ! 100901
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_OUT_OF_MEMORY        = 100902 ! 100902
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_INVALID_OPERATION    = 100904 ! 100904
!  NurbsDisplay values
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_OUTLINE_POLYGON      = 100240 ! 100240
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_OUTLINE_PATCH        = 100241 ! 100241
!  NurbsCallback values
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR          = 100103 ! 100103
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_ERROR                = 100103 ! 100103
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_BEGIN          = 100164 ! 100164
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_BEGIN_EXT      = 100164 ! 100164
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_VERTEX         = 100165 ! 100165
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_VERTEX_EXT     = 100165 ! 100165
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_NORMAL         = 100166 ! 100166
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_NORMAL_EXT     = 100166 ! 100166
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_COLOR          = 100167 ! 100167
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_COLOR_EXT      = 100167 ! 100167
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_TEXTURE_COORD  = 100168 ! 100168
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_TEX_COORD_EXT  = 100168 ! 100168
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_END            = 100169 ! 100169
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_END_EXT        = 100169 ! 100169
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_BEGIN_DATA     = 100170 ! 100170
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_BEGIN_DATA_EXT = 100170 ! 100170
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_VERTEX_DATA    = 100171 ! 100171
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_VERTEX_DATA_EXT = 100171 ! 100171
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_NORMAL_DATA    = 100172 ! 100172
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_NORMAL_DATA_EXT = 100172 ! 100172
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_COLOR_DATA     = 100173 ! 100173
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_COLOR_DATA_EXT = 100173 ! 100173
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_TEXTURE_COORD_DATA = 100174 ! 100174
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_TEX_COORD_DATA_EXT = 100174 ! 100174
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_END_DATA       = 100175 ! 100175
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_END_DATA_EXT   = 100175 ! 100175
!  NurbsError values
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR1         = 100251 ! 100251
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR2         = 100252 ! 100252
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR3         = 100253 ! 100253
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR4         = 100254 ! 100254
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR5         = 100255 ! 100255
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR6         = 100256 ! 100256
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR7         = 100257 ! 100257
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR8         = 100258 ! 100258
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR9         = 100259 ! 100259
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR10        = 100260 ! 100260
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR11        = 100261 ! 100261
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR12        = 100262 ! 100262
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR13        = 100263 ! 100263
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR14        = 100264 ! 100264
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR15        = 100265 ! 100265
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR16        = 100266 ! 100266
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR17        = 100267 ! 100267
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR18        = 100268 ! 100268
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR19        = 100269 ! 100269
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR20        = 100270 ! 100270
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR21        = 100271 ! 100271
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR22        = 100272 ! 100272
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR23        = 100273 ! 100273
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR24        = 100274 ! 100274
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR25        = 100275 ! 100275
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR26        = 100276 ! 100276
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR27        = 100277 ! 100277
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR28        = 100278 ! 100278
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR29        = 100279 ! 100279
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR30        = 100280 ! 100280
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR31        = 100281 ! 100281
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR32        = 100282 ! 100282
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR33        = 100283 ! 100283
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR34        = 100284 ! 100284
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR35        = 100285 ! 100285
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR36        = 100286 ! 100286
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_ERROR37        = 100287 ! 100287
!  NurbsProperty values
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_AUTO_LOAD_MATRIX     = 100200 ! 100200
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_CULLING              = 100201 ! 100201
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_SAMPLING_TOLERANCE   = 100203 ! 100203
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_DISPLAY_MODE         = 100204 ! 100204
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_PARAMETRIC_TOLERANCE = 100202 ! 100202
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_SAMPLING_METHOD      = 100205 ! 100205
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_U_STEP               = 100206 ! 100206
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_V_STEP               = 100207 ! 100207
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_MODE           = 100160 ! 100160
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_MODE_EXT       = 100160 ! 100160
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_TESSELLATOR    = 100161 ! 100161
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_TESSELLATOR_EXT = 100161 ! 100161
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_RENDERER       = 100162 ! 100162
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NURBS_RENDERER_EXT   = 100162 ! 100162
!  NurbsSampling values
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_OBJECT_PARAMETRIC_ERROR = 100208 ! 100208
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_OBJECT_PARAMETRIC_ERROR_EXT = 100208 ! 100208
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_OBJECT_PATH_LENGTH   = 100209 ! 100209
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_OBJECT_PATH_LENGTH_EXT = 100209 ! 100209
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_PATH_LENGTH          = 100215 ! 100215
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_PARAMETRIC_ERROR     = 100216 ! 100216
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_DOMAIN_DISTANCE      = 100217 ! 100217
!  NurbsTrim values
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_MAP1_TRIM_2          = 100210 ! 100210
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_MAP1_TRIM_3          = 100211 ! 100211
!  QuadricDrawStyle values
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_POINT                = 100010 ! 100010
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_LINE                 = 100011 ! 100011
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_FILL                 = 100012 ! 100012
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_SILHOUETTE           = 100013 ! 100013
!  QuadricCallback values
!  QuadricNormal values
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_SMOOTH               = 100000 ! 100000
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_FLAT                 = 100001 ! 100001
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_NONE                 = 100002 ! 100002
!  QuadricOrientation values
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_OUTSIDE              = 100020 ! 100020
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_INSIDE               = 100021 ! 100021
!  TessCallback values
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_BEGIN           = 100100 ! 100100
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_BEGIN                = 100100 ! 100100
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_VERTEX          = 100101 ! 100101
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_VERTEX               = 100101 ! 100101
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_END             = 100102 ! 100102
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_END                  = 100102 ! 100102
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_ERROR           = 100103 ! 100103
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_EDGE_FLAG       = 100104 ! 100104
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_EDGE_FLAG            = 100104 ! 100104
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_COMBINE         = 100105 ! 100105
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_BEGIN_DATA      = 100106 ! 100106
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_VERTEX_DATA     = 100107 ! 100107
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_END_DATA        = 100108 ! 100108
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_ERROR_DATA      = 100109 ! 100109
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_EDGE_FLAG_DATA  = 100110 ! 100110
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_COMBINE_DATA    = 100111 ! 100111
!  TessContour values
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_CW                   = 100120 ! 100120
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_CCW                  = 100121 ! 100121
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_INTERIOR             = 100122 ! 100122
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_EXTERIOR             = 100123 ! 100123
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_UNKNOWN              = 100124 ! 100124
!  TessProperty values
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_WINDING_RULE    = 100140 ! 100140
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_BOUNDARY_ONLY   = 100141 ! 100141
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_TOLERANCE       = 100142 ! 100142
!  TessError values
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_ERROR1          = 100151 ! 100151
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_ERROR2          = 100152 ! 100152
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_ERROR3          = 100153 ! 100153
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_ERROR4          = 100154 ! 100154
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_ERROR5          = 100155 ! 100155
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_ERROR6          = 100156 ! 100156
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_ERROR7          = 100157 ! 100157
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_ERROR8          = 100158 ! 100158
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_MISSING_BEGIN_POLYGON = 100151 ! 100151
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_MISSING_BEGIN_CONTOUR = 100152 ! 100152
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_MISSING_END_POLYGON = 100153 ! 100153
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_MISSING_END_CONTOUR = 100154 ! 100154
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_COORD_TOO_LARGE = 100155 ! 100155
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_NEED_COMBINE_CALLBACK = 100156 ! 100156
!  TessWinding values
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_WINDING_ODD     = 100130 ! 100130
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_WINDING_NONZERO = 100131 ! 100131
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_WINDING_POSITIVE = 100132 ! 100132
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_WINDING_NEGATIVE = 100133 ! 100133
INTEGER(GLenum), PARAMETER, PUBLIC :: GLU_TESS_WINDING_ABS_GEQ_TWO = 100134 ! 100134


!  GLU  commands

!  BeginCurve(nurb)
!  return  void
!  param   nurb		NurbsObj in value
PUBLIC :: gluBeginCurve
INTERFACE
SUBROUTINE gluBeginCurve(nurb) BIND(C,NAME="gluBeginCurve")
IMPORT
TYPE(C_PTR), VALUE :: nurb
END SUBROUTINE gluBeginCurve

END INTERFACE

!  BeginPolygon(tess)
!  return  void
!  param   tess		TesselatorObj in value
PUBLIC :: gluBeginPolygon
INTERFACE
SUBROUTINE gluBeginPolygon(tess) BIND(C,NAME="gluBeginPolygon")
IMPORT
TYPE(C_PTR), VALUE :: tess
END SUBROUTINE gluBeginPolygon

END INTERFACE

!  BeginSurface(nurb)
!  return  void
!  param   nurb		NurbsObj in value
PUBLIC :: gluBeginSurface
INTERFACE
SUBROUTINE gluBeginSurface(nurb) BIND(C,NAME="gluBeginSurface")
IMPORT
TYPE(C_PTR), VALUE :: nurb
END SUBROUTINE gluBeginSurface

END INTERFACE

!  BeginTrim(nurb)
!  return  void
!  param   nurb		NurbsObj in value
PUBLIC :: gluBeginTrim
INTERFACE
SUBROUTINE gluBeginTrim(nurb) BIND(C,NAME="gluBeginTrim")
IMPORT
TYPE(C_PTR), VALUE :: nurb
END SUBROUTINE gluBeginTrim

END INTERFACE

!  Build1DMipmapLevels(target, internalFormat, width, format, type, level, base, max, data)
!  return  Int32
!  param   target		TextureTarget in value
!  param   internalFormat	Int32 in value
!  param   width		SizeI in value
!  param   format		PixelFormat in value
!  param   type		PixelType in value
!  param   level		Int32 in value
!  param   base		Int32 in value
!  param   max		Int32 in value
!  param   data		void in reference
PUBLIC :: gluBuild1DMipmapLevels
INTERFACE
FUNCTION gluBuild1DMipmapLevels(target, internalFormat, width, format, type,  &
    level, base, max, data) BIND(C,NAME="gluBuild1DMipmapLevels")
IMPORT
INTEGER(GLint) :: gluBuild1DMipmapLevels
INTEGER(GLenum), VALUE :: target, format, type
INTEGER(GLint), VALUE :: internalFormat, level, base, max
INTEGER(GLsizei), VALUE :: width
TYPE(C_PTR), VALUE :: data
END FUNCTION gluBuild1DMipmapLevels

END INTERFACE

!  Build1DMipmaps(target, internalFormat, width, format, type, data)
!  return  Int32
!  param   target		TextureTarget in value
!  param   internalFormat	Int32 in value
!  param   width		SizeI in value
!  param   format		PixelFormat in value
!  param   type		PixelType in value
!  param   data		void in reference
PUBLIC :: gluBuild1DMipmaps
INTERFACE
FUNCTION gluBuild1DMipmaps(target, internalFormat, width, format, type,       &
    data) BIND(C,NAME="gluBuild1DMipmaps")
IMPORT
INTEGER(GLint) :: gluBuild1DMipmaps
INTEGER(GLenum), VALUE :: target, format, type
INTEGER(GLint), VALUE :: internalFormat
INTEGER(GLsizei), VALUE :: width
TYPE(C_PTR), VALUE :: data
END FUNCTION gluBuild1DMipmaps

END INTERFACE

!  Build2DMipmapLevels(target, internalFormat, width, height, format, type, level, base, max, data)
!  return  Int32
!  param   target		TextureTarget in value
!  param   internalFormat	Int32 in value
!  param   width		SizeI in value
!  param   height		SizeI in value
!  param   format		PixelFormat in value
!  param   type		PixelType in value
!  param   level		Int32 in value
!  param   base		Int32 in value
!  param   max		Int32 in value
!  param   data		void in reference
PUBLIC :: gluBuild2DMipmapLevels
INTERFACE
FUNCTION gluBuild2DMipmapLevels(target, internalFormat, width, height,        &
    format, type, level, base, max, data)                                         &
    BIND(C,NAME="gluBuild2DMipmapLevels")
IMPORT
INTEGER(GLint) :: gluBuild2DMipmapLevels
INTEGER(GLenum), VALUE :: target, format, type
INTEGER(GLint), VALUE :: internalFormat, level, base, max
INTEGER(GLsizei), VALUE :: width, height
TYPE(C_PTR), VALUE :: data
END FUNCTION gluBuild2DMipmapLevels

END INTERFACE

!  Build2DMipmaps(target, internalFormat, width, height, format, type, data)
!  return  Int32
!  param   target		TextureTarget in value
!  param   internalFormat	Int32 in value
!  param   width		SizeI in value
!  param   height		SizeI in value
!  param   format		PixelFormat in value
!  param   type		PixelType in value
!  param   data		void in reference
PUBLIC :: gluBuild2DMipmaps
INTERFACE
FUNCTION gluBuild2DMipmaps(target, internalFormat, width, height, format,     &
    type, data) BIND(C,NAME="gluBuild2DMipmaps")
IMPORT
INTEGER(GLint) :: gluBuild2DMipmaps
INTEGER(GLenum), VALUE :: target, format, type
INTEGER(GLint), VALUE :: internalFormat
INTEGER(GLsizei), VALUE :: width, height
TYPE(C_PTR), VALUE :: data
END FUNCTION gluBuild2DMipmaps

END INTERFACE

!  Build3DMipmapLevels(target, internalFormat, width, height, depth, format, type, level, base, max, data)
!  return  Int32
!  param   target		TextureTarget in value
!  param   internalFormat	Int32 in value
!  param   width		SizeI in value
!  param   height		SizeI in value
!  param   depth		SizeI in value
!  param   format		PixelFormat in value
!  param   type		PixelType in value
!  param   level		Int32 in value
!  param   base		Int32 in value
!  param   max		Int32 in value
!  param   data		void in reference
PUBLIC :: gluBuild3DMipmapLevels
INTERFACE
FUNCTION gluBuild3DMipmapLevels(target, internalFormat, width, height,        &
    depth, format, type, level, base, max, data)                                  &
    BIND(C,NAME="gluBuild3DMipmapLevels")
IMPORT
INTEGER(GLint) :: gluBuild3DMipmapLevels
INTEGER(GLenum), VALUE :: target, format, type
INTEGER(GLint), VALUE :: internalFormat, level, base, max
INTEGER(GLsizei), VALUE :: width, height, depth
TYPE(C_PTR), VALUE :: data
END FUNCTION gluBuild3DMipmapLevels

END INTERFACE

!  Build3DMipmaps(target, internalFormat, width, height, depth, format, type, data)
!  return  Int32
!  param   target		TextureTarget in value
!  param   internalFormat	Int32 in value
!  param   width		SizeI in value
!  param   height		SizeI in value
!  param   depth		SizeI in value
!  param   format		PixelFormat in value
!  param   type		PixelType in value
!  param   data		void in reference
PUBLIC :: gluBuild3DMipmaps
INTERFACE
FUNCTION gluBuild3DMipmaps(target, internalFormat, width, height, depth,      &
    format, type, data) BIND(C,NAME="gluBuild3DMipmaps")
IMPORT
INTEGER(GLint) :: gluBuild3DMipmaps
INTEGER(GLenum), VALUE :: target, format, type
INTEGER(GLint), VALUE :: internalFormat
INTEGER(GLsizei), VALUE :: width, height, depth
TYPE(C_PTR), VALUE :: data
END FUNCTION gluBuild3DMipmaps

END INTERFACE

!  CheckExtension(extName, extString)
!  return  Boolean
!  param   extName		UInt8 in array [COMPSIZE()]
!  param   extString	UInt8 in array [COMPSIZE()]
PUBLIC :: gluCheckExtension
INTERFACE
FUNCTION gluCheckExtension(extName, extString)                                &
    BIND(C,NAME="gluCheckExtension")
IMPORT
INTEGER(GLboolean) :: gluCheckExtension
INTEGER(GLbyte), DIMENSION(*), INTENT(IN) :: extName, extString
END FUNCTION gluCheckExtension

END INTERFACE

!  Cylinder(quad, base, top, height, slices, stacks)
!  return  void
!  param   quad		QuadricObj in value
!  param   base		Float64 in value
!  param   top		Float64 in value
!  param   height		Float64 in value
!  param   slices		Int32 in value
!  param   stacks		Int32 in value
PUBLIC :: gluCylinder
INTERFACE
SUBROUTINE gluCylinder(quad, base, top, height, slices, stacks)               &
    BIND(C,NAME="gluCylinder")
IMPORT
INTEGER(GLint), VALUE :: slices, stacks
REAL(GLdouble), VALUE :: base, top, height
TYPE(C_PTR), VALUE :: quad
END SUBROUTINE gluCylinder

END INTERFACE

!  DeleteNurbsRenderer(nurb)
!  return  void
!  param   nurb		NurbsObj in value
PUBLIC :: gluDeleteNurbsRenderer
INTERFACE
SUBROUTINE gluDeleteNurbsRenderer(nurb)                                       &
    BIND(C,NAME="gluDeleteNurbsRenderer")
IMPORT
TYPE(C_PTR), VALUE :: nurb
END SUBROUTINE gluDeleteNurbsRenderer

END INTERFACE

!  DeleteQuadric(quad)
!  return  void
!  param   quad		QuadricObj in value
PUBLIC :: gluDeleteQuadric
INTERFACE
SUBROUTINE gluDeleteQuadric(quad) BIND(C,NAME="gluDeleteQuadric")
IMPORT
TYPE(C_PTR), VALUE :: quad
END SUBROUTINE gluDeleteQuadric

END INTERFACE

!  DeleteTess(tess)
!  return  void
!  param   tess		TesselatorObj in value
PUBLIC :: gluDeleteTess
INTERFACE
SUBROUTINE gluDeleteTess(tess) BIND(C,NAME="gluDeleteTess")
IMPORT
TYPE(C_PTR), VALUE :: tess
END SUBROUTINE gluDeleteTess

END INTERFACE

!  Disk(quad, inner, outer, slices, loops)
!  return  void
!  param   quad		QuadricObj in value
!  param   inner		Float64 in value
!  param   outer		Float64 in value
!  param   slices		Int32 in value
!  param   loops		Int32 in value
PUBLIC :: gluDisk
INTERFACE
SUBROUTINE gluDisk(quad, inner, outer, slices, loops)                         &
    BIND(C,NAME="gluDisk")
IMPORT
INTEGER(GLint), VALUE :: slices, loops
REAL(GLdouble), VALUE :: inner, outer
TYPE(C_PTR), VALUE :: quad
END SUBROUTINE gluDisk

END INTERFACE

!  EndCurve(nurb)
!  return  void
!  param   nurb		NurbsObj in value
PUBLIC :: gluEndCurve
INTERFACE
SUBROUTINE gluEndCurve(nurb) BIND(C,NAME="gluEndCurve")
IMPORT
TYPE(C_PTR), VALUE :: nurb
END SUBROUTINE gluEndCurve

END INTERFACE

!  EndPolygon(tess)
!  return  void
!  param   tess		TesselatorObj in value
PUBLIC :: gluEndPolygon
INTERFACE
SUBROUTINE gluEndPolygon(tess) BIND(C,NAME="gluEndPolygon")
IMPORT
TYPE(C_PTR), VALUE :: tess
END SUBROUTINE gluEndPolygon

END INTERFACE

!  EndSurface(nurb)
!  return  void
!  param   nurb		NurbsObj in value
PUBLIC :: gluEndSurface
INTERFACE
SUBROUTINE gluEndSurface(nurb) BIND(C,NAME="gluEndSurface")
IMPORT
TYPE(C_PTR), VALUE :: nurb
END SUBROUTINE gluEndSurface

END INTERFACE

!  EndTrim(nurb)
!  return  void
!  param   nurb		NurbsObj in value
PUBLIC :: gluEndTrim
INTERFACE
SUBROUTINE gluEndTrim(nurb) BIND(C,NAME="gluEndTrim")
IMPORT
TYPE(C_PTR), VALUE :: nurb
END SUBROUTINE gluEndTrim

END INTERFACE

!  ErrorString(error)
!  return  String
!  param   error		ErrorCode in value
PUBLIC :: gluErrorString
INTERFACE
FUNCTION gluErrorString(error) BIND(C,NAME="gluErrorString")
IMPORT
TYPE(C_PTR) :: gluErrorString
INTEGER(GLenum), VALUE :: error
END FUNCTION gluErrorString

END INTERFACE

!  GetString(name)
!  return  String
!  param   name		StringName in value
PUBLIC :: gluGetString
INTERFACE
FUNCTION gluGetString(name) BIND(C,NAME="gluGetString")
IMPORT
TYPE(C_PTR) :: gluGetString
INTEGER(GLenum), VALUE :: name
END FUNCTION gluGetString

END INTERFACE

!  GetNurbsProperty(nurb, property, data)
!  return  void
!  param   nurb		NurbsObj in value
!  param   property	NurbsProperty in value
!  param   data		Float32Pointer out value
PUBLIC :: gluGetNurbsProperty
INTERFACE
SUBROUTINE gluGetNurbsProperty(nurb, property, data)                          &
    BIND(C,NAME="gluGetNurbsProperty")
IMPORT
INTEGER(GLenum), VALUE :: property
TYPE(C_PTR), VALUE :: nurb
REAL(GLfloat), DIMENSION(*), INTENT(OUT) :: data
END SUBROUTINE gluGetNurbsProperty

END INTERFACE

!  GetTessProperty(tess, which, data)
!  return  void
!  param   tess		TesselatorObj in value
!  param   which		TessProperty in value
!  param   data		Float64Pointer out value
PUBLIC :: gluGetTessProperty
INTERFACE
SUBROUTINE gluGetTessProperty(tess, which, data)                              &
    BIND(C,NAME="gluGetTessProperty")
IMPORT
INTEGER(GLenum), VALUE :: which
TYPE(C_PTR), VALUE :: tess
REAL(GLdouble), DIMENSION(*), INTENT(OUT) :: data
END SUBROUTINE gluGetTessProperty

END INTERFACE

!  LoadSamplingMatrices(nurb, model, perspective, view)
!  return  void
!  param   nurb		NurbsObj in value
!  param   model		Float32 in array [16]
!  param   perspective	Float32 in array [16]
!  param   view		Int32 in array [4]
PUBLIC :: gluLoadSamplingMatrices
INTERFACE
SUBROUTINE gluLoadSamplingMatrices(nurb, model, perspective, view)            &
    BIND(C,NAME="gluLoadSamplingMatrices")
IMPORT
TYPE(C_PTR), VALUE :: nurb
INTEGER(GLint), DIMENSION(4), INTENT(IN) :: view
REAL(GLfloat), DIMENSION(16), INTENT(IN) :: model, perspective
END SUBROUTINE gluLoadSamplingMatrices

END INTERFACE

!  LookAt(eyeX, eyeY, eyeZ, centerX, centerY, centerZ, upX, upY, upZ)
!  return  void
!  param   eyeX		Float64 in value
!  param   eyeY		Float64 in value
!  param   eyeZ		Float64 in value
!  param   centerX		Float64 in value
!  param   centerY		Float64 in value
!  param   centerZ		Float64 in value
!  param   upX		Float64 in value
!  param   upY		Float64 in value
!  param   upZ		Float64 in value
PUBLIC :: gluLookAt
INTERFACE
SUBROUTINE gluLookAt(eyeX, eyeY, eyeZ, centerX, centerY, centerZ, upX, upY,   &
    upZ) BIND(C,NAME="gluLookAt")
IMPORT
REAL(GLdouble), VALUE :: eyeX, eyeY, eyeZ, centerX, centerY, centerZ, upX, upY, upZ
END SUBROUTINE gluLookAt

END INTERFACE

!  NewNurbsRenderer()
!  return  NurbsObj
PUBLIC :: gluNewNurbsRenderer
INTERFACE
FUNCTION gluNewNurbsRenderer() BIND(C,NAME="gluNewNurbsRenderer")
IMPORT
TYPE(C_PTR) :: gluNewNurbsRenderer
END FUNCTION gluNewNurbsRenderer

END INTERFACE

!  NewQuadric()
!  return  QuadricObj
PUBLIC :: gluNewQuadric
INTERFACE
FUNCTION gluNewQuadric() BIND(C,NAME="gluNewQuadric")
IMPORT
TYPE(C_PTR) :: gluNewQuadric
END FUNCTION gluNewQuadric

END INTERFACE

!  NewTess()
!  return  TesselatorObj
PUBLIC :: gluNewTess
INTERFACE
FUNCTION gluNewTess() BIND(C,NAME="gluNewTess")
IMPORT
TYPE(C_PTR) :: gluNewTess
END FUNCTION gluNewTess

END INTERFACE

!  NextContour(tess, type)
!  return  void
!  param   tess		TesselatorObj in value
!  param   type		TessContour in value
PUBLIC :: gluNextContour
INTERFACE
SUBROUTINE gluNextContour(tess, type) BIND(C,NAME="gluNextContour")
IMPORT
INTEGER(GLenum), VALUE :: type
TYPE(C_PTR), VALUE :: tess
END SUBROUTINE gluNextContour

END INTERFACE

!  NurbsCallback(nurb, which, CallBackFunc)
!  return  void
!  param   nurb		NurbsObj in value
!  param   which		NurbsCallback in value
!  param   CallBackFunc	FunctionPointer in value
PUBLIC :: gluNurbsCallback
INTERFACE
SUBROUTINE gluNurbsCallback(nurb, which, CallBackFunc)                        &
    BIND(C,NAME="gluNurbsCallback")
IMPORT
INTEGER(GLenum), VALUE :: which
TYPE(C_FUNPTR), VALUE :: CallBackFunc
TYPE(C_PTR), VALUE :: nurb
END SUBROUTINE gluNurbsCallback

END INTERFACE

!  NurbsCallbackData(nurb, userData)
!  return  void
!  param   nurb		NurbsObj in value
!  param   userData	VoidPointer in value
PUBLIC :: gluNurbsCallbackData
INTERFACE
SUBROUTINE gluNurbsCallbackData(nurb, userData)                               &
    BIND(C,NAME="gluNurbsCallbackData")
IMPORT
TYPE(C_PTR), VALUE :: nurb, userData
END SUBROUTINE gluNurbsCallbackData

END INTERFACE

!  NurbsCallbackDataEXT(nurb, userData)
!  return  void
!  param   nurb		NurbsObj in value
!  param   userData	VoidPointer in value
PUBLIC :: gluNurbsCallbackDataEXT
INTERFACE
SUBROUTINE gluNurbsCallbackDataEXT(nurb, userData)                            &
    BIND(C,NAME="gluNurbsCallbackDataEXT")
IMPORT
TYPE(C_PTR), VALUE :: nurb, userData
END SUBROUTINE gluNurbsCallbackDataEXT

END INTERFACE

!  NurbsCurve(nurb, knotCount, knots, stride, control, order, type)
!  return  void
!  param   nurb		NurbsObj in value
!  param   knotCount	Int32 in value
!  param   knots		Float32 out reference
!  param   stride		Int32 in value
!  param   control		Float32 out reference
!  param   order		Int32 in value
!  param   type		MapTarget in value
PUBLIC :: gluNurbsCurve
INTERFACE
SUBROUTINE gluNurbsCurve(nurb, knotCount, knots, stride, control, order,      &
    type) BIND(C,NAME="gluNurbsCurve")
IMPORT
INTEGER(GLenum), VALUE :: type
INTEGER(GLint), VALUE :: knotCount, stride, order
TYPE(C_PTR), VALUE :: nurb
REAL(GLfloat), DIMENSION(*), INTENT(OUT) :: knots, control
END SUBROUTINE gluNurbsCurve

END INTERFACE

!  NurbsProperty(nurb, property, value)
!  return  void
!  param   nurb		NurbsObj in value
!  param   property	NurbsProperty in value
!  param   value		Float32 in value
PUBLIC :: gluNurbsProperty
INTERFACE
SUBROUTINE gluNurbsProperty(nurb, property, value)                            &
    BIND(C,NAME="gluNurbsProperty")
IMPORT
INTEGER(GLenum), VALUE :: property
REAL(GLfloat), VALUE :: value
TYPE(C_PTR), VALUE :: nurb
END SUBROUTINE gluNurbsProperty

END INTERFACE

!  NurbsSurface(nurb, sKnotCount, sKnots, tKnotCount, tKnots, sStride, tStride, control, sOrder, tOrder, type)
!  return  void
!  param   nurb		NurbsObj in value
!  param   sKnotCount	Int32 in value
!  param   sKnots		Float32Pointer in value
!  param   tKnotCount	Int32 in value
!  param   tKnots		Float32Pointer in value
!  param   sStride		Int32 in value
!  param   tStride		Int32 in value
!  param   control		Float32Pointer in value
!  param   sOrder		Int32 in value
!  param   tOrder		Int32 in value
!  param   type		MapTarget in value
PUBLIC :: gluNurbsSurface
INTERFACE
SUBROUTINE gluNurbsSurface(nurb, sKnotCount, sKnots, tKnotCount, tKnots,      &
    sStride, tStride, control, sOrder, tOrder, type)                              &
    BIND(C,NAME="gluNurbsSurface")
IMPORT
INTEGER(GLenum), VALUE :: type
INTEGER(GLint), VALUE :: sKnotCount, tKnotCount, sStride, tStride, sOrder, tOrder
TYPE(C_PTR), VALUE :: nurb
REAL(GLfloat), DIMENSION(*), INTENT(IN) :: sKnots, tKnots, control
END SUBROUTINE gluNurbsSurface

END INTERFACE

!  Ortho2D(left, right, bottom, top)
!  return  void
!  param   left		Float64 in value
!  param   right		Float64 in value
!  param   bottom		Float64 in value
!  param   top		Float64 in value
PUBLIC :: gluOrtho2D
INTERFACE
SUBROUTINE gluOrtho2D(left, right, bottom, top) BIND(C,NAME="gluOrtho2D")
IMPORT
REAL(GLdouble), VALUE :: left, right, bottom, top
END SUBROUTINE gluOrtho2D

END INTERFACE

!  PartialDisk(quad, inner, outer, slices, loops, start, sweep)
!  return  void
!  param   quad		QuadricObj in value
!  param   inner		Float64 in value
!  param   outer		Float64 in value
!  param   slices		Int32 in value
!  param   loops		Int32 in value
!  param   start		Float64 in value
!  param   sweep		Float64 in value
PUBLIC :: gluPartialDisk
INTERFACE
SUBROUTINE gluPartialDisk(quad, inner, outer, slices, loops, start, sweep)    &
    BIND(C,NAME="gluPartialDisk")
IMPORT
INTEGER(GLint), VALUE :: slices, loops
REAL(GLdouble), VALUE :: inner, outer, start, sweep
TYPE(C_PTR), VALUE :: quad
END SUBROUTINE gluPartialDisk

END INTERFACE

!  Perspective(fovy, aspect, zNear, zFar)
!  return  void
!  param   fovy		Float64 in value
!  param   aspect		Float64 in value
!  param   zNear		Float64 in value
!  param   zFar		Float64 in value
PUBLIC :: gluPerspective
INTERFACE
SUBROUTINE gluPerspective(fovy, aspect, zNear, zFar)                          &
    BIND(C,NAME="gluPerspective")
IMPORT
REAL(GLdouble), VALUE :: fovy, aspect, zNear, zFar
END SUBROUTINE gluPerspective

END INTERFACE

!  PickMatrix(x, y, delX, delY, viewport)
!  return  void
!  param   x		Float64 in value
!  param   y		Float64 in value
!  param   delX		Float64 in value
!  param   delY		Float64 in value
!  param   viewport	Int32 out array [4]
PUBLIC :: gluPickMatrix
INTERFACE
SUBROUTINE gluPickMatrix(x, y, delX, delY, viewport)                          &
    BIND(C,NAME="gluPickMatrix")
IMPORT
REAL(GLdouble), VALUE :: x, y, delX, delY
INTEGER(GLint), DIMENSION(4), INTENT(OUT) :: viewport
END SUBROUTINE gluPickMatrix

END INTERFACE

!  Project(objX, objY, objZ, model, proj, view, winX, winY, winZ)
!  return  Int32
!  param   objX		Float64 in value
!  param   objY		Float64 in value
!  param   objZ		Float64 in value
!  param   model		Float64 in array [16]
!  param   proj		Float64 in array [16]
!  param   view		Int32 in array [4]
!  param   winX		Float64Pointer in value
!  param   winY		Float64Pointer in value
!  param   winZ		Float64Pointer in value
PUBLIC :: gluProject
INTERFACE
FUNCTION gluProject(objX, objY, objZ, model, proj, view, winX, winY, winZ)    &
    BIND(C,NAME="gluProject")
IMPORT
INTEGER(GLint) :: gluProject
REAL(GLdouble), VALUE :: objX, objY, objZ
INTEGER(GLint), DIMENSION(4), INTENT(IN) :: view
REAL(GLdouble), DIMENSION(*), INTENT(IN) :: winX, winY, winZ
REAL(GLdouble), DIMENSION(16), INTENT(IN) :: model, proj
END FUNCTION gluProject

END INTERFACE

!  PwlCurve(nurb, count, data, stride, type)
!  return  void
!  param   nurb		NurbsObj in value
!  param   count		Int32 in value
!  param   data		Float32Pointer in value
!  param   stride		Int32 in value
!  param   type		NurbsTrim in value
PUBLIC :: gluPwlCurve
INTERFACE
SUBROUTINE gluPwlCurve(nurb, count, data, stride, type)                       &
    BIND(C,NAME="gluPwlCurve")
IMPORT
INTEGER(GLenum), VALUE :: type
INTEGER(GLint), VALUE :: count, stride
TYPE(C_PTR), VALUE :: nurb
REAL(GLfloat), DIMENSION(*), INTENT(IN) :: data
END SUBROUTINE gluPwlCurve

END INTERFACE

!  QuadricCallback(quad, which, CallBackFunc)
!  return  void
!  param   quad		QuadricObj in value
!  param   which		QuadricCallback in value
!  param   CallBackFunc	FunctionPointer in value
PUBLIC :: gluQuadricCallback
INTERFACE
SUBROUTINE gluQuadricCallback(quad, which, CallBackFunc)                      &
    BIND(C,NAME="gluQuadricCallback")
IMPORT
INTEGER(GLenum), VALUE :: which
TYPE(C_FUNPTR), VALUE :: CallBackFunc
TYPE(C_PTR), VALUE :: quad
END SUBROUTINE gluQuadricCallback

END INTERFACE

!  QuadricDrawStyle(quad, draw)
!  return  void
!  param   quad		QuadricObj in value
!  param   draw		QuadricDrawStyle in value
PUBLIC :: gluQuadricDrawStyle
INTERFACE
SUBROUTINE gluQuadricDrawStyle(quad, draw)                                    &
    BIND(C,NAME="gluQuadricDrawStyle")
IMPORT
INTEGER(GLenum), VALUE :: draw
TYPE(C_PTR), VALUE :: quad
END SUBROUTINE gluQuadricDrawStyle

END INTERFACE

!  QuadricNormals(quad, normal)
!  return  void
!  param   quad		QuadricObj in value
!  param   normal		QuadricNormal in value
PUBLIC :: gluQuadricNormals
INTERFACE
SUBROUTINE gluQuadricNormals(quad, normal) BIND(C,NAME="gluQuadricNormals")
IMPORT
INTEGER(GLenum), VALUE :: normal
TYPE(C_PTR), VALUE :: quad
END SUBROUTINE gluQuadricNormals

END INTERFACE

!  QuadricOrientation(quad, orientation)
!  return  void
!  param   quad		QuadricObj in value
!  param   orientation	QuadricOrientation in value
PUBLIC :: gluQuadricOrientation
INTERFACE
SUBROUTINE gluQuadricOrientation(quad, orientation)                           &
    BIND(C,NAME="gluQuadricOrientation")
IMPORT
INTEGER(GLenum), VALUE :: orientation
TYPE(C_PTR), VALUE :: quad
END SUBROUTINE gluQuadricOrientation

END INTERFACE

!  QuadricTexture(quad, texture)
!  return  void
!  param   quad		QuadricObj in value
!  param   texture		Boolean in value
PUBLIC :: gluQuadricTexture
INTERFACE
SUBROUTINE gluQuadricTexture(quad, texture)                                   &
    BIND(C,NAME="gluQuadricTexture")
IMPORT
INTEGER(GLboolean), VALUE :: texture
TYPE(C_PTR), VALUE :: quad
END SUBROUTINE gluQuadricTexture

END INTERFACE

!  ScaleImage(format, wIn, hIn, typeIn, dataIn, wOut, hOut, typeOut, dataOut)
!  return  Int32
!  param   format		PixelFormat in value
!  param   wIn		SizeI in value
!  param   hIn		SizeI in value
!  param   typeIn		PixelType in value
!  param   dataIn		void in reference
!  param   wOut		SizeI in value
!  param   hOut		SizeI in value
!  param   typeOut		PixelType in value
!  param   dataOut		VoidPointer out value
PUBLIC :: gluScaleImage
INTERFACE
FUNCTION gluScaleImage(format, wIn, hIn, typeIn, dataIn, wOut, hOut,          &
    typeOut, dataOut) BIND(C,NAME="gluScaleImage")
IMPORT
INTEGER(GLint) :: gluScaleImage
INTEGER(GLenum), VALUE :: format, typeIn, typeOut
INTEGER(GLsizei), VALUE :: wIn, hIn, wOut, hOut
TYPE(C_PTR), VALUE :: dataIn, dataOut
END FUNCTION gluScaleImage

END INTERFACE

!  Sphere(quad, radius, slices, stacks)
!  return  void
!  param   quad		QuadricObj in value
!  param   radius		Float64 in value
!  param   slices		Int32 in value
!  param   stacks		Int32 in value
PUBLIC :: gluSphere
INTERFACE
SUBROUTINE gluSphere(quad, radius, slices, stacks) BIND(C,NAME="gluSphere")
IMPORT
INTEGER(GLint), VALUE :: slices, stacks
REAL(GLdouble), VALUE :: radius
TYPE(C_PTR), VALUE :: quad
END SUBROUTINE gluSphere

END INTERFACE

!  TessBeginContour(tess)
!  return  void
!  param   tess		TesselatorObj in value
PUBLIC :: gluTessBeginContour
INTERFACE
SUBROUTINE gluTessBeginContour(tess) BIND(C,NAME="gluTessBeginContour")
IMPORT
TYPE(C_PTR), VALUE :: tess
END SUBROUTINE gluTessBeginContour

END INTERFACE

!  TessBeginPolygon(tess, data)
!  return  void
!  param   tess		TesselatorObj in value
!  param   data		VoidPointer in value
PUBLIC :: gluTessBeginPolygon
INTERFACE
SUBROUTINE gluTessBeginPolygon(tess, data)                                    &
    BIND(C,NAME="gluTessBeginPolygon")
IMPORT
TYPE(C_PTR), VALUE :: tess, data
END SUBROUTINE gluTessBeginPolygon

END INTERFACE

!  TessCallback(tess, which, CallBackFunc)
!  return  void
!  param   tess		TesselatorObj in value
!  param   which		TessCallback in value
!  param   CallBackFunc	FunctionPointer in value
PUBLIC :: gluTessCallback
INTERFACE
SUBROUTINE gluTessCallback(tess, which, CallBackFunc)                         &
    BIND(C,NAME="gluTessCallback")
IMPORT
INTEGER(GLenum), VALUE :: which
TYPE(C_FUNPTR), VALUE :: CallBackFunc
TYPE(C_PTR), VALUE :: tess
END SUBROUTINE gluTessCallback

END INTERFACE

!  TessEndContour(tess)
!  return  void
!  param   tess		TesselatorObj in value
PUBLIC :: gluTessEndContour
INTERFACE
SUBROUTINE gluTessEndContour(tess) BIND(C,NAME="gluTessEndContour")
IMPORT
TYPE(C_PTR), VALUE :: tess
END SUBROUTINE gluTessEndContour

END INTERFACE

!  TessEndPolygon(tess)
!  return  void
!  param   tess		TesselatorObj in value
PUBLIC :: gluTessEndPolygon
INTERFACE
SUBROUTINE gluTessEndPolygon(tess) BIND(C,NAME="gluTessEndPolygon")
IMPORT
TYPE(C_PTR), VALUE :: tess
END SUBROUTINE gluTessEndPolygon

END INTERFACE

!  TessNormal(tess, valueX, valueY, valueZ)
!  return  void
!  param   tess		TesselatorObj in value
!  param   valueX		Float64 in value
!  param   valueY		Float64 in value
!  param   valueZ		Float64 in value
PUBLIC :: gluTessNormal
INTERFACE
SUBROUTINE gluTessNormal(tess, valueX, valueY, valueZ)                        &
    BIND(C,NAME="gluTessNormal")
IMPORT
REAL(GLdouble), VALUE :: valueX, valueY, valueZ
TYPE(C_PTR), VALUE :: tess
END SUBROUTINE gluTessNormal

END INTERFACE

!  TessProperty(tess, which, data)
!  return  void
!  param   tess		TesselatorObj in value
!  param   which		TessProperty in value
!  param   data		Float64 in value
PUBLIC :: gluTessProperty
INTERFACE
SUBROUTINE gluTessProperty(tess, which, data)                                 &
    BIND(C,NAME="gluTessProperty")
IMPORT
INTEGER(GLenum), VALUE :: which
REAL(GLdouble), VALUE :: data
TYPE(C_PTR), VALUE :: tess
END SUBROUTINE gluTessProperty

END INTERFACE

!  TessVertex(tess, location, data)
!  return  void
!  param   tess		TesselatorObj in value
!  param   location	Float64 out array [3]
!  param   data		VoidPointer in value
PUBLIC :: gluTessVertex
INTERFACE
SUBROUTINE gluTessVertex(tess, location, data) BIND(C,NAME="gluTessVertex")
IMPORT
TYPE(C_PTR), VALUE :: tess, data
REAL(GLdouble), DIMENSION(3), INTENT(OUT) :: location
END SUBROUTINE gluTessVertex

END INTERFACE

!  UnProject(winX, winY, winZ, model, proj, view, objX, objY, objZ)
!  return  Int32
!  param   winX		Float64 in value
!  param   winY		Float64 in value
!  param   winZ		Float64 in value
!  param   model		Float64 in array [16]
!  param   proj		Float64 in array [16]
!  param   view		Int32 in array [4]
!  param   objX		Float64Pointer in value
!  param   objY		Float64Pointer in value
!  param   objZ		Float64Pointer in value
PUBLIC :: gluUnProject
INTERFACE
FUNCTION gluUnProject(winX, winY, winZ, model, proj, view, objX, objY, objZ)  &
    BIND(C,NAME="gluUnProject")
IMPORT
INTEGER(GLint) :: gluUnProject
REAL(GLdouble), VALUE :: winX, winY, winZ
INTEGER(GLint), DIMENSION(4), INTENT(IN) :: view
REAL(GLdouble), DIMENSION(*), INTENT(IN) :: objX, objY, objZ
REAL(GLdouble), DIMENSION(16), INTENT(IN) :: model, proj
END FUNCTION gluUnProject

END INTERFACE

!  UnProject4(winX, winY, winZ, clipW, model, proj, view, near, far, objX, objY, objZ, objW)
!  return  Int32
!  param   winX		Float64 in value
!  param   winY		Float64 in value
!  param   winZ		Float64 in value
!  param   clipW		Float64 in value
!  param   model		Float64 in array [16]
!  param   proj		Float64 in array [16]
!  param   view		Int32 in array [4]
!  param   near		Float64 in value
!  param   far		Float64 in value
!  param   objX		Float64Pointer in value
!  param   objY		Float64Pointer in value
!  param   objZ		Float64Pointer in value
PUBLIC :: gluUnProject4
INTERFACE
FUNCTION gluUnProject4(winX, winY, winZ, clipW, model, proj, view, near,      &
    far, objX, objY, objZ, objW) BIND(C,NAME="gluUnProject4")
IMPORT
INTEGER(GLint) :: gluUnProject4
REAL(GLdouble), VALUE :: winX, winY, winZ, clipW, near, far
INTEGER(GLint), DIMENSION(4), INTENT(IN) :: view
REAL(GLdouble), DIMENSION(*), INTENT(IN) :: objX, objY, objZ, objW
REAL(GLdouble), DIMENSION(16), INTENT(IN) :: model, proj
END FUNCTION gluUnProject4

END INTERFACE


END MODULE  OpenGL_glu
