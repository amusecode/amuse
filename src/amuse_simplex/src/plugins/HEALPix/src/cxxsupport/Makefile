MPICXX ?= mpicxx

VPATH	= $(INCDIR)
HEADERS=cxxutils.h arr.h paramfile.h \
	fitshandle.h message_error.h vec3.h \
	lsconstants.h rotmatrix.h pointing.h \
	planck_rng.h geom_utils.h simparams.h \
	tga_image.h xcomplex.h trafos.h \
	fftpack_support.h datatypes.h \
	openmp_support.h

#include $(PARAMFILE)
INCDIR	= ../../include
LIBDIR	= ../../lib

INCLUDE = -I$(INCDIR)

all: libcxxsupport.a

%.o : %.cc
	$(MPICXX) $(CXXCFLAGS) -o $@ -c $< $(INCLUDE)

clean:
	$(RM) *.o *.a

MESSAGE_ERROR_H := message_error.h
LSCONSTANTS_H := lsconstants.h
DATATYPES_H := datatypes.h $(MESSAGE_ERROR_H)
CXXUTILS_H := cxxutils.h $(MESSAGE_ERROR_H) $(LSCONSTANTS_H)
SIMPARAMS_H := simparams.h $(CXXUTILS_H)
PARAMFILE_H := paramfile.h $(SIMPARAMS_H) $(CXXUTILS_H)
ARR_H := arr.h $(CXXUTILS_H)
FITSHANDLE_H := fitshandle.h fitsio.h $(ARR_H) $(DATATYPES_H)
VEC3_H := vec3.h
ROTMATRIX_H := rotmatrix.h $(CXXUTILS_H) $(VEC3_H)
POINTING_H := pointing.h $(VEC3_H) $(CXXUTILS_H)
PLANCK_RNG_H := planck_rng.h $(CXXUTILS_H)
GEOM_UTILS_H := geom_utils.h $(CXXUTILS_H) $(VEC3_H)
TGA_IMAGE_H := tga_image.h $(ARR_H)
XCOMPLEX_H := xcomplex.h
TRAFOS_H := trafos.h $(ROTMATRIX_H) $(GEOM_UTILS_H) $(POINTING_H)
FFTPACK_SUPPORT_H := fftpack_support.h $(ARR_H) $(XCOMPLEX_H)

SUPPORT_OBJ= cxxutils.o fitshandle.o rotmatrix.o simparams.o tga_image.o \
	trafos.o

cxxutils.o: $(CXXUTILS_H)
fitshandle.o: $(FITSHANDLE_H) $(CXXUTILS_H)
rotmatrix.o: $(ROTMATRIX_H) $(VEC3_H) $(LSCONSTANTS_H)
simparams.o: $(SIMPARAMS_H) $(FITSHANDLE_H)
tga_image.o: $(TGA_IMAGE_H) font_data.inc
trafos.o: $(TRAFOS_H) $(LSCONSTANTS_H)

libcxxsupport.a: $(SUPPORT_OBJ) $(HEADERS)
	$(ARCREATE) libcxxsupport.a $(SUPPORT_OBJ)
