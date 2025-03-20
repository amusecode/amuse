MPICXX ?= mpicxx

SPHERE_HEADERS= alm.h alm_fitsio.h \
	alm_powspec_tools.h powspec.h \
	powspec_fitsio.h \
	ylmgen.h alm_map_tools.h

HEADERS= healpix_base.h healpix_map.h \
	healpix_map_fitsio.h alm_healpix_tools.h \
	healpix_data_io.h healpix_base2.h \
	$(SPHERE_HEADERS)

SPHERE_OBJ= alm_fitsio.o powspec_fitsio.o alm_powspec_tools.o powspec.o \
	alm_map_tools.o

HEALPIX_OBJ= healpix_base.o healpix_map.o healpix_map_fitsio.o \
	alm_healpix_tools.o healpix_data_io.o healpix_base2.o $(SPHERE_OBJ)

#include $(PARAMFILE)
INCDIR	= ../include
LIBDIR	= ../lib

INCLUDE = -I$(INCDIR)
ARCREATE= ar crv
export ARCREATE


all: libhealpix_cxx.a

libcxxsupport.a:
	@cd cxxsupport ; \
	make all
	cp -f cxxsupport/libcxxsupport.a $(LIBDIR)/libcxxsupport.a
	cp -f cxxsupport/libcxxsupport.a libcxxsupport.a

libfftpack.a:
	@cd libfftpack ; \
	make all
	cp -f libfftpack/libfftpack.a $(LIBDIR)/libfftpack.a
	cp -f libfftpack/libfftpack.a libfftpack.a

libhealpix_cxx.a: $(HEALPIX_OBJ) $(HEADERS)
	$(ARCREATE) libhealpix_cxx.a $(HEALPIX_OBJ)
	cp -f libhealpix_cxx.a $(LIBDIR)/libhealpix_cxx.a

%.o : %.cc
	$(MPICXX) $(CXXCFLAGS) -o $@ -c $< $(INCLUDE)

clean:
	$(RM) *.o *.a $(LIBDIR)/*.a
	@cd cxxsupport ; \
	make clean
	@cd libfftpack ; \
	make clean

healpix_base.o: healpix_base.h libcxxsupport.a
healpix_base2.o: healpix_base.h healpix_base2.h libcxxsupport.a
healpix_map.o: healpix_base.h healpix_map.h libcxxsupport.a
healpix_map_fitsio.o: healpix_map_fitsio.h healpix_map.h healpix_base.h \
	libcxxsupport.a
alm_fitsio.o: alm_fitsio.h alm.h \
	libcxxsupport.a
healpix_data_io.o: healpix_data_io.h libcxxsupport.a
powspec_fitsio.o: powspec.h libcxxsupport.a
powspec.o: powspec.h libcxxsupport.a
alm_healpix_tools.o: healpix_base.h healpix_map.h alm.h \
	ylmgen.h alm_healpix_tools.h alm_map_tools.h \
	libfftpack.a libcxxsupport.a
alm_map_tools.o: alm.h ylmgen.h alm_map_tools.h libfftpack.a libcxxsupport.a
alm_powspec_tools.o: powspec.h alm.h alm_powspec_tools.h libcxxsupport.a

