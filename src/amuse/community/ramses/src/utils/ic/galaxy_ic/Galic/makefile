EXEC   = bin/galic_MakeDiskGalaxy

SRCS   = main.c  globvars.c \
	halo.c haloset.c escape.c \
	disk.c diskset.c \
	bulge.c bulgeset.c \
	gas.c  gasset.c \
	save.c init.c structure.c \
	misc.c toomre.c force.c \
	forcetree.c \
	nrsrc/bessi0.c   nrsrc/bessj1.c   \
	nrsrc/erff.c     nrsrc/polint.c   nrsrc/spline.c \
	nrsrc/bessi1.c   nrsrc/bessk0.c   nrsrc/gasdev.c  \
	nrsrc/qromb.c    nrsrc/splint.c  \
	nrsrc/bessj0.c   nrsrc/bessk1.c   \
	nrsrc/nrutil.c   nrsrc/ran1.c     nrsrc/trapzd.c \
	nrsrc/gammp.c    nrsrc/gcf.c      nrsrc/gser.c  \
        nrsrc/gammln.c   nrsrc/indexx.c   nrsrc/zriddr.c

OBJS   = $(SRCS:.c=.o)
INCL   = prototypes.h globvars.h nrsrc/nrutil.h nrsrc/nr.h


CFLAGS = -O2 

LIBS   =  -lm

CC     =  gcc

all: $(EXEC)
	#cd part_ic_converter ; make ; cd ../

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(EPATH)$(EXEC)  $(LIBS)

$(OBJS): $(INCL) 

clean: 
	rm *.o nrsrc/*.o

cleanall: clean
	cd part_ic_converter/ ; make clean ; cd ../
