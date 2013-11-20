from f90 import rd_data

class RamsesData:

    def __init__(self,idir=1,type=1,xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1,lmax=5):
        """ Read RAMSES data and return a cube """

        #Read data in file
        dirname='output_'+str_suffix(idir)

        #Cube size parameter
        n=2**lmax
        nx=int((xmax-xmin)*n)
        ny=int((ymax-ymin)*n)
        nz=int((zmax-zmin)*n)

        self.xmin=xmin ; self.xmax=xmax
        self.ymin=ymin ; self.ymax=ymax
        self.zmin=zmin ; self.zmax=zmax

        self.cube,self.time=rd_data.amr2cube(dirname,type,xmin,xmax,\
                                       ymin,ymax,zmin,zmax,lmax,nx,ny,nz)

def str_suffix(n,length=5):
    fewzero=''
    for i in range(length-len(str(n))):
        fewzero=fewzero+'0'
    return fewzero+str(n)
