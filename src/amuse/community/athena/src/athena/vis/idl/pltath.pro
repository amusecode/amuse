; FILE pltath.pro
;
; PURPOSE: contains IDL procedures to read and make plots from Athena dumps and
;   outputs.  Contains the following procedures:
;
; PRO readbin,filename:        reads .bin file 'filename'
; PRO sod_plot,filename:       plots analytic Sod solution over numerical
; PRO four_plot,filename:      plots d,P,Vx,P/d
; PRO nine_plot,filename,flag: use flag=0 for points, flag=1 for line
; PRO flines,nlev              plot 2D field lines
; PRO readvtk,filename:        reads .vtk file 'filename'
;
COMMON SHARE1,nx,ny,nz,nvar,nscalars
COMMON SHARE2,x,y,z
COMMON SHARE3,time,dt,gamm1,isocs
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz,s,phi

;-------------------------------------------------------------------------------
; Procedure READBIN: Reads ATHENA binary dumps
;
PRO readbin,filename
COMMON SHARE1,nx,ny,nz,nvar,nscalars
COMMON SHARE2,x,y,z
COMMON SHARE3,time,dt,gamm1,isocs
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz,s,phi
openr,1,filename
;
; Read coordsys
;
readu,1,coordsys
;
; Read number of zones and variables
;
ndata=LONARR(7)
readu,1,ndata
nx=ndata[0]
ny=ndata[1]
nz=ndata[2]
nvar=ndata[3]
nscalars=ndata[4]
ngrav=ndata[5]
npart=ndata[6]
;
; Read (gamma-1) and isothermal sound speed
;
dat=fltarr(2)
readu,1,dat
gamm1=dat[0]
isocs=dat[1]
;
; Read time,dt
;
readu,1,dat
time=dat[0]
dt=dat[1]
;
; Read grid coordinates
;
x=fltarr(nx)
readu,1,x
y=fltarr(ny)
readu,1,y
z=fltarr(nz)
readu,1,z
;
;  Allocate arrays.  Note arrays allocated even if not used.
;
d =fltarr(nx,ny,nz)
e =fltarr(nx,ny,nz)
vx=fltarr(nx,ny,nz)
vy=fltarr(nx,ny,nz)
vz=fltarr(nx,ny,nz)
bx=fltarr(nx,ny,nz)
by=fltarr(nx,ny,nz)
bz=fltarr(nx,ny,nz)
phi=fltarr(nx,ny,nz)
nscal = 1
IF nscalars gt 1 THEN nscal = nscalars
s = fltarr(nx,ny,nz,nscal)
;
; Read data.
; nvar=4 means isothermal hydro.  nvar=5 means adiabatic hydro
; nvar=7 means isothermal MHD.    nvar=8 means adiabatic mhd
;
readu,1,d
readu,1,vx
readu,1,vy
readu,1,vz
IF (nvar-nscalars) eq 5 OR (nvar-nscalars) eq 8 THEN readu,1,e
IF (nvar-nscalars) eq 7 OR (nvar-nscalars) eq 8 THEN BEGIN
  readu,1,bx
  readu,1,by
  readu,1,bz
ENDIF
IF (nscalars) gt 0 THEN readu,1,s
IF (ngrav) gt 0 THEN readu,1,phi
;
; compute velocities and pressure
;
vx=vx/d
vy=vy/d
vz=vz/d
IF gamm1 NE 0 THEN p=gamm1*(e-0.5*d*(vx^2+vy^2+vz^2)-0.5*(bx^2+by^2+bz^2))
IF gamm1 EQ 0 THEN p=isocs*isocs*d
;
close,1
END

;-------------------------------------------------------------------------------
; Procedure FOUR_PLOT: plots d,P,Vx,P/d
;
PRO four_plot,filename
COMMON SHARE1,nx,ny,nz,nvar,nscalars
COMMON SHARE2,x,y,z
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz,s,phi
;
readbin,filename
!P.MULTI=[0,2,2,0,0]
dmin=MIN(d)
dmax=MAX(d)
plot,x,d,YTITLE='D',XTITLE='X',psym=6,symsize=.4,YRANGE=[dmin,dmax],XSTYLE=1
dmin=MIN(p)
dmax=MAX(p)
plot,x,p,YTITLE='P',XTITLE='X',psym=6,symsize=.4,YRANGE=[dmin,dmax],XSTYLE=1
dmin=MIN(vx)
dmax=MAX(vx)
plot,x,vx,YTITLE='Vx',XTITLE='X',psym=6,symsize=.4,YRANGE=[dmin,dmax],XSTYLE=1
dmin=MIN(p/d)
dmax=MAX(p/d)
plot,x,p/d,YTITLE='P/D',XTITLE='X',psym=6,symsize=.4,YRANGE=[dmin,dmax],XSTYLE=1
!P.MULTI=0
END

;-------------------------------------------------------------------------------
; Procedure NINE_PLOT: plots d,P,E,Vx,Vy,Vz,By,Bz,Phi - same plot as in RJ
;   Use flag=0 to plot points, flag=1 to plot line
;
PRO nine_plot,filename,flag
COMMON SHARE1,nx,ny,nz,nvar,nscalars
COMMON SHARE2,x,y,z
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz,s,phi
;
readbin,filename
!P.MULTI=[0,3,3,0,0]
dmin=MIN(d)
dmax=MAX(d)
IF (flag EQ 0) THEN plot,x,d,psym=6,symsize=.4,YTITLE='D',XTITLE='X',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,d,YTITLE='D',XTITLE='X',YRANGE=[dmin,dmax]
dmin=MIN(p)
dmax=MAX(p)
IF (flag EQ 0) THEN plot,x,p,psym=6,symsize=.4,YTITLE='P',XTITLE='X',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,p,YTITLE='P',XTITLE='X',YRANGE=[dmin,dmax]
dmin=MIN(e)
dmax=MAX(e)
IF (flag EQ 0) THEN plot,x,e,psym=6,symsize=.4,YTITLE='E',XTITLE='X',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,e,YTITLE='E',XTITLE='X',YRANGE=[dmin,dmax]
dmin=MIN(vx)
dmax=MAX(vx)
IF (flag EQ 0) THEN plot,x,vx,psym=6,symsize=.4,YTITLE='Vx',XTITLE='X',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,vx,YTITLE='Vx',XTITLE='X',YRANGE=[dmin,dmax]
dmin=MIN(vy)
dmax=MAX(vy)
IF (flag EQ 0) THEN plot,x,vy,psym=6,symsize=.4,YTITLE='Vy',XTITLE='X',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,vy,YTITLE='Vy',XTITLE='X',YRANGE=[dmin,dmax]
dmin=MIN(vz)
dmax=MAX(vz)
IF (flag EQ 0) THEN plot,x,vz,psym=6,symsize=.4,YTITLE='Vz',XTITLE='X',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,vz,YTITLE='Vz',XTITLE='X',YRANGE=[dmin,dmax]
dmin=MIN(by)
dmax=MAX(by)
IF (flag EQ 0) THEN plot,x,by,psym=6,symsize=.4,YTITLE='By',XTITLE='X',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,by,YTITLE='By',XTITLE='X',YRANGE=[dmin,dmax]
dmin=MIN(bz)
dmax=MAX(bz)
IF (flag EQ 0) THEN plot,x,bz,psym=6,symsize=.4,YTITLE='Bz',XTITLE='X',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,bz,YTITLE='Bz',XTITLE='X',YRANGE=[dmin,dmax]
phi = 180*(atan(bz/(by+1.0e-10))/3.1415927)
dmin=MIN(phi)
dmax=MAX(phi)
IF (flag EQ 0) THEN plot,x,phi,psym=6,symsize=.4,YTITLE='PHI',XTITLE='X',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,phi,YTITLE='PHI',XTITLE='X',YRANGE=[dmin,dmax]
!P.MULTI=0
END

;------------------------------------------------------------------------------
; Procedure SOD_PLOT: plots analytic solution for Sod shocktube over numerical
;
PRO sod_plot,filename
COMMON SHARE1,nx,ny,nz,nvar,nscalars
COMMON SHARE2,x,y,z
COMMON SHARE3,time,dt,gamm1,isocs
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz,s,phi
;
readbin,filename
vs = 1.7522
vc = 0.92745
vf = 0.07027
vh = 1.1832
xint = x[0] + (x[nx-1]-x[0])/2.0
xs = xint+vs*time
xc = xint+vc*time
xf = xint-vf*time
xh = xint-vh*time
dsod=FLTARR(500)
vsod=FLTARR(500)
esod=FLTARR(500)
psod=FLTARR(500)
xsod=FLTARR(500)
xsod[0] = x[0]
FOR I=1,499 DO xsod[I]=xsod[I-1] + (x[nx-1]-x[0])/500.
FOR I=0,499 DO BEGIN
  IF xsod[I] GT xs THEN BEGIN
    dsod[I] = 0.125
    psod[I] = 0.1
    vsod[I] = 0.0
    esod[I] = 2.0
  ENDIF
  IF xsod[I] GT xc AND xsod[I] LT xs THEN BEGIN
    dsod[I] = 0.26557
    psod[I] = 0.30313
    vsod[I] = 0.92745
    esod[I] = 2.8535
  ENDIF
  IF xsod[I] GT xf AND xsod[I] LT xc THEN BEGIN
    dsod[I] = 0.42632
    psod[I] = 0.30313
    vsod[I] = 0.92745
    esod[I] = 1.7776
  ENDIF
  IF xsod[I] GT xh AND xsod[I] LT xf THEN BEGIN
    vsod[I] = 0.92745*(xsod[I]-xh)/(xf-xh)
    dsod[I] = 0.42632*(1.0+0.20046*(0.92745-vsod[I]))^5
    psod[I] = 0.30313*(1.0+0.20046*(0.92745-vsod[I]))^7
    esod[I] = psod[I]/(0.4*dsod[I])
  ENDIF
  IF xsod[I] LT xh THEN BEGIN
    dsod[I] = 1.0
    psod[I] = 1.0
    vsod[I] = 0.0
    esod[I] = 2.5
  ENDIF
ENDFOR
;
!P.MULTI=[0,2,2,0,0]
dmin=MIN(d)
dmax=MAX(d)
plot,x,d,YTITLE='D',XTITLE='X',psym=6,symsize=.4,YRANGE=[dmin,dmax],XSTYLE=1
oplot,xsod,dsod
dmin=MIN(p)
dmax=MAX(p)
plot,x,p,YTITLE='P',XTITLE='X',psym=6,symsize=.4,YRANGE=[dmin,dmax],XSTYLE=1
oplot,xsod,psod
dmin=MIN(vx)
dmax=MAX(vx)
plot,x,vx,YTITLE='Vx',XTITLE='X',psym=6,symsize=.4,YRANGE=[dmin,dmax],XSTYLE=1
oplot,xsod,vsod
dmin=MIN(p/d)
dmax=MAX(p/d)
plot,x,p/d,YTITLE='P/D',XTITLE='X',psym=6,symsize=.4,YRANGE=[dmin,dmax],XSTYLE=1
oplot,xsod,psod/dsod
!P.MULTI=0
END
;------------------------------------------------------------------------------
; Procedure FLINES:  2D plot of field lines
;
PRO flines,nlev
COMMON SHARE1,nx,ny,nz,nvar,nscalars
COMMON SHARE2,x,y,z
COMMON SHARE3,time,dt,gamm1,isocs
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz,s,phi
vecpot=fltarr(nx,ny)
dx = x[1]-x[0]
dy = y[1]-y[0]
vecpot[0,0] = 0.0
FOR J=1,ny-1 DO vecpot[0,J] = vecpot[0,j-1] + dy*bx[0,j]
FOR I=1,nx-1 DO vecpot[I,*] = vecpot[i-1,*] - dx*by[i,*]
contour,vecpot,x,y,nlevels=nlev,/ISOTROPIC,XSTYLE=1,YSTYLE=1
END

;------------------------------------------------------------------------------
; Procedure READVTK: Reads ATHENA VTK files
;
PRO readvtk,filename,pfact
COMMON SHARE1,nx,ny,nz,nvar,nscalars
COMMON SHARE2,x,y,z
COMMON SHARE3,time,dt,gamm1,isocs
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz,s,phi
; -----------------------------------------
; pfact = gamm1 for adiabatic
; pfact = isocs for isothermal
; leave out parameter to use existing value
; -----------------------------------------
;
; Read header information, which is assumed
; to be (roughly) in the following form:
;
; # vtk DataFile Version 3.0
; Really cool Athena data at time = 0.000000e+00
; BINARY
; DATASET STRUCTURED_POINTS
; DIMENSIONS 257 257 513
; ORIGIN -5.000000e-02 -5.000000e-02 -1.000000e-01
; SPACING 3.906250e-04 3.906250e-04 3.906250e-04
; CELL_DATA 33554432
; SCALARS density float
; LOOKUP_TABLE default
; (array of dim[nx,ny,nz])
; VECTORS velocity float
; (array of dim[3,nx,ny,nz])
; SCALARS total_energy float
; LOOKUP_TABLE default
; (array of dim[nx,ny,nz])
; VECTORS cell-centered-B float
; (array of dim[3,nx,ny,nz])
;
; There are differences in the VTK files output from
; different versions of Athena and also join_vtk_dump!

string = ' '
string_array=STRARR(8)
ndata = LONARR(3)

openr,1,filename
; read line 1 (do nothing)
readf,1,string
readf,1,string
string_array=STRSPLIT(string,' ',count=cnt,/EXTRACT)
if (cnt eq 8) then begin
  reads,string_array[7],time
  print,"Time:", time
end
; read lines 3,4 (do nothing)
readf,1,string
readf,1,string
; read line 5, get dimensions
readf,1,string
string_array=STRSPLIT(string,' ',/EXTRACT)
reads,string_array[1],nxs
reads,string_array[2],nys
reads,string_array[3],nzs
nx = LONG(nxs)-1
ny = LONG(nys)-1
nz = LONG(nzs)-1
print,"nx,ny,nz:",nx,ny,nz
; read line 6, get origin
readf,1,string
string_array=STRSPLIT(string,' ',/EXTRACT)
reads,string_array[1],x0
reads,string_array[2],y0
reads,string_array[3],z0
print,"x0,y0,z0:",x0,y0,z0
; read line 7, get grid spacing
readf,1,string
string_array=STRSPLIT(string,' ',/EXTRACT)
reads,string_array[1],dx
reads,string_array[2],dy
reads,string_array[3],dz
print,"dx,dy,dz:",dx,dy,dz
; read line 8, do nothing
readf,1,string

nvar = 0
isothermal = 1
mhd = 0
while (not eof(1)) do begin
  readf,1,string
  block = matchsechead(string)
  if ((block le 0) and (not eof(1))) then begin
    readf,1,string
    block = matchsechead(string)
  end
  if ((block le 0) and (not eof(1))) then begin
    print,"Unexpected file format"
    stop
  end else begin
    print,string
    if (block lt 20) then begin
      ; SCALARS block
      if (block eq 11) then begin
        readscalarblock,1,d
        nvar = nvar + 1
      end else if (block eq 12) then begin
        readscalarblock,1,e
        isothermal = 0
        nvar = nvar + 1
      end else begin
        readscalarblock,1
        print,"Data from unrecognized SCALARS block not stored."
      end
    end else if (block lt 30) then begin
      ; VECTORS block
      if (block eq 21) then begin
        readvectblock,1,vx,vy,vz
        nvar = nvar + 3
      end else if (block eq 22) then begin
        readvectblock,1,bx,by,bz
        mhd = 1
        nvar = nvar + 3
      end else begin
        readvectblock,1
        print,"Data from unrecognized VECTORS block not stored."
      end
    end else begin
      print,"Unrecognized block type!"
      stop
    end
  end
endwhile

close,1

if (isothermal eq 0) then begin
  if (mhd eq 0) then begin
    print,"Assuming adiabatic hydro"
  end else begin
    print,"Assuming adiabatic MHD"
  end
end else begin
  if (mhd eq 0) then begin
    print,"Assuming isothermal hydro"
  end else begin
    print,"Assuming isothermal MHD"
  end
end

if (isothermal eq 0) then begin
  if (N_Elements(pfact) eq 1) then gamm1 = pfact
  if (N_Elements(gamm1) eq 0) then begin
    print,"Pressure not set because gamm1 undefined!"
  end else begin
    print,"Pressure set assuming gamm1=",gamm1
    if (mhd eq 0) then begin
      p=gamm1*(e-0.5*d*(vx^2+vy^2+vz^2))
    end else begin
      p=gamm1*(e-0.5*d*(vx^2+vy^2+vz^2)-0.5*(bx^2+by^2+bz^2))
    end
  end
end else begin
  if (N_Elements(pfact) eq 1) then isocs = pfact
  if (N_Elements(isocs) eq 0) then begin
    print,"Pressure not set because isocs undefined!"
  end else begin
    print,"Pressure set assuming isocs=",isocs
    p=isocs*isocs*d
  end
end

END
;------------------------------------------------------------------------------
FUNCTION matchsechead,string

if (string eq "SCALARS density float") then begin
  return,11
end else if (string eq "VECTORS velocity float") then begin
  return,21
end else if (string eq "SCALARS total_energy float") then begin
  return,12
end else if (string eq "VECTORS cell_centered_B float") then begin
  return,22
end else if (strpos(string, "SCALARS") ne -1) then begin
  return,10
end else if (strpos(string, "VECTORS") ne -1) then begin
  return,20
end else if (string eq "") then begin
  return,0
end else begin
  return,-1
end

END
;------------------------------------------------------------------------------
PRO readvectblock,fp,vectx,vecty,vectz
COMMON SHARE1,nx,ny,nz,nvar,nscalars
; VECTORS block

var=fltarr(3,nx,ny,nz)
readu,fp,var
var=SWAP_ENDIAN(temporary(var),/SWAP_IF_LITTLE_ENDIAN)

vectx=temporary(reform(var[0,*,*,*]))
vecty=temporary(reform(var[1,*,*,*]))
vectz=temporary(reform(var[2,*,*,*]))

END
;------------------------------------------------------------------------------
PRO readscalarblock,fp,var
COMMON SHARE1,nx,ny,nz,nvar,nscalars
; SCALARS block

string = ' '
readf,fp,string ; "LOOKUP_TABLE default"
var=fltarr(nx,ny,nz)
readu,fp,var
var=SWAP_ENDIAN(temporary(var),/SWAP_IF_LITTLE_ENDIAN)

END
;------------------------------------------------------------------------------
