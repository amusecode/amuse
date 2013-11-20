;+
; NAME:
;	MAKE_PHOTON
;
; PURPOSE:
;       This procedure computes the number of photons per second
;       received by a detector and computes a realization.
;
; CATEGORY:
;	Input/Output.
;
; CALLING SEQUENCE:
;       MAKE_PHOTON $
;              , amrfile=amrfile, hydrofile=hydrofile $
;              , xcenter=xcenter, ycenter=zcenter, zcenter=zcenter $
;              , scale_factor=scale_factor, levelmax=levelmax $
;              , first_bin=first_bin, last_bin=last_bin $
;              , degrade=degrade
;
; INPUTS
;       None
;
; OPTIONAL INPUTS:
;       AMRFILE:   string array the AMR filename to be read. If not
;                  specified, a PICKFILE widget is launched.
;
;       HYDROFILE: string array the HYDRO filename to be read. If not
;                  specified, a PICKFILE widget is launched.
;
;       XCENTER:   set the x coordinate of the center. Default: true
;                  center of AMR grid. 
;
;       YCENTER:   set the y coordinate of the center. Default: true
;                  center of AMR grid. 
;
;       ZCENTER:   set the z coordinate of the center. Default: true
;                  center of AMR grid. 
;
;       LEVELMAX:  if set, specifies the maximum level to be shown.
;                  Default: the levelmax of amr file.
;
;       FIRST_BIN: integer. first bin of energy. Default:0
;
;       LAST_BIN:  integer. last bin of energy. Default:4094
;
;       DEGRADE:  integer. number of old bin per new bin. Default:1
;
;       SCALE_FACTOR: float. number to rescale hydro variables.
;
; OUTPUTS:
;	NONE
;
; COMMON BLOCKS:
;       EMISSIVITY_TYPE
;       MEWE
;
; EXAMPLE:
;       make_photon,xcenter=2.,ycenter=2.,zcenter=2.,lev=3  $
;            ,first_bin=1,last_bin=36,degrade=11,scale_factor=0.6
;
; MODIFICATION HISTORY:
; 	Written by:	Matthias Gonzalez, 26/09/2001.
;         	     	Comments and header by Matthias Gonzalez.
;       20-11-01:       Routine name changed by R. Teyssier
;       20-11-01:       Some keywords added by R. Teyssier
;-
;*******************************************************************
;*******************************************************************
;*******************************************************************
pro make_photon $
           , amrfile=amrfile, hydrofile=hydrofile $
           , xcenter=xcenter,ycenter=ycenter,zcenter=zcenter $
           , levelmax=levelmax $
           , first_bin=first_bin, last_bin=last_bin $
           , degrade=degrade $
           , scale_factor=scale_factor

common emissivity_type, e_type
common mewe,ok_read,spec,ener,bin_number,delta_hnu,t_min,t_max,delta_t

;*******************************************************************
;  PARAMETRES DU PROGRAMME EVENTUELLEMENT CHANGEABLES
;*******************************************************************

; Parametres spectraux
filespectre='/home/storage/teyssier/Spec_DB_image_60kev.fits'
nphoton_max=2500000L
random_seed=123L
random_seed2=256L
surface_ccd=6000.             ; in cm2
t_obs=30.d3                   ; observation time in s
z_amas=0.2                    ; redshift of the cluster

; Sous-cube en unites grille "coarse"
xr=[1,3]
yr=[1,3]
zr=[1,3]

; Parametres physiques
omegam=0.3          ; omega matiere
omegav=0.7          ; omega vide
hubble=0.7          ; constante de hubble
rhoc=1.88d-29       ; h2 g cm-3
alpha=0.28          ; fraction of He in mass

;*******************************************************************
;  FIN DES PARAMETRES CHANGEABLES
;*******************************************************************

if not keyword_set(degrade)   then degrade=1
if not keyword_set(first_bin) then first_bin=0
if not keyword_set(last_bin)  then last_bin=-1
bin_path=degrade
if not keyword_set(xcenter) then xcenter=(xr[1]+xr[0])/2.0d0
if not keyword_set(ycenter) then ycenter=(yr[1]+yr[0])/2.0d0
if not keyword_set(zcenter) then zcenter=(zr[1]+zr[0])/2.0d0
x0in=xcenter
y0in=ycenter
z0in=zcenter

mu=4./(8.-5.*alpha) ; mean molecular weight for ntot
mue=4./(4.-2.*alpha); mean molecular weight for ne
muh=1./(1.-alpha)   ; mean molecular weight for nh
; Physical constants
mh=1.67d-24         ; proton mass in g
me=9.11d-28         ; electron mass in g 
kb=1.38d-16         ; Boltzmann constant in erg K-1
sigmat=0.66d-24     ; Thomson cross section en cm2
c=3d10              ; Light speed in cm s-1

; Reading AMR and HYDRO data
if not keyword_set(amrfile) then begin
    amrfile=DIALOG_PICKFILE(/READ,filter='*amr*',get_path=path,title='Please select the AMR file')
endif
if not keyword_set(hydrofile) then begin
    hydrofile=DIALOG_PICKFILE(/READ,filter='*hydro*',path=path,title='Please select the HYDRO file')
endif
fileamr=amrfile
filehydro=hydrofile
print,''
rd_amr,a,file=fileamr
print,''       
rd_hydro,a,h,file=filehydro
print,''       

;Rescaling
if keyword_set(scale_factor)  then begin
    print,'Rescaling AMR and HYDRO variables to new box size'
    a.boxlen=a.boxlen*scale_factor
    print,'New box size (h-1 Mpc)=',a.boxlen
    print,' '
    (*h.levelh[0]).u(*,*,*,1)=(*h.levelh[0]).u(*,*,*,1)*scale_factor
    (*h.levelh[0]).u(*,*,*,2)=(*h.levelh[0]).u(*,*,*,2)*scale_factor
    (*h.levelh[0]).u(*,*,*,3)=(*h.levelh[0]).u(*,*,*,3)*scale_factor
    (*h.levelh[0]).u(*,*,*,4)=(*h.levelh[0]).u(*,*,*,4)*scale_factor^2
    (*h.levelh[0]).u(*,*,*,5)=(*h.levelh[0]).u(*,*,*,5)*scale_factor^2
    for i=1,a.nlevelmax-1 do begin
        (*h.levelh[i]).u(*,*,1)=(*h.levelh[i]).u(*,*,1)*scale_factor
        (*h.levelh[i]).u(*,*,2)=(*h.levelh[i]).u(*,*,2)*scale_factor
        (*h.levelh[i]).u(*,*,3)=(*h.levelh[i]).u(*,*,3)*scale_factor
        (*h.levelh[i]).u(*,*,4)=(*h.levelh[i]).u(*,*,4)*scale_factor^2
        (*h.levelh[i]).u(*,*,5)=(*h.levelh[i]).u(*,*,5)*scale_factor^2
    endfor
endif

; Parametres de la simulation
lbox=a.boxlen      ; longueur de la grille "coarse" en h-1 Mpc
ngridcoarse=a.nx   ; nombre de cellules "coarse" dans une direction

; Read Mewe data
READ_MEWE2,filespectre,bin_path,first_bin,last_bin,/old

; Other parameters
gamma=h.gamma
aexp=a.aexp
if not keyword_set(level) then level=a.nlevelmax-1

; Build a list of all AMR cells
print,' '
amr2cell,a,h,cell,xr=xr,yr=yr,zr=zr,levelmax=level

n   = long  (cell.n)
dx  = double(cell.dx)
x   = double(cell.x)
y   = double(cell.y)
z   = double(cell.z)
p   = double(cell.var(*,5))
rho = double(cell.var(*,0))

; Cleaning memory
print,' '
del_amr,a
del_hydro,h
cell=0   

t = mu*mh/kb*p/rho*1.d14    ; Temperature en K
t = t/11604.5d0/1000.d0     ; Temperature en keV
print,' '
print,'    minmax rho    =',minmax(rho)
print,'    minmax T (keV)=',minmax(t)

indinf=where(t lt t_min)
indsup=where(t gt t_max)
if indinf(0) ne -1 then t(indinf)=t_min
if indsup(0) ne -1 then t(indsup)=t_max
indkev=fix((t-t_min)/delta_t)
print,'New minmax T (keV)=',minmax(t)

dist=1
print,' '
make_dist,z_amas,omega_m=omegam,omega_v=omegav,save=dist
print,'distance angulaire  (h-1 Mpc)=',dist.d_diam
print,'distance luminosite (h-1 Mpc)=',dist.d_lum

dxMpc=dx*lbox/ngridcoarse
dxcm =3.087d24*dxMpc
lMpc =double(xr(1)-xr(0))*lbox/ngridcoarse  ; Longueur de projection en h-1Mpc
lcm  =3.087d24*lMpc                         ; Longueur de projection en h-1cm
domegapix=(dxMpc/dist.d_lum)^2/(4.0d0*!DPi) ; Cell solid angle

; Allocate arrays for photon list
nphoton_tot=0L
dx_photon=DBLARR(nphoton_max)
x_photon =DBLARR(nphoton_max)
y_photon =DBLARR(nphoton_max)
z_photon =DBLARR(nphoton_max)
e_photon =DBLARR(nphoton_max)

scale_emissivity=aexp*dxcm/(mue*muh)*(omegam*rhoc/aexp^3/mh)^2
scale_emissivity=scale_emissivity*hubble^3
scale_emissivity=scale_emissivity*domegapix
scale_emissivity=scale_emissivity*surface_ccd
scale_emissivity=scale_emissivity*t_obs

rcell = sqrt((x-x0in)^2+(y-y0in)^2+(z-z0in)^2)

for i=first_bin,last_bin do begin

    bin_number=i
    print,''
    print,'Energy bin number',bin_number*bin_path
    print,'Computing emissivity'
  
    ; Emissivity in this bin
    emtot = rho^2*double(spec(indkev,bin_number)) ;    photons
    emtot = scale_emissivity*emtot

    ; Generate Poisson deviate
    ind=where(emtot gt 1.d-7, nok)

    print,'Keep ',nok,' out of ',n,' cells'

    print,'Minmax photons s-1 in the image =',minmax(emtot)/t_obs
    print,'Total  photons s-1 in the image =',double(total(emtot))/t_obs
    
    for icell=0L,nok-1L do begin

        nphotons=LONG(randomu(random_seed,POISSON=emtot(ind(icell)),/double))

        if (nphoton_tot+nphotons) gt nphoton_max then begin
            print,'Number of photons too large'
            print,'Increase nphoton_max'
            return
        endif

        if nphotons gt 0L then begin
            dx_photon(nphoton_tot:nphoton_tot+nphotons-1L)=dx(ind(icell))
            x_photon (nphoton_tot:nphoton_tot+nphotons-1L)=x (ind(icell))
            y_photon (nphoton_tot:nphoton_tot+nphotons-1L)=y (ind(icell))
            z_photon (nphoton_tot:nphoton_tot+nphotons-1L)=z (ind(icell))
            e_photon (nphoton_tot:nphoton_tot+nphotons-1L)=ener(bin_number)
            nphoton_tot=nphoton_tot+nphotons
        endif

    endfor
    print,'Number of photons in the CCD so far    =',nphoton_tot
    
endfor
print,'Final number of photons =',nphoton_tot
dx_photon=dx_photon(0:nphoton_tot-1)
x_photon =x_photon (0:nphoton_tot-1)
y_photon =y_photon (0:nphoton_tot-1)
z_photon =z_photon (0:nphoton_tot-1)
e_photon =e_photon (0:nphoton_tot-1)

rrr=randomu(random_seed2,nphoton_tot)-0.5d0
x_photon=x_photon+dx_photon*rrr
rrr=randomu(random_seed2,nphoton_tot)-0.5d0
y_photon=y_photon+dx_photon*rrr
rrr=randomu(random_seed2,nphoton_tot)-0.5d0
z_photon=z_photon+dx_photon*rrr
rrr=randomu(random_seed2,nphoton_tot)-0.5d0
e_photon=e_photon+(bin_path*delta_hnu/1.d3)*rrr

x_photon =(x_photon-2.0d0)/2.0d0*lMpc/dist.d_diam*180.d0/!DPI*60.d0
y_photon =(y_photon-2.0d0)/2.0d0*lMpc/dist.d_diam*180.d0/!DPI*60.d0
z_photon =(z_photon-2.0d0)/2.0d0*lMpc/dist.d_diam*180.d0/!DPI*60.d0

print,minmax(x_photon)
print,minmax(y_photon)
print,minmax(z_photon)
print,minmax(e_photon)

; Write photon list to file
mwrfits,x_photon,'photon_x_coor.fits',/create
mwrfits,y_photon,'photon_y_coor.fits',/create
mwrfits,z_photon,'photon_z_coor.fits',/create
mwrfits,e_photon,'photon_energy.fits',/create
spawn,'ls -als photon*.fits'

; Cleaning memory before exiting
e_photon=0
x_photon=0
y_photon=0
z_photon=0
dx_photon=0
rrr=0

end
