;+
; NAME:
;	MAKE_DIST
;
; PURPOSE:
;	This procedure computes all the cosmological distances
;       given a redshift.
;
; CATEGORY:
;	Input/Output.
;
; CALLING SEQUENCE:
;       MAKE_DIST, z,omega_m=omega_m,omega_r=omega_r
;                  ,omega_v=omega_v,save=save
;
; INPUTS
;       z: redshift
;
; OPTIONAL INPUTS:
;       OMEGA_M: matter density. Default: 0
;
;       OMEGA_R: radiative density. Default: 0
;
;       OMEGA_V: cosmological constant. Default: 0
;
;       SAVE: if set, record all lengths in structure save
;
; OUTPUTS:
;	NONE
;
; COMMON BLOCKS:
;       None.
;
; EXAMPLE:
;       make_dist,2,omega_m=0.4,omega_v=0.2,omega_r=0.2,save=save
;
; MODIFICATION HISTORY:
; 	Written by:	Matthias Gonzalez, 08/10/2001.
;	     	Comments and header by Matthias Gonzalez.
;-
pro make_dist, z,omega_m=omega_m,omega_r=omega_r $
               ,omega_v=omega_v,save=save

if  N_PARAMS() ne 1 then begin
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'make_dist'
    return
endif

if z eq 0 then begin
    print,'everything equal 0...'
    return
endif

if not keyword_set(omega_m) then omega_m=0.
if not keyword_set(omega_v) then omega_v=0.
if not keyword_set(omega_r) then omega_r=0.

h0     =100.           ; Hubble constant in h km s-1 Mpc-1

omega_k=1.-omega_m-omega_v-omega_r

if omega_k eq 0 then begin
    k=0
endif else if omega_k lt 0 then begin
    k=1
endif else begin
    k=-1
endelse

c=3.e5                  ; ligth speed in km.s-1 
nb_echan=10000.         ; number of points for the integration

if k ne 0 then r0=sqrt(-k*(c/h0)^2/omega_k)

z_array=findgen(nb_echan)/(nb_echan-1.)*z
z_array=z_array(sort(z_array))
f=1./sqrt(omega_v+omega_k*(1.+z_array)^2+omega_m*(1.+z_array)^3+omega_r*(1.+z_array)^4)
d=c/h0*f
temps_array=1./(1.+z_array)/h0*f

lookback=int_tabulated(z_array,temps_array) ; in h-1 s km-1 Mpc
lookback=lookback*3.09e24/1.e5                 ; in h-1 s
lookback=lookback/3600./24./365./10.^9         ; in h-1 Gyr

d_comobile=int_tabulated(z_array,d)
d_propre=d_comobile/(1.+z)

a_array=(findgen(nb_echan)+0.5)/(nb_echan-0.5)*1./(1.+z)
z_array=1./a_array-1.
f=1./sqrt(omega_v+omega_k*(1.+z_array)^2+omega_m*(1.+z_array)^3+omega_r*(1.+z_array)^4)
temps_array=(1.+z_array)/h0*f
temps=int_tabulated(a_array,temps_array)
temps=temps*3.09e24/1.e5                 ; in h-1 s
temps=temps/3600./24./365./10.^9         ; in h-1 Gyr


if k eq 0 then begin
    d_diam=1./(1.+z)*d_comobile
    d_lum =(1.+z)   *d_comobile
endif else if k eq 1 then begin
    d_diam=1./(1.+z)*sin(d_comobile/r0)*r0
    d_lum =(1.+z)   *sin(d_comobile/r0)*r0
endif else if k eq -1 then begin
    d_diam=1./(1.+z)*sinh(d_comobile/r0)*r0
    d_lum =(1.+z)   *sinh(d_comobile/r0)*r0
endif

if not keyword_set(save) then begin
    print,'k =',k
    print,'z =',z
    print,'Look back time            (h-1 Gyr)=',lookback*z/abs(z)
    print,'Age                       (h-1 Gyr)=',temps
    print,'Comoving         distance (h-1 Mpc)=',d_comobile
    print,'Proper           distance (h-1 Mpc)=',d_propre
    print,'Angular-Diameter distance (h-1 Mpc)=',d_diam
    print,'Luminosity       distance (h-1 Mpc)=',d_lum
endif else begin
    print,'writing all distances in structure... '
    save={d_lum:d_lum,d_diam:d_diam,d_propre:d_propre $
          ,d_comobile:d_comobile,temps:temps,lookback:lookback*z/abs(z)}
endelse

end
