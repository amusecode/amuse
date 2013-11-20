;+
; NAME:
;       TUBE
;
; PURPOSE:
;       This procedure reads 1D profiles from a RAMSES LOG ASCII file
;       for the Sod test (shock tube test) and plot the results versus
;       the analytical solution.
;
; CATEGORY:
;       Plotting hydro test case.
;
; CALLING SEQUENCE:
;       TUBE,file=file,ana=ana,xr=xr,symsize=symsize
;
; OPTIONAL INPUTS:
;       file:   if set, input the scalar string containing the name of
;       the file to be read. Otherwise, a PICKFILE widget is launched. 
;
;       ana:    if set, input the scalar string containing the name of
;       the file to be read for the analytical solution.
;       Default: '$HOME/ramses/idl/tubeana.dat'
;       The file should be an ASCII file with N lines containing in 
;       each column: zone number, position, velocity, density,
;       pressure and internal energy
;
;       xr:     the X axis range. Default: [0,1]
;
;       symsize: the size of the square symbol to plot the numerical
;       solution. Default=1.
;
; OUTPUTS:
;       None.
;;
; EXAMPLE:
;       To plot shock tube profiles from a RAMSES LOG ASCII file, type:
;
;               TUBE,file='tube.log',ana='mytubeana.dat'
;
; MODIFICATION HISTORY:
;       Written by:     Romain Teyssier, 01/01/2000.
;                       e-mail: Romain.Teyssier@cea.fr
;       Fevrier, 2001:  Comments and header added by Romain Teyssier.
;-
;###################################################
;###################################################
;###################################################
pro tube,file=file,ana=ana,xr=xr,symsize=symsize

if not keyword_set(file) then file=pickfile(/READ)
if not keyword_set(ana) then ana='$HOME/ramses/idl/tubeana.dat'
if not keyword_set(symsize) then symsize=1.
if not keyword_set(xr) then xr=[0.,1.]

ss=symsize
tek_color
rd_1d,output,file=file
ntime=n_elements(output)
readcol,ana,ia,ra,ua,da,pa,ea,/silent

!p.multi=[0,2,2]
plot,(*output[ntime-1]).x,(*output[ntime-1]).u $
  ,color=1,xr=xr,symsize=ss,xtitle='!17x',ytitle='u',psym=6
oplot,ra,ua,color=2

!p.multi=[1,2,2]
plot,(*output[ntime-1]).x,(*output[ntime-1]).d $
  ,color=1,xr=xr,symsize=ss,xtitle='x',ytitle='!7q!17',psym=6
oplot,ra,da,color=2

!p.multi=[2,2,2]
plot,(*output[ntime-1]).x,(*output[ntime-1]).p $
  ,color=1,xr=xr,symsize=ss,xtitle='x',ytitle='P',psym=6
oplot,ra,pa,color=2

!p.multi=[3,2,2]
plot,(*output[ntime-1]).x,(*output[ntime-1]).l $
  ,xr=xr,psym=10,color=1,symsize=ss,xtitle='x',ytitle='Level'

end
;###################################################
;###################################################
;###################################################
