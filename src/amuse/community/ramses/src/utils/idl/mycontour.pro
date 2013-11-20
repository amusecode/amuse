;+
; NAME:
;	MYCONTOUR
;
; PURPOSE:
;	This procedure performs a contour plot using the FILL method.
;	The color table is plotted on the right side with corresonding
;	data values. The input image can have different X and Y
;	elements, but the X range and Y range of the resulting contour
;	plot are always equal (defining a square region).
;
; CATEGORY:
;	Plotting routines.
;
; CALLING SEQUENCE:
;       MYCONTOUR, Image, X, Y, LBOX=lbox, TITLE=title, CLT=clt,
;                         UNIT=unit, LOG=log,
;                         NOERASE=noerase, NCONTOUR=ncontour
;
; INPUTS:
;       Image = A 2D rectangular image.
;
; OPTIONAL INPUTS:
;
;       X:     If set, defines the X axis mesh points
;     
;       Y:     If set, defines the Y axis mesh points
;
;       LOG:   If set, contours use a logarithmic color coding.
;              Default: linear color coding.
;
;       LBOX:   If set, defines the size of the square region. The X
;       and Y box length are supposed to be equal. Default: 1.
;
;       UNIT:   Character string. The name of the unit length in the
;       plot. Used to define the X title and Y title. Default:
;       '(arbitrary units)'. 
;
;       TITLE:  Title of the plot. Default: 'No title'.
;       NOERASE: standard IDL option.
;
;       NCONTOUR:  Number of contours. Default: 10.
;
; OUTPUTS:
;       None.
;       
; COMMON BLOCKS:
;       None.
;
; EXAMPLE:
;       Generate an image:
;                   
;                   image = DIST(400)
; 
;       Use MYCONTOUR, specifying the box units and the title:
; 
;                   mycontour, image, lbox=100., unit='(cm)'
;
; MODIFICATION HISTORY:
; 	Written by:	Romain Teyssier, 01/01/2000.
;                       e-mail: Romain.Teyssier@cea.fr
;	Fevrier, 2001:	Comments and header added by Romain Teyssier.
;-
;###################################################
;###################################################
;###################################################
pro mycontour,image,x,y,lbox=lbox,title=title,clt=clt $
              ,unit=unit,log=log, verbose=verbose, table=table $
              ,noerase=noerase,ncontour=ncontour $
              ,minval=minval, maxval=maxval $
              ,xtitle=xtitle, ytitle=ytitle, isotropic=isotropic, reverse=reverse

IF N_PARAMS() NE 1 AND N_PARAMS() NE 3 THEN BEGIN
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'mycontour'
    RETURN
ENDIF

if not keyword_set(title) then title='No title'
if not keyword_set(lbox) then lbox=1.
if not keyword_set(unit) then unit='(arbitrary units)'
if not keyword_set(xtitle) then xtitle='!17x '+unit
if not keyword_set(ytitle) then ytitle='!17z '+unit

; Reset possible 3D transformation
t3d,/reset

; Size of image
ext=double(image)
ss=size(ext)
nx=ss(1)
ny=ss(2)
IF N_PARAMS() EQ 1 THEN BEGIN
    x=(0.5+FINDGEN(nx))/float(nx)*lbox-lbox/2.
    y=(0.5+FINDGEN(ny))/float(ny)*lbox-lbox/2.
ENDIF
maxi=max(ext)
mini=min(ext)
if keyword_set(log) then begin
    ind_plus=where(ext gt 0.d0,n_plus)
    if(n_plus gt 0)then begin
        mini=min(ext(ind_plus))
    endif
    ind_minus=where(ext le mini,n_minus)
    if(n_minus gt 0)then begin
        ext(ind_minus)=mini+0.d0*ext(ind_minus)
    endif
endif
if keyword_set(verbose) then begin
    print,nx,ny    ,format="('Image size  : nx=',I4,' ny=',I4)"
    print,mini,maxi,format="('Image values: min=',E10.3,' max=',E10.3)"
endif

; Set min and max
if keyword_set(minval) then begin
    ind_min=where(ext lt minval, n_min)
    if(n_min gt 0)then ext(ind_min)=minval
    mini=minval
    ind_min=0
endif
if keyword_set(maxval) then begin
    ind_max=where(ext gt maxval, n_max)
    if(n_max gt 0)then ext(ind_max)=maxval
    maxi=maxval
endif

; Number of contours
if not keyword_set(ncontour) then ncontour=10
ncol=ncontour+1
ncolm1=ncol-1

; Log scaling
if keyword_set(log) then begin
    
    pas=alog10(maxi/mini)/ncolm1
    niv=mini*10^(pas*findgen(ncol))
endif else begin
    pas=(maxi-mini)/ncolm1
    niv=findgen(ncol)*pas+mini
endelse

; Colors
if keyword_set(clt) then loadct,clt
n_colors=MIN([!d.n_colors,256])
tvlct,255,255,255,n_colors-1
tvlct,0,0,0,0
cmax=n_colors-20.

clr=(findgen(ncol)+1.)/float(ncol)*cmax+15
if keyword_set(reverse) then begin
   print,'Don''t forget to .r reverse'
   clr2=reverse(clr)
endif else begin
   clr2=clr
endelse

clr=clr2
; Color table versus values
a=niv
c=fltarr(3,ncol)
b=findgen(3)
for j=0,2 do c(j,*)=a(*)

; Multiple windows
nxscale=MAX([!p.multi[1],1L])
nyscale=MAX([!p.multi[2],1L])
xscale=1./float(nxscale)
yscale=1./float(nyscale)
if !p.multi[0] gt 0 then noerase=1
indwin=(!p.multi[0] + (nxscale*nyscale)-1) mod (nxscale*nyscale)
joffset=indwin / nxscale
ioffset=indwin - nxscale*joffset
sizewin=[xscale,yscale,xscale,yscale]
offswin=[ioffset,joffset,ioffset,joffset]*sizewin
!p.multi[0]=(!p.multi[0]-1) mod (nxscale*nyscale)

; Plot color table
if keyword_set(table) then begin
    posleg=[0.9,0.1,0.95,0.9]*sizewin+offswin
    contour,/follow,/fill,levels=niv,c,b,a $
      ,title='!17Value/Color' $
      ,xr=[0,b[1]],yr=[min(a),max(a)],/xs,/ys,ticklen=0. $
      ,xticks=1,xtickn=[' ',' '],c_colors=clr $
      ,position=posleg,ylog=log,noerase=noerase
    
; Plot contours
    pos=[0.1,0.1,0.8,0.8]*sizewin+offswin
endif else begin
    max_x=0.9
    max_y=0.9
    if keyword_set(isotopic) then begin
        if(nx gt ny)then begin
            max_y=0.1+(max_x-0.1)*double(ny)/double(nx)
        endif else begin
            max_x=0.1+(max_y-0.1)*double(nx)/double(ny)
        endelse
    endif
    pos=[0.1,0.1,max_x,max_y]*sizewin+offswin
endelse
contour,/follow,/fill,ext,x,y,levels=niv $
  ,xr=[min(x),max(x)],yr=[min(y),max(y)],/xs,/ys,c_colors=clr $
  ,position=pos,title='!17'+title $
  ,xtitle=xtitle,ytitle=ytitle,/noerase

end
