;##################################################
;##                 ROUTINE mytv.pro             ##
;## For problems, contact Romain.Teyssier@cea.fr ##
;##################################################
pro mytv,im,vmax=vmax,vmin=vmin,save=save,noadjust=noadjust
if not keyword_set(vmax) then vmax=max(im)
if not keyword_set(vmin) then vmin=min(im)
vmax1=vmax
vmin1=vmin
if keyword_set(noadjust) then begin
image=bytscl(im,min=vmin1,max=vmax1)
endif else begin
image=bytscl(im,top=MIN([!d.table_size,256])-20,min=vmin1,max=vmax1)+byte(15)
endelse
ssim=size(im)
xeff=ssim(1)
yeff=ssim(2)

if not keyword_set(save) then begin
    xmax=!d.x_size
    ymax=!d.y_size
endif else begin
    ssz=size(save)
    xmax=ssz(1)
    ymax=ssz(2)
endelse
xscale=double(xmax)/double(xeff)
yscale=double(ymax)/double(yeff)
lscale=min([xscale,yscale])
;imtv=CONGRID(image,lscale*xeff,lscale*yeff)
imtv=REBIN(image,lscale*xeff,lscale*yeff,/sample)

if not keyword_set(save) then begin
    tv,imtv
endif else begin
    save=imtv
endelse
end
