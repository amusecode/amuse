;+
; NAME:
;	PROJVOL3D
;
; PURPOSE:
;       This function returns a two dimensional image that is the
;       projection of a 3-D volume of data onto a plane (similar to an
;       X-ray). The returned image is a translucent rendering of the
;       volume (the highest data values within the volume show up as the
;       brightest regions in the returned image). Depth queing and
;       opacity may be used to affect the image. The volume is
;       projected using a 4x4 matrix, so any type of projection may
;       be used including perspective. Typically the system viewing
;       matrix (!P.T) is used as the 4x4 matrix.
;
;       This routine is largely inspired from the IDL routine
;       PROJECT_VOL. 
;
; CATEGORY:
;	Volume Rendering.
;
; CALLING SEQUENCE:
;       Image = PROJVOL3D(Vol, NX_ray, NY_ray, NZ_plane)
;
; INPUTS:
;       Vol:        A structure containing the three dimensional
;                   volume of data to project. This structure can be
;                   generated using function MESHIJK.
;       NX_ray:     The number of rays to project along the X dimension
;                   of the image. (The returned image will have the
;                   dimensions NX_ray by NY_ray).
;                   Data type : Long.
;       NY_ray:     The number of rays to project along the Y dimension
;                   of the image.
;                   Data type : Long.
;       NZ_plane:   The number of samples to take along each ray.
;                   Higher values for NX_ray, NY_ray, and NZ_plane
;                   increase the image resolution as well as execution time.
;                   Data type : Long.
;
; KEYWORD PARAMETERS:
;       SHOW_PROGRESS:
;                   If set, the function plots to screen the image
;                   obtained during the ray-tracing process, every
;                   SHOW_PROGRESS planes. A window should be open for
;                   this option to work.
;       AVERAGE:    If set, then the average intensity method of projection
;                   is used.   The default is a maximum intensity projection.
;       OPACITY:    A 3-D array with the same size and dimensions as Vol.
;                   This array specifies the opacity of each cell in the
;                   volume. Opaque values of 0 allow all light to
;                   pass through. Opaque values are cumulative.
;                   For example, if a ray eminates from a data value of 50,
;                   and then passes through 10 opaque cells (each with a
;                   data value of 0 and an opacity value of 5) then that
;                   ray would be completely blocked out (the cell with the
;                   data value of 50 would be invisible on the returned
;                   image).   The default is no opacity.
;                   Data type : Any 3-D array except string or structure
;                               (usually the same type as Vol).
;       XSIZE:      The X size of the image to return.   Congrid is used to
;                   resize the final image to be XSIZE by YSIZE.   The default
;                   is the X size of the current window.   If there is
;                   no current window then the default is NX_ray.
;                   Data type: Int or Long.
;       YSIZE:      The Y size of the image to return.   Congrid is used to
;                   resize the final image to be XSIZE by YSIZE.   The default
;                   is the Y size of the current window.   If there is
;                   no current window then the default is NY_ray.
;                   Data type: Int or Long.
;
; OUTPUTS:
;       This function returns the projected volume as a two dimensional
;       array. The dimensions of the returned array are XSIZE by YSIZE.
;
; EXAMPLE:
;       Read some AMR grid data:
;
;               rd_amr,a
;
;       Read some hydro data:
;
;               rd_hydro,a,h
;
;       Use "PP_AMR3D" with the NODATA keyword to set up a viewing
;       projection without plotting the mesh: 
;
;               pp_amr3d,a,/col,ax=30,ay=30,scale=2,persp=2,/nodata
;
;       Use "MESHIJK" to extract a cartesian grid covering the chosen
;       volume, and for the chosen variable:
;
;               vol = meshijk(a,h,xr=[0,32],yr=[0,32],zr=[0,32],lev=3)
;               
;       Render the volume of data using "PROJVOL3D":
;       
;               image = projvol3d(vol,256,256,256)
;
;       Plot the image:
;
;               tvscl,image
;
; MODIFICATION HISTORY:
;       ROUTINE PROJECT_VOL written by Daniel Carr 01/09/1992 
;
;       ROUTINE PROJVOL3D written by Romain Teyssier 17/09/2001
;                                    Romain.Teyssier@cea.fr
;
;-
FUNCTION My_Tv_Loc,img,back,ind_back,xsize,ysize
n_colors=MIN([!d.n_colors,256])
img2 = Congrid(img, xsize, ysize, /Interp, /Minus_One)
tmp_img = 0*BYTSCL(img2)
tmp_img(ind_back)=BYTSCL(alog10(img2(ind_back)+1d-5),TOP=n_colors-15)
ind_scale=where(tmp_img > 0, n_non_zero)
if n_non_zero gt 0 then tmp_img(ind_scale)=tmp_img(ind_scale)+13
tv, back + tmp_img
img2=0
tmp_img=0
ind_scale=0
return,1
END

FUNCTION ProjVol3D, vol, NX_ray, NY_ray, NZ_plane $
         , Opacity=opacity, Average=average $
         , Xsize=xsize, Ysize=ysize $
         , show_progress=show_progress

; *** Test inputs.
size_vol = Size(vol.data)
vol_type = size_vol[size_vol[0]+1]
xr=vol.xr
yr=vol.yr
zr=vol.zr
x_sample = Long(NX_ray[0])
y_sample = Long(NY_ray[0])
z_sample = Long(NZ_plane[0])
xy_sample = x_sample * y_sample
zf_sample = Float(z_sample)
zf_sample_m1 = zf_sample - 1.0
z_sample_m1 = z_sample - 1L
z_max = Float(z_sample - 1L)

; Compute final image size
IF (n_elements(xsize) LE 0L) THEN BEGIN
   IF (!D.Window GE 0L) THEN xsize = !D.X_Size ELSE xsize = x_sample
ENDIF

IF (n_elements(ysize) LE 0L) THEN BEGIN
   IF (!D.Window GE 0L) THEN ysize = !D.Y_Size ELSE ysize = y_sample
ENDIF

; Check for opacity array
block_out = 0B
IF (N_Elements(opacity) GT 0L) THEN BEGIN
   IF (N_Elements(opacity) EQ N_Elements(vol.data)) THEN BEGIN
      opacity = Reform(opacity, size_vol[1], size_vol[2], size_vol[3])
      block_out = 1B
   ENDIF ELSE BEGIN
      Print, 'Opacity array must be the same size as volume data array.'
      RETURN, Fltarr(xsize, ysize)
   ENDELSE
ENDIF

; Invert 3D rotation matrix
trans = !P.T
trans = Invert(trans, status)
IF (status NE 0) THEN BEGIN
   Print, 'Unable to invert transformation matrix.'
   Print, 'You must set the 3D view using "PP_AMR3D" first.'
   RETURN, Fltarr(xsize, ysize)
ENDIF

; Data boundaries
xmin=xr[0] & xmax=xr[1]
ymin=yr[0] & ymax=yr[1]
zmin=zr[0] & zmax=zr[1]

; Minimize number of planar cuts
; Define corners
x_corners=[xmin,xmin,xmin,xmax,xmin,xmax,xmax,xmax] 
y_corners=[ymin,ymin,ymax,ymin,ymax,ymin,ymax,ymax]
z_corners=[zmin,zmax,zmin,zmin,zmax,zmax,zmin,zmax]

; Convert to normalized coordinates
xc_norm=!x.s[0]+!x.s[1]*x_corners
yc_norm=!y.s[0]+!y.s[1]*y_corners
zc_norm=!z.s[0]+!z.s[1]*z_corners

; Apply rotation matrix
index=[[xc_norm[*]],[yc_norm[*]],[zc_norm[*]],[1,1,1,1,1,1,1,1]] # !p.t
zcut_min=MIN(index(*,2)/index(*,3))
zcut_max=MAX(index(*,2)/index(*,3))

; Save background of image
if keyword_set(show_progress) then begin
    IF (!D.Window GE 0L) then begin
        back=tvrd()
        ind_back=where(back eq 0)
    endif else begin
        show_progress=0B
    endelse
endif

; Get min and max of data
vol_max=MAX(vol.data)
vol_min=MIN(vol.data)

; Ray tracing
print,'Starting ray tracing...'
img = Fltarr(x_sample,y_sample)

for i=0L,z_sample-1L do begin

    zcut=zcut_min+float(i)/Float(z_sample)*(zcut_max-zcut_min)
    
; Generate plane in normalized coordinates
    x_rays_norm=(Findgen(x_sample)+0.5)/Float(x_sample) # Replicate(1.0, y_sample)
    y_rays_norm=Replicate(1.0,x_sample)  # (Findgen(y_sample)+0.5)/Float(y_sample)
    z_rays_norm=Replicate(zcut, xy_sample)
    
; Rotate plane using inverse rotation matrix
    index=[[x_rays_norm[*]],[y_rays_norm[*]],[z_rays_norm[*]],[replicate(1.0, xy_sample)]] # trans
    x_rays_rot = index[*, 0] / index[*, 3]
    y_rays_rot = index[*, 1] / index[*, 3]
    z_rays_rot = index[*, 2] / index[*, 3]
    
; Convert to data coordinates
    x_rays_data=(x_rays_rot-!x.s[0])/!x.s[1]
    y_rays_data=(y_rays_rot-!y.s[0])/!y.s[1]
    z_rays_data=(z_rays_rot-!z.s[0])/!z.s[1]
    
; Convert to index in data cube
    indx = (x_rays_data-xmin)/(xmax-xmin)*Float(size_vol[1])
    indy = (y_rays_data-ymin)/(ymax-ymin)*Float(size_vol[2])
    indz = (z_rays_data-zmin)/(zmax-zmin)*Float(size_vol[3])
    
    ind = where (x_rays_data ge xmin and x_rays_data lt xmax and $
                 y_rays_data ge ymin and y_rays_data lt ymax and $
                 z_rays_data ge zmin and z_rays_data lt zmax, n_in_cube)
    
; Opacity
    if (block_out) then begin
        subim=0.*x_rays_data
        if n_in_cube gt 0 then begin    
            subim(ind)=opacity(fix(indx(ind)),fix(indy(ind)),fix(indz(ind)))
        endif
        img = img - (subim < img)
    endif

; Emissivity
    subim=0.*x_rays_data
    if n_in_cube gt 0 then begin
        subim(ind)=vol.data(fix(indx(ind)),fix(indy(ind)),fix(indz(ind)))
    endif
    if keyword_set(average) then begin
        img = img + subim
    endif else begin
        img = img > (subim > img)
    endelse

; Show temporary image
    if keyword_set(show_progress) then begin
        if i mod FIX(show_progress) eq 0 then begin
            ok=my_tv_loc(img,back,ind_back,xsize,ysize)
        endif
    endif

; Write progress counter
    print,float(i)/float(z_sample)*100.,format='(F5.2,"% complete")'

endfor

if keyword_set(show_progress) then begin
    ok=my_tv_loc(img,back,ind_back,xsize,ysize)
endif

IF ((xsize NE x_sample) OR (ysize NE y_sample)) THEN $
  img = Congrid(img, xsize, ysize, /Interp, /Minus_One)

RETURN, img
END
