pro convert_ieee,file=file,type=type,float=float

if not keyword_set(type) then type='part'

if type eq 'part' then begin

    if not keyword_set(file) then file=DIALOG_PICKFILE(/READ,filter='*part*')
    print,'Reading file ',file
    ndimp=0L
    npart=0L
    openr,1,file,/f77_unformatted,swap_endian=1
    readu,1,ndimp
    readu,1,npart
    print,ndimp,npart,format='("ndim=",I1," npart=",I8)'    
    if not keyword_set(float)then begin
; If file is in REAL*8
        xp=dblarr(npart,ndimp)
        vp=dblarr(npart,ndimp)
        mp=dblarr(npart)
        lp=lonarr(npart)
    endif else begin
; If file is in REAL*4
        xp=fltarr(npart,ndimp)
        vp=fltarr(npart,ndimp)
        mp=fltarr(npart)
        lp=intarr(npart)
    endelse
    readu,1,xp
    readu,1,vp
    readu,1,mp
    readu,1,lp
    close,1
    
    fileout=file+'.ieee'
    print,'Writing file ',fileout
    openw,1,fileout,/f77_unformatted
    writeu,1,ndimp
    writeu,1,npart
    writeu,1,xp
    writeu,1,vp
    writeu,1,mp
    writeu,1,lp
    close,1

; Free memory
    xp=0.
    vp=0.
    mp=0.
    lp=0.

endif

end
