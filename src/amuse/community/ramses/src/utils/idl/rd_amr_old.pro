;+
; NAME:
;	RD_AMR
;
; PURPOSE:
;	This procedure reads the mesh structure from a RAMSES AMR file.
;
; CATEGORY:
;	Input/Output.
;
; CALLING SEQUENCE:
;	RD_AMR, Grid, FILE=file, SWAP=swap, NCPU=ncpu, ICPU=icpu,
;	VERBOSE=verbose 
;
; OPTIONAL INPUTS:
;	FILE:   if set, input the scalar string containing the name of
;	        the file to be read. Otherwise, a PICKFILE widget is
;	        launched.  
;
;       SWAP:   if set, reverse the bit ordering (Little Endian versus
;               Big Endian)
;
;       ICPU:   first cpu file to be read. Default: 1.
;
;       NCPU:   number of cpu files to read, starting from
;               icpu. Default: all files from icpu to ncpu_max.  
;
; OUTPUTS:
;	Grid:   store the AMR tree in structure Grid.
;
; COMMON BLOCKS:
;       None.
;
; EXAMPLE:
;       To read on a SGI architecture a RAMSES AMR file created on a
;       COMPAQ Workstation, type:
;
;	        RD_AMR, Grid, file='amr_00001.out',/swap
;
;       If the file was generated on the same IEEE system, just type:
;
;               RD_AMR, Grid, file='amr_00001.out'
;
; MODIFICATION HISTORY:
; 	Written by:	Romain Teyssier, 01/01/2000.
;                       e-mail: Romain.Teyssier@cea.fr
;	Fevrier, 2001:	Comments and header added by Romain Teyssier.
;-
pro rd_amr, grid, file=file, swap=swap, ncpu=ncpu, icpu=icpu, verbose=verbose, nout=nout

IF N_PARAMS() NE 1 THEN BEGIN
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'rd_amr'
    RETURN
ENDIF

if not keyword_set(icpu) then icpu=1

suffix=getcarnum(icpu)
if not keyword_set(file) and not keyword_set(nout) then begin
    key='*amr*.out'+suffix(icpu-1)
    file=DIALOG_PICKFILE(/READ,filter=key)
endif
if keyword_set(nout) then begin
    suffnout=getcarnum(nout)
    file='output_'+suffnout(nout-1)+'/amr_'+suffnout(nout-1)+'.out'
endif
if not keyword_set(file) then return
base_offset=strpos(file,'.out')+4
file_base=strmid(file,0,base_offset)

; Free memory associated to grid
del_amr,grid

ncpu_run=0L & ndim=0L & nx=0L & ny=0L & nz=0L
nlevelmax=0L & ngridmax=0L & nstep=0L 

print,'Reading file ',trim(file_base)
file=trim(file_base+suffix(icpu-1))
openr,1,file,/f77_unformatted,swap_endian=swap
readu,1,ncpu_run
readu,1,ndim
readu,1,nx,ny,nz
readu,1,nlevelmax
readu,1,ngridmax
readu,1,nstep
readu,1,boxlen
readu,1,t,aexp,hexp
readu,1,omega_m,omega_l,omega_k,omega_b
readu,1,unit_l,unit_d,unit_t
close,1

print,'ncpu      =',ncpu_run
print,'ndim      =',ndim
print,'nlevelmax =',nlevelmax
print,'nstep     =',nstep
print,'boxlen    =',boxlen
print,'time      =',t
if(hexp gt 0.0)then begin ; detect cosmo run
print,'aexp      =',aexp
print,'omega_m   =',omega_m
print,'omega_l   =',omega_l
print,'omega_k   =',omega_k
print,'omega_b   =',omega_b
print,'unit_l    =',unit_l,' cm'
print,'unit_d    =',unit_d,' g/cc'
print,'unit_t    =',unit_t,' s'
endif

if not keyword_set(ncpu) then ncpu=ncpu_run-icpu+1

suffix=getcarnum(ncpu_run)

ncell=2L^ndim
ngrid=lonarr(nlevelmax,ncpu)
level=PTRARR(nlevelmax,ncpu)

; Read AMR grids
ngridtot=0L
for jcpu=0,ncpu-1 do begin

    file=trim(file_base+suffix(jcpu+icpu-1))
    if keyword_set(verbose) then print,'Reading file ',trim(file)
    openr,1,file,/f77_unformatted,swap_endian=swap
    readu,1,ncpu_run
    readu,1,ndim
    readu,1,nx,ny,nz
    readu,1,nlevelmax
    readu,1,ngridmax
    readu,1,nstep
    readu,1,boxlen
    readu,1,t,aexp,hexp
    readu,1,omega_m,omega_l,omega_k,omega_b
    readu,1,unit_l,unit_d,unit_t

    nlevel=0L & ilevel=0L & ng=0L
    for i=0L,nlevelmax-1L do begin
        readu,1,ilevel
        readu,1,ng
        ngrid(i,jcpu)=ng
        if(ngrid(i,jcpu) gt 0)then begin
            ngridtot=ngridtot+ngrid(i,jcpu)
            if keyword_set(verbose) then print,ilevel,ngrid(i,jcpu) $
              ,format='("Level ",i2," has ",i6," grids")'
            nlevel=nlevel+1L
            ; Define level structure
            mesh2={ilevel:ilevel,nc:ng,xg:fltarr(ng,ndim),son:lonarr(ng,ncell)}
            xx=fltarr(ngrid(i,jcpu))
            for idim=0,ndim-1 do begin
                readu,1,xx
                mesh2.xg(*,idim)=xx
            endfor
            xx=0.
            pc=ptr_new(mesh2)
            level(i,jcpu)=pc
        endif
    endfor    
    for i=0L,nlevelmax-1L do begin
        if(ngrid(i,jcpu) gt 0)then begin
            iig=lonarr(ngrid(i,jcpu))
            readu,1,iig
            for icell=0,ncell-1 do begin
                readu,1,iig
                (*level[i,jcpu]).son(*,icell)=iig
            endfor
            iig=0
        endif
    endfor
    close,1
endfor

ngrid=ngrid[0:nlevelmax-1,0:ncpu-1]
level=level[0:nlevelmax-1,0:ncpu-1]

grid={ncpu:ncpu,ndim:ndim,time:t,aexp:aexp,nlevelmax:nlevelmax,boxlen:boxlen $
      ,ngridtot:ngridtot,ngrid:ngrid,level:level,unit_l:unit_l,unit_d:unit_d,unit_t:unit_t}

mesh2=0

return

bad_luck:  print,'I/O Error, exiting...'
           close,1
           mesh2=0
           return

end
;###################################################
;###################################################
;###################################################
