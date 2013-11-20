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
pro rd_amr, grid, file=file, swap=swap, icpu=icpu, verbose=verbose, nout=nout

IF N_PARAMS() NE 1 THEN BEGIN
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'rd_amr'
    RETURN
ENDIF

if not keyword_set(icpu) then icpu=0
if icpu eq 0 then jcpu=1 else jcpu=icpu

suffix=getcarnum(jcpu)
if not keyword_set(file) and not keyword_set(nout) then begin
    key='*amr*.out'+suffix(jcpu-1)
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

; Initialize header variables
ncpu_run=0L & ndim=0L & nx=0L & ny=0L & nz=0L
nlevelmax=0L & ngridmax=0L & nboundary=0L & ngridactual=0L & nstep=0L 
noutput=0L & boxlen=0.0d0 & t=0.0d0
iout=0L & ifout=0L
aexp=0.0d0 & hexp=0.0d0 & aexp_old=0.0d0 & epot_tot_int=0.0d0
epot_tot_old=0.0d0
const=0.0d0 & mass_tot_0=0.0d0 & rho_tot=0.0d0
omega_m=0.0d0 & omega_l=0.0d0 & omega_k=0.0d0 & omega_b=0.0d0 & h0=0.0d0
aexp_ini=0.0d0 & mass_sph=0.0d0
headf=0L & tailf=0L & numbf=0L & used_mem=0L & used_mem_tot=0L

; Read first file to get header
print,'Reading file ',trim(file_base)
file=trim(file_base+suffix(jcpu-1))
openr,1,file,/f77_unformatted,swap_endian=swap
readu,1,ncpu_run
readu,1,ndim
readu,1,nx,ny,nz
ncoarse=nx*ny*nz
readu,1,nlevelmax
readu,1,ngridmax
readu,1,nboundary
readu,1,ngridactual
readu,1,boxlen
readu,1,noutput,iout,ifout
tout=dblarr(noutput)
aout=dblarr(noutput)
readu,1,tout
readu,1,aout
readu,1,t
dtold=dblarr(nlevelmax)
dtnew=dblarr(nlevelmax)
readu,1,dtold
readu,1,dtnew
readu,1,nstep,nstep_coarse
readu,1,const,mass_tot_0,rho_tot
readu,1,omega_m,omega_l,omega_k,omega_b,h0,aexp_ini
readu,1,aexp,hexp,aexp_old,epot_tot_int,epot_tot_old
readu,1,mass_sph
close,1

; Write header to screen
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
endif
if nboundary eq 0 then begin
    print,"Periodic boundary conditions"
endif

; Allocate arrays
ncpu=ncpu_run
icpumin=1L & icpumax=ncpu & listmax=ncpu
if icpu gt 0 then begin
    icpumin=icpu & icpumax=icpu & listmax=ncpu+nboundary
endif
suffix=getcarnum(ncpu_run)
ncell=2L^ndim
ngrid=LONARR(nlevelmax,listmax)
level=PTRARR(nlevelmax,listmax)
headl=LONARR(ncpu,nlevelmax)
taill=LONARR(ncpu,nlevelmax)
numbl=LONARR(ncpu,nlevelmax)
if(nboundary gt 0)then begin
    headb=LONARR(nboundary,nlevelmax)
    tailb=LONARR(nboundary,nlevelmax)
    numbb=LONARR(nboundary,nlevelmax)
endif
numbtot=LONARR(10,nlevelmax)
bound_key=DBLARR(ncpu+1)
son=LONARR(ncoarse)
flag1=LONARR(ncoarse)
cpu_map=LONARR(ncoarse)
xbound=[0d0,0d0,0d0]
ordering='                       '
; Read AMR grids

; Loop over cpu files
ngridtot=0L
list=0L
for jcpu=icpumin,icpumax do begin

    file=trim(file_base+suffix(jcpu-1))
    if keyword_set(verbose) then print,'Reading file ',trim(file)
    openr,1,file,/f77_unformatted,swap_endian=swap
    
; Read header
    readu,1,ncpu_run
    readu,1,ndim
    readu,1,nx,ny,nz
    readu,1,nlevelmax
    readu,1,ngridmax
    readu,1,nboundary
    readu,1,ngridactual
    readu,1,boxlen
    readu,1,noutput,iout,ifout
    readu,1,tout
    readu,1,aout
    readu,1,t
    readu,1,dtold
    readu,1,dtnew
    readu,1,nstep,nstep_coarse
    readu,1,const,mass_tot_0,rho_tot
    readu,1,omega_m,omega_l,omega_k,omega_b,h0,aexp_ini
    readu,1,aexp,hexp,aexp_old,epot_tot_int,epot_tot_old
    readu,1,mass_sph
    readu,1,headl
    readu,1,taill
    readu,1,numbl
    readu,1,numbot
    if nboundary gt 0 then begin
        readu,1,headb
        readu,1,tailb
        readu,1,numbb
        xbound=[double(nx/2),double(ny/2),double(nz/2)]
    endif
    readu,1,headf,tailf,numbf,used_mem,used_mem_tot
    readu,1,ordering
    if(strcompress(ordering) eq 'bisection ')then begin
        readu,1
        readu,1
        readu,1
        readu,1
        readu,1
    endif else begin
        readu,1,bound_key
    endelse
    readu,1,son
    readu,1,flag1
    readu,1,cpu_map
    
; Read fine levels
    nlevel=0L & ilevel=0L & ng=0L
    kcpumin=1L & kcpumax=nboundary+ncpu
    for ilevel=0L,nlevelmax-1L do begin
        for kcpu=kcpumin,kcpumax do begin
            if(kcpu le ncpu)then begin
                ng=numbl(kcpu-1,ilevel)
            endif else begin
                ng=numbb(kcpu-ncpu-1,ilevel)
            endelse
            if icpu eq 0 then begin
                if (kcpu eq jcpu) then ngrid(ilevel,kcpu-1)=ng
            endif else begin
                ngrid(ilevel,kcpu-1)=ng
            endelse
            if(ng gt 0)then begin
                ngridtot=ngridtot+ng
                if keyword_set(verbose) then begin
                    print,ilevel+1,ng,kcpu $
                      ,format='("Level ",i2," has ",i6," grids in proc",i6)'
                endif
                nlevel=nlevel+1L
                mesh2={ilevel:ilevel,nc:ng,xg:dblarr(ng,ndim),son:lonarr(ng,ncell)}
                ii=lonarr(ng)   ; Define level structure
                xx=dblarr(ng)
                readu,1,ii      ; Read grid index
                readu,1,ii      ; Read next index
                readu,1,ii      ; Read prev index
                for idim=0,ndim-1 do begin
                    readu,1,xx  ; Read grid center
                    mesh2.xg(*,idim)=xx-xbound(idim)
                endfor
                readu,1,ii      ; Read father index
                for idim=0,2*ndim-1 do begin
                    readu,1,ii  ; Read nbor index
                endfor
                for idim=0,2^ndim-1 do begin
                    readu,1,ii  ; Read son index
                    mesh2.son(*,idim)=ii
                endfor
                for idim=0,2^ndim-1 do begin
                    readu,1,ii  ; Read cpu map
                endfor
                for idim=0,2^ndim-1 do begin
                    readu,1,ii  ; Read refinement map
                endfor
                xx=0d0
                if icpu eq 0 then begin
                    if(kcpu eq jcpu)then begin
                        pc=ptr_new(mesh2)
                        level(ilevel,jcpu-1)=pc
                    endif
                endif else begin
                    pc=ptr_new(mesh2)
                    level(ilevel,kcpu-1)=pc
                endelse
                list=list+1L
            endif
        endfor
    endfor    
    close,1
endfor 
ngrid=ngrid[0:nlevelmax-1,0L:listmax-1L]
level=level[0:nlevelmax-1,0L:listmax-1L]

grid={ncpu:listmax,ndim:ndim,time:t,aexp:aexp,nlevelmax:nlevelmax,boxlen:boxlen $
      ,ngridtot:ngridtot,ngrid:ngrid,level:level}

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
