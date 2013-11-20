;+
; NAME:
;	RD_LOG
;
; PURPOSE:
;	This procedure reads particles from a RAMSES LOG ASCII file.
;
; CATEGORY:
;	Input/Output.
;
; CALLING SEQUENCE:
;	RD_HYDRO, Log, FILE=file, NMAX=nmax, LMAX=lmax
;
; OPTIONAL INPUTS:
;	FILE:   if set, input the scalar string containing the name of
;	the file to be read. Otherwise, a PICKFILE widget is launched. 
;
;	NMAX:  if set, the maximum number of COARSE time steps to be
;	read from the file. Default: 1000.
;
;       LMAX:   if set, the maximum number of levels to be read from
;       the file. Default: 10.
;	
; OUTPUTS:
;	Log: structure containing the control variables of the RAMSES
;	run. 
;
; EXAMPLE:
;       To read a RAMSES LOG ASCII file, type:
;
;	        RD_LOG,log,file='run12.log'
;
; MODIFICATION HISTORY:
; 	Written by:	Romain Teyssier, 01/01/2000.
;                       e-mail: Romain.Teyssier@cea.fr
;	Fevrier, 2001:	Comments and header added by Romain Teyssier.
;-
pro rd_log,log,file=file,nmax=nmax,lmax=lmax,star=star

IF N_PARAMS() NE 1 THEN BEGIN
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'rd_log'
    RETURN
ENDIF

if not keyword_set(file) then file=pickfile(/READ)
if not keyword_set(nmax) then nmax=1000L
if not keyword_set(lmax) then lmax=20

print,'Reading ',TRIM(file)
openr,1,file
st=''
jpos=-1
while not eof(1) and jpos eq -1 do begin
    readf,1,st
    jpos=strpos(st,'nproc =')
endwhile
ncpu=long(strmid(st,jpos+8,100))
close,1

openr,1,file
st=''
i=0
econs=dblarr(nmax)
epot=dblarr(nmax)
ekin=dblarr(nmax)
eint=dblarr(nmax)
emag=dblarr(nmax)
aexp=dblarr(nmax)
t   =dblarr(nmax)
is  =lonarr(nmax)
n   =lonarr(nmax)
nl  =lonarr(nmax,lmax)
dt  =dblarr(nmax)
mstar=dblarr(nmax)
mem =dblarr(nmax)
memp =dblarr(nmax)

j=0L
ilevmax=0
iipos=-1
while not eof(1) and j lt nmax do begin
    jjpos=-1
    while iipos eq -1 and jjpos eq -1 and not eof(1) do begin
        readf,1,st
        iipos=strpos(st,'Mesh')
        jjpos=strpos(st,'Main')
    endwhile
    if not eof(1) then begin
        nnn=0
        nnnl=lonarr(lmax)
        ilev=0
        if(jjpos eq -1)then begin
            ipos=-1
            while ipos eq -1 and not eof(1) do begin
                readf,1,st
                ipos=strpos(st,'econs')
                jpos=strpos(st,'has')
                if jpos ne -1 then begin
                    ilev=ilev+1
                    nnnl(ilev)=long(strmid(st,jpos+3,100))
                    nnn=nnn+long(strmid(st,jpos+3,100))
                endif
            endwhile
        endif else begin
            ipos=strpos(st,'econs')
        endelse
        ilevmax=max([ilevmax,ilev])
        eee=double(strmid(st,ipos+6,100))

        ipos=strpos(st,'step')
        nstep=long(strmid(st,ipos+6,100))
        ipos=strpos(st,'epot')
        www=double(strmid(st,ipos+6,100))
        ipos=strpos(st,'ekin')
        kkk=double(strmid(st,ipos+6,100))
        ipos=strpos(st,'eint')
        iii=0.
        if(ipos ne -1) then iii=double(strmid(st,ipos+6,100))
        ipos=strpos(st,'emag')
        mag=0.
        if(ipos ne -1) then mag=double(strmid(st,ipos+6,100))

        readf,1,st
        ipos=strpos(st,' a=')
        aaa=double(strmid(st,ipos+3,100))
        ipos=strpos(st,' dt=')
        ttt=double(strmid(st,ipos+4,100))
        ipos=strpos(st,' t=')
        uuu=double(strmid(st,ipos+3,100))
        ipos=strpos(st,' mem=')
        mmm=double(strmid(st,ipos+5,100))
        ipos=strpos(st,' mem=')
        mmp=double(strmid(st,ipos+11,100))
        econs[j]=eee
        ekin [j]=kkk
        eint [j]=iii
        emag [j]=mag
        epot [j]=www
        aexp [j]=aaa
        n    [j]=nnn
        is   [j]=nstep
        for il=0,lmax-1 do begin
            nl   [j,il]=nnnl[il] 
        endfor
        dt   [j]=ttt
        t    [j]=uuu
        mem  [j]=mmm
        memp [j]=mmp
        ipos=-1
        while ipos eq -1 and jpos eq -1 and not eof(1) do begin
            readf,1,st
            ipos=strpos(st,'Mesh')
            jpos=strpos(st,'New star')
        endwhile
        if jpos ne -1 then begin
            kpos=strpos(st,'Mass=')
            mstar[j]=double(strmid(st,kpos+5,100))
        endif
        iipos=ipos
        j=j+1L
    endif
endwhile
close,1
ntime=j
tmean=eint/0.13*1d10*1.66d-24/1.38d-16*2./3./(aexp+1.0d-15)^2
mem=mem/100.
memp=memp/100.
log={ncpu:ncpu,ntime:ntime, $
     t:t[0:ntime-1],dt:dt[0:ntime-1],aexp:aexp[0:ntime-1],nstep:is[0:ntime-1], $
     noct:n[0:ntime-1],noctl:nl[0:ntime-1,0:ilevmax],ekin:ekin[0:ntime-1],emag:emag[0:ntime-1], $
     eint:eint[0:ntime-1],epot:epot[0:ntime-1],econs:econs[0:ntime-1],tmean:tmean[0:ntime-1],mstar:mstar[0:ntime-1],mem:mem[0:ntime-1],memp:memp[0:ntime-1]}

end
