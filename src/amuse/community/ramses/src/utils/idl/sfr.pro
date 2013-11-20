pro sfr,s $
        ,smoot=smoot,hist=hist,disp=disp,time=time, nbin=nbin $
        ,opl=opl,color=color,psym=psym,linestyle=linestyle,charsize=charsize $
        ,title=title,xtitle=xtitle,ytitle=ytitle,boite=boite,ncoarse=ncoarse,csfr=csfr,fstars=fstars,info=info,ylog=ylog,thick=thick $
        ,xr=xr,yr=yr,adapt=adapt $
        ,gal=gal,mvirspec=mvirspec,sfrout=sfrout,zout=zout,noplot=noplot,oldhist=oldhist

;###################################################
;TO SAVE IN A FILE
;set_plot,'PS'& device, file='~/SFR/rapport/testc/sfrjoliceq0.ps',/portrait, /color, bits=8 ,xs=20, ys=18,xoffset=0.5,yoffset=5.
;sfr,s,boite=10,ncoarse=128,csfr=0.1,/fstars,charsize=1.4
;device,/close&set_plot,'x'
;##################################################
;+
; NAME:
;       SFR
;
; PURPOSE:
;       This procedure computes the comoving star formation rate
;       from the structure s created by rd_star.pro
;        
; CATEGORY:
;       Stars analyze routine.
;
; CALLING SEQUENCE:
;         sfr,s 
;        ,smoot=smoot,hist=hist,disp=disp,time=time, nbin=nbin 
;        ,opl=opl,color=color,psym=psym,linestyle=linestyle,charsize=charsize 
;        ,title=title,xtitle=xtitle,ytitle=ytitle,boite=boite,ncoarse=ncoarse,csfr=csfr,fstars=fstars,info=info 
;        ,xr=xr,yr=yr

; INPUTS:
;       S: structure containing star particle masses and star
;       particule formation expansion factor
;
;  OPTIONAL INPUTS:       
;       SMOOT (Default): Smooth the star formation rate with a
;       constant number of particles in each of the nbin bin
;       
;       HIST: Histogram mode with a bin redshift. Use this mode for
;       the too discontinuous star formation!!!
;       
;       DISP: To have the rough SFR
;       
;       TIME: To have SFR vs time (and not vs redshift (default))
;   
;       NBIN: Number of bins 
;
;
;  GRAPHIC OPTIONS: 
;
;       OPL: To overplot SFR
;       
;       COLOR: Tek_color color of the smooth or hist plot (default red)
;
;       PSYM, LINESTYLE,CHARSIZE,XR,YR,TITLE,XTITLE,YTITLE: Plot
;       options. Default: psym=-1, linestyle=0., charsize=1.2,
;       xr=[0,12] and yr=[0.001,1] and xtitle='z'for the redshift mode,  
;       xr=[0,15] and yr=[0.,0.3] and xtitle='t(Gyr)' for the time mode, 
;       title='Comoving SFR', ytitle='SFR(M!d!9n!17!n/Mpc!u3!n/yr)'
;
;       BOITE: Size of the box in Mpc for the title
;       NCOARSE: (Number of cells for the coarse grid)^(1./3.) for the
;       title
;       CSFR: Normalisation of the Kennicutt SFR  for the title
;       FSTARS: Fraction of stars at the final redshift for the title
;       INFO: Other informations for the title
;
;       GAL: SFR in Msun/yr. You need to give the value of BOITE and
;       adjust XR and YR
;       MVIRSPEC: SFR in Msun/Gyr/Msun_gas(=0.13*Mvir). You need to give Mvir in Msun
;       and adjust XR and YR
;
; OUTPUTS:
;       None.
;       
; COMMON BLOCKS:
;       None.
;
; EXAMPLE:
;      To overplot on a graph, the smooth sfr vs time   
;              sfr,s,/time,/opl,color=4
;      To plot rough sfr and histogram sfr vs redshift  with a beautiful title
;              sfr,s,/disp,/hist,nbin=200,boite=10,ncoarse=128,csfr=0.1,/fstars,charsize=1.
;
; MODIFICATION HISTORY:
;       Written by:     Yann Rasera, 15/04/2003.
;                       e-mail: yann.rasera@cea.fr
;       Fevrier, 2001:  Comments and header added by Yann Rasera
;-
;###################################################
;###################################################
;###################################################

if s.npart eq 0 then return
;parametres cosmo
omega_m=0.3
omega_lambda=0.7
unsurH0=9.8/0.7
;on a choisi h=0.7   

;tableaux
deltam=fltarr(100000)
deltaz=fltarr(100000)
dtsurdz=fltarr(100000)
z=fltarr(100000)

;initialisation
deltamtmp=0.
ind=sort(s.ap)
as=s.ap(ind)
mp=s.mp(ind)
aptmp=as(0)
j=0L

;histogramme en masse des redshifts pour avoir deltaM/deltaZ
for i=0L,long(n_elements(as)-1) do begin
    if (as(i) gt aptmp*0.99999 and as(i) lt aptmp*1.00001 ) then begin
        deltamtmp=deltamtmp+mp(i) 
    endif else begin
        z(j)=1./aptmp-1.
        if (j gt 0) then deltaz(j)=z(j-1)-z(j) else deltaz(j)=1000. ;valeur grande ne servant a rien       
        dtsurdz(j)=unsurH0/(1.+z(j))*1./sqrt( (1.+z(j))^2*(1.+omega_m*z(j) )-z(j)*(2.+z(j))*omega_lambda )
        deltam(j)=deltamtmp
        
        deltamtmp=mp(i)
        aptmp=as(i)
        j=j+1
    endelse
endfor
z(j)=1./aptmp-1.
deltaz(j)=z(j-1)-z(j)        
dtsurdz(j)=unsurH0/(1.+z(j))*1./sqrt( (1.+z(j))^2*(1.+omega_m*z(j) )-z(j)*(2.+z(j))*omega_lambda )
deltam(j)=deltamtmp

deltam=deltam(0:j)
deltaz=deltaz(0:j)
dtsurdz=dtsurdz(0:j)
z=z(0:j)
 
;conversion en temps
deltat=deltaz*dtsurdz ;approximation deltat<<age Univers
deltat=deltat*1.d9
deltam=deltam*4.08d10 ;(h=0.7)
;differentes unites possibles
if keyword_set(mvirspec) then gal=1
if keyword_set(gal) then begin
    if not keyword_set(boite) then begin
        print,'Attention: taille de boite par defaut: 10 h^-1 Mpc!'
        boite=10.
    endif
    deltam=deltam*(boite/0.7)^3 ;pour avoir sfr en Msun/yr!
    ;boite=0. 
    if keyword_set(mvirspec) then begin 
        deltam=deltam/(0.13*mvirspec)*1d9 ;pour avoir sfr en Msun/Gyr/Msun
    endif
endif
SFRcom=deltam/deltat
dist=1.
make_dist,z(0),omega_m=0.3,omega_v=0.7,omega_r=0.,save=dist
deltat(0)=dist.temps*1.d9/0.7   ;h0=0.7
dist=0.
 
;calcul de l'age de l'Univers
t=fltarr(n_elements(sfrcom))
for i=0,n_elements(sfrcom)-1 do begin 
t(i)=total(deltat(0:i))
endfor
t=t/1.d9


;MOYENNES

if (keyword_set(nbin)) then begin
    nbin=nbin
    npas=nbin 
endif else begin 
    nbin=50
    npas=200 
endelse



;Valeur par defaut
if not keyword_set(disp) and not keyword_set(hist)  and not keyword_set(oldhist) then smoot=1



;moyenne a N constant pour avoir nbin=nbin (pour les cas tres discontinu la
;mieux vaut utiliser sfr3,s,npas!!!!!)
if keyword_set(smoot) then begin
     largeur=n_elements(sfrcom)/nbin
;    sfrmean=smooth(sfrcom,largeur,/edge_truncate) 
     
     if (largeur lt 3) then begin 
         print, 'trop peu d''etoiles dans chaque bin, ce nombre est fixe a 3'
         largeur=3
     endif

    sfrmean=smooth(deltam, largeur,/edge_truncate)/smooth(deltat, largeur,/edge_truncate)
endif

;histogramme a pas constant en z en utilisant directement 
;psym=10 sans reflechir pour les cas tres discontinus
if keyword_set(hist) then begin
    ttot=fltarr(npas)&ttot(0)=13.45
    zmin=min(z)
    zmax=max(z)+5.
    ztab=bindec(zmin,zmax,npas) ;le plus 5 c'est pour etre bien avant debut
    deltat=(zmax-zmin)/npas*dtsurdz;si npas est grand!

    sfrpas=fltarr(npas)
    zpas=fltarr(npas)
    for i=0,npas-1 do begin
        ind=where(z ge ztab(i) and z lt ztab(i+1),nok)
        zpas(i)=(ztab(i+1)+ztab(i))/2.
        if (nok gt 0) then begin
            zpas(i)=mean(z(ind))
            sfrpas(i)=total(deltam(ind))/deltat(ind(0))/1.d9    
            ttot(i)=mean(t(ind))
        endif
        if (ttot(i) lt 0.0001 and zpas(i) lt z(0)) then begin
            ttot(i)=ttot(i-1)-20./npas*unsurH0/(1.+zpas(i))*1./sqrt( (1.+zpas(i))^2*(1.+omega_m*zpas(i) )-zpas(i)*(2.+zpas(i))*omega_lambda )
        endif
    endfor  
endif

if keyword_set(oldhist) then begin
    ttot=fltarr(npas)&ttot(0)=13.45 ;age de l'Univers   
    ztab=findgen(npas+1)*20./npas  ;redshift de depart 20!
    tpas=20./npas*dtsurdz ;si npas est grand!
    sfrpas=fltarr(npas)
    zpas=fltarr(npas)
    for i=0,npas-1 do begin
        ind=where(z ge ztab(i) and z lt ztab(i+1),nok)
        zpas(i)=(ztab(i+1)+ztab(i))/2.
        if (nok gt 0) then begin
            zpas(i)=mean(z(ind))
            sfrpas(i)=total(deltam(ind))/tpas(ind(0))/1.d9    
            ttot(i)=mean(t(ind))
        endif
        if (ttot(i) lt 0.0001 and zpas(i) lt z(0)) then begin
            ttot(i)=ttot(i-1)-20./npas*unsurH0/(1.+zpas(i))*1./sqrt( (1.+zpas(i))^2*(1.+omega_m*zpas(i) )-zpas(i)*(2.+zpas(i))*omega_lambda )
        endif
    endfor  
endif






;PLOT
fstar=total(s.mp,/double)/0.13
!p.background=255 ;!!!!!!!!!!!!yes!!!
!p.color=0
tek_color

;style
if not keyword_set(color) then color=0
if not keyword_set(linestyle) then linestyle=0
if not keyword_set(charsize) then charsize=1.2
if not keyword_set(psym) then psym=0

;Titre
if not keyword_set(title) then begin
    title='Comoving SFR'
    if keyword_set(boite) then title= title+'('+string(boite,format='(I3)')+'Mpc'
    if keyword_set(ncoarse) then title=title+'/'+string(ncoarse,format='(I3)')+'!u3!n'
    if keyword_set(csfr) then title=title+'/c='+string(csfr,format='(F3.1)')
    if keyword_set(info) then title=title+info
    if keyword_set(fstars) then title= title+'/f!d*!n='+string(fstar,format='(F4.2)')+')'
endif 

;Abscisse et ordonnee
if not keyword_set(ytitle) then ytitle='SFR(M!d!9n!17!n/Mpc!u3!n/yr)'
if keyword_set(gal) then begin
ytitle='SFR(M!d!9n!17!n/yr)'
if keyword_set(mvirspec) then ytitle='SFR(M!d!9n!17!n/Gyr/M!d!9n!17!n)'
endif

if keyword_set(time) then begin
    if not keyword_set(xr) then xr=[0.,15.]
    if not keyword_set(yr) then yr=[0.,0.3]
    if not keyword_set(xtitle) then xtitle='t(Gyr)'
    if keyword_set(adapt) then begin
        xr=[0.9*min(t),1.1*max(t)]
        yr=[0.8*min(sfrcom),1.2*max(sfrcom)]
    endif

endif else begin
ylog=1
if not keyword_set(xr) then xr=[0.,12.]
if not keyword_set(yr) then yr=[0.001,1.]
if not keyword_set(xtitle) then xtitle='z'
if keyword_set(adapt) then begin
    xr=[0.9*min(z),1.1*max(z)]
    yr=[0.9*min(sfrcom),1.1*max(sfrcom)]
endif

endelse

;definition des plots
if not keyword_set(noplot) then begin
    if not(keyword_set(opl)) then plot,xr,yr,xr=xr,yr=yr,/xs,/ys,title=title,xtitle=xtitle,ytitle=ytitle,psym=3,charsize=charsize,ylog=ylog,thick=thick
endif
;oplots

if keyword_set(smoot) then begin
    sfrout=sfrmean
    zout=z
    if keyword_set(time) then x=t  else  x=z
    y=sfrmean
    if not keyword_set(noplot) then oplot,x,y,psym=psym,linestyle=linestyle,color=color,thick=thick
endif

if keyword_set(disp) then begin
    sfrout=sfrcom
    zout=z
    if keyword_set(time) then x=t else x=z
    y=sfrcom
    if not keyword_set(noplot) then oplot,x,y,psym=3,thick=thick
endif

if keyword_set(hist) or keyword_set(oldhist) then begin
    if not keyword_set(psym) then psym=10
    zout=zpas
    sfrout=sfrpas
    if keyword_set(time) then x=ttot else x=zpas
    y=sfrpas
    if not keyword_set(noplot) then oplot,x,y,psym=psym,linestyle=linestyle,color=color,thick=thick
endif

;histogramme brut
;redshift
;plot,[z(0),z(0),z(1)],[1.d-4,sfrcom(1),sfrcom(1)],/ylog,xr=[0.,12.],yr=[0.001,1.],color=0,title='Comoving SFR (10Mpc/64!u3!n/c=0.1/f!d*!n='+string(fstar,format='(F4.2)')+')',ytitle='SFR(M!d!9n!17!n/Mpc!u3!n/yr)',xtitle='z'
;for i=1,n_elements(z)-2 do oplot,[z(i),z(i),z(i+1)],[sfrcom(i),sfrcom(i+1),sfrcom(i+1)]
;temps
;plot,[t(0),t(0),t(1)],[0.,sfrcom(1),sfrcom(1)],title='Comoving SFR (10Mpc/64!u3!n/c=0.1/f!d*!n='+string(fstar,format='(F4.2)')+')',ytitle='SFR(M!d!9n!17!n/Mpc!u3!n/yr)',xtitle='t(Gyr)',psym=10,color=0,xr=[0,15],yr=[0,0.3]   
;for i=1,n_elements(z)-2 do oplot,[t(i),t(i),t(i+1)],[sfrcom(i),sfrcom(i+1),sfrcom(i+1)]   
  
end
    








