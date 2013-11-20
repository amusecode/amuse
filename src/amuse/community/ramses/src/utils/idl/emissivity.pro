PRO EMISSIVITY, Var, Emi

common emissivity_type, e_type
common mewe,ok_read,spec,ener,bin_number,delta_hnu,t_min,t_max,delta_t

alpha=0.28                      ; fraction of He in mass
mu=4./(8.-5.*alpha)
mh=1.67d-24                     ; proton mass in g
kb=1.38d-16                     ; Boltzmann constant in erg.K-1

if e_type eq 0 then begin
    d = Var(*,0)
    p = Var(*,5) 
    t = mu*mh/kb*p/d*1.d14    ; Temperature en K
    t = t/11604.5/1000.       ; Temperature en keV
    Emi = d^2*1.2d-24*sqrt(t)
endif else if e_type eq 1 then begin
    d = Var(*,0)
    p = Var(*,5) 
    t = mu*mh/kb*p/d*1.d14    ; Temperature en K
    t = t/11604.5d0/1000.d0   ; Temperature en keV
    Emi = d^2*1.2d-24*sqrt(t)*t
endif else if e_type eq 2 then begin
    d = Var(*,0)
    p = Var(*,5) 
    t = mu*mh/kb*p/d*1.d14    ; Temperature en K
    t = t/11604.5d0/1000.d0   ; Temperature en keV

    indinf=where(t lt t_min)
    indsup=where(t gt t_max)
    if indinf(0) ne -1 then t(indinf)=t_min
    if indsup(0) ne -1 then t(indsup)=t_max
    
    indkev=fix((t-t_min)/delta_t)

    Emi = d^2*double(spec(indkev,bin_number))
endif else if e_type ge 3 and e_type le 5 then begin
    d = Var(*,0)
    u = Var(*,e_type-2) 
    Emi = d*u
endif else begin
    print,'Invalid emissivity type'
    print,'e_type=',e_type
    stop
endelse

end

