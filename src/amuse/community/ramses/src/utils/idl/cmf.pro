;+
; Name: cmf.pro
; routine to produce clump mass functions from ramses/clumpfinder autput
; Andreas Bleuler 2011
;
; calling sequence:
; 
; cmf,nbins
; 
; with nbins the desired number of bins
;
; clump textfiles (clump_masses.txt from clumpfinder output) can be selected from a window
; multiple files can be selected at the same time
;-
pro cmf,bins

IF N_PARAMS() LT 1 THEN BEGIN
    DOC_LIBRARY,'cmf'
    RETURN
ENDIF

bins=1.0*bins

round=1
answer=1
while ((round LE 5) and (answer EQ 1)) do begin

    file = Dialog_Pickfile(/Read, Title='select clump_mass files', /MULTIPLE_FILES,FILTER = '*.txt')
    print,file
    print,size(file)
    

;get total number of clumps   
    n_tot=0
    n=0
    for i=0,(size(file))(1)-1 do begin 
        openr,10,file(i)
        readf,10,n
        nn=file_lines(file(i))-1
        if (nn ne n) then print,'attention, problem with old bug in clump_merger routine. only clumps with peak to saddle ratio of > 2 where printed, even if a different value was chosen in the code'
        print,'number of clumps = ',nn
        n_tot=n_tot+nn
        m=fltarr(n)
        readf,10,m
        if (i eq 0) then m_tot=m else m_tot=[m_tot,m]
        close,10
    end 
    
;construct and plot histogram
    m_log=alog10(m_tot)
    ma=max(m_log)+0.001
    mi=-4.-0.001
;mi=min(m_log)-0.001
    h=histogram(m_log,nbins=bins,max=ma,min=mi)
    h=alog10(h)
    x=(0.5+indgen(bins))*(ma-mi)/(1.*(bins-1.)) + mi
    if (round EQ 1)then begin
        loadct,0
        window,0,title='CMF',xs=1024,ys=768
        plot,x,h,psym=10,color=255,/nodata,xtitle='log (M/M_sol)',ytitle='log N',title='clump mass function from simulations compared to Kroupa IMF'
        loadct,4
    endif
    oplot,x,h,psym=10,color=250-round*25,linestyle=round-1,thick=2
    
    ;define Kroupa IMF 
    if (round EQ 1)then begin
        y_max=max(h)
        x_max=alog10(0.08)
        x_0=-2.
        y_0=y_max+(x_0-x_max)*0.7
        x_1=alog10(0.5)
        y_1=y_max+(x_1-x_max)*(-0.3)
        x_2=1.
        y_2=y_1+(x_2-x_1)*(-1.3)
        x2=[x_0,x_max,x_1,x_2]
        y2=[y_0,y_max,y_1,y_2]
        oplot,x2,y2,color=160,linestyle=0,thick=1
    endif
    
;another cmf if desired
    answer=0
    read,'do you want to add another cmf? (type 1 if yes) ',answer
    round=round+1
endwhile
   
if (round GT 5)then print,'you should not overload your plot anyway '

;save picture if wanted
answer=0
read,'do you want to save plot to file? (type 1 if yes) ',answer

if (answer eq 1) then begin
    filename=''
    read,'filename: ',filename
    image=tvrd()
    tv,image
    write_jpeg,filename,image
endif


end
;###################################################
;###################################################
;###################################################
