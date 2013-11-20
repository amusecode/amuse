pro rd_all,a,h,nout=nout

if not keyword_set(nout) then nout=1

rd_amr,a,nout=nout
rd_hydro,h,nout=nout


end
