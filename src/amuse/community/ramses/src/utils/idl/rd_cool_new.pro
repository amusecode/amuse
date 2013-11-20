pro rd_cool,cooling,file=file,nout=nout

IF N_PARAMS() NE 1 THEN BEGIN
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'rd_cool'
    RETURN
ENDIF

if not keyword_set(file) and not keyword_set(nout) then begin
    file=DIALOG_PICKFILE(/READ,filter='*cool*.*')
endif
if keyword_set(nout) then begin
    suffnout=getcarnum(nout)
    file='output_'+suffnout(nout-1)+'/cooling_'+suffnout(nout-1)+'.out'
endif
if not keyword_set(file) then return

n1=161L
n2=101L
openr,1,file,/f77_unformatted
readu,1,n1,n2
print,n1,n2
spec=dblarr(n1,n2,6)
metal=dblarr(n1,n2)
cool=dblarr(n1,n2)
heat=dblarr(n1,n2)
cool_com=dblarr(n1,n2)
heat_com=dblarr(n1,n2)
metal_prime=dblarr(n1,n2)
cool_prime=dblarr(n1,n2)
heat_prime=dblarr(n1,n2)
cool_com_prime=dblarr(n1,n2)
heat_com_prime=dblarr(n1,n2)
mu  =dblarr(n1,n2)
n   =dblarr(n1)
T   =dblarr(n2)
Teq =dblarr(n1)
readu,1,n
readu,1,T
readu,1,cool
readu,1,heat
readu,1,cool_com
readu,1,heat_com
readu,1,metal
readu,1,cool_prime
readu,1,heat_prime
readu,1,cool_com_prime
readu,1,heat_com_prime
readu,1,metal_prime
readu,1,mu
close,1

cooling={n1:n1,n2:n2,n:n,t:t,teq:teq,cool:cool,heat:heat,metal:metal,cool_com:cool_com,heat_com:heat_com,cool_prime:cool_prime,heat_prime:heat_prime,cool_com_prime:cool_com_prime,heat_com_prime:heat_com_prime,metal_prime:metal_prime,mu:mu}

end
