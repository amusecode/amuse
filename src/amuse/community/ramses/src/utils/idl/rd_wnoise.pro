pro rd_wnoise,d,file=file


n1=0L
n2=n1
n3=n1
iseed=0L
openr, 1,file,/f77_unf

readu,1,n1,n2,n3,iseed
print,n1,n2,n3,iseed

din_plane=fltarr(n1,n2)
d=fltarr(n1,n2,n3)
for i3=0,n3-1 do begin

    print,i3

    readu,1,din_plane
    d(*,*,i3)=din_plane
endfor




end
