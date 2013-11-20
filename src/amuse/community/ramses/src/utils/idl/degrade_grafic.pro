pro degrade_grafic,dir_in,dir_out


n1=0L
n2=n1
n3=n1
xoff1=0.0
xoff2=xoff1
xoff3=xoff1
astart=0.0
omega_m=0.0
omega_l=0.0
h0=0.0

openr, 1,dir_in +'/ic_d',/f77_unf
openr, 2,dir_in +'/ic_u',/f77_unf
openr, 3,dir_in +'/ic_v',/f77_unf
openr, 4,dir_in +'/ic_w',/f77_unf
openr, 5,dir_in +'/ic_p',/f77_unf
openw,11,dir_out+'/ic_d',/f77_unf
openw,12,dir_out+'/ic_u',/f77_unf
openw,13,dir_out+'/ic_v',/f77_unf
openw,14,dir_out+'/ic_w',/f77_unf
openw,15,dir_out+'/ic_p',/f77_unf

readu,1,n1,n2,n3,dxini,xoff1,xoff2,xoff3,astart,omega_m,omega_l,h0
readu,2,n1,n2,n3,dxini,xoff1,xoff2,xoff3,astart,omega_m,omega_l,h0
readu,3,n1,n2,n3,dxini,xoff1,xoff2,xoff3,astart,omega_m,omega_l,h0
readu,4,n1,n2,n3,dxini,xoff1,xoff2,xoff3,astart,omega_m,omega_l,h0
readu,5,n1,n2,n3,dxini,xoff1,xoff2,xoff3,astart,omega_m,omega_l,h0

print,n1,n2,n3,dxini

n1_out=n1/2
n2_out=n2/2
n3_out=n3/2
dxini_out=2.0*dxini

writeu,11,n1_out,n2_out,n3_out,dxini_out,xoff1,xoff2,xoff3,astart,omega_m,omega_l,h0
writeu,12,n1_out,n2_out,n3_out,dxini_out,xoff1,xoff2,xoff3,astart,omega_m,omega_l,h0
writeu,13,n1_out,n2_out,n3_out,dxini_out,xoff1,xoff2,xoff3,astart,omega_m,omega_l,h0
writeu,14,n1_out,n2_out,n3_out,dxini_out,xoff1,xoff2,xoff3,astart,omega_m,omega_l,h0
writeu,15,n1_out,n2_out,n3_out,dxini_out,xoff1,xoff2,xoff3,astart,omega_m,omega_l,h0

din_plane_1=fltarr(n1,n2)
din_plane_2=fltarr(n1,n2)
uin_plane_1=fltarr(n1,n2)
uin_plane_2=fltarr(n1,n2)
vin_plane_1=fltarr(n1,n2)
vin_plane_2=fltarr(n1,n2)
win_plane_1=fltarr(n1,n2)
win_plane_2=fltarr(n1,n2)
pin_plane_1=fltarr(n1,n2)
pin_plane_2=fltarr(n1,n2)

dout_plane=fltarr(n1_out,n2_out)
uout_plane=fltarr(n1_out,n2_out)
vout_plane=fltarr(n1_out,n2_out)
wout_plane=fltarr(n1_out,n2_out)
pout_plane=fltarr(n1_out,n2_out)

for i3=0,n3-1,2 do begin

    print,i3

    readu,1,din_plane_1
    readu,1,din_plane_2
    readu,2,uin_plane_1
    readu,2,uin_plane_2
    readu,3,vin_plane_1
    readu,3,vin_plane_2
    readu,4,win_plane_1
    readu,4,win_plane_2
    readu,5,pin_plane_1
    readu,5,pin_plane_2

    uin_plane_1=uin_plane_1*din_plane_1
    uin_plane_2=uin_plane_2*din_plane_2
    vin_plane_1=vin_plane_1*din_plane_1
    vin_plane_2=vin_plane_2*din_plane_2
    win_plane_1=win_plane_1*din_plane_1
    win_plane_2=win_plane_2*din_plane_2


    dout_plane=0.5*(rebin(din_plane_1,n1_out,n2_out)+rebin(din_plane_2,n1_out,n2_out))
    uout_plane=0.5*(rebin(uin_plane_1,n1_out,n2_out)+rebin(uin_plane_2,n1_out,n2_out))/dout_plane
    vout_plane=0.5*(rebin(vin_plane_1,n1_out,n2_out)+rebin(vin_plane_2,n1_out,n2_out))/dout_plane
    wout_plane=0.5*(rebin(win_plane_1,n1_out,n2_out)+rebin(win_plane_2,n1_out,n2_out))/dout_plane
    pout_plane=0.5*(rebin(pin_plane_1,n1_out,n2_out)+rebin(pin_plane_2,n1_out,n2_out))

    writeu,11,dout_plane
    writeu,12,uout_plane
    writeu,13,vout_plane
    writeu,14,wout_plane
    writeu,15,pout_plane

endfor
close,1
close,2
close,3
close,4
close,5
close,11
close,12
close,13
close,14
close,15


end
