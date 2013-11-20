
pro pall,nin,typin

car=getcarnum(nin)
charin=car(nin-1)

window,0,xs=512,ys=512
window,1,xs=512,ys=512
window,2,xs=512,ys=512
window,3,xs=512,ys=512
window,4,xs=512,ys=512
window,5,xs=512,ys=512

loadct,33

rd_im,dx,file='dens_'+charin+'_'+typin+'_dirx.map'
wset,0 & mycontour,dx,/log,/tab,ncont=200
rd_im,tx,file='temp_'+charin+'_'+typin+'_dirx.map'
wset,1 & mycontour,tx,/log,/tab,ncont=200
rd_im,zx,file='met_'+charin+'_'+typin+'_dirx.map'
wset,2 & mycontour,zx,/log,/tab,ncont=200
rd_im,bx,file='blast_'+charin+'_'+typin+'_dirx.map'
wset,3 & mycontour,bx,/log,/tab,ncont=200,min=1d-4
rd_im,sx,file='star_'+charin+'_'+typin+'_dirx.map'
wset,4 & mycontour,sx,/log,/tab,ncont=200,min=1d-16
rd_im,cx,file='dark_'+charin+'_'+typin+'_dirx.map'
wset,5 & mycontour,cx,/log,/tab,ncont=200,min=100.
read,itest

rd_im,dy,file='dens_'+charin+'_'+typin+'_diry.map'
wset,0 & mycontour,dy,/log,/tab,ncont=200
rd_im,ty,file='temp_'+charin+'_'+typin+'_diry.map'
wset,1 & mycontour,ty,/log,/tab,ncont=200
rd_im,zy,file='met_'+charin+'_'+typin+'_diry.map'
wset,2 & mycontour,zy,/log,/tab,ncont=200
rd_im,by,file='blast_'+charin+'_'+typin+'_diry.map'
wset,3 & mycontour,by,/log,/tab,ncont=200,min=1d-4
rd_im,sy,file='star_'+charin+'_'+typin+'_diry.map'
wset,4 & mycontour,sy,/log,/tab,ncont=200,min=1d-16
rd_im,cy,file='dark_'+charin+'_'+typin+'_diry.map'
wset,5 & mycontour,cy,/log,/tab,ncont=200,min=100.
read,itest

rd_im,dz,file='dens_'+charin+'_'+typin+'_dirz.map'
wset,0 & mycontour,dz,/log,/tab,ncont=200
rd_im,tz,file='temp_'+charin+'_'+typin+'_dirz.map'
wset,1 & mycontour,tz,/log,/tab,ncont=200
rd_im,zz,file='met_'+charin+'_'+typin+'_dirz.map'
wset,2 & mycontour,zz,/log,/tab,ncont=200
rd_im,bz,file='blast_'+charin+'_'+typin+'_dirz.map'
wset,3 & mycontour,bz,/log,/tab,ncont=200,min=1d-4
rd_im,sz,file='star_'+charin+'_'+typin+'_dirz.map'
wset,4 & mycontour,sz,/log,/tab,ncont=200,min=1d-16
rd_im,cz,file='dark_'+charin+'_'+typin+'_dirz.map'
wset,5 & mycontour,cz,/log,/tab,ncont=200,min=100.

end
