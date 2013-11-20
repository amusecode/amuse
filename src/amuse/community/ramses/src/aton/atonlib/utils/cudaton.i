// Yorick program for reading the binary files produced by aton_debug_dump_.

struct snapaton{
  int ncells;
  int nsource;
  float time;
  pointer egy,flx,xion,temp,dens,src,srcpos;
  float aexp;
}



func plotx(s,n=)
{
  if(is_void(n)) n=(s.ncells-4)/2;
  plfc,1-(*s.xion)(,,n),span(0,1,s.ncells-4)(-:1:s.ncells-4,),span(0,1,s.ncells-4)(,-:1:s.ncells-4),levs=spanl(1e-6,1,100);
}

func plott(s,n=)
{
  if(is_void(n)) n=(s.ncells-4)/2;
  plfc,(*s.temp)(,,n),span(0,1,s.ncells-4)(-:1:s.ncells-4,),span(0,1,s.ncells-4)(,-:1:s.ncells-4),levs=spanl(1e2,1e5,100);
}

func readsnap(name,nbnd=)
{
  if(is_void(nbnd)) nbnd=1;
  
  ff=open(name,"rb");
  adress=0;
  ncells=array(int);_read,ff,adress,ncells;adress+=sizeof(ncells);
  nsource=array(int);_read,ff,adress,nsource;adress+=sizeof(nsource);
  time=array(float);_read,ff,adress,time;adress+=sizeof(time);
  egy=array(float,ncells*ncells*ncells);_read,ff,adress,egy;adress+=sizeof(egy);
  flx=array(float,3*ncells*ncells*ncells);_read,ff,adress,flx;adress+=sizeof(flx);
  xion=array(float,ncells*ncells*ncells);_read,ff,adress,xion;adress+=sizeof(xion);
  temp=array(float,ncells*ncells*ncells);_read,ff,adress,temp;adress+=sizeof(temp);
  if(nsource==0)
    {
      src=array(float,ncells*ncells*ncells);_read,ff,adress,src;adress+=sizeof(src);
    }
  else
    {
      src=array(float,nsource);_read,ff,adress,src;adress+=sizeof(src);
      srcpos=array(int,nsource*3);_read,ff,adress,srcpos;adress+=sizeof(srcpos);
    }
  dens=array(float,ncells*ncells*ncells);_read,ff,adress,dens;adress+=sizeof(dens);
  aexp=array(float);_read,ff,adress,aexp;adress+=sizeof(aexp);
  close,ff;

  egy=reform(egy,[3,ncells,ncells,ncells]);
  egy=egy(nbnd+1:ncells-nbnd,nbnd+1:ncells-nbnd,nbnd+1:ncells-nbnd);
  flx=reform(flx,[4,ncells,ncells,ncells,3]);
  flx=flx(nbnd+1:ncells-nbnd,nbnd+1:ncells-nbnd,nbnd+1:ncells-nbnd,);
  xion=reform(xion,[3,ncells,ncells,ncells]);
  xion=xion(nbnd+1:ncells-nbnd,nbnd+1:ncells-nbnd,nbnd+1:ncells-nbnd);
  temp=reform(temp,[3,ncells,ncells,ncells]);
  temp=temp(nbnd+1:ncells-nbnd,nbnd+1:ncells-nbnd,nbnd+1:ncells-nbnd);
  dens=reform(dens,[3,ncells,ncells,ncells]);
  dens=dens(nbnd+1:ncells-nbnd,nbnd+1:ncells-nbnd,nbnd+1:ncells-nbnd);
  if(nsource==0)
    {
      src=reform(src,[3,ncells,ncells,ncells]);
      src=src(nbnd+1:ncells-nbnd,nbnd+1:ncells-nbnd,nbnd+1:ncells-nbnd);
    }
  //error;
  s=snapaton();
  s.ncells=ncells;
  s.nsource=nsource;
  s.time=time;
  s.aexp=aexp;
  s.egy=&egy;
  s.flx=&flx;
  s.xion=&xion;
  s.temp=&temp;
  s.src=&src;
  s.dens=&dens;
  if(nsource!=0)  s.srcpos=&srcpos;
  return s;
}

func readsnappave(name,nbnd=,npave=,keep=,getflux=)
{
  if(is_void(nbnd)) nbnd=1;
  if(is_void(npave)) npave=2;
  if(is_void(getflux)) getflux=0;
  if(is_void(keep)) keep=0;
  
  ff=open(name,"rb");
  adress=0;
  ncells=array(int);_read,ff,adress,ncells;adress+=sizeof(ncells);
  nsource=array(int);_read,ff,adress,nsource;adress+=sizeof(nsource);
  time=array(float);_read,ff,adress,time;adress+=sizeof(time);
  egy=array(float,ncells*ncells*(npave*(ncells-2*nbnd)+2*nbnd));_read,ff,adress,egy;adress+=sizeof(egy);
  if(getflux)
    {
      flx=array(float,3*ncells*ncells*(npave*(ncells-2*nbnd)+2*nbnd));_read,ff,adress,flx;adress+=sizeof(flx);
    }
  xion=array(float,ncells*ncells*(npave*(ncells-2*nbnd)+2*nbnd));_read,ff,adress,xion;adress+=sizeof(xion);
  temp=array(float,ncells*ncells*(npave*(ncells-2*nbnd)+2*nbnd));_read,ff,adress,temp;adress+=sizeof(temp);

  if(nsource==0)
    {
      src=array(float,ncells*ncells*(npave*(ncells-2*nbnd)+2*nbnd));_read,ff,adress,src;adress+=sizeof(src);
    }
  else
    {
      src=array(float,nsource);_read,ff,adress,src;adress+=sizeof(src);
      srcpos=array(int,nsource*3);_read,ff,adress,srcpos;adress+=sizeof(srcpos);
    }

  dens=array(float,ncells*ncells*(npave*(ncells-2*nbnd)+2*nbnd));_read,ff,adress,dens;adress+=sizeof(dens);
  aexp=array(float);_read,ff,adress,aexp;adress+=sizeof(aexp);
  close,ff;
  
  egy=reform(egy,[3,ncells,ncells,(npave*(ncells-2*nbnd)+2*nbnd)]);
  egy=egy(nbnd+1-keep:ncells-nbnd+keep,nbnd+1-keep:ncells-nbnd+keep,nbnd+1-keep:(npave*(ncells-2*nbnd)+2*nbnd)-nbnd+keep);
  if(getflux)
    {
      flx=reform(flx,[4,ncells,ncells,(npave*(ncells-2*nbnd)+2*nbnd),3]);
      //      if(is_void(keep))flx=flx(nbnd+1:ncells-nbnd,nbnd+1:ncells-nbnd,nbnd+1:(npave*(ncells-2*nbnd)+2*nbnd)-nbnd,);
      flx=flx(nbnd+1-keep:ncells-nbnd+keep,nbnd+1-keep:ncells-nbnd+keep,nbnd+1-keep:(npave*(ncells-2*nbnd)+2*nbnd)-nbnd+keep,);
    }
  
  xion=reform(xion,[3,ncells,ncells,(npave*(ncells-2*nbnd)+2*nbnd)]);
  xion=xion(nbnd+1-keep:ncells-nbnd+keep,nbnd+1-keep:ncells-nbnd+keep,nbnd+1-keep:(npave*(ncells-2*nbnd)+2*nbnd)-nbnd+keep);
  //  if(is_void(keep))xion=xion(nbnd+1:ncells-nbnd,nbnd+1:ncells-nbnd,nbnd+1:(npave*(ncells-2*nbnd)+2*nbnd)-nbnd);

  temp=reform(temp,[3,ncells,ncells,(npave*(ncells-2*nbnd)+2*nbnd)]);
  temp=temp(nbnd+1-keep:ncells-nbnd+keep,nbnd+1-keep:ncells-nbnd+keep,nbnd+1-keep:(npave*(ncells-2*nbnd)+2*nbnd)-nbnd+keep);

  //if(is_void(keep))  temp=temp(nbnd+1:ncells-nbnd,nbnd+1:ncells-nbnd,nbnd+1:(npave*(ncells-2*nbnd)+2*nbnd)-nbnd);

  dens=reform(dens,[3,ncells,ncells,(npave*(ncells-2*nbnd)+2*nbnd)]);
  dens=dens(nbnd+1-keep:ncells-nbnd+keep,nbnd+1-keep:ncells-nbnd+keep,nbnd+1-keep:(npave*(ncells-2*nbnd)+2*nbnd)-nbnd+keep);
  if(is_void(keep))dens=dens(nbnd+1:ncells-nbnd,nbnd+1:ncells-nbnd,nbnd+1:(npave*(ncells-2*nbnd)+2*nbnd)-nbnd);

  if(nsource==0)
    {
      src=reform(src,[3,ncells,ncells,(npave*(ncells-2*nbnd)+2*nbnd)]);
      //if(is_void(keep))src=src(nbnd+1:ncells-nbnd,nbnd+1:ncells-nbnd,nbnd+1:(npave*(ncells-2*nbnd)+2*nbnd)-nbnd);
      src=reform(src,[3,ncells,ncells,(npave*(ncells-2*nbnd)+2*nbnd)]);

    }

  //error;
  s=snapaton();
  s.ncells=ncells;
  s.nsource=nsource;
  s.time=time;
  s.aexp=aexp;
  s.egy=&egy;
  if(getflux) s.flx=&flx;
  s.xion=&xion;
  s.temp=&temp;
  s.src=&src;
  s.dens=&dens;
  if(nsource!=0)  s.srcpos=&srcpos;
  return s;
}


func draw_frame(i)
{
  plfc,1-axion(,,i),indgen(128)(-:1:128,),indgen(128)(,-:1:128),levs=spanl(1e-5,1,50);
  plt,swrite(format="snap #%d",i),0.2,0.5,color="white",height=20.;
  return (i!=dimsof(axion)(0));
}

#if 0
ll=exec("ls -d runs/cosmo.*");
for(i=1;i<=numberof(ll)-1;i++)
{
  if(mod(i,5)==0) i;
  ff=open(ll(i),"rb");
  adress=0;
  ncells=array(int);_read,ff,adress,ncells;adress+=sizeof(ncells);
  if(i==1) axion=array(double,ncells-4,ncells-4,numberof(ll));
  time=array(float);_read,ff,adress,time;adress+=sizeof(time);
  egy=array(float,ncells*ncells*ncells);_read,ff,adress,egy;adress+=sizeof(egy);
  flx=array(float,3*ncells*ncells*ncells);_read,ff,adress,flx;adress+=sizeof(flx);
  xion=array(float,ncells*ncells*ncells);_read,ff,adress,xion;adress+=sizeof(xion);
  temp=array(float,ncells*ncells*ncells);_read,ff,adress,temp;adress+=sizeof(temp);
  close,ff;
  egy=reform(egy,[3,ncells,ncells,ncells]);
  egy=egy(3:ncells-2,3:ncells-2,3:ncells-2);
  flx=reform(flx,[4,3,ncells,ncells,ncells]);
  flx=flx(,3:ncells-2,3:ncells-2,3:ncells-2);
  xion=reform(xion,[3,ncells,ncells,ncells]);
  xion=xion(3:ncells-2,3:ncells-2,3:ncells-2);
  temp=reform(temp,[3,ncells,ncells,ncells]);
  temp=temp(3:ncells-2,3:ncells-2,3:ncells-2);
  axion(,,i)=xion(,,1);
}

