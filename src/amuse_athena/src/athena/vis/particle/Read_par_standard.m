clear all;
figlab = 5;
figure(figlab);
fname=''; % Input the file name here

 fid=fopen(['../../bin/',fname,'.lis'],'rb');
 
 % Read the coordinate limits
 coorlim = fread(fid,12,'float');
 x1l = coorlim(1);     x1u = coorlim(2);
 x2l = coorlim(3);     x2u = coorlim(4);
 x3l = coorlim(5);     x3u = coorlim(6);
 x1dl = coorlim(7);    x1du = coorlim(8);
 x2dl = coorlim(9);    x2du = coorlim(10);
 x3dl = coorlim(11);   x3du = coorlim(12);
 
 % Read the time
 time = fread(fid,1,'float');
 dt   = fread(fid,1,'float');
 
 % Read the particle number
 n =  fread(fid,1,'int64');
 
 % Read all the particle information
 for i=1:n
     parinfo = fread(fid,8,'float');
     x1(i) = parinfo(1);
     x2(i) = parinfo(2);
     x3(i) = parinfo(3);
     v1(i) = parinfo(4);
     v2(i) = parinfo(5);
     v3(i) = parinfo(6);
     rad(i) = parinfo(7);
     mas(i) = parinfo(8);
     pid(i) = fread(fid,1,'int64');
     cpuid(i) = fread(fid,1,'int32');
 end
 
 plot(x1,x2,'.');
 %axis([x1l x1u x2l x2u]);
 %view(0,0);
 grid on;
