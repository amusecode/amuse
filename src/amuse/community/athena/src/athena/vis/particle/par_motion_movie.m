clear all;
figlab = 2;
figM = figure(figlab); hold off;
fname=''; % file name
fn1=1;    % first figure index
fn2=200;  % lasr figure index

for i=fn1:fn2
 if (i<10)
     numlab = ['000',num2str(i)];
     else if (i<100)
         numlab = ['00',num2str(i)];
         else if (i<1000)
             numlab = ['0',num2str(i)];
             else
                 numlab = num2str(i);
             end
         end
 end
 
 fid=fopen(['../../bin/',fname,'.',numlab,'.lis'],'rb');
 
 % Read the coordinate limits
 coorlim = fread(fid,12,'float');
 x1l = coorlim(1); x1u = coorlim(2); x2l = coorlim(3); x2u = coorlim(4); x3l = coorlim(5); x3u = coorlim(6);
 x1dl = coorlim(7); x1du = coorlim(8); x2dl = coorlim(9); x2du = coorlim(10); x3dl = coorlim(11); x3du = coorlim(12);
 
 % Read the time
 time = fread(fid,1,'float');
 dt   = fread(fid,1,'float');
 
 % Read the particle number
 n =  fread(fid,1,'int64')
 
 % Read all the particle information
 for j=1:n
     parinfo = fread(fid,8,'float');
     x1(j) = parinfo(1); x2(j) = parinfo(2); x3(j) = parinfo(3);
     v1(j) = parinfo(4); v2(j) = parinfo(5); v3(j) = parinfo(6);
     rad(j) = parinfo(7); mas(j) = parinfo(8);
     pid(j) = fread(fid,1,'int64');
     cpuid(j) = fread(fid,1,'int32');
 end
 
 plot(x1,x2,'.'); hold off;
 axis([x1l x1u x2l x2u]);
 getframe(figM);
 %view(0,0);
 grid on;
 
 fclose(fid);
end
