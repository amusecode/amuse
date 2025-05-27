function animateppm(filename,nimages)
%
[basename,remain] = strtok(filename,'.')
[start,extension] = strtok(remain,'.')
nstart=str2num(start);
for k=nstart:(nstart+nimages)
    filename = strcat(basename,'.',sprintf('%04d',k),extension);
    [A,map]=imread(filename);
    image(A),colormap(map)
    axis image;
    M = getframe;
end;
