function displayppm(filename)
%
[A,map]=imread(filename);
image(A),colormap(map);
axis image;
