%ATH_PCOLOR    2D Pseudocolor plot
% 
%   STATUS = ATH_PCOLOR(X,Y,C) creates a pseudocolor plot of the 2D data
%   contained in C.  Convert coordinates if necessary.  X and Y must be
%   meshgrid structures containing the coordinates of the edges of the grid
%   zones.
%
%   AUTHOR:  Aaron Skinner
%   LAST MODIFIED:  2/1/2010
function status = ath_pcolor(X,Y,C);

status = 0;

% CHECK COMPATIBILITY OF ARGUMENTS
[nx1,nx2] = size(X);
[ny1,ny2] = size(Y);
[nc1,nc2] = size(C);
if ~(nx1==ny1 && nx2==ny2)
    status = -1;
    fprintf(2,'[ath_pcolor]:  X and Y must be the same size!\n');
    return;
end;
if ~((nc1==nx1 || nc1==nx1-1) && (nc2==nx2 || nc2==nx2-1))
    status = -1;
    fprintf(2,'[ath_pcolor]:  C is not of compatible size!\n');
    return;
end;
    
% MAKE PSEUDOCOLOR PLOT
fvc = surf2patch(X,Y,zeros(nx1,nx2),C);
fvc.facevertexcdata = reshape(C,(nx1-1)*(nx2-1),1);
fvc.vertices = fvc.vertices(:,1:2);
patch(fvc);
shading flat;
view(2);
return;