%ATH_XYZ_TO_IJK    Convert from xyz to ijk coordinates
% 
%   [I,J,K] = ATH_XYZ_TO_IJK(GRID,X,Y,Z) converts the coordinate(s) of the
%   spatial (but not necessarily geometric) center(s) (X,Y,Z) to the grid
%   coordinate(s) (I,J,K), like an inverse of Athena's cc_pos().
%
%   AUTHOR:  Aaron Skinner
%   LAST MODIFIED:  2/1/2010
function [i,j,k,onfacei,onfacej,onfacek] = ath_xyz_to_ijk(Grid,x,y,z)

tol = 1E-6;
onfacei = 0;
onfacej = 0;
onfacek = 0;

i = 1; j = 1; k = 1;
if (Grid.dx1>0) 
    i = (x-Grid.x1min)./Grid.dx1;
    if (abs(i-round(i)) < tol)
        onfacei = 1;
        i = round(i);
    else
        i = ceil(i);
    end;
    i = max(i,1);  i = min(i,Grid.nx1);
end;
if (Grid.dx2>0)
%     j = round((y-Grid.x2min)./Grid.dx2-0.5)+1;
    j = (y-Grid.x2min)./Grid.dx2;
    if (abs(j-round(j)) < tol)
        onfacej = 1;
        j = round(j);
    else
        j = ceil(j);
    end;
    j = max(j,1);  j = min(j,Grid.nx2);
end;
if (Grid.dx3>0)
%     k = round((z-Grid.x3min)./Grid.dx3-0.5)+1;
    k = (z-Grid.x3min)./Grid.dx3;
    if (abs(k-round(k)) < tol)
        onfacek = 1;
        k = round(k);
    else
        k = ceil(k);
    end;
    k = max(k,1);  k = min(k,Grid.nx3);
end;