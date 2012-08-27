%ATH_SLICE_XYZ    Slice resampling using xyz coordinates
%
%   [X,Y,Z,STATUS] = ATH_SLICE_XYZ(GRID,VAR,X,Y,Z,INTERP) creates a 2D
%   slice of a given variable VAR using spatial coordinates.  Exactly two
%   of X, Y, Z must be vectors of coordinates, the remaining one being the
%   fixed value through which to take the slice.  If it is determined that
%   this fixed value lies between grid zones (e.g. for an even grid), the
%   interpolation flag INTERP may be set to 1 and the average of the two
%   adjacent slices will be returned.  Note that in the case of
%   (anti)-symmetry for an even grid, this interpolation is exact.  
% 
%   See also ATH_SLICE_IJK
% 
%   AUTHOR:  Aaron Skinner
%   LAST MODIFIED:  2/1/2010
function [X,Y,Z,status] = ath_slice_xyz(Grid,var,x,y,z,interp)

% TRANSFORM IN CYLINDRICAL CASE?
transform = true;

X = [];  Y = [];  Z = [];
status = 0;

% DETERMINE WHICH DIRECTION TO SLICE
lenx = length(x);  leny = length(y);  lenz = length(z);
if (lenx==1 && leny>1 && lenz>1) 
    dir = 1;
elseif (lenx>1 && leny==1 && lenz>1) 
    dir = 2;
elseif (lenx>1 && leny>1 && lenz==1) 
    dir = 3;
else
    status = -1;
    fprintf(2,'[ath_slice_xyz]:  Exactly two of x,y,z must be vectors!\n');
    return;
end;

[i,j,k,onfacei,onfacej,onfacek] = ath_xyz_to_ijk(Grid,x,y,z);

% COMPUTE SLICE
switch (dir)
    case 1
        [X,Y] = meshgrid(y,z);
        Z = squeeze(var(i,:,:))';
        if (onfacei && interp)
            Z = 0.5*(Z+squeeze(var(i+1,:,:))');
        end;
    case 2
        [X,Y] = meshgrid(x,z);
        Z = squeeze(var(:,j,:))';
        if (onfacej && interp)
            Z = 0.5*(Z+squeeze(var(:,j+1,:))');
        end;
    case 3
        if ((Grid.coordsys == -2) && transform) % CYLINDRICAL
            [R,PHI] = meshgrid(x,y);
            [X,Y] = pol2cart(PHI,R);
        else
            [X,Y] = meshgrid(x,y);
        end;
        Z = squeeze(var(:,:,k))';
        if (onfacek && interp)
            Z = 0.5*(Z+squeeze(var(:,:,k+1))');
        end;
end;

return;