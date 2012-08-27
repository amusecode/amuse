%ATH_LINEOUT_XYZ    Lineout resampling using xyz coordinates
% 
%   [X,Y,STATUS] = ATH_LINEOUT_XYZ(GRID,VAR,X,Y,Z,INTERP) creates two 1D
%   vectors X and Y containing a sampling of the variable VAR using
%   coordinates from the GRID metadata structure.  This function is
%   basically a wrapper for ath_lineout_ijk, but has the added option to
%   set the interpolation flag INTERP to 1 so that if a requested spatial
%   coordinate is closer to a zone edge than a zone center, ATH_LINEOUT_XYZ
%   will return an average of the data on either side of the edge.  This is
%   useful for grids with an even number of zones in a given direction.
% 
%   E.g. Read in the 3D density data array contained in my_data_file and
%   produce a lineout plot along the entire x1-direction, where x2=x3=0.0
%   are fixed. Note that X will contain the zone center coordinates along
%   the x1-direction and Y will contain the density sampled at (x1,0.0,0.0)
%   for all values of x1 in my_grid (i.e. in the my_data array). Since the
%   interpolation is set to 1, if x2=0.0 or x3=0.0 are the coordinates of a
%   zone edge, then ath_lineout_xyz will return an average of the data at
%   both adjacent zone centers.
%
%       [time,dt,my_data,status] = ath_getvar(my_grid,my_data_file,'d');
%       x = my_grid.x1zones;
%       y = 0.0;
%       z = 0.0;
%       [X,Y,status] = ath_lineout_xyz(my_grid,my_data,x,y,z,1);
%       plot(X,Y);
%
%   See also ATH_LINEOUT_IJK
%
%   AUTHOR:  Aaron Skinner
%   LAST MODIFIED:  2/1/2010
function [X,Y,status] = ath_lineout_xyz(Grid,var,x,y,z,interp)

[i,j,k,onfacei,onfacej,onfacek] = xyz_to_ijk(Grid,x,y,z);

[X,Y,status] = lineout_ijk(Grid,var,i,j,k);

if (interp && onfacei)
    [X,Y2,status] = lineout_ijk(Grid,var,i+1,j,k);
    Y = 0.5*(Y+Y2);
end;
if (interp && onfacej)
    [X,Y2,status] = lineout_ijk(Grid,var,i,j+1,k);
    Y = 0.5*(Y+Y2);
end;
if (interp && onfacek)
    [X,Y2,status] = lineout_ijk(Grid,var,i,j,k+1);
    Y = 0.5*(Y+Y2);
end;

if (status == -1)
    printf(2,'[ath_lineout_xyz]:  ath_lineout_ijk returned in error!');
end;

return;