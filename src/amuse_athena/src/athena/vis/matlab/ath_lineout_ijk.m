%ATH_LINEOUT_IJK    Lineout resampling using ijk coordinates
% 
%   [X,Y,STATUS] = ATH_LINEOUT_IJK(GRID,VAR,I,J,K) creates two 1D vectors X
%   and Y containing a sampling of the variable VAR using coordinates from
%   the GRID metadata structure.  Exactly one of I, J, or K must be a
%   vector of indices, the remaining two being the fixed indices through
%   which to take the lineout sample.  The vector X contains the spatial
%   coordinates in the indicated direction, and Y contains the resampled
%   data.  Typically, this command will follow a call to the ath_getvar
%   command and will precede a call to MATLAB's plot command.
% 
%   E.g. Read in the 3D density data array contained in my_data_file and
%   produce a lineout plot along the entire i-direction, where j=2 and k=3
%   are fixed. Note that X will contain the zone center coordinates along
%   the x1-direction and Y will contain the density sampled at (i,2,3) for
%   all values of i in my_grid (i.e. in the my_data array).
%
%       [time,dt,my_data,status] = ath_getvar(my_grid,my_data_file,'d');
%       i = 1:my_grid.nx1;
%       j = 2;
%       k = 3;
%       [X,Y,status] = ath_lineout_ijk(my_grid,my_data,i,j,k);
%       plot(X,Y);
%
%   See also ATH_LINEOUT_XYZ
%
%   AUTHOR:  Aaron Skinner
%   LAST MODIFIED:  2/1/2010
function [X,Y,status] = ath_lineout_ijk(Grid,var,i,j,k)

X = [];  Y = [];
status = 0;

[x,y,z] = ath_ijk_to_xyz(Grid,i,j,k);

lenx = length(x);
leny = length(y);
lenz = length(z);

if (lenx>1 && leny==1 && lenz==1)
    X = x;
    Y = squeeze(var(:,j,k));
elseif (lenx==1 && leny>1 && lenz==1)
    X = y;
    Y = squeeze(var(i,:,k));
elseif (lenx==1 && leny==1 && lenz>1)
    X = z;
    Y = squeeze(var(i,j,:));
else
    status = -1;
    fprintf(2, ...
        '[ath_lineout_ijk]:  Exactly one of i,j,k must be a vector!\n');
    return;
end;

return;