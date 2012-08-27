%ATH_SLICE_IJK    Slice resampling using ijk coordinates
%
%   [X,Y,Z,STATUS] = ATH_SLICE_IJK(GRID,VAR,I,J,K) creates a 2D slice of a
%   given variable VAR using spatial coordinates.  Exactly two of I, J, K
%   must be vectors of grid indices, the remaining one being the fixed
%   index value through which to take the slice.  This function is
%   basically a wrapper for the function ath_slice_xyz.
% 
%   See also ATH_SLICE_XYZ
% 
%   AUTHOR:  Aaron Skinner
%   LAST MODIFIED:  2/1/2010
function [X,Y,Z,status] = ath_slice_ijk(Grid,var,i,j,k)

[x,y,z] = ath_ijk_to_xyz(Grid,i,j,k);

[X,Y,Z,status = ath_slice_xyz(Grid,var,x,y,z);

return;