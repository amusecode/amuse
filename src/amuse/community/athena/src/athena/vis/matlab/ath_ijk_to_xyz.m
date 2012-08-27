%ATH_IJK_TO_XYZ    Convert from ijk to xyz coordinates
% 
%   [X,Y,Z] = ATH_IJK_TO_XYZ(GRID,I,J,K) converts the grid coordinate(s)
%   (I,J,K) to the coordinate(s) of the corresponding spatial (but not
%   necessarily geometric) center(s) (X,Y,Z), like Athena's cc_pos().
%
%   AUTHOR:  Aaron Skinner
%   LAST MODIFIED:  2/1/2010
function [x,y,z] = ath_ijk_to_xyz(Grid,i,j,k)

x = Grid.x1min + (i-1+0.5)*Grid.dx1;
y = Grid.x2min + (j-1+0.5)*Grid.dx2;
z = Grid.x3min + (k-1+0.5)*Grid.dx3;