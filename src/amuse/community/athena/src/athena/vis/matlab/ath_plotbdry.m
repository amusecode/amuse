%ATH_PLOTBDRY    Plot the 2D boundary of given grid.
%
%   AUTHOR:  Aaron Skinner
%   LAST MODIFIED:  2/1/2010
function ath_plotbdry(grid)

X = grid.x1nodes;
Y = grid.x2nodes;
Xbdry = [X,max(X)*ones(size(Y)),fliplr(X),min(X)*ones(size(Y))];
Ybdry = [min(Y)*ones(size(X)),Y,max(Y)*ones(size(X)),fliplr(Y)];

if (grid.coordsys == -2)
    [Xbdry,Ybdry] = pol2cart(Ybdry,Xbdry);
end;

plot(Xbdry,Ybdry,'k');