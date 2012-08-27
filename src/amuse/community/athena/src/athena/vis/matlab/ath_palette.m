%ATH_PALETTE    Useful colormaps
%
%   MAP = ATH_PALETTE(NAME,SIZE) creates a colormap MAP of size SIZE-by-3,
%   where NAME is a string representing the desired colormap.
%
%   Current NAME options include:
%       'xray'      - An inverted grayscale colormap 
%       'hot'       - The standard 'hot' colormap as used in VisIt
%       'jh_colors' - John Hawley's famous red-gray-blue colormap (returns  
%                       a map of size floor(SIZE/8) )
%
%   AUTHOR:  Aaron Skinner
%   LAST MODIFIED:  2/1/2010
function map = ath_palette(name,size)

switch (name)
    case 'xray'
        map = flipud(gray(size));
    case 'hot'
        map = hsv2rgb([linspace(2/3,0,size)',ones(size,1),ones(size,1)]);
    case 'jh_colors'  % JOHN HAWLEY'S FAMOUS RED-GRAY-BLUE COLORMAP
        n = floor(size/8);
        x = linspace(0,.5,n+1)';
        x = x(1:n,1);
        o = ones(n,1);
        z = zeros(n,1);
        map =  [z     , z     , 0.5+x ;
                z     , x     , o     ;
                z     , 0.5+x , o     ;
                x     , 1.0-x , 1.0-x ;
                0.5+x , 0.5+x , 0.5-x ;
                o     , 1.0-x , z     ;
                o     , 0.5-x , z     ;
                1.0-x , z     , z      ];
    otherwise
        fprintf(2,'[ath_palette]:  %s is not a known colormap!\n',name);
        map = jet(size);
end;

return;