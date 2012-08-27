%ATH_TEXT    Plot a text string using map-direction strings.
%
%   ATH_TEXT(AX,POS,TEXTSTRING,FONTSIZE) plots TEXTSTRING on the axes AX at
%   the position indicated by the map-direction string POS.
%
%   Currently, the supported strings are
%       'north'
%       'northeast'
%       'northwest'
%
%   AUTHOR:  Aaron Skinner
%   LAST MODIFIED:  2/1/2010
function ath_text(ax,pos,textstring,fontsize)

a = 0.02;
margin = 2;
dx = a*(ax(2) - ax(1));
dy = a*(ax(4) - ax(3));
xmin = ax(1) + dx;
xmax = ax(2) - dx;
ymin = ax(3) + dy; 
ymax = ax(4) - dy;
xavg = 0.5*(xmin + xmax);
yavg = 0.5*(ymin + ymax);

switch pos
    case 'north'
        x = xavg;
        y = ymax;
        ha = 'center';
        va = 'top';
    case 'northeast'
        x = xmax;
        y = ymax;
        ha = 'right';
        va = 'top';
    case 'northwest'
        x = xmin;
        y = ymax;
        ha = 'left';
        va = 'top';
    otherwise
        x = xavg;
        y = ymax;
        ha = 'center';
        va = 'top';        
end;

text(x,y,0,textstring,...
    'HorizontalAlignment',ha,'VerticalAlignment',va,...
    'FontSize',fontsize,'BackgroundColor','w','Margin',margin);
        
