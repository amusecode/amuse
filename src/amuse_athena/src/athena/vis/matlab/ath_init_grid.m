%ATH_INIT_GRID    Grid initialization
%
%   [GRID,STATUS] = ATH_INIT_GRID(FILENAME) initializes a GRID metadata
%   structure by opening the .bin file specified by the string FILENAME
%   (any cycle will do) and reading the header information only.  Some
%   extra variables are computed for clarity.
%
%   Currently, the elements of the GRID metadata structure are
%       coordsys        - coordinate system (-1=Cartesian, -2=cylindrical)
%       nx1,nx2,nx3     - number of zones in the x1-, x2-, or x3-direction
%       dx1,dx2,dx3     - grid spacing in the x1-, x2-, or x3-direction
%       x1min,x1max     - grid spatial extents in the x1-direction
%       x2min,x2max     - grid spatial extents in the x2-direction
%       x3min,x3max     - grid spatial extents in the x3-direction
%       ndim            - number of (non-singleton) spatial dimensions
%       nvar            - number of stored variables
%       nscalars        - number of stored passive-advection scalars
%       gamma_1         - (gamma-1) for adiabatic gas
%       iso_csound      - isothermal sound speed
%       gravity         - self-gravity (1=enabled, 0=disabled)
%       adiabatic       - equation of state (1=adiabatic, 0=isothermal)
%       mhd             - gas type (1=MHD, 2=hydro only)
%       x1zones         - vector of zone center coords in x1-direction
%       x2zones         - vector of zone center coords in x2-direction
%       x3zones         - vector of zone center coords in x3-direction
%       x1nodes         - vector of edge coords in x1-direction
%       x2nodes         - vector of edge coords in x2-direction
%       x3nodes         - vector of edge coords in x3-direction
%
%   AUTHOR:  Aaron Skinner
%   LAST MODIFIED:  2/1/2010
function [Grid,status] = ath_init_grid(filename)

status = 0;

% PARSE FILENAME AND TEST TO SEE IF .bin
[path,basename,step,ext] = ath_parse_filename(filename);
if (~strcmp(ext,'.bin'))
    fprintf(2,'[ath_init_grid]:  %s is not a .bin file!\n', filename);
    status = -1;
    return;
end;

% OPEN FILE FOR READING
[fid, message] = fopen(filename,'r');
if (fid==-1)
    fprintf(2,'[ath_init_grid]:  %s could not be opened!\n', filename);
    fprintf(2,'%s', message);
    status = -1;
    return;
end;

% READ COORDINATE SYSTEM INFORMATION
coordsys = fread(fid,1,'int');

% READ NUMBER OF DATA POINTS, VARIABLES, SCALARS, ETC.
dat = fread(fid,7,'int');
nx1      = dat(1);
nx2      = dat(2);
nx3      = dat(3);
nvar     = dat(4);
nscalars = dat(5);
ifgrav   = dat(6);
ifpart   = dat(7);
ndim     = (nx1 > 1) + (nx2 > 1) + (nx3 > 1);
if (~(ndim==1 || ndim==2 || ndim==3))
    fprintf(2,'[ath_init_grid]:  %d is an invalid dimension!\n', ndims);
    status = -1;
    return;
end;
if (~((nvar==4) || (nvar==5) || (nvar==7) || (nvar==8)))
    fprintf(2,...
        '[ath_init_grid]:  %d is an invalid number of variables!\n', nvar);
    status = -1;
    return;
end;

% READ (Gamma-1), ISOTHERMAL SOUND SPEED, TIME, AND dt
dat         = fread(fid,2,'float');
gamma_1     = dat(1);
iso_csound  = dat(2);
time_offset = ftell(fid);  % GET POSITION OF time, dt
dat         = fread(fid,2,'float');
time        = dat(1);  % READ IN, BUT NOT USED
dt          = dat(2);  % READ IN, BUT NOT USED

% READ X1,X2,X3 COORDINATES
x1zones = fread(fid,nx1,'float');
x2zones = fread(fid,nx2,'float');
x3zones = fread(fid,nx3,'float');
data_offset = ftell(fid);

% CLOSE FILE
status = fclose(fid);
if (status == -1)
    fprintf(2,'[ath_init_grid]:  %s could not be closed!\n', filename);
end;

% COMPUTE SOME DERIVED QUANTITIES
dx1 = 0; dx2 = 0; dx3 = 0;
if (nx1>1) 
    dx1 = (max(x1zones)-min(x1zones))/(nx1-1);
end;
if (nx2>1) 
    dx2 = (max(x2zones)-min(x2zones))/(nx2-1);
end;
if (nx3>1) 
    dx3 = (max(x3zones)-min(x3zones))/(nx3-1);
end;
% X1MIN, X1MAX, ETC. ARE THE ABSOLUTE LIMITS OF THE GRID
x1min = min(x1zones) - 0.5*dx1;  x1max = max(x1zones) + 0.5*dx1;
x2min = min(x2zones) - 0.5*dx2;  x2max = max(x2zones) + 0.5*dx2;
x3min = min(x3zones) - 0.5*dx3;  x3max = max(x3zones) + 0.5*dx3;
x1nodes = linspace(x1min,x1max,nx1+1);
x2nodes = linspace(x2min,x2max,nx2+1);
x3nodes = linspace(x3min,x3max,nx3+1);

% INITIALIZE GRID STRUCTURE
Grid.coordsys    = coordsys;
Grid.nx1         = nx1;
Grid.nx2         = nx2;
Grid.nx3         = nx3;
Grid.dx1         = dx1;
Grid.dx2         = dx2;
Grid.dx3         = dx3;
Grid.x1min       = x1min;
Grid.x1max       = x1max;
Grid.x2min       = x2min;
Grid.x2max       = x2max;
Grid.x3min       = x3min;
Grid.x3max       = x3max;
Grid.ndim        = ndim;
Grid.nvar        = nvar;
Grid.nscalars    = nscalars;
Grid.gamma_1     = gamma_1;
Grid.iso_csound  = iso_csound;
Grid.gravity     = (ifgrav == 1);
Grid.particles   = (ifpart == 1);
Grid.adiabatic   = (nvar==5 || nvar==8);
Grid.mhd         = (nvar==7 || nvar==8);
Grid.x1zones     = x1zones;
Grid.x2zones     = x2zones;
Grid.x3zones     = x3zones;
Grid.x1nodes     = x1nodes;
Grid.x2nodes     = x2nodes;
Grid.x3nodes     = x3nodes;
Grid.time_offset = time_offset;
Grid.data_offset = data_offset;

return;