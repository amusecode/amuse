%ATH_GETVAR    Reads a desired variable from a given file.
%
%   [TIME,DT,VAR,STATUS] = ATH_GETVAR(GRID,FILENAME,VARNAME) returns the
%   current simulation time, TIME, current time-step, DT, and the data
%   structure VAR containing the desired variable data in an array reshaped
%   according to the dimensions specified in the metadata structure GRID
%   read from the .bin file located at FILENAME.
% 
%   Current options for the variable string VARNAME are
%       'd'             - density
%       'M1','M2','M3'  - momentum components
%       'E'             - total energy
%       'B1','B2','B3'  - magnetic field components
%       'V1','V2','V3'  - velocity components
%       'M','V'         - momentum/velocity magnitude
%       'Emag'          - magnetic energy density (pressure)
%       'Ekin'          - kinetic energy density
%       'SpecEkin'      - specific kinetic energy density
%       'Eint'          - internal energy density
%       'P'             - pressure
%       'Mach'          - Mach number (V/c_s)
%       'A2','A3'       - Magnetic vector potential components
%
%   AUTHOR:  Aaron Skinner
%   LAST MODIFIED:  2/1/2010
function [time,dt,var,status] = ath_getvar(Grid,filename,varname)

status = 0;
time = 0.0;
dt = 0.0;
var = [];

switch(varname)
    case {'d','M1','M2','M3','E','B1','B2','B3'}
        [time,dt,var,status] = ath_readbin(Grid,filename,varname);
        return;
    case {'V1','V2','V3'}
        if (strcmp(varname,'V1'))
            [time,dt,var,status] = ath_readbin(Grid,filename,'M1');
        elseif (strcmp(varname,'V2'))
            [time,dt,var,status] = ath_readbin(Grid,filename,'M2');
        elseif (strcmp(varname,'V3'))
            [time,dt,var,status] = ath_readbin(Grid,filename,'M3');
        end;
        tmp = var;
        [time,dt,var,status] = ath_readbin(Grid,filename,'d');
        var = tmp./d;
        return;
    case 'Emag'
        if (Grid.mhd)
            [time,dt,var,status] = ath_readbin(Grid,filename,'B1');
            tmp = var.^2;
            [time,dt,var,status] = ath_readbin(Grid,filename,'B2');
            tmp = tmp + var.^2;
            [time,dt,var,status] = ath_readbin(Grid,filename,'B3');
            tmp = 0.5*(tmp + var.^2);
            var = tmp;
            return;
        end;
    case {'Ekin','SpecEkin'}
        [time,dt,var,status] = ath_readbin(Grid,filename,'M1');
        tmp = var.^2;
        [time,dt,var,status] = ath_readbin(Grid,filename,'M2');
        tmp = tmp + var.^2;
        [time,dt,var,status] = ath_readbin(Grid,filename,'M3');
        tmp = 0.5*(tmp + var.^2);
        [time,dt,var,status] = ath_readbin(Grid,filename,'d');
        tmp = tmp./var;
        if (strcmp(varname,'SpecEkin'))
            tmp = tmp./var;
        end;
        var = tmp;
        return;
    case {'Eint','P'}
        if (Grid.adiabatic)
            [time,dt,var,status] = ath_readbin(Grid,filename,'E');
            tmp = var;
            [time,dt,var,status] = ath_getvar(Grid,filename,'Ekin');
            tmp = tmp - var;
            [time,dt,var,status] = ath_getvar(Grid,filename,'Emag');
            tmp = tmp - var;
            if (strcmp(varname,'P'))
                tmp = tmp*Grid.gamma_1;
            end;
            var = tmp;
            return;
        end;
    case {'M','V'}
        [time,dt,var,status] = ath_getvar(Grid,filename,'Ekin');
        tmp = 2.0*var;
        [time,dt,var,status] = ath_readbin(Grid,filename,'d');
        if (strcmp(varname,'M'))
            var = sqrt(tmp.*var);
            return;
        end;
        var = sqrt(tmp./var);
        return;
    case 'Mach'
        if (Grid.adiabatic)
            [time,dt,var,status] = ath_getvar(Grid,filename,'Ekin');
            tmp = 2.0*var/Grid.gamma_1;
            [time,dt,var,status] = ath_getvar(Grid,filename,'P');
            var = sqrt(tmp./var);
            return;
        end;
    case 'A2'
        if (Grid.mhd)
            [time,dt,B1,status] = ath_readbin(Grid,filename,'B1');
            [time,dt,B3,status] = ath_readbin(Grid,filename,'B3');
            nx1 = Grid.nx1;
            nx2 = Grid.nx2;
            nx3 = Grid.nx3;
            dx1 = Grid.dx1;
            dx3 = Grid.dx3;
            x1 = Grid.x1zones;
            var = zeros(nx1,nx2,nx3);
            for j = 1:nx2
                for i = 2:nx1
                    x1i = 0.5*(x1(i) + x1(i-1));
                    B3i = (x1(i)*B3(i,j,1) + x1(i-1)*B3(i-1,j,1))/(2*x1i);
                    var(i,j,1) = (x1(i-1)*var(i-1,j,1)+x1i*dx1*B3i)/x1(i);
                end;
                for k = 2:nx3
                    for i = 1:nx1
                        B1i = 0.5*(B1(i,j,k) + B1(i,j,k-1));
                        var(i,j,k) = var(i,j,k-1) - B1i*dx3;
                    end;
                end;
            end;
            return;
        end;
    case 'A3'
        if (Grid.mhd)
            [time,dt,B1,status] = ath_readbin(Grid,filename,'B1');
            [time,dt,B2,status] = ath_readbin(Grid,filename,'B2');
            nx1 = Grid.nx1;
            nx2 = Grid.nx2;
            nx3 = Grid.nx3;
            dx1 = Grid.dx1;
            dx2 = Grid.dx2;
            x1 = Grid.x1zones;
            var = zeros(nx1,nx2,nx3);
            for k = 1:nx3

%                 % INTEGRATE X1, THEN X2
%                 for i = 2:nx1
%                     B2i = 0.5*(B2(i-1,1,k) + B2(i,1,k));
%                     var(i,1,k) = var(i-1,1,k) - B2i*dx1;
%                 end;
%                 for j = 2:nx2
%                     for i = 1:nx1
%                         B1i = 0.5*(B1(i,j-1,k) + B1(i,j,k));
%                         if (Grid.coordsys==-1)  % CARTESIAN
%                             var(i,j,k) = var(i,j-1,k) + B1i*dx2;
%                         elseif (Grid.coordsys==-2)  % CYLINDRICAL
%                             var(i,j,k) = var(i,j-1,k) + B1i*x1(i)*dx2;
%                         end;
%                     end;
%                 end;

                % INTEGRATE X2, THEN X1
                for j = 2:nx2
                    B1i = 0.5*(B1(1,j-1,k) + B1(1,j,k));
                    if (Grid.coordsys==-1)  % CARTESIAN
                        var(1,j,k) = var(1,j-1,k) + B1i*dx2;
                    elseif (Grid.coordsys==-2)  % CYLINDRICAL
                        var(1,j,k) = var(1,j-1,k) + B1i*x1(1)*dx2;
                    end;
                end;
                for i = 2:nx1
                    for j = 1:nx2
                        B2i = 0.5*(B2(i-1,j,k) + B2(i,j,k));
                        var(i,j,k) = var(i-1,j,k) - B2i*dx1;
                    end;
                end;
                
                
            end;
            
            return;
        end;
end;

status = -1;
fprintf('[ath_getvar]:  %s is not a valid variable!\n',varname);
return;