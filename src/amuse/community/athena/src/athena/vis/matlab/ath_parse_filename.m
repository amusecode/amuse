%ATH_PARSE_FILENAME    Break up a full-path filename into its component
%parts to check the extension, make it more readable, and extract the step
%number.  
%
%   [PATH,BASENAME,STEP,EXT] = ATH_PARSE_FILENAME(FILENAME)
%
%   E.g. If FILENAME='/home/Blast.0000.bin', then PATH='/home',
%   BASENAME='Blast', STEP=0, and EXT='.bin'.
%
%   AUTHOR:  Aaron Skinner
%   LAST MODIFIED:  2/1/2010
function [path,basename,step,ext] = ath_parse_filename(filename)

[path, file, ext, versn] = fileparts(filename);
[basename,remain] = strtok(file,'.');
while 1
    [token,remain] = strtok(remain,'.');
    if (isempty(remain))
        step = sscanf(token,'%d');
        break;
    end;
    basename = strcat(basename,token,'.');
end;
