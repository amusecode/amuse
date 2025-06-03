%ATH_CONSTRUCT_FILENAME    Construct a full-path filename from its
%component parts.
%
%   FILENAME = ATH_CONSTRUCT_FILENAME(PATH,BASENAME,STEP) 
%
%   E.g. for the file '/home/Blast.0000.bin', PATH is '/home', BASENAME is
%   'Blast', and STEP is 0.  Files from parallel runs must be previously
%   joined together.
%
%   AUTHOR:  Aaron Skinner
%   LAST MODIFIED:  2/1/2010
function filename = ath_construct_filename(path,basename,step)

filename = strcat(path,'/',basename,'.',sprintf('%04d',step),'.bin');