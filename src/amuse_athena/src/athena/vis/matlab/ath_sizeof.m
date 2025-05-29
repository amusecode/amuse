%ATH_SIZEOF    Returns the size, in bytes, of a given data type.
% 
%   NBTYES = ATH_SIZEOF(DATATYPE)
%
%   AUTHOR:  Aaron Skinner
%   LAST MODIFIED:  2/1/2010
function nbytes = ath_sizeof(datatype)

try
    z = zeros(1,datatype);
catch
    error('[ath_sizeof]:  Unsupported data type!');
end

w = whos('z');
nbytes = w.bytes;