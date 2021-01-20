function val = check_mosek()
% A wrapper for mosekdiag, checking if mosek works correctly.
%
% Usage: 
%    ``check_mosek()``
%
% Output: 
%    ``1`` if working correctly, ``0`` otherwise.
%
try
    str = evalc('mosekdiag');
    if contains(str, 'mosekopt works OK.') || contains(str, 'mosekopt is working correctly')
       val = 1;
    else
       val = 0;
    end
catch
    val = 0; 
end
end