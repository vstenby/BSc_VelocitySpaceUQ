function val = check_mosek()
% Calls the mosekdiag but silently.
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