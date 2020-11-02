function varargout = welford(Mk, Sk, count, s, varargin)
%Count is the number of samples seen so far.

switch s
    case 'update'
        x = varargin{1};
        
        Mprev = Mk;
        Sprev = Sk;
        count = count + 1;
        
        Mk = Mprev + (x - Mprev)/count; 
        Sk = Sprev + (x - Mprev).*(x - Mk);
        
        varargout{1} = Mk;
        varargout{2} = Sk;
        varargout{3} = count;
        
    case 'finalize'
        
        xmean = Mk; %mean
        xstd = (Sk/(count)).^(1/2);

        varargout{1} = xmean;
        varargout{2} = xstd;
end
end