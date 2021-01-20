function evals = varargin_to_eval(input,validvars)
%If valid vares are passed, then we should check whether or not
%each variable is valid.

if nargin == 1
   check_valid = 0;
else
    check_valid = 1;
end

if mod(length(input),2) == 1, error('Odd number of optional inputs.'), end

vars  = {input{1:2:end}};
nvars = length(vars);
idx = 1:2:2*nvars;
evals = cell(nvars,1);

%Check if each of the variables are a valid variable.
for i=1:nvars
   var = char(vars{i}); 
   var = strtrim(var); %In case the user writes a variable with spaces.
   if check_valid
     if ~any(strcmpi(var,validvars)), error('%s is not a valid input.', var), end
   end
   vars{i} = var;
   evals{i} = sprintf('%s = varargin{%d+1};',var,idx(i));
end
end