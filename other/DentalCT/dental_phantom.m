function [p, info, mask] = dental_phantom(varargin)
%   dental_phantom makes a dental phantom with the specified values in the
%   vector. 
%   
%   Usage: 
%   [p, info, mask] = dental_phantom()      default 256 x 256
%   [p, info, mask] = dental_phantom(N)     default N x N
%   [p, info, mask] = dental_phantom(info)  custom dental phantom
%   
%   info is a struct containing the parameters of the phantom.
%
%   info.n : The dimension of the output image (n x n)
%
%   info.I : Vector specifying the intensities of the different areas
%     info.I(1): soft tissue
%     info.I(2): teeth
%     info.I(3): spine
%     info.I(4): bone marrow
%
%   info.geom : Vector specifying the geometric properties of the phantom.
%       All values in info.geom should be between 0 and 1.

switch length(varargin)
    case 0
        %The default values are used. 
        info.n = 256;
        info.geom = [0.6, 0.35, 0.3, 0.3, 0.25, 0.15, 0.05, 0.05, 0.2, 0.2, 0.18, 0.35, 0.3, 0.05, 0.04];
        info.I = [0.18, 0.8, 0.65, 0];
    case 1
        %User input is used for generating the dental phantom.
        info = varargin{1};
        if length(info) == 1
           %Then we assume that N is given, and we set the remaining info
           %to the default values.
           n = info; clear info;
           info.n = n;
           info.geom = [0.6, 0.35, 0.3, 0.3, 0.25, 0.15, 0.05, 0.05, 0.2, 0.2, 0.18, 0.35, 0.3, 0.05, 0.04];
           info.I = [0.18, 0.8, 0.65, 0];
        end
    otherwise
        error('Too many inputs or wrong number of inputs')
end
[n, I, geom_mat] = extract_values(info);
g = info.geom;
p = zeros([n, n]);
mask = zeros([n, n]);
%g(14) er tooth radius for de små tænder.
 

xax =  ((0:n-1)-(n-1)/2 ) / ((n-1)/2); 
x = repmat(xax, n, 1);  % x coordinates, the y coordinates are rot90(xg)
                

idx = x <= g(1) & -g(1) <= x & rot90(x) <= g(2) & -g(2) <= rot90(x);
  
p(idx) = info.I(1);
mask(idx) = 1;

for k=1:6
    idx = [];
    maskidx = [];
    asq = geom_mat(k,2)^2; % a^2
    bsq = geom_mat(k,3)^2; % b^2
    y0 = geom_mat(k,5);    % y offset
    y=rot90(x)-y0;  
    if k==1||k==2
        idx=find((x.^2)./asq + (y.^2)./bsq <= 1 & y < 0); 
        if k==1, maskidx = idx; end
    elseif k==3||k==4
        idx=find((x.^2)./asq + (y.^2)./bsq <= 1 & y > 0);
        if k==3, maskidx = idx; end
    else
        idx=find((x.^2)./asq + (y.^2)./bsq <= 1);
    end
    p(idx) = geom_mat(k, 1);
    mask(maskidx) = 1;
end

%Square part of the mouth
idx = x <= g(12) & -g(12) <= x & rot90(x) <= -g(2)+g(11) & -g(2)-g(11) <= rot90(x);  
p(idx) = 0;


%Adding the corner teeth.
%Right corner teeth:
idx = sqrt(((g(12)-g(14)) - x).^2 + (-g(2)+g(11)-g(14) - rot90(x)).^2) <= g(14); p(idx) = I(2);

%idx = sqrt(((g(12)-g(14)) - x).^2 + (-g(2)+g(11)-3.5*g(14) - rot90(x)).^2) <= g(14); 
%p(idx) = 1;
%Left corner teeth:
idx = sqrt(((-g(12)+g(14)) - x).^2 + (-g(2)+g(11)-g(14) - rot90(x)).^2) <= g(14); p(idx) = I(2);
idx = sqrt(((-g(12)+g(14)) - x).^2 + (-g(2)+g(11)-3*g(14) - rot90(x)).^2) <= g(14); p(idx) = I(2);

%Adding teeth in the outer mouth.
r = @(u) [(g(12) - g(15)) * cos(u); (g(13) - g(15))*sin(u) - g(2)-g(11)];

teeth = linspace(0, -1*pi, 10);
for i=1:length(teeth)
    c = r(teeth(i));
    idx = sqrt((c(1) - x).^2 + (c(2) - rot90(x)).^2) <= g(15);
    p(idx) = I(2);
end

mask = logical(mask);
%Auxilliary functions

function [n, I, geom_mat] = extract_values(info)
%geom_mat is a matrix for drawing the ellipses.
n = info.n;
I = info.I;
g = info.geom;
geom_mat = zeros([2 5]);
geom_mat(1,:) = [I(1), g(1), 2*g(3), 0, -g(2)];                  %Cheek, soft tissue
geom_mat(2,:) = [0,    g(12), g(13), 0, -g(2)-g(11)];            %The front of the mouth, background
geom_mat(3,:) = [I(1), g(1), 2*g(4), 0,  g(2)];                  %Back of the head, soft tissue
geom_mat(4,:) = [0,    g(9),  g(10), 0, -g(2)+g(11)];            %The back of the throat, background
geom_mat(5,:) = [I(3), g(5),   g(6), 0,  g(2)+0.5*g(4)];         %The spine
geom_mat(6,:) = [I(4), g(7),   g(8), 0,  g(2)+0.5*g(4)];         %The bone marrow
end





end