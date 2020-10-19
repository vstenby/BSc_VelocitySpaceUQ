function A = transferMatrix(vpara, vperp, phi, u)
%Constructing the transfer matrix A.
M = length(phi)*length(u);
N = numel(vpara); 
A = zeros(M,N);

du = u(2)-u(1);
if du <= 0, error('du should be positive'), end

if ismatrix(vpara)
    %Here, we assume that vpara is a grid.
    dvpara = diff(vpara(1,1:2));
else
    %Throw an error.
    error('Throw an error');
end

if ismatrix(vperp)
    %Here, we assume that vperp is a grid.
    dvperp = diff(vperp(1:2,1));
else
    error('Throw an error');
end

if dvperp == 0 || dvpara == 0, error('dvperp or dvpara is 0'), end

A = [];
row_idx=1; %Index for filling out the rows of A.

%if any(phi==0), error('phi shouldnt be zero'); end

for phi_k=phi
    %Calculate beforehand for simplicity.
    cos_phik = cos(phi_k/180*pi);
    sin_phik = sin(phi_k/180*pi);
    for i=1:length(u)   
        u_i = u(i);
        %Calculate the i'th row of A^{(k)}
        aik = (1/pi).*(acos((u_i-du/2 - vpara.*cos_phik)./(vperp.*sin_phik)) - ...
                       acos((u_i+du/2 - vpara.*cos_phik)./(vperp.*sin_phik)));
                  
        aik = real(aik);
        %Consider rewriting to increase speed.
        A(row_idx,:) = aik(:);
        row_idx = row_idx + 1;
    end 
end
%Multiply A with the factors.
A = 1/du * dvpara * dvperp * A;
end