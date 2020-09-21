function A = transfer_matrix(vpara, vperp, phivec, uvec, dures)
%Constructing the transfer matrix A.
M = length(phivec)*length(uvec);
N = numel(vpara); 
A = zeros(M,N);

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
for phi=phivec
    for u=uvec
        
        gamma1=acos((u-dures./2-cos(phi/180*pi).*vpara)./(sin(phi/180*pi).*vperp));
        gamma2=acos((u+dures./2-cos(phi/180*pi).*vpara)./(sin(phi/180*pi).*vperp));
        wv=real((gamma1-gamma2))/pi/dures*dvpara*dvperp;
       
        A(row_idx,:) = wv(:);
        
        row_idx = row_idx + 1;
    end 
end
end