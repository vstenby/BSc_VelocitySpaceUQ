function [xalpha, alpha] = TikhonovNonNeg(transfer_matrix,G,vpara,vperp,L)


ax=[-10.5 10.5 -0.5 20.5];
n = size(transfer_matrix,2);

L0 = eye(n);

[L1vpara,L1vperp] = gradient_v_space_matrix(vpara,vperp,'custom');
%[L1vpara,L1vperp] = gradient_v_space_matrix(vpara,vperp,'null');


scaling_factor = 1/max(max(transfer_matrix));
%scaling_factor = 1;
transfer_matrix = transfer_matrix*scaling_factor;

switch L
    case 0
        H = L0'*L0;
        %alpha = logspace(-5,5,100)';
        alpha=[0.0033 0.01 0.033 0.1 0.33]'.^2;
        norm_operator = 1;
        
    case 1
        H = L1vperp'*L1vperp + L1vpara'*L1vpara;
        %synth
  %      alpha=([4e2 8e2 16e2 32e2 64e2 128e2]'*1).^2;
  %5 spectra
        alpha=([1e1 1e2 1e3 1e4]').^2;
  
 
        
    otherwise
        error('L must be 0 or 1 (0th or 1st order penalty function)')
end

xalpha = zeros(n,length(alpha));


for i = 1:length(alpha)
   i
   switch L
       case 0
           GalphaL0=zeros(size(L0,2),1);
           WalphaL=double([transfer_matrix; sqrt(alpha(i))*L0]);
           GalphaL=double([G; GalphaL0]);
           xalpha(:,i) = lsqnonneg(WalphaL,GalphaL);
           %xalpha(:,i) = (transfer_matrix'*transfer_matrix + alpha(i)*H)\transfer_matrix'*G;
       case 1 
           GalphaL1=zeros(2*size(L1vpara,2),1);
           GalphaL0=zeros(size(L0,2),1);
           WalphaL=double([transfer_matrix; sqrt(alpha(i))*L1vpara; sqrt(alpha(i))*L1vperp; 2e-6*sqrt(alpha(i))*L0]);
           GalphaL=double([G; GalphaL1; GalphaL0]);
           xalpha(:,i) = lsqnonneg(WalphaL,GalphaL);
   %xalpha(:,i) = (transfer_matrix'*transfer_matrix + alpha(i)*H)\transfer_matrix'*G;
   end
end  

xalpha = xalpha*scaling_factor;

