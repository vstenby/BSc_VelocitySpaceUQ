function [Delta_v_par,Delta_v_perp] = gradient_v_space_matrix(v_par,v_perp,scheme)

v_par_dim = length(v_par);
v_perp_dim = length(v_perp);

dv_par = v_par(2)-v_par(1);
dv_perp = v_perp(2)-v_perp(1);

if strcmpi(scheme,'central')
    
    Delta_v_par = diag([-2*ones(1,v_perp_dim) zeros(1,(v_par_dim-2)*v_perp_dim) 2*ones(1,v_perp_dim)]);
    Delta_v_par = Delta_v_par + diag([2*ones(1,v_perp_dim) ones(1,(v_par_dim-2)*v_perp_dim)],v_perp_dim) - diag([ones(1,(v_par_dim-2)*v_perp_dim) 2*ones(1,v_perp_dim)],-v_perp_dim);
    
    Delta_v_par = Delta_v_par/(2*dv_par);
    
        
    Dv_perp_c0 = zeros(v_par_dim*v_perp_dim,1);
    Dv_perp_c1 = zeros(v_par_dim*v_perp_dim,1);
    Dv_perp_c2 = zeros(v_par_dim*v_perp_dim,1);
    
    for i = 1:v_par_dim
        Dv_perp_c0((i-1)*v_perp_dim+1:i*v_perp_dim) = [-2 zeros(1,v_perp_dim-2) 2];
        Dv_perp_c1((i-1)*v_perp_dim+1:i*v_perp_dim) = [2 ones(1,v_perp_dim-2) 0];
        Dv_perp_c2((i-1)*v_perp_dim+1:i*v_perp_dim) = [-ones(1,v_perp_dim-2) -2 0];
    end
    
    Dv_perp_c1 = Dv_perp_c1(1:end-1);
    Dv_perp_c2 = Dv_perp_c2(1:end-1);
    
    Delta_v_perp = diag(Dv_perp_c0) + diag(Dv_perp_c1,1) + diag(Dv_perp_c2,-1);
    
    Delta_v_perp = Delta_v_perp/(2*dv_perp);
    
elseif strcmpi(scheme,'backward')
    
    Delta_v_par = diag([-1*ones(1,v_perp_dim) ones(1,(v_par_dim-1)*v_perp_dim)])...
        + diag([ones(1,v_perp_dim) zeros(1,(v_par_dim-2)*v_perp_dim)],v_perp_dim)...
        + diag(-1*ones(1,(v_par_dim-1)*v_perp_dim),-v_perp_dim);
    
    Delta_v_par = Delta_v_par/dv_par;
    
       
    Dv_perp_c0 = zeros(v_par_dim*v_perp_dim,1);
    Dv_perp_c1 = zeros(v_par_dim*v_perp_dim,1);
    Dv_perp_c2 = zeros(v_par_dim*v_perp_dim,1);
    
    for i = 1:v_par_dim
        Dv_perp_c0((i-1)*v_perp_dim+1:i*v_perp_dim) = [-1 ones(1,v_perp_dim-1)];
        Dv_perp_c1((i-1)*v_perp_dim+1:i*v_perp_dim) = [1 zeros(1,v_perp_dim-1)];
        Dv_perp_c2((i-1)*v_perp_dim+1:i*v_perp_dim) = [-ones(1,v_perp_dim-1) 0];
    end
    
    Dv_perp_c1 = Dv_perp_c1(1:end-1);
    Dv_perp_c2 = Dv_perp_c2(1:end-1);
    
    Delta_v_perp = diag(Dv_perp_c0) + diag(Dv_perp_c1,1) + diag(Dv_perp_c2,-1);
    
    Delta_v_perp = Delta_v_perp/dv_perp;
    
elseif strcmpi(scheme,'custom')
    
    Delta_v_par = diag(ones(1,v_perp_dim*v_par_dim))...
        + diag(-1*ones(1,(v_par_dim-1)*v_perp_dim),-v_perp_dim);
        
    Delta_v_par(1:v_perp_dim,:) = 0;
    
    Delta_v_par = Delta_v_par/dv_par;
    
        
    Dv_perp_c0 = zeros(v_par_dim*v_perp_dim,1);
    Dv_perp_c1 = zeros(v_par_dim*v_perp_dim,1);
        
    for i = 1:v_par_dim
        Dv_perp_c0((i-1)*v_perp_dim+1:i*v_perp_dim) = [-1*ones(1,v_perp_dim-1) 0];
        Dv_perp_c1((i-1)*v_perp_dim+1:i*v_perp_dim) = [ones(1,v_perp_dim-1) 0];
        %Dp_c0((i-1)*p_dim+1:i*p_dim) = [-1*ones(1,p_dim-1) -1];
        %Dp_c1((i-1)*p_dim+1:i*p_dim) = [ones(1,p_dim-1) 1];
    end
    
    Dv_perp_c1 = Dv_perp_c1(1:end-1);
    
    Delta_v_perp = diag(Dv_perp_c0) + diag(Dv_perp_c1,1);
    
    Delta_v_perp = Delta_v_perp/dv_perp;
    
    
    elseif strcmpi(scheme,'null')
    
    nulls2D=ones(v_perp_dim,v_par_dim);

for j=1:length(v_perp)
    if v_perp(j)>13.5e6
    nulls2D(j,:)=1+(v_perp(j)-13.5e6)/13.5e6*100;
    end
end
if 1
figure;  imagesc(v_par,v_perp,nulls2D);  set(gca,'ydir','normal');
colorbar;
colormap(flipud(hot));
set(gca,'Fontsize',20)
% xlabel('E [keV]')
% ylabel('p [-]')
% saveas(gcf,['C:\Users\msal\Documents\MATLAB\tomo\results\32323\lambda2D.png'])
% saveas(gcf,['C:\Users\msal\Documents\MATLAB\tomo\results\32323\lambda2D.eps'],'epsc')
end
nulls=reshape(nulls2D,v_par_dim*v_perp_dim,1);    
        
    Delta_v_par = diag(ones(1,v_perp_dim*v_par_dim))...
        + diag(-1*ones(1,(v_par_dim-1)*v_perp_dim),-v_perp_dim);
        
    Delta_v_par(1:v_perp_dim,:) = 0;
    
    Delta_v_par = Delta_v_par/dv_par;
    for i=1:length(nulls)
        Delta_v__par(:,i)=Delta_v_par(:,i)*nulls(i);
    end
        
    Dv_perp_c0 = zeros(v_par_dim*v_perp_dim,1);
    Dv_perp_c1 = zeros(v_par_dim*v_perp_dim,1);
        
    for i = 1:v_par_dim
        Dv_perp_c0((i-1)*v_perp_dim+1:i*v_perp_dim) = [-1*ones(1,v_perp_dim-1) 0];
        Dv_perp_c1((i-1)*v_perp_dim+1:i*v_perp_dim) = [ones(1,v_perp_dim-1) 0];
        %Dp_c0((i-1)*p_dim+1:i*p_dim) = [-1*ones(1,p_dim-1) -1];
        %Dp_c1((i-1)*p_dim+1:i*p_dim) = [ones(1,p_dim-1) 1];
    end
    
    Dv_perp_c1 = Dv_perp_c1(1:end-1);
    
    Delta_v_perp = diag(Dv_perp_c0) + diag(Dv_perp_c1,1);
    
    Delta_v_perp = Delta_v_perp/dv_perp;
    for i=1:length(nulls)
        Delta_v_perp(i,:)=Delta_v_perp(i,:)*nulls(i);
    end
end