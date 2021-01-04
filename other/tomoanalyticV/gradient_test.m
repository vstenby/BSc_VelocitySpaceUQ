clear, clc, close all

%Parameters for the (vpara,vperp)-grid.
vparamin=-4e6;
vparamax=4e6;
vperpmin=1e4;
vperpmax=4e6;

vparadim = 40; vperpdim = 20;

v_par = linspace(vparamin, vparamax, vparadim);
v_perp = linspace(vperpmin, vperpmax, vperpdim);

v_par_dim = length(v_par);
v_perp_dim = length(v_perp);

dv_par = v_par(2)-v_par(1);
dv_perp = v_perp(2)-v_perp(1);

Delta_v_par = diag(ones(1,v_perp_dim*v_par_dim))...
    + diag(-1*ones(1,(v_par_dim-1)*v_perp_dim),-v_perp_dim);

Delta_v_par(1:v_perp_dim,:) = 0;

Delta_v_par = Delta_v_par/dv_par;

%%
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
