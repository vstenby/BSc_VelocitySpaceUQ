clear, clc, close all


s = linspace(-5,0,5)';
t = linspace(5,10,4)'; 

disp(size(s))
disp(size(t))


f = @(x,y) 2*x - 1*y;

[T,S] = meshgrid(t,s);

X = f(T,S);

disp(size(X))

%Find the derivative
hs = s(2)-s(1); ht = t(2)-t(1);

m = size(s,1);
n = size(t,1);

%If we want to pad it with zeros.
D = @(d) [spdiags([zeros(d-1,1) -ones(d-1,1) ones(d-1,1)],-1:1,d-1,d) ; zeros(1,d)];
I = @(d) speye(d);

test1 = 1/hs * D(m) * X; 

%test1 = test1(:);
test2 = kron(I(n),1/hs*D(m))*X(:);

test3 = 1/ht * X * D(n)'; %test3 = test3(:);
test4 = kron(1/ht*D(n), I(m))*X(:);


myL1vpar = kron(1/ht*D(n), I(m)); myL1vpar = full(myL1vpar);
myL1vperp = kron(I(n),1/hs*D(m)); myL1vperp = full(myL1vperp);

[L1vpar, L1vperp] = gradient_v_space_matrix(t,s,'custom');


%%
mu = zeros(20,1);
L = [myL1vpar ; myL1vperp];
delta = 0.1;





R = mvnrnd(mu,Sigma)

