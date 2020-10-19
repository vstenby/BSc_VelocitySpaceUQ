%Simulation test.

function test(n)
    x = randn(n);
    save(strcat(num2str(n),'x.mat'),'x')
end
