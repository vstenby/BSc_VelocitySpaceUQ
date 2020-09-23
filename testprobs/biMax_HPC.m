%% Program to run on HPC computer.
clear, clc, close all

idx = 1;
while 1
    run('biMax_UQ')
    dataname = strcat('test',num2str(idx),'.mat');
    save(dataname,'lamsamp','delsamp','xsamp')
    idx = idx + 1;
end