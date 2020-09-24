%% Program to run on HPC computer.
%clear, clc, close all

while 1
    run('biMax_UQ') 
    n = find_newest_name();
    dataname = strcat('sim',num2str(n),'.mat');
    save(dataname,'A','S','x_sim','del_sim','lam_sim','alph_sim')
end

function n = find_newest_name()
%Finds the next sim number in the folder.
sims = dir('sim*.mat');
sims = {sims.name};
n = 0;
for i=1:length(sims)
    temp_str = sims{i};
    temp_str = strrep(temp_str,'sim','');
    temp_str = strrep(temp_str,'.mat','');
    num = str2double(temp_str);
    if num > n
        n = num;
    end
end
n = n+1;
end