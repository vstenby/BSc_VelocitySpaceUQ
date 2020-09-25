function n = nextsimnum(folder)
%Finds the next sim number in the folder.
sims = dir(strcat('./',folder,'/sim*.mat'));
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