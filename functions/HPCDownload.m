function HPCDownload(folder,destination,username)
% Downloads from DTU's HPC clusters using scp.
%
%

cmd = sprintf('scp -r %s@transfer.gbar.dtu.dk:~/Desktop/DTU/BSc/BSc_VelocitySpaceUQ/simulation/%s %s',username,folder,destination);
system(cmd);
end