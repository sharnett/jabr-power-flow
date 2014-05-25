casefile = 'case5_renumber_tree.m';
opt = mpoption('OUT_ALL', 0, 'VERBOSE', 0);
%%
% computes matpower solution and prints out a python-style list of 
% (voltage magnitude, voltage angle) tuples. angles are in degrees
mpc = loadcase(casefile);
r = runpf(mpc, opt);
n = size(r.bus, 1);

stringify = @(i) sprintf('(%f, %f)', r.bus(i, 8), r.bus(i, 9));
fprintf('[');
for i=1:n-1, fprintf('%s, ', stringify(i)); end;
fprintf('%s]\n', stringify(n));