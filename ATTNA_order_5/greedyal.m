function [core,D,rse]=greedyal(core,I,D,X)
D_old = D;
% the network has 10 edges
for i = 1:10
    D(i) = 1;
    [core_new]=change_r(core,I,D);
    rse(i) = calculate_rse(X,core_new,I,5);
    D(i) = D_old(i);
end