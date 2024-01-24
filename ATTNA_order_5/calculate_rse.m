function [rse,core] = calculate_rse(X,core,I,i)

for iter = 1: i
% update core
    [core]=update_core_5_1(X,core,I);

% calculate new X
    n1 = ncon({core{1},core{2},core{3}},{[-1,1,3,-4,-5],[-2,2,-6,-7,1],[-3,-8,-9,3,2]});
    n2 = ncon({core{4},core{5}},{[-1,1,-3,-4,-5],[-2,-6,-7,-8,1]});
    X_new = ncon({n1,n2},{[-1,-2,-3,1,2,3,4,5,6],[-4,-5,1,3,5,2,4,6]});   
end

leq1=X(:)-X_new(:);
rse=sqrt(sum(leq1.^2)/sum(X(:).^2));
