% update TN with 5 core tensors
function core=update_core_5_1(X,core,I)
N = ndims(X);
for n = 1:N
    core_size{n} = size(core{n});
end

%% update core{i}
n1 = zeros([size(core{1},1),size(core{2},1),size(core{1},3:5),size(core{2},2:4)]);
n2 = zeros([size(core{4},1),size(core{5},1),size(core{4},3:5),size(core{5},2:4)]);
ccore = zeros([size(n1,5),size(n1,8),size(n2,3),size(n2,6),size(n1,1:2),size(n2,1:2)]);

for n = 1:N
    matx = reshape(permute(X,[n,n+1:N 1:n-1]),I(n),[]);
    core = circshift(core,-(n-1));
%     ccore = ncon({core{2}, core{3}, core{4}, core{5}}, ...
%         {[-5,1,2,3,-1],[-6,4,5,-2,1],[-7,6,-3,2,4],[-8,-4,3,5,6]}, [1,2,3,4,5,6]);
    n1 = ncon({core{2},core{3}},{[-1,1,-3,-4,-5],[-2,-6,-7,-8,1]});
    n2 = ncon({core{4},core{5}},{[-1,1,-3,-4,-5],[-2,-6,-7,-8,1]});
    ccore = ncon({n1,n2},{[-5,-6,1,2,-1,3,4,-2],[-7,-8,-3,1,3,-4,2,4]});

    mat_ccore = reshape(ccore,[],prod(I)/I(n));
    
    right = mat_ccore'*pinv(mat_ccore*mat_ccore');
    mat_core_new = matx * right;
    core = circshift(core,(n-1));
    core{n} = reshape(mat_core_new,core_size{n});
end
end
