%% update tensor, rank remain the same
function [X_new,rank,G] = new_ten(X,G,rank,sX,ATTNA_para)
X = reshape(X,ATTNA_para.rX);

% calculate X by network contraction
n1 = ncon({G{1},G{2},G{3}},{[-1,1,3,-4,-5],[-2,2,-6,-7,1],[-3,-8,-9,3,2]});
n2 = ncon({G{4},G{5}},{[-1,1,-3,-4,-5],[-2,-6,-7,-8,1]});
X_new = ncon({n1,n2},{[-1,-2,-3,1,2,3,4,5,6],[-4,-5,1,3,5,2,4,6]});

%% main procedure
for iter = 1: ATTNA_para.new_ten_iter_max
    
    X_last = X_new;
% update G (5 core tensors)
    [G]=update_core_5_1(X,G,ATTNA_para.rX);

% calculate new X by network contraction
    n1 = ncon({G{1},G{2},G{3}},{[-1,1,3,-4,-5],[-2,2,-6,-7,1],[-3,-8,-9,3,2]});
    n2 = ncon({G{4},G{5}},{[-1,1,-3,-4,-5],[-2,-6,-7,-8,1]});
    X_new = ncon({n1,n2},{[-1,-2,-3,1,2,3,4,5,6],[-4,-5,1,3,5,2,4,6]});
   
   leq = X_new(:) - X_last(:);
   err = sqrt(sum(leq.^2)/sum(X_last(:).^2));
   % fprintf('iter_inner = %d, dif = %.8f\n', iter,err);
   
 
 % check convergence
   if err<ATTNA_para.new_ten_tol
        break;
   end
   
end

X_new = reshape(X_new,sX); 