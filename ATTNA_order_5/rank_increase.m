function [G,rank,RSE,total_size,X_new]=rank_increase(X,G,rank,ATTNA_para)

n1 = ncon({G{1},G{2},G{3}},{[-1,1,3,-4,-5],[-2,2,-6,-7,1],[-3,-8,-9,3,2]});
n2 = ncon({G{4},G{5}},{[-1,1,-3,-4,-5],[-2,-6,-7,-8,1]});
X_new = ncon({n1,n2},{[-1,-2,-3,1,2,3,4,5,6],[-4,-5,1,3,5,2,4,6]});


tic;
%% main procedure
for iter = 1: ATTNA_para.ri_iter_max
    
    X_last = X_new;
% update G
    [G]=update_core_5_1(X,G,ATTNA_para.rX);

% calculate new X
    n1 = ncon({G{1},G{2},G{3}},{[-1,1,3,-4,-5],[-2,2,-6,-7,1],[-3,-8,-9,3,2]});
    n2 = ncon({G{4},G{5}},{[-1,1,-3,-4,-5],[-2,-6,-7,-8,1]});
    X_new = ncon({n1,n2},{[-1,-2,-3,1,2,3,4,5,6],[-4,-5,1,3,5,2,4,6]});
   
   leq = X_new(:) - X_last(:);
   err = sqrt(sum(leq.^2)/sum(X_last(:).^2));
   %fprintf('iter_inner = %d, dif = %.8f\n', iter,err);
   
 % rank increase  
    leq1=X(:)-X_new(:);
    RSE=sqrt(sum(leq1.^2)/sum(X(:).^2));
    %fprintf('iter_inner = %d, dif = %.8f\n', iter,RSE);
    
    if mod(iter,ATTNA_para.ri_update_iter) == 0  &&  RSE > ATTNA_para.ri_RSE_threshold
        [G,rank]=greedyal_add(X,G,rank,ATTNA_para); 
        n1 = ncon({G{1},G{2},G{3}},{[-1,1,3,-4,-5],[-2,2,-6,-7,1],[-3,-8,-9,3,2]});
        n2 = ncon({G{4},G{5}},{[-1,1,-3,-4,-5],[-2,-6,-7,-8,1]});
        X_new = ncon({n1,n2},{[-1,-2,-3,1,2,3,4,5,6],[-4,-5,1,3,5,2,4,6]});
    elseif RSE <= ATTNA_para.ri_RSE_threshold
        break;    
    end
 
 % check convergence
%    if err<tol
%         break;
%    end
   
end

leq1=X(:)-X_new(:);
RSE=sqrt(sum(leq1.^2)/sum(X(:).^2));

tr_size=storage_size(G);

total_size=tr_size;