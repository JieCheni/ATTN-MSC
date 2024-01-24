%% prune the rank of tensor
function [X_new,rank,G] = prune_ten(X,G,rank,sX,ATTNA_para)
X = reshape(X,ATTNA_para.rX);

%% Initialization

n1 = ncon({G{1},G{2},G{3}},{[-1,1,3,-4,-5],[-2,2,-6,-7,1],[-3,-8,-9,3,2]});
n2 = ncon({G{4},G{5}},{[-1,1,-3,-4,-5],[-2,-6,-7,-8,1]});
X_new = ncon({n1,n2},{[-1,-2,-3,1,2,3,4,5,6],[-4,-5,1,3,5,2,4,6]});


tic;
%% main procedure
for iter = 1: ATTNA_para.prue_ten_iter_max
    
    X_last = X_new;
% update G
    [G]=update_core_5_1(X,G,ATTNA_para.rX);

% calculate new X
    n1 = ncon({G{1},G{2},G{3}},{[-1,1,3,-4,-5],[-2,2,-6,-7,1],[-3,-8,-9,3,2]});
    n2 = ncon({G{4},G{5}},{[-1,1,-3,-4,-5],[-2,-6,-7,-8,1]});
    X_new = ncon({n1,n2},{[-1,-2,-3,1,2,3,4,5,6],[-4,-5,1,3,5,2,4,6]});

   
   leq = X_new(:) - X_last(:);
   err = sqrt(sum(leq.^2)/sum(X_last(:).^2));
   % fprintf('iter = %d, dif = %.8f\n', iter,err);
   if err<ATTNA_para.prue_ten_tol
        break;
   end
   
end

leq1=X(:)-X_new(:);
RSE=sqrt(sum(leq1.^2)/sum(X(:).^2));
% fprintf('iter_inner = %d, rse = %.8f\n', iter,RSE);
ground = RSE;

% calculate each rse when a edge is cut 
[G,rank,rse]=greedyal(G,ATTNA_para.rX,rank,X);

% the network has 10 edges，calculate the change of rse
for i = 1:10
    change(i) = rse(i)-ground;
end
[~,index2] = sort(change,'ascend');
for i = 1:9   
    ratio_change(i+1) = (change(index2(i+1)) - change(index2(i)))/change(index2(1));
end
ratio_change(1) = 0;
%% solution 1
% [~,a]=max(ratio_change(1:7));
% for i = 1:10
%     if i <= a-1
%         rank(index2(i)) = 1;
%     end
% end

%% solution 2
% for i = 1:10      
%       if ratio_change(i)< ratio_change (i+1) && ratio_change (i+1) > ratio_change(i+2)
%          for jj = 1: i
%              rank(index2(jj)) = 1;
%          end
%          break;
%       end
% end
% [~,index1] = sort(rse,'ascend');
% for i = 1:9
%     if rse(index(i+1)) ~= 1
%         ratio_rse(i+1) = (rse(index1(i+1)) - rse(index1(i)))/rse(index1(i));
%     end
% end
% ratio_rse(1) = 0;

%% solution 3

% rate = 0.4;
% len = length(index);
% for n = 1:(len-rate*len)
%     rank(index(n)) = 2;
% end
% for n = (len-rate*len)+1:len
%     rank(index(n)) = 1;
% end

%% solution 4
if min(change)<ATTNA_para.prune_theshold
    [~,idx] = min(change);
    rank(idx) = 1;
end
%%
[G]=change_r(G,ATTNA_para.rX,rank);

[G,rank,RSE,total_size,X_new]=rank_increase(X,G,rank,ATTNA_para);



%fprintf('RSE = %.8f, total_size = %d\n',RSE,total_size);
X_new = reshape(X_new,sX);

%% greed search， cut one edge at one time
function [core,rank,rse]=greedyal(core,I,rank,X)
D_old = rank;
% the network has 10 edges
for i = 1:10
    rank(i) = 1;
    [core_new]=change_r(core,I,rank);
    rse(i) = calculate_rse(X,core_new,I,5);
    rank(i) = D_old(i);
end


%% update core
function [core]=change_r(core,I,RA)
    core_new{1} = core{1}(1:I(1),1:RA(1),1:RA(2),1:RA(3),1:RA(4));   
    core{1}=core_new{1};
    
    core_new{2} = core{2}(1:I(2),1:RA(5),1:RA(6),1:RA(7),1:RA(1));   
    core{2}=core_new{2};
    
    core_new{3} = core{3}(1:I(3),1:RA(8),1:RA(9),1:RA(2),1:RA(5));  
    core{3}=core_new{3};
    
    core_new{4} = core{4}(1:I(4),1:RA(10),1:RA(3),1:RA(6),1:RA(8));    
    core{4}=core_new{4};
    
    core_new{5} = core{5}(1:I(5),1:RA(4),1:RA(7),1:RA(9),1:RA(10));    
    core{5}=core_new{5};
