% ri_RSE_threshold_list = [0.5];
lambda_list = 0.5;
ri_RSE_threshold_list = [0.5];
% lambda_list = [0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1,10,100,1000];
for iii = 1:length(ri_RSE_threshold_list)
    for jjj = 1:length(lambda_list)
%% For convinience, we assume the order of the tensor is always 3;
clearvars -except iii jjj ri_RSE_threshold_list lambda_list
close all; 
% clc;
addpath('datasets','ATTNA_order_5','ConfusionMatrices');
addpath('ClusteringMeasure', 'LRR', 'Nuclear_norm_l21_Algorithm');
rand('seed',30);

%% Load Yale，Note: each column is an sample
load('yale.mat');
cls_num = length(unique(gt));
X{1} = X2; X{2} =X1; X{3} = X3;
clear X1 X2 X3
V = length(X); % views number

for v=1:V
    [X{v}]=NormalizeData(X{v});
     %X{v} = zscore(X{v},1);
end

N = size(X{1},2); % total samples number

for v=1:V
    Z{v} = zeros(N,N); %Z{2} = zeros(N,N);
    W{v} = zeros(N,N);
    S{v} = zeros(N,N);
    E{v} = zeros(size(X{v},1),N); %E{2} = zeros(size(X{v},1),N);
    Y{v} = zeros(size(X{v},1),N); %Y{2} = zeros(size(X{v},1),N);
end

sX = [N, N, V];
%% Parameter
% usually turn ir,lambda

% reshape X to 5th-order tensor
rX = [N/cls_num,cls_num,N/cls_num,cls_num,V];
% rx = [];

        ri_RSE_threshold = ri_RSE_threshold_list(iii);
        lambda = lambda_list(jjj);
        fprintf("ri_RSE_threshold = %f, lambda = %f\n",ri_RSE_threshold,lambda);



% initial rank,usually choose 2 or 3
ir = 2;
rank=ir*ones(1,10);
G{1} = rand(rX(1),rank(1),rank(2),rank(3),rank(4));
G{2} = rand(rX(2),rank(5),rank(6),rank(7),rank(1));
G{3} = rand(rX(3),rank(8),rank(9),rank(2),rank(5));
G{4} = rand(rX(4),rank(10),rank(3),rank(6),rank(8));
G{5} = rand(rX(5),rank(4),rank(7),rank(9),rank(10));

%%
Isconverg = 0;
epson = 1e-7;

% trade off between 2,1-norm of E and F-norm of Z minus phi
% lambda = 0.5; %0.5 best

% set Max_iter of
Max_iter = 50;
iter = 0;

mu = 10e-5; max_mu = 10e10; pho_mu = 2;
rho = 10e-5; max_rho = 10e12; pho_rho = 2;

%% ATTNA parameter
ir = 2;
rank=ir*ones(1,10);
G{1} = rand(rX(1),rank(1),rank(2),rank(3),rank(4));
G{2} = rand(rX(2),rank(5),rank(6),rank(7),rank(1));
G{3} = rand(rX(3),rank(8),rank(9),rank(2),rank(5));
G{4} = rand(rX(4),rank(10),rank(3),rank(6),rank(8));
G{5} = rand(rX(5),rank(4),rank(7),rank(9),rank(10));

%% ATTNA parameter
ATTNA_para.ir = 2;
% reshape X to 5th-order tensor
ATTNA_para.rX = [N/cls_num,cls_num,N/cls_num,cls_num,V];

% new_ten.m parameter
ATTNA_para.new_ten_iter_max = 10;
ATTNA_para.new_ten_tol=1e-3;

% prue_ten.m parameter
ATTNA_para.prune_theshold = 0.2;
ATTNA_para.prue_gap = 4;
ATTNA_para.prue_ten_iter_max = 50;
ATTNA_para.prue_ten_tol=1e-3;

%rank_increase parameter
ATTNA_para.ri_RSE_threshold = ri_RSE_threshold;%epsilon,0.5
ATTNA_para.ri_step = 1;
ATTNA_para.ri_update_iter = 5;
ATTNA_para.ri_iter_max = 50;
ATTNA_para.ri_tol=1e-3;

tic;
t1 = cputime;
while(Isconverg == 0)
    % fprintf('----processing iter %d--------\n', iter+1);
    for v=1:V
        %1 update Z^v
        
        tmp = (X{v}'*Y{v} + mu*X{v}'*X{v} - mu*X{v}'*E{v} - W{v})./rho +  S{v};
        Z{v}=inv(eye(N,N)+ (mu/rho)*X{v}'*X{v})*tmp;
        
        %2 update E^v
        D = [];
        for ii = 1:V
            D = [D;X{ii}-X{ii}*Z{ii}+Y{ii}/mu];
        end
        [Econcat] = solve_l1l2(D,lambda/mu);       
        for ii = 1:v % reshape Econcat into tensor
            idx1=1;idx2=0;
            for jj = 1:ii-1
                idx1 = idx1+size(X{jj},1);
            end
            for jj = 1:ii
                idx2 = idx2+size(X{jj},1);
            end
            E{ii} = Econcat(idx1:idx2,:);
        end

        %3 update Yk
        Y{v} = Y{v} + mu*(X{v}-X{v}*Z{v}-E{v});
    end
    
    %4 update S by ATTNA
    Z_tensor = cat(3, Z{:,:});
    W_tensor = cat(3, W{:,:});
    z = Z_tensor(:);
    w = W_tensor(:);
    
    [S_tensor,rank,G] = ATTNA(Z_tensor + 1/rho*W_tensor,G,rank,iter,sX,ATTNA_para);
    s = S_tensor(:);

    %5 update W
    w = w + rho*(z - s);

    %% coverge condition
    Isconverg = 1;
    for v=1:V
        T(v)=norm(X{v}-X{v}*Z{v}-E{v},inf);
        if (norm(X{v}-X{v}*Z{v}-E{v},inf)>epson)
            history.norm_Z = norm(X{v}-X{v}*Z{v}-E{v},inf);
            % fprintf('    norm_Z %7.10f    ', history.norm_Z);
            Isconverg = 0;
        end
        
        S{v} = S_tensor(:,:,v);
        W_tensor = reshape(w, sX);
        W{v} = W_tensor(:,:,v);
        
        Ti(v)=norm(Z{v}-S{v},inf);
        if (norm(Z{v}-S{v},inf)>epson)
            history.norm_Z_G = norm(Z{v}-S{v},inf);
            % fprintf('norm_Z_G %7.10f    \n', history.norm_Z_G);
            Isconverg = 0;
        end
    end
    
   if iter>0
        Tm(iter)=max(T);
        Tim(iter)=max(Ti); 
   end
    
    if (iter>Max_iter)
        Isconverg  = 1;
    end
    iter = iter + 1;
    mu = min(mu*pho_mu, max_mu);
    rho = min(rho*pho_rho, max_rho);
end
fprintf('----processing iter %d--------\n', iter+1);
t2 = cputime;
total_time = t2-t1;

M = 0;
for v=1:V
    M = M + abs(Z{v})+abs(Z{v}');
end
figure(1); imagesc(M);colorbar;
% S_bar = CLR(M, cls_num, 0, 0 );
% figure(2); imagesc(S_bar);



g=1:1:(iter-1);
figure();
plot(g,Tm,'r',g,Tim,'b','LineWidth',2);
legend('Reconstruction','Match');
xlabel('Iteration');
ylabel('Error');

% C = SpectralClustering(M,cls_num);% 
% [A nmi avgent] = compute_nmi(gt,C);
% %C = SpectralClustering(abs(Z{1})+abs(Z{1}'),cls_num);
% %[A nmi avgent] = compute_nmi(gt,C)
% % C = SpectralClustering(abs(Z{2})+abs(Z{2}'),cls_num);
% % [A nmi avgent] = compute_nmi(gt,C)
% % C = SpectralClustering(abs(Z{3})+abs(Z{3}'),cls_num);
% % [A nmi avgent] = compute_nmi(gt,C)
% ACC = Accuracy(C,double(gt));
% [f,p,r] = compute_f(gt,C);
% [AR,RI,MI,HI]=RandIndex(gt,C);
% toc;
% %save('my_new_COIL20MV_res.mat','M','ACC','nmi','AR','f','p','r');
 for i=1:10
        C = SpectralClustering(M,cls_num);% C = kmeans(U,numClust,'EmptyAction','drop');
        [Fi(i),Pi(i),Ri(i)] = compute_f(gt,C);
        [A nmii(i) avgenti(i)] = compute_nmi(gt,C);    
        ACCi(i) = Accuracy(C,double(gt));
        if (min(gt)==0)
            [ARi(i),RIi(i),MIi(i),HIi(i)]=RandIndex(gt+1,C);
        else
            [ARi(i),RIi(i),MIi(i),HIi(i)]=RandIndex(gt,C);
        end       
 end
    
%  figure();
% predict_label = Accuracy_map(C,double(gt));
% [confusion_matrix]=compute_confusion_matrix(predict_label,yalenum_in_class,yalename_class); 

    F(1) = mean(Fi); F(2) = std(Fi);
    P(1) = mean(Pi); P(2) = std(Pi);
    R(1) = mean(Ri); R(2) = std(Ri);
    nmi(1) = mean(nmii); nmi(2) = std(nmii);
    avgent(1) = mean(avgenti); avgent(2) = std(avgenti);
    AR(1) = mean(ARi); AR(2) = std(ARi);
    ACC(1)=mean(ACCi); ACC(2)=std(ACCi);
    Time = toc;
    fprintf('F: %.3f(%.3f)\n', F(1), std(Fi));
    fprintf('P: %.3f(%.3f)\n', P(1), std(Pi));    
    fprintf('R: %.3f(%.3f)\n', R(1), std(Ri));
    fprintf('nmi:%.3f(%.3f)\n', nmi(1), std(nmii));
%     fprintf('avgent: %f(%f)\n', avgent(1), std(avgenti));
    fprintf('AR: %.3f(%.3f)\n', AR(1), std(ARi));
    fprintf('ACC: %.3f(%.3f)\n',  ACC(1), std(ACCi));
    dlmwrite('ATTN_yale.txt',[ri_RSE_threshold lambda F(1) F(2) P(1) P(2) R(1) R(2) nmi(1), nmi(2) AR(1), AR(2) ACC(1) ACC(2) Time],'-append','delimiter','\t','newline','pc');

    end
end