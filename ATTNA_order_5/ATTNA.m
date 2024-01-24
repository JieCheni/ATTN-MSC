function [S_tensor,rank,G] = ATTNA(X,G,rank,iter,sX,ATTNA_para)
    if (iter == ATTNA_para.prue_gap)
        % prune the rank of tensor and increase the rank
        [S_tensor,rank,G] = prune_ten(X,G,rank,sX,ATTNA_para);
    else
        % update tensor
        [S_tensor,rank,G] = new_ten(X,G,rank,sX,ATTNA_para);
    end