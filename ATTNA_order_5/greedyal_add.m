%% greedy algorithm, add rank
function [G_last,rank_last]=greedyal_add(X,G,rank,ATTNA_para)
calculate_rse_iter = 3;% 
rank_new = rank;
    % check whether the edge is cut
    if rank(1) > 1
        rank_new(1) = rank(1) + ATTNA_para.ri_step;
        [G_new]=change_r_add(G,ATTNA_para.rX,rank_new,rank);
        [rse_old,G_last] = calculate_rse(X,G_new,ATTNA_para.rX,calculate_rse_iter);
        rank_last = rank_new;
        rank_new(1) = rank(1);
    else
        rse_old = 1; 
        G_last = G;
    end


for i = 2:10
    if rank(i) > 1
        rank_new(i) = rank(i) + ATTNA_para.ri_step;
        [G_new]=change_r_add(G,ATTNA_para.rX,rank_new,rank);
        [rse_new,G_nn] = calculate_rse(X,G_new,ATTNA_para.rX,calculate_rse_iter);
        if rse_new < rse_old
            rse_old = rse_new;
            G_last = G_nn;
            rank_last = rank_new;
        end
        rank_new(i) = rank(i);
    end
end
