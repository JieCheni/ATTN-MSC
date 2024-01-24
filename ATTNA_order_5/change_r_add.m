function [core]=change_r_add(core,I,rank_new,rank)
%  rand('seed',100);
%% ¸üĞÂ core
    core_new{1} = rand(I(1),rank_new(1),rank_new(2),rank_new(3),rank_new(4));
    core_new{1}(1:I(1),1:rank(1),1:rank(2),1:rank(3),1:rank(4)) = core{1};   
    core{1}=core_new{1};
    
    core_new{2} = rand(I(2),rank_new(5),rank_new(6),rank_new(7),rank_new(1));
    core_new{2}(1:I(2),1:rank(5),1:rank(6),1:rank(7),1:rank(1)) = core{2};   
    core{2}=core_new{2};
    
    core_new{3} = rand(I(3),rank_new(8),rank_new(9),rank_new(2),rank_new(5));
    core_new{3}(1:I(3),1:rank(8),1:rank(9),1:rank(2),1:rank(5)) = core{3};  
    core{3}=core_new{3};
    
    core_new{4} = rand(I(4),rank_new(10),rank_new(3),rank_new(6),rank_new(8));
    core_new{4}(1:I(4),1:rank(10),1:rank(3),1:rank(6),1:rank(8)) = core{4};    
    core{4}=core_new{4};
    
    core_new{5} = rand(I(5),rank_new(4),rank_new(7),rank_new(9),rank_new(10));
    core_new{5}(1:I(5),1:rank(4),1:rank(7),1:rank(9),1:rank(10)) = core{5};    
    core{5}=core_new{5};
end