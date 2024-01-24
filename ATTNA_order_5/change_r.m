function [core]=change_r(core,I,RA)

%% update core
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
end