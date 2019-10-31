%function to add velocity info to dataset at subject level

function Block = addvelsubj(Block)

for b = 1:size(Block.pos,1)
    for c = 1:size(Block.pos,2)
            
        %if the velocity is already calculated, or if there is no valid
        %position data to calculate the velocity, move on.
        if isfield(Block,'vel') && (size(Block.vel,1) >= b && size(Block.vel,2) >= c) && ~isempty(Block.vel{b,c})
            continue;
        elseif isempty(Block.pos{b,c})
            Block.vel{b,c} = [];
            continue;
        end
        
        dt = mean(diff(Block.time{b,c}));
        
        
        for d = 1:size(Block.pos{b,c},2)
            for e = 1:size(Block.pos{b,c},3)
                
                Block.vel{b,c}(:,d,e) = gradient(sgolayfilt(Block.pos{b,c}(:,d,e),2,min([19,2*floor((size(Block.pos{b,c},1)-2)/2)+1])),dt);
            end
        end
    end
end


