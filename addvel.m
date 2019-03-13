%function to add velocity info to dataset

function Grp = addvel(Grp)

for a = 1:length(Grp)
    
    for b = 1:size(Grp(a).pos,1)
        for c = 1:size(Grp(a).pos,2)
            
            %if the velocity is already calculated, or if there is no valid
            %position data to calculate the velocity, move on.
            if isfield(Grp(a),'vel') && (size(Grp(a).vel,1) >= b && size(Grp(a).vel,2) >= c) && ~isempty(Grp(a).vel{b,c})
                continue;
            elseif isempty(Grp(a).pos{b,c})
                Grp(a).vel{b,c} = [];
                continue;
            end
            
            for d = 1:size(Grp(a).pos{b,c},2)
                for e = 1:size(Grp(a).pos{b,c},3)
                
                    Grp(a).vel{b,c}(:,d,e) = gradient(sgolayfilt(Grp(a).pos{b,c}(:,d,e),2,min([19,2*floor((size(Grp(a).pos{b,c},1)-2)/2)+1])));
                end
            end
        end
    end
end

