%merge two processed data files together

fprintf('\nSelect primary data file.\n');
[fname1,pname1] = uigetfile('*.mat','Select primary data file.');

fprintf('\nSelect data file to merge in.\n');
[fname2,pname2] = uigetfile('*.mat','Select data file to merge in.');

load(fullfile(pname1,fname1));

tmp = load(fullfile(pname2,fname2));

%%

for blk = 1:length(BlockData)
    if isempty(tmp.BlockData{blk})
        continue;
    elseif isempty(BlockData{blk}) && ~isempty(tmp.BlockData{blk})
        BlockData{blk} = tmp.BlockData{blk};
    else %~isempty(BlockData{blk}) && ~isempty(tmp.BlockData{blk})
        
        %for now we will just merge the Grp field, assuming that all other
        %fields are the same (BlockID, BlockName, Items)
        
        for a = 1:3  %there are at most 3 groups!
            if length(tmp.BlockData{blk}.Grp) < a || isempty(tmp.BlockData{blk}.Grp(a)) || isempty(tmp.BlockData{blk}.Grp(a).Group)
                %this group doesn't exist the to-be-merged data set, so skip it
                continue;
            elseif isempty(BlockData{blk}.Grp(a)) || isempty(BlockData{blk}.Grp(a).Group)
                %this is the easy case, we there is only data in the
                %to-be-merged data set to we can just copy all the fields
                %in
                f = fieldnames(tmp.BlockData{blk}.Grp(a));
                for b = 1:length(f)
                    BlockData{blk}.Grp(a).(f{b}) = tmp.BlockData{blk}.Grp(a).(f{b});
                end
            else
                %this is the hard case where there is data in each block;
                %we need to see what overlaps!
                sub_orig = BlockData{blk}.Grp(a).Subjects;
                sub_new = tmp.BlockData{blk}.Grp(a).Subjects;
                f = fieldnames(tmp.BlockData{blk}.Grp(a));
                
                for b = 1:length(sub_new)
                    sub_ind = strcmpi(sub_orig,sub_new{b});
                    if any(sub_ind)
                        %both data sets have the same subject in them
                        todo = input(sprintf('Duplicate entries exist for S%s. (r)eplace, (d)uplicate, (s)kip: ',sub_new{b}),'s');
                        if strcmpi(todo,'s')
                            continue;
                        elseif strcmpi(todo,'r')  %replace entry
                            sub_ind = find(sub_ind == 1);
                            if length(sub_ind) > 1
                                sub_rep = input(sprintf('%d duplicate entries exist. Which entry do you want to replace? ',length(sub_ind)));
                                if ~isempty(sub_rep)
                                    sub_ind = sub_ind(sub_rep(1));
                                else
                                    sub_ind = sub_ind(1);
                                end
                            end
                            for c = 1:length(f)
                                if strcmpi(f{c},'Group') || isempty(tmp.BlockData{blk}.Grp(a).(f{c}))
                                    continue;
                                elseif any(size(BlockData{blk}.Grp(a).(f{c})) == 1)
                                    BlockData{blk}.Grp(a).(f{c})(sub_ind) = tmp.BlockData{blk}.Grp(a).(f{c})(b);
                                else
                                    BlockData{blk}.Grp(a).(f{c})(sub_ind,:) = tmp.BlockData{blk}.Grp(a).(f{c})(b,:);
                                end
                            end
                        else %strcmpi(todo,'d')  %duplicate entry
                            for c = 1:length(f)
                                if strcmpi(f{c},'Group') || isempty(tmp.BlockData{blk}.Grp(a).(f{c}))
                                    continue;
                                elseif any(size(BlockData{blk}.Grp(a).(f{c})) == 1)
                                    BlockData{blk}.Grp(a).(f{c})(end+1) = tmp.BlockData{blk}.Grp(a).(f{c})(b);
                                else
                                    BlockData{blk}.Grp(a).(f{c})(end+1,:) = tmp.BlockData{blk}.Grp(a).(f{c})(b,:);
                                end
                            end
                        end
                        
                    else
                        %this is a new subject, we can just append
                        for c = 1:length(f)
                            if strcmpi(f{c},'Group') || isempty(tmp.BlockData{blk}.Grp(a).(f{c}))
                                continue;
                            elseif any(size(BlockData{blk}.Grp(a).(f{c})) == 1)
                                BlockData{blk}.Grp(a).(f{c})(end+1) = tmp.BlockData{blk}.Grp(a).(f{c})(b);
                            else
                                BlockData{blk}.Grp(a).(f{c})(end+1,:) = tmp.BlockData{blk}.Grp(a).(f{c})(b,:);
                            end
                        end
                    end
                end %end for subjects in the to-be-merged data structure
            
            end %end how to merge the data
        end %end for the 3 possible data groups
    end %end is the block empty or needs to be merged
end %end for all blocks

%%
tosave = input('Replace existing data file? (y/n)','s');
if strcmpi(tosave,'y')
    save(fullfile(pname1,fname1),'BlockData');
else
    uisave('BlockData','compiled_data.mat');
end

