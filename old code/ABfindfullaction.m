%In this section, we pre-process each trajectory by automatically marking
%  the "full" action (from lift-off to return) based on a velocity
%  criterion, and asking the user to manually verify the marks. This will
%  allow us to perform a faster analysis on the kinematic data without
%  needing to identify the "core" portion of the action.
%You will be prompted to select which blocks, items, and subjects you wish
%  to process. This information can be entered in two ways: either by
%  indexing into the BlockData cell appropriately, or by typing the name.
%  For example, you can ask for:
%      Block ID/Name: 
%          1, 2, 4
%          meaningful_named, meaningless
%      Item ID/Name:
%           4, 8
%           comb, nailclippers
%      Subject (Patient or Control):
%           2, 4, 6
%           S1000, S2000
%
%Refer to the instructions in the markdataGUI file for detailed information
%  regarding how to view and/or mark points in the data. Note that it is
%  necessary to close the markdataGUI window when you are done
%  viewing/marking data from one participant/trial in order to move on to
%  the next one. Be sure that there are either only 0 or 2 marks on each
%  trial. Anything else may throw an error later in post-processing.
%At the conclusion of this section, the data are resaved. It is possible to
%  modify the various loops by hand to analyze different subsets of data,
%  but the best way is to use the input dialog box. The load/resave
%  commands at the start and end of this section will help to ensure that
%  your previous work is preserved.  

clear all;
close all;


fprintf('\nSelect data file to process.\n');
[fname,fpath] = uigetfile('*.mat','Select data file to process.');

load(fullfile(fpath,fname));



%% ***********************************************************************

%Hard code: Which trials to analyze? Identify block numbers and item
%           numbers to analyze (write in list format, e.g., "1,2,4", or
%           set to "-1" to do all existing blocks/items). Also do the same
%           for controls and patients (i.e., to choose to analyze only a
%           subset of the existing patients and controls in the current
%           datafile. For controls and patients, if you wish to analyze NO
%           files, enter "0". In all cases, the index of the desired
%           block/item/subject should be noted, NOT the subject number. If
%           you wish to enter subjectID numbers, enter them preceeded by an
%           "S", e.g. S1000, S2000, S3000".
BlockNumbers = -1;
ItemNumbers = -1;
CtlNumbers = -1;
PatNumbers = -1;
ClearMarks = -1;

inputVals = inputdlg({'Block Numbers/Names:','Item Numbers/Names:','Control Numbers: ','Patient Numbers: ','Clear Existing Marks (1 = Yes, 0 = No)'},'Input analysis parameters.',1,{'-1','-1','-1 or S000','-1 or S000','0'},'on');

%parse the block number/names
BlockNumbers = str2num(inputVals{1}); %see if the blocks were entered numerically
if isempty(BlockNumbers) %entered as block names, not indexes
    BlockNumbers = inputVals{1};
    BlockNumbers = strrep(BlockNumbers,' ','');
    tmpdelim = strfind(BlockNumbers,',');
    tmpdelim = [0 tmpdelim length(BlockNumbers)+1];
    cn = [];
    for a = 1:length(tmpdelim)-1
        cn(a) = decodetrials(BlockNumbers(tmpdelim(a)+1:tmpdelim(a+1)-1));
    end
    BlockNumbers = cn; %save the decoded block numbers
end
if isempty(BlockNumbers) || BlockNumbers == -1 %if we asked for all of them or invalid blocks, search for any valid blocks in the data file
    cn = [];
    for blk = 1:length(BlockData)
        if isfield(BlockData{blk},'Grp') && isfield(BlockData{blk}.Grp(1),'Subjects') && (~isempty(BlockData{blk}.Grp(1).Subjects) || ~isempty(BlockData{blk}.Grp(2).Subjects) || ~isempty(BlockData{blk}.Grp(3).Subjects))
            %the block has valid data in it
            cn = [cn blk];
        end
    end
    BlockNumbers = cn;
end
if isempty(BlockNumbers) %no valid blocks were found!
    fprintf('\n\nERROR: No valid blocks found.\n\n');
    return;
end

%parse the item number/names
ItemNumbers = str2num(inputVals{2}); %see if the items were entered numerically
if isempty(ItemNumbers) %entered as item names, not indexes
    ItemNumbers = inputVals{2};
    ItemNumbers = strrep(ItemNumbers,' ','');
    tmpdelim = strfind(ItemNumbers,',');
    tmpdelim = [0 tmpdelim length(ItemNumbers)+1];
    cn = [];
    for a = 1:length(tmpdelim)-1
        %we need to know the block number! We will assume it's the first
        %one in the BlockNumbers list.
        cn(a) = decodetrials(BlockNumbers(1),ItemNumbers(tmpdelim(a)+1:tmpdelim(a+1)-1));
    end
    ItemNumbers = cn;
end
if isempty(ItemNumbers) %if no item numbers could be decoded, error out!
    fprintf('\n\nERROR: No valid items found.\n\n');
    return;
end

%parse the control subjects
CtlNumbers = inputVals{3};
if isempty(CtlNumbers) || (~isempty(str2num(CtlNumbers)) && any(str2num(CtlNumbers) == 0))
    CtlNumbers = 0;
elseif strcmpi(CtlNumbers,'-1 or S000') %asked to search for all valid options
    CtlNumbers = []; 
elseif ~contains(lower(CtlNumbers),'s') %subject indices entered numerically
    CtlNumbers = str2num(CtlNumbers);
else %subject ID numbers entered
    tmps = strfind(lower(CtlNumbers),'s');
    tmpdelim = strfind(CtlNumbers,',');
    if length(tmps) > length(tmpdelim)
        tmpdelim(end+1) = length(CtlNumbers)+1;
    end
    cn = {};
    for a = 1:length(tmps)
        cn{a} = CtlNumbers(tmps(a)+1:tmpdelim(a)-1);
    end
    CtlNumbers = cn;
end

%parse the patient numbers
PatNumbers = inputVals{4};
if isempty(PatNumbers) || (~isempty(str2num(PatNumbers)) && any(str2num(PatNumbers) == 0))
    PatNumbers = 0;
elseif strcmpi(PatNumbers,'-1 or S000') %asked to search for all valid options
    PatNumbers = [];
elseif ~contains(lower(PatNumbers),'s') %subject indices entered numerically
    PatNumbers = str2num(PatNumbers);
else %subject ID numbers entered
    tmps = strfind(lower(PatNumbers),'s');
    tmpdelim = strfind(PatNumbers,',');
    if length(tmps) > length(tmpdelim)
        tmpdelim(end+1) = length(PatNumbers)+1;
    end
    cn = {};
    for a = 1:length(tmps)
        cn{a} = PatNumbers(tmps(a)+1:tmpdelim(a)-1);
    end
    PatNumbers = cn;
end

ClearMarks = str2num(inputVals{5});
if isempty(ClearMarks)
    ClearMarks = 0;
elseif ClearMarks > 0
    ClearMarks = 1;
end


%% ***********************************************************************


domarkmodel = 1;
showmarkmodel = 0;

%parse which blocks to analyze, or do them all
%if isempty(BlockNumbers) || BlockNumbers == -1
%    BlockNumbers = [1:length(BlockData)];
%end

fprintf('\n\n');

loopcont = 1;

for blk = BlockNumbers %specify which blocks to mark

    if isempty(BlockData{blk})
        fprintf('No data available for block %d, block will be skipped.\n',blk);
        continue;
    end
    
    %add velocity to dataset
    BlockData{blk}.Grp = addvel(BlockData{blk}.Grp);
    
    %parse what the loop indices below should be based on the inputs above.
    if isempty(ItemNumbers) || ItemNumbers == -1  %figure out which items to do
        if ~isempty(BlockData{blk}.Grp(1).pos)
            ItemNumbers = [1:size(BlockData{blk}.Grp(1).pos,2)];
        elseif ~isempty(BlockData{blk}.Grp(2).pos)
            ItemNumbers = [1:size(BlockData{blk}.Grp(2).pos,2)];
        elseif ~isempty(BlockData{blk}.Grp(3).pos)
            ItemNumbers = [1:size(BlockData{blk}.Grp(3).pos,2)];
        else
            fprintf('\nNo items found for block %d; block is skipped.\n',blk);
            continue;
        end
    end
    if (length(BlockData{blk}.Grp) >= 2 && ~isempty(BlockData{blk}.Grp(2).Subjects))
        if ~iscell(CtlNumbers) && any(CtlNumbers == 0)
            CtlNumbers = [];
        elseif isempty(CtlNumbers) || ~iscell(CtlNumbers)
            if isempty(CtlNumbers) || all(CtlNumbers == -1)
                CtlNumbers = [1:length(BlockData{blk}.Grp(2).Subjects)];
            end
            CtlInds = CtlNumbers;
        else
            sinds = [];
            for a = 1:length(CtlNumbers)
                
                tmpstrcmp = strcmpi(BlockData{blk}.Grp(2).Subjects, CtlNumbers{a});
                if isempty(tmpstrcmp)
                    fprintf('  Control S%s not found.\n',CtlNumbers{a});
                elseif sum(tmpstrcmp) > 1
                    fprintf('  %d multiple entires found for Control S%s.\n',sum(tmpstrcmp), CtlNumbers{a});
                    repn = input('    Enter which repetition(s) would you like to analyze: ','s');
                    repn = str2num(repn);
                    sreps = find(tmpstrcmp == 1);
                    sinds = [sinds sreps(repn)];
                else
                    sinds = [sinds find(tmpstrcmp == 1)];
                end
            end
            CtlInds = sinds;
        end
    else
        CtlInds = [];  %no data exists, so we will skip it
    end
    if (length(BlockData{blk}.Grp) >= 3 && ~isempty(BlockData{blk}.Grp(3).Subjects))
        if ~iscell(PatNumbers) && any(PatNumbers == 0)
            PatNumbers = [];
        elseif isempty(PatNumbers) || ~iscell(PatNumbers)
            if isempty(PatNumbers) || all(PatNumbers == -1)
                PatNumbers = [1:length(BlockData{blk}.Grp(3).Subjects)];
            end
            PatInds = PatNumbers;
        else
            sinds = [];
            for a = 1:length(PatNumbers)
                tmpstrcmp = strcmpi(BlockData{blk}.Grp(3).Subjects, PatNumbers{a});
                
                if isempty(tmpstrcmp)
                    fprintf('  Patient S%s not found.\n',PatNumbers{a});
                elseif sum(tmpstrcmp) > 1
                    fprintf('  %d multiple entires found for Patient S%s.\n',sum(tmpstrcmp), PatNumbers{a});
                    repn = input('    Enter which repetition(s) would you like to analyze: ','s');
                    repn = str2num(repn);
                    sreps = find(tmpstrcmp == 1);
                    sinds = [sinds sreps(repn)];
                else
                    sinds = [sinds find(tmpstrcmp == 1)];
                end
            end
            PatInds = sinds;
        end
    else
        PatInds = [];   %no data exists, so we will skip it
    end
    
    
    if isfield(BlockData{blk},'Grp') && isfield(BlockData{blk}.Grp(1),'inds') && ~isempty(BlockData{blk}.Grp(1).inds)
        domarkmodel = input('Remark Model? (Y = 1, N = 0): ');
    end
    if (~domarkmodel) && (isfield(BlockData{blk},'Grp') && isfield(BlockData{blk}.Grp(1),'inds')) && ~isempty(BlockData{blk}.Grp(1).inds)
        showmarkmodel = input('Show Model? (Y = 1, N = 0): ');
    end
    

    for c = ItemNumbers %process all the movements in this block, or set to fixed values to process only a subset of trials

        %loops to control what data get processed
        if (domarkmodel) && isfield(BlockData{blk},'Grp') && isfield(BlockData{blk}.Grp(1),'pos') && (size(BlockData{blk}.Grp(1).pos,2) >= c)
            
            if ~isfield(BlockData{blk}.Grp(1),'inds') || (length(BlockData{blk}.Grp(1).inds)<c) || ~isfield(BlockData{blk}.Grp(1).inds(1,c),'Full') || isempty(BlockData{blk}.Grp(1).inds(1,c).Full) || all(isnan(BlockData{blk}.Grp(1).inds(1,c).Full))
                inds{1} = [];
            else
                inds{1} = BlockData{blk}.Grp(1).inds(1,c).Full;
            end
            %if ~isfield(BlockData{blk}.Grp(1),'inds') || (length(BlockData{blk}.Grp(1).inds)<c) || ~isfield(BlockData{blk}.Grp(1).inds(1,c),'SubAction') || isempty(BlockData{blk}.Grp(1).inds(1,c).SubAction) || all(isnan(BlockData{blk}.Grp(1).inds(1,c).SubAction))
            inds{2} = [];
            %else
            %    inds{2} = BlockData{blk}.Grp(1).inds(1,c).SubAction;
            %end
            if isempty(inds{1}) && isempty(inds{2}) || ClearMarks == 1
                inds = addmarks(BlockData{blk}.Grp(1).vel{1,c},'Nchan',2,'vthreshMin',8,'throwmid');
            end
            
            inds = markdataGUI(BlockData{blk}.Grp(1).pos{1,c},'title',sprintf('Model: %s',BlockData{blk}.Items{c}),'ang',BlockData{blk}.Grp(1).ang{1,c},'mark',inds,'full');
            if isempty(inds{1})
                BlockData{blk}.Grp(1).inds(1,c).Full = [NaN NaN];
            else
                BlockData{blk}.Grp(1).inds(1,c).Full = inds{1};
            end
            %if isempty(inds{2})
            %    BlockData{blk}.Grp(1).inds(1,c).SubAction = [];
            %else
            %    BlockData{blk}.Grp(1).inds(1,c).SubAction = inds{2};
            %end
            
        elseif (showmarkmodel) && isfield(BlockData{blk},'Grp') && isfield(BlockData{blk}.Grp(1),'pos') && (size(BlockData{blk}.Grp(1).pos,2) >= c)
            
            if ~isfield(BlockData{blk}.Grp(1),'inds') || (length(BlockData{blk}.Grp(1).inds)<c) || ~isfield(BlockData{blk}.Grp(1).inds(1,c),'Full') || isempty(BlockData{blk}.Grp(1).inds(1,c).Full) || all(isnan(BlockData{blk}.Grp(1).inds(1,c).Full))
                inds{1} = [];
            else
                inds{1} = BlockData{blk}.Grp(1).inds(1,c).Full;
            end
            %if ~isfield(BlockData{blk}.Grp(1),'inds') || (length(BlockData{blk}.Grp(1).inds)<c) || ~isfield(BlockData{blk}.Grp(1).inds(1,c),'SubAction') || isempty(BlockData{blk}.Grp(1).inds(1,c).SubAction) || all(isnan(BlockData{blk}.Grp(1).inds(1,c).SubAction))
            inds{2} = [];
            %else
            %    inds{2} = BlockData{blk}.Grp(1).inds(1,c).SubAction;
            %end
            
            if isempty(inds{1}) && isempty(inds{2})
                inds = addmarks(BlockData{blk}.Grp(1).vel{1,c},'Nchan',2,'vthreshMin',8,'throwmid');
            end

            tmp = markdataGUI(BlockData{blk}.Grp(1).pos{1,c},'title',sprintf('Model: %s',BlockData{blk}.Items{c}),'ang',BlockData{blk}.Grp(1).ang{1,c},'mark',inds,'full');
            
            if ~isequal(tmp,inds)
                fprintf('\n\nChange of markers detected.\n')
                updateMdl = input(' Do you want to update? (Y = 1, N = 0): ');
                if updateMdl == 1
                    if isempty(tmp{1})
                        BlockData{blk}.Grp(1).inds(1,c).Full = [NaN NaN];
                    else
                        BlockData{blk}.Grp(1).inds(1,c).Full = tmp{1};
                    end
                    %if isempty(tmp{2})
                    %    BlockData{blk}.Grp(1).inds(1,c).SubAction = [];
                    %else
                    %    BlockData{blk}.Grp(1).inds(1,c).SubAction = tmp{2};
                    %end
                end
            end
            
        end %end if Mdl
        
        
        if ~isempty(CtlInds) && (size(BlockData{blk}.Grp(2).pos,2) >= c)
            
            for b = CtlInds %process all specified subjects for that movement
                
                %fprintf('  Subject: %d\n',BlockData{blk}.Grp(2).subj(b));
                if isempty(BlockData{blk}.Grp(2).pos{b,c})
                    fprintf('Item %d not available for Subject %s - Rep %d',b,BlockData{blk}.Grp(2).Subjects{b},BlockData{blk}.Grp(2).SubjRep(b));
                    continue;
                end
                
                if ~isfield(BlockData{blk}.Grp(2),'inds') || (size(BlockData{blk}.Grp(2).inds,1) < b || size(BlockData{blk}.Grp(2).inds,2) < c) || isempty(BlockData{blk}.Grp(2).inds(b,c)) || ~isfield(BlockData{blk}.Grp(2).inds(b,c),'Full') || isempty(BlockData{blk}.Grp(2).inds(b,c).Full) || all(isnan(BlockData{blk}.Grp(2).inds(b,c).Full))
                    inds{1} = [];
                else
                    inds{1} = BlockData{blk}.Grp(2).inds(b,c).Full;
                end
                %if ~isfield(BlockData{blk}.Grp(2),'inds') || (size(BlockData{blk}.Grp(2).inds,1) < b || size(BlockData{blk}.Grp(2).inds,2) < c) || isempty(BlockData{blk}.Grp(2).inds(b,c)) || ~isfield(BlockData{blk}.Grp(2).inds(b,c),'SubAction') || isempty(BlockData{blk}.Grp(2).inds(b,c).SubAction) || all(isnan(BlockData{blk}.Grp(2).inds(b,c).SubAction))
                inds{2} = [];
                %else
                %    inds{2} = BlockData{blk}.Grp(2).inds(b,c).SubAction;
                %end
                
                if isempty(inds{1}) && isempty(inds{2}) || ClearMarks == 1
                    inds = addmarks(BlockData{blk}.Grp(2).vel{b,c},'Nchan',2,'vthreshMin',8,'throwmid');
                end

                inds = markdataGUI(BlockData{blk}.Grp(2).pos{b,c},'title',sprintf('S%s-Rep%d: %s',BlockData{blk}.Grp(2).Subjects{b},BlockData{blk}.Grp(2).SubjRep(b),BlockData{blk}.Items{c}),'ang',BlockData{blk}.Grp(2).ang{b,c},'mark',inds,'full');
                
                if isempty(inds{1})
                    BlockData{blk}.Grp(2).inds(b,c).Full = [NaN NaN];
                else
                    BlockData{blk}.Grp(2).inds(b,c).Full = inds{1};
                end
                %if isempty(inds{2})
                %    BlockData{blk}.Grp(2).inds(b,c).SubAction = [];
                %else
                %    BlockData{blk}.Grp(2).inds(b,c).SubAction = inds{2};
                %end
                
            end
        end %end if Ctl
                
        if ~isempty(PatInds) && (size(BlockData{blk}.Grp(3).pos,2) >= c)
            
            for b = PatInds %process all the specified subjects, for that movement
                
                %fprintf('  Subject: %d\n',BlockData{blk}.Grp(3).subj(b));
                
                if isempty(BlockData{blk}.Grp(3).pos{b,c})
                    fprintf('Item %d not available for Subject %s - Rep %d',b,BlockData{blk}.Grp(3).Subjects{b},BlockData{blk}.Grp(3).SubjRep(b));
                    continue;
                end

                
                if ~isfield(BlockData{blk}.Grp(3),'inds') || (size(BlockData{blk}.Grp(3).inds,1) < b || size(BlockData{blk}.Grp(3).inds,2) < c) || isempty(BlockData{blk}.Grp(3).inds(b,c)) || ~isfield(BlockData{blk}.Grp(3).inds(b,c),'Full') || isempty(BlockData{blk}.Grp(3).inds(b,c).Full) || all(isnan(BlockData{blk}.Grp(3).inds(b,c).Full))
                    inds{1} = [];
                else
                    inds{1} = BlockData{blk}.Grp(3).inds(b,c).Full;
                end
                %if ~isfield(BlockData{blk}.Grp(3),'inds') || (size(BlockData{blk}.Grp(3).inds,1) < b || size(BlockData{blk}.Grp(3).inds,2) < c) || isempty(BlockData{blk}.Grp(3).inds(b,c)) || ~isfield(BlockData{blk}.Grp(3).inds(b,c),'SubAction') || isempty(BlockData{blk}.Grp(3).inds(b,c).SubAction) || all(isnan(BlockData{blk}.Grp(3).inds(b,c).SubAction))
                inds{2} = [];
                %else
                %    inds{2} = BlockData{blk}.Grp(3).inds(b,c).SubAction;
                %end
                
                if isempty(inds{1}) && isempty(inds{2}) || ClearMarks == 1
                    inds = addmarks(BlockData{blk}.Grp(3).vel{b,c},'Nchan',2,'vthreshMin',8,'throwmid');
                end
                
                inds = markdataGUI(BlockData{blk}.Grp(3).pos{b,c},'title',sprintf('S%s-Rep%d: %s',BlockData{blk}.Grp(3).Subjects{b},BlockData{blk}.Grp(3).SubjRep(b),BlockData{blk}.Items{c}),'ang',BlockData{blk}.Grp(3).ang{b,c},'mark',inds,'full');
                if isempty(inds{1})
                    BlockData{blk}.Grp(3).inds(b,c).Full = [NaN NaN];
                else
                    BlockData{blk}.Grp(3).inds(b,c).Full = inds{1};
                end
                %if isempty(inds{2})
                %    BlockData{blk}.Grp(3).inds(b,c).SubAction = [];
                %else
                %    BlockData{blk}.Grp(3).inds(b,c).SubAction = inds{2};
                %end
                
            end
        end %end if Pat
        
        
        loopcont = input('Continue to the next item? (1 = yes, 0 = no/quit): ');
        if loopcont == 0
            break;
        end

    end  %end for c

    if loopcont == 0
        break;
    end


end  %end for blk


%%


%resave the data to ensure that your marks have been recorded.
tosave = input('Save file (note, this will overwrite your existing datafile)? (Y = 1, N = 0): ');
if (tosave)
    
    %save('./GroupData_Austin_training.mat');
    save(fullfile(fpath,fname),'BlockData');
    
end




