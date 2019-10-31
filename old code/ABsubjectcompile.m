%do a simple data compile across all subjects

function [] = ABsubjectcompile(varargin)

if nargin == 0
%     fprintf('\nSelect data subject folder(s) to be compiled.\n');
%     SUBJpaths = uigetmultidir('/Users/mrri/Desktop/Kinematic Apraxia Battery','Select data folder(s) to be compiled.');
    %the data have been labeled as control or patient so we can auto-sort it!
    
    fprintf('\nSelect data subject file(s) to be compiled.\n');
    [SUBJpaths,SUBJpathdir] = uigetfile('/Users/mrri/Desktop/Kinematic Apraxia Battery/*.mat','Select data file(s) to be compiled.','Multiselect','on');
    
else
    
    SUBJpaths = varargin{1};
    
end

if isempty(SUBJpaths) || (isnumeric(SUBJpaths) && SUBJpaths == 0)
    fprintf('\n\nNo files selected for analysis.\n\n')
    return;
end

if nargin < 2
    fprintf('\nSelect data file to append to (or cancel to create a new data file).\n');
    [datafname,datafpath] = uigetfile('*.mat','Select data file to append to.');
else
    datafname = varargin{2};
    datafpath = '';
end

if ~isempty(datafname) && all(datafname ~= 0)
    
    %load existing data file
    
    load(fullfile(datafpath,datafname));
    
else

    %no existing data file to append to, initialize everything to zero
    
    BlockData = cell(10,1);
    
end

%%
fname = '*_data.mat';


for a = 1:length(SUBJpaths)

%     % this code below handles the case when multiple folders are identified and we have to find the data file in the folder.    
%     tmp = sort([strfind(SUBJpaths{a},'/') strfind(SUBJpaths{a},'\')]);
% %     
% %     %tmpfname = regexp(SUBJpaths{a}(tmp(end-1)+1:end),'\d*','Match');
% %     tmpfname = strfind(SUBJpaths{a}(tmp(end-1)+1:end),'AB');
% %     if isempty(tmpfname)
% %         tmpfname = strfind(SUBJpaths{a}(tmp(end-1)+1:end),'CS');
% %     end
% %     if ~isempty(tmpfname)
% %         tmpfname = SUBJpaths{a}(tmp(end-1)+tmpfname:end);
% %         tmpfname = tmpfname{1};
% %         tmpfname = strrep(fname,'*',tmpfname);
% %     else  %we assume this is the model, which is subject 0
% %         tmpfname = strrep(fname,'*','0');
% %     end
% 
%     %we will assume that the first data.mat file in the folder is our data
%     %file
%     tmpdir = dir(fullfile(SUBJpaths{a},'*data.mat'));
%     if ~isempty(tmpdir)
%         tmpfname = tmpdir(1).name;
%     else
%         fprintf('No valid data file exists in folder .\\%s\n',SUBJpaths{a}(tmp(end-1)+1:end));
%         continue;
%     end
%     
%     tmp = load(fullfile(SUBJpaths{a},tmpfname));

    %if we are choosing data files instead of folders, we just have to load the file
    tmp = load(fullfile(SUBJpathdir,SUBJpaths{a}));

    fprintf('Processing data for %s\n',SUBJpaths{a});
    
    if strcmpi(tmp.Group,'model')
        
        for b = 1:length(tmp.BlockData)
            if isempty(tmp.BlockData{b})
                continue;
            end
            
            for c = 1:size(tmp.BlockData{b}.Data,1)
                for c2 = 1:size(tmp.BlockData{b}.Data,2)
                    tmppos = [];
                    tmprotang = [];
                    tmpang = [];
                    for d = 1:length(tmp.BlockData{b}.Data(c,c2).m)
                        tmppos = cat(3,tmppos,tmp.BlockData{b}.Data(c,c2).m(d).pos);
                        tmprotang = cat(4,tmprotang,tmp.BlockData{b}.Data(c,c2).m(d).rotang);
                        tmpang = cat(3,tmpang,tmp.BlockData{b}.Data(c,c2).m(d).angs);
                    end
                    BlockData{b}.Grp(1).Group = 'Model';
                    BlockData{b}.Grp(1).pos{1,c,c2} = tmppos;
                    BlockData{b}.Grp(1).ja{1,c,c2} = tmp.BlockData{b}.Data(c,c2).ja;
                    BlockData{b}.Grp(1).rotang{1,c,c2} = tmprotang;
                    BlockData{b}.Grp(1).ang{1,c,c2} = tmpang;
                    BlockData{b}.Grp(1).time{1,c,c2} = tmp.BlockData{b}.Data(c,c2).time;
                end
            end
            
            BlockData{b}.BlockID = tmp.BlockData{b}.BlockID;
            BlockData{b}.BlockName = tmp.BlockData{b}.BlockName;
            BlockData{b}.Items = tmp.BlockData{b}.Items;
            if size(tmp.BlockData{b}.Data,2) > 1
                BlockData{b}.Grp(1).VisFdbk = {'yes','no'};  %we will hard code this for now...
            else
                BlockData{b}.Grp(1).VisFdbk = {'yes'};
            end
            BlockData{b}.Grp(1).Subjects{1} = '0';
            BlockData{b}.Grp(1).SubjN = 0;
            BlockData{b}.Grp(1).SubjRep = 1;
            
            
        end
        
    elseif strcmpi(tmp.Group,'control')
        
        for b = 1:length(tmp.BlockData)
            if isempty(tmp.BlockData{b})
                continue;
            end
            
            if ~isfield(BlockData{b}, 'Grp') || length(BlockData{b}.Grp) < 2 || isempty(BlockData{b}.Grp(2)) || ~isfield(BlockData{b}.Grp(2),'Subjects') || isempty(BlockData{b}.Grp(2).Subjects)
                %first subject in this group
                ictl = 1;
                repn = 1;
                
            elseif all(~strcmpi(BlockData{b}.Grp(2).Subjects, tmp.BlockData{b}.SubjID))
                %add a new subject to this group
                ictl = length(BlockData{b}.Grp(2).Subjects)+1;
                repn = 1;
                
            elseif any(strcmpi(BlockData{b}.Grp(2).Subjects, tmp.BlockData{b}.SubjID))
                %subject exists in this group, decide what to do:
                dorep = input(sprintf('Subject %s already exists. (r)eplace, (n)ew entry, or (s)kip: ',tmp.BlockData{b}.SubjID),'s');
                if strcmpi(dorep,'r') %overwrite
                    duplicates = find(strcmpi(BlockData{b}.Grp(2).Subjects, tmp.BlockData{b}.SubjID) == 1); %find( == tmp.BlockData{b}.SubjID);
                    if length(duplicates) > 1
                        whichentry = input(sprintf(' %d entries exist for %s. Select which to overwrite: ',length(duplicates),tmp.BlockData{b}.SubjID));
                        ictl = duplicates(whichentry);
                        repn = BlockData{b}.Grp(2).SubjRep(ictl);
                    else
                    	ictl = duplicates;
                        fprintf('  Existing entry of Subject %s overwritten.\n',tmp.BlockData{b}.SubjID);
                        repn = 1;
                    end
                elseif strcmpi(dorep,'n') %add a new entry
                    ictl = length(BlockData{b}.Grp(2).Subjects)+1;
                    fprintf('  Subject %s duplicate entry created.\n',tmp.BlockData{b}.SubjID);
                    duplicates = find(strcmpi(BlockData{b}.Grp(2).Subjects, tmp.BlockData{b}.SubjID) == 1); %find( == tmp.BlockData{b}.SubjID);
                    repn = length(duplicates)+1;
                else 
                    fprintf('  Subject %s skipped.\n',tmp.BlockData{b}.SubjID);
                    continue;
                end
            
            end
            
            %originalfields = fieldnames(tmp.BlockData{b}.Data);
            for c = 1:size(tmp.BlockData{b}.Data,1)
                for c2 = 1:size(tmp.BlockData{b}.Data,2)
                    tmppos = [];
                    tmprotang = [];
                    tmpang = [];
                    for d = 1:length(tmp.BlockData{b}.Data(c,c2).m)
                        tmppos = cat(3,tmppos,tmp.BlockData{b}.Data(c,c2).m(d).pos);
                        tmprotang = cat(4,tmprotang,tmp.BlockData{b}.Data(c,c2).m(d).rotang);
                        tmpang = cat(3,tmpang,tmp.BlockData{b}.Data(c,c2).m(d).angs);
                    end
                    BlockData{b}.Grp(2).Group = 'Control';
                    BlockData{b}.Grp(2).pos{ictl,c,c2} = tmppos;
                    BlockData{b}.Grp(2).ja{ictl,c,c2} = tmp.BlockData{b}.Data(c,c2).ja;
                    BlockData{b}.Grp(2).rotang{ictl,c,c2} = tmprotang;
                    BlockData{b}.Grp(2).ang{ictl,c,c2} = tmpang;
                    BlockData{b}.Grp(2).time{ictl,c,c2} = tmp.BlockData{b}.Data(c,c2).time;
                end
            end
            
            BlockData{b}.BlockID = tmp.BlockData{b}.BlockID;
            BlockData{b}.BlockName = tmp.BlockData{b}.BlockName;
            BlockData{b}.Items = tmp.BlockData{b}.Items;
            if size(tmp.BlockData{b}.Data,2) > 1
                BlockData{b}.Grp(2).VisFdbk = {'yes','no'};  %we will hard code this for now...
            else
                BlockData{b}.Grp(2).VisFdbk = {'yes'};
            end
            BlockData{b}.Grp(2).Subjects{ictl} = tmp.BlockData{b}.SubjID;      %record the subject ID to the list of subjects
            BlockData{b}.Grp(2).SubjN(ictl) = str2double(tmp.BlockData{b}.SubjID); 
            BlockData{b}.Grp(2).SubjRep(ictl) = repn;
        end
        
    elseif strcmpi(tmp.Group,'patient')
        
        for b = 1:length(tmp.BlockData)
            if isempty(tmp.BlockData{b})
                continue;
            end
            
            if ~isfield(BlockData{b}, 'Grp') || length(BlockData{b}.Grp) < 3 ||  isempty(BlockData{b}.Grp(3)) || ~isfield(BlockData{b}.Grp(3),'Subjects') || isempty(BlockData{b}.Grp(3).Subjects)
                %first subject in this group
                ipat = 1;
                repn = 1;
                
            elseif all(~strcmpi(BlockData{b}.Grp(3).Subjects, tmp.BlockData{b}.SubjID))
                %add a new subject to this group
                ipat = length(BlockData{b}.Grp(3).Subjects)+1;
                repn = 1;
                
            elseif any(strcmpi(BlockData{b}.Grp(3).Subjects, tmp.BlockData{b}.SubjID))
                %subject exists in this group, decide what to do:
                dorep = input(sprintf('Subject %s already exists. (r)eplace, (n)ew entry, or (s)kip: ',tmp.BlockData{b}.SubjID),'s');
                if strcmpi(dorep,'r') %overwrite
                    duplicates = find(strcmpi(BlockData{b}.Grp(3).Subjects, tmp.BlockData{b}.SubjID) == 1); %find( == tmp.BlockData{b}.SubjID);
                    if length(duplicates) > 1
                        whichentry = input(sprintf(' %d entries exist for %s. Select which to overwrite: ',length(duplicates),tmp.BlockData{b}.SubjID));
                        ipat = duplicates(whichentry);
                        repn = BlockData{b}.Grp(3).SubjRep(ictl);
                    else
                    	ipat = duplicates;
                        fprintf('  Existing entry of Subject %s overwritten.\n',tmp.BlockData{b}.SubjID);
                        repn = 1;
                    end
                elseif strcmpi(dorep,'n') %add a new entry
                    ipat = length(BlockData{b}.Grp(3).Subjects)+1;
                    fprintf('  Subject %s duplicate entry created.\n',tmp.BlockData{b}.SubjID);
                    duplicates = find(strcmpi(BlockData{b}.Grp(3).Subjects, tmp.BlockData{b}.SubjID) == 1);
                    repn = length(duplicates)+1;
                    
                else 
                    fprintf('  Subject %s skipped.\n',tmp.BlockData{b}.SubjID);
                    continue;
                end
            
            end
            
            %originalfields = fieldnames(tmp.BlockData{b}.Data);
            for c = 1:size(tmp.BlockData{b}.Data,1)
                for c2 = 1:size(tmp.BlockData{b}.Data,2)
                    tmppos = [];
                    tmprotang = [];
                    tmpang = [];
                    for d = 1:length(tmp.BlockData{b}.Data(c,c2).m)
                        tmppos = cat(3,tmppos,tmp.BlockData{b}.Data(c,c2).m(d).pos);
                        tmprotang = cat(4,tmprotang,tmp.BlockData{b}.Data(c,c2).m(d).rotang);
                        tmpang = cat(3,tmpang,tmp.BlockData{b}.Data(c,c2).m(d).angs);
                    end
                    BlockData{b}.Grp(3).Group = 'Patients';
                    BlockData{b}.Grp(3).pos{ipat,c,c2} = tmppos;
                    BlockData{b}.Grp(3).ja{ipat,c,c2} = tmp.BlockData{b}.Data(c,c2).ja;
                    BlockData{b}.Grp(3).rotang{ipat,c,c2} = tmprotang;
                    BlockData{b}.Grp(3).ang{ipat,c,c2} = tmpang;
                    if isfield(tmp.BlockData{b}.Data(c,c2),'time')
                        BlockData{b}.Grp(3).time{ipat,c,c2} = tmp.BlockData{b}.Data(c,c2).time;
                    elseif isfield(tmp.BlockData{b}.Data(c,c2),'Time')
                        BlockData{b}.Grp(3).time{ipat,c,c2} = tmp.BlockData{b}.Data(c,c2).Time;
                    end
                end
            end
            
            BlockData{b}.BlockID = tmp.BlockData{b}.BlockID;
            BlockData{b}.BlockName = tmp.BlockData{b}.BlockName;
            BlockData{b}.Items = tmp.BlockData{b}.Items;
            if size(tmp.BlockData{b}.Data,2) > 1
                BlockData{b}.Grp(3).VisFdbk = {'yes','no'};  %we will hard code this for now...
            else
                BlockData{b}.Grp(3).VisFdbk = {'yes'};
            end
            BlockData{b}.Grp(3).Subjects{ipat} = tmp.BlockData{b}.SubjID;      %record the subject ID to the list of subjects
            BlockData{b}.Grp(3).SubjN(ipat) = str2double(tmp.BlockData{b}.SubjID); 
            BlockData{b}.Grp(3).SubjRep(ipat) = repn;
        end
        
        
        
%         for b = 1:length(tmp.BlockData)
%             if isempty(tmp.BlockData{b})
%                 continue;
%             end
%             
%             overwriteflag = 1;
%             
%             if ~isfield(BlockData{b}, 'Patient') || isempty(BlockData{b}.Patient)
%                 %first subject in this group
%                 ipat = 1;
%                 
%             elseif all(BlockData{b}.Patient ~= tmp.BlockData{b}.SubjID)
%                 %add a new subject to this group
%                 iclt = length(BlockData{b}.Controls)+1;
%                 
%             elseif any(BlockData{b}.Patient == tmp.BlockData{b}.SubjID)
%                 %subject exists in this group, decide what to do:
%                 dorep = input(sprintf('Subject %d already exists. (o)verwrite fully, (r)eplace data but preserve markings, (n)ew entry, or (s)kip: ',subjID),'s');
%                 if strcmpi(dorep,'o') %overwrite
%                     duplicates = find(BlockData{b}.Patient == tmp.BlockData{b}.SubjID);
%                     if length(duplicates) > 1
%                         whichentry = input(sprintf(' %d entries exist for %d. Select which to overwrite: ',length(duplicates),subjID));
%                         ipat = duplicates(whichentry);
%                     else
%                     	ipat = duplicates;
%                         fprintf('  Existing entry of Subject %d overwritten.\n',subjID);
%                     end
%                 elseif strcmpi(dorep,'r') %replace
%                     ipat = find(BlockData{b}.Patient == tmp.BlockData{b}.SubjID);
%                     fprintf('  Subject %d replaced.\n',subjID);
%                     overwriteflag = 0;
%                 elseif strcmpi(dorep,'n') %add a new entry
%                     iclt = length(BlockData{b}.Patient)+1;
%                     fprintf('  Subject %d duplicate entry created.\n',subjID);
%                 else 
%                     fprintf('  Subject %d skipped.\n',subjID);
%                     continue;
%                 end
%             
%             end
%             
%             if overwriteflag == 0
%                 originalfields = fieldnames(tmp.BlockData{b}.Data);
%                 for c = 1:length(tmp.BlockData{b}.Data)
%                     for d = 1:length(originalfields)
%                         BlockData{b}.Pat(ipat).Data(c).(originalfields{d}) = tmp.BlockData{b}.Data(c).(originalfields{d});
%                     end
%                 end
%             else
%                 BlockData{b}.Pat(ipat).Data = tmp.BlockData{b}.Data;    %transfer subject data into a single structure
%             end
%             
%             BlockData{b}.BlockID = tmp.BlockData{b}.BlockID;
%             BlockData{b}.BlockName = tmp.BlockData{b}.BlockName;
%             BlockData{b}.Items = tmp.BlockData{b}.Items;
%             BlockData{b}.Pat(ipat).SubjID = tmp.BlockData{b}.SubjID;    %Record the Subject ID
%             BlockData{b}.Patient(ipat) = tmp.BlockData{b}.SubjID;      %add subject ID to the list of subjects
%         end
        
        
    end
    
end


%%

clear tmp* fname;

if ~isempty(datafname) && all(datafname ~= 0)
    save(fullfile(datafpath,datafname),'BlockData','-v7.3');
else
    %uisave('BlockData','compiled_data.mat');
    fprintf('\n\nChoose path and file name to save as...\n\n');
    [savefile,savepath] = uiputfile('*.mat','Save as','compiled_data.mat');
    save(fullfile(savepath,savefile),'BlockData','-v7.3');
end

fprintf('\nFile saved!\n\n');

