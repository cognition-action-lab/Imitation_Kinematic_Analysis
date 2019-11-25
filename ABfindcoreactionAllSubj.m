%This code takes a set of processed data files at the block level (from one
%  or more subjects, all of the same experiment block), temporarily
%  combines them to support data marking, then breaks them up again
%  afterward. 
%We assume that ALL data, including the model, live in this folder. The
%  Model data should be named as "Model" for identification purposes. Note,
%  this may require rearranging the data so that the files coming out of
%  ABblockcompile get put into a unified folder, not back into the parent
%  subject folder!
%We will assume 1 entry per subject (for now). We can modify this later as
%  the need arises.
%
%In this section, we pre-process each trajectory by manually marking the
%  "core" portion of each the movement.
%You will be prompted to select which subjects and items you wish
%  to process. You will be prompted to select these from a list.
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

%choose which folder to process
fprintf('\nSelect data folder to process.\n');
%[fname,fpath] = uigetfile('*.mat','Select data file to process.');
fpath = uigetdir('*.*','Select data folder to process');

%extract the subjects in the folder
fnamelist = cellstr(ls(fpath));
[~,idx]=sort(fnamelist); %ensure the files are sorted alphabetically
fnamelist=fnamelist(idx);

%subgroup by type: model, control, patient
% imdl = 0;
% ipat = 0;
% ictl = 0;
% Model = {};
% Control = {};
% Patient = {};
isubj = 0;
Subj = {};
for a = 1:length(fnamelist)
%     if contains(fnamelist{a},'model','IgnoreCase',true)
%         imdl = imdl+1;
%         fnames.model{imdl} = fnamelist{a};
%         Model{imdl} = 'Model';
%     elseif contains(fnamelist{a},'cb','IgnoreCase',true)
%         ictl = ictl+1;
%         fnames.model{ictl} = fnamelist{a};
%         itmp = strfind(lower(fnamelist{a}),'cb');
%         isep = strfind(fnamelist{a},'_');
%         Control{ictl} = fnamelist{a}(itmp:isep(1)-1);
%     elseif contains(fnamelist{a},'ab','IgnoreCase',true)
%         ipat = ipat+1;
%         fnames.model{ipat} = fnamelist{a};
%         itmp = strfind(lower(fnamelist{a}),'ab');
%         isep = strfind(fnamelist{a},'_');
%         Patient{ipat} = fnamelist{a}(itmp:isep(1)-1);
    

    if contains(fnamelist{a},'model','IgnoreCase',true)
        isubj = isubj+1;
        fnames{isubj} = fnamelist{a};
        Subj{isubj} = 'Model';
    elseif contains(fnamelist{a},'cb','IgnoreCase',true)
        itmp = strfind(lower(fnamelist{a}),'cb');
        isep = strfind(fnamelist{a},'_');
        subjid = fnamelist{a}(itmp:isep(1)-1);
        if ~any(contains(Subj,subjid))        
            isubj = isubj+1;
            fnames{isubj} = fnamelist{a};
            Subj{isubj} = subjid;
        end
    elseif contains(fnamelist{a},'ab','IgnoreCase',true)
        itmp = strfind(lower(fnamelist{a}),'ab');
        isep = strfind(fnamelist{a},'_');
        subjid = fnamelist{a}(itmp:isep(1)-1);
        if ~any(contains(Subj,subjid))        
            isubj = isubj+1;
            fnames{isubj} = fnamelist{a};
            Subj{isubj} = subjid;
        end        
    elseif contains(fnamelist{a},'.mat','IgnoreCase',true)
        fprintf('\nFile %s not recognized.',fnamelist{a});
    end
end

%choose which of the existing subjects we actually want to process
isubj = listdlg('ListString',Subj,'PromptString',{'Choose subjects to process.','','Select multiple options by holding CNTL or Shift',''});
if isempty(isubj)
    error('Error: No subjects selected to be analyzed.')
end

Subj = Subj(isubj);
fnames = fnames(isubj);

%if the model is in the list, put it first
imdl = contains(Subj,'model','IgnoreCase',true);
if any(imdl)
    tmpSubj{1} = Subj{imdl};
    tmpSubj = [tmpSubj Subj(~imdl)];
    Subj = tmpSubj;
end

%next, choose what items to process. We'll do that by loading the first
%file and checking what items are there; we will assume that all the other
%files contain the same items!
temp = load(fullfile(fpath,fnames{1}));

items = temp.BlockData.Items;
iitems = listdlg('ListString',items,'PromptString',{'Choose items to process.','','Select multiple options by holding CNTL or Shift',''});
if isempty(iitems)
    error('Error: No items selected to be analyzed.')
end


%% step through the data; we will do this as a just-in-time load to save memory
%if the file has been marked previously it will have a _markdata appended
%to it in the file name. So we will look for existing files that contain
%the subject name and the _markdata tag; if it exists we will
%preferentially load that. Otherwise we will load the raw data file and
%do some additional pre-processing.

%we will do this repetitively for every item!

for b = iitems
    
    for a = 1:length(Subj)
        
        clear temp BlockData tmpsavfile;
        
        tmpfname = cellstr(ls(fpath));
        ifile = contains(tmpfname,Subj{a},'IgnoreCase',true);
        imark = contains(tmpfname,'_markdata','IgnoreCase',true);
        
        if any(ifile & imark)
            ifile = find(ifile & imark);
            checkoverwrite = 1;
            tmpsavfile = tmpfname{ifile};
            
            temp = load(fullfile(fpath,tmpfname{ifile}));
            BlockData = temp.BlockData;
        elseif any(ifile)
            ifile = find(ifile);
            tmpsavfile = tmpfname{ifile};
            tmpsavfile = strrep(tmpsavfile,'_data','_markdata');
            checkoverwrite = 0;
            
            temp = load(fullfile(fpath,tmpfname{ifile}));
            BlockData = temp.BlockData;
            
            %do some additional pre-processing
            for c = 1:size(BlockData.Data,1)
                for c2 = 1:size(BlockData.Data,2)
                    tmppos = [];
                    tmprotang = [];
                    tmpang = [];
                    for d = 1:length(BlockData.Data(c,c2).m)
                        tmppos = cat(3,tmppos,BlockData.Data(c,c2).m(d).pos);
                        tmprotang = cat(4,tmprotang,BlockData.Data(c,c2).m(d).rotang);
                        tmpang = cat(3,tmpang,BlockData.Data(c,c2).m(d).angs);
                    end
                    BlockData.pos{c,c2} = tmppos; %*scalefactor;
                    BlockData.ja{c,c2} = BlockData.Data(c,c2).ja;
                    BlockData.rotang{c,c2} = tmprotang;
                    BlockData.ang{c,c2} = tmpang;
                    if isfield(BlockData.Data(c,c2),'time')
                        BlockData.time{c,c2} = BlockData.Data(c,c2).time;
                    elseif isfield(BlockData.Data(c,c2),'TrackerTime')
                        BlockData.time{c,c2} = BlockData.Data(c,c2).TrackerTime;
                    elseif isfield(BlockData.Data(c,c2),'Time')
                        BlockData.time{c,c2} = BlockData.Data(c,c2).Time;
                    end
                end
            end
            
            %remove the bulky Data field now that we've extracted pos and angle info
            BlockData = rmfield(BlockData,'Data');
            
            %add velocity to dataset
            BlockData = addvelsubj(BlockData);
    
        else
            continue;
        end
        
        
        
        %mark the item
        scaletype = 'scalemeters';
        
        if isempty(BlockData.pos{b})
            BlockData.inds(b).Full = [NaN NaN];
            BlockData.inds(b).SubAction = [NaN NaN];
            
            if checkoverwrite == 1
                if input(sprintf('Overwrite %s existing marked data file? (Y = 1; N = 0): ',Subj{a}))
                    save(fullfile(fpath,tmpsavfile),'BlockData');
                else
                    fprintf('\n\nData from subject %s, item %s not saved!\n  All changes to this subject/item have been lost.\n',Subj{a},BlockData.Items{b});
                end
            else
                save(fullfile(fpath,tmpsavfile),'BlockData');
            end
            continue;
        end
                
        if ~isfield(BlockData,'inds') || (length(BlockData.inds) < b) || isempty(BlockData.inds(b)) || ~isfield(BlockData.inds(b),'Full') || isempty(BlockData.inds(b).Full) || all(isnan(BlockData.inds(b).Full))
            inds{1} = [];
        else
            inds{1} = BlockData.inds(b).Full;
        end
        if ~isfield(BlockData,'inds') || (length(BlockData.inds) < b) || isempty(BlockData.inds(b)) || ~isfield(BlockData.inds(b),'SubAction') || isempty(BlockData.inds(b).SubAction) || all(isnan(BlockData.inds(b).SubAction))
            inds{2} = [];
        else
            inds{2} = BlockData.inds(b).SubAction;
        end
        
        if isempty(inds{1}) && isemty(inds{2})
            inds = addmarks(BlockData.vel{b}(:,:,1:6),'Nchan',3,'vthresh',0.2,'vthreshMin',0.1);
        end
        
        inds = markdataGUI(BlockData.pos{b},'title',sprintf('Item %d: %s',b,BlockData.Items{b}),'ang',BlockData.ang{b},'mark',inds,scaletype);
                
        if isempty(inds{1})
            BlockData.inds(b).Full = [NaN NaN];
        else
            BlockData.inds(b).Full = inds{1};
        end
        if isempty(inds{2})
            BlockData.inds(b).SubAction = [NaN NaN];
        else
            BlockData.inds(b).SubAction = inds{2};
        end
                
        %save the marked data for this subject
        if checkoverwrite == 1
            if input(sprintf('Overwrite %s existing marked data file? (Y = 1; N = 0): ',Subj{a}))
                save(fullfile(fpath,tmpsavfile),'BlockData');
            else
                fprintf('\n\nData from subject %s, item %s not saved!\n  All changes to this subject/item have been lost.\n',Subj{a},BlockData.Items{b});
            end
        else
            save(fullfile(fpath,tmpsavfile),'BlockData');
        end
        
    end %end for all subjects
    
    loopcont = input('Continue to the next item? (1 = yes, 0 = no/quit): ');
    if loopcont == 0
        break;
    end
    
end %end for all items


%the data are saved DURING the loop above, so we don't need to add in any
%code to save the data again here.


