%In this section, we pre-process each trajectory by automatically marking
%  the "full" action (from lift-off to return) based on a velocity
%  criterion, and asking the user to manually verify the marks. This will
%  allow us to perform a faster analysis on the kinematic data without
%  needing to identify the "core" portion of the action. We will do this on
%  a per-subject basis, processing the files coming directly out of
%  ABblockcompile.
%Refer to the instructions in the markdataGUI file for detailed information
%  regarding how to view and/or mark points in the data. Note that it is
%  necessary to close the markdataGUI window when you are done
%  viewing/marking data from one trial in order to move on to
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


ClearMarks = 0;


%% ***********************************************************************

%do some additional pre-processing

if ~contains(fname,'markdata')
    
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
    
end


%% ***********************************************************************

%the data file contains data from only one block

fprintf('\n\n');

% loopcont = 1;

    %set the appropriate scale factor so the data get converted to meters
    % --> ABLoad should already take care of this!
%     if blk == 10
%         scalefactor = 1; %note, this is to convert meters to inches (1/0.0254)
%         scaletype = 'scalemeters';
%     else
%         scalefactor = (0.0254/1);  %note, this is to convert inches to meters (0.0254/1)
%         scaletype = 'scaleinches';
%     end
    
    
scaletype = 'scalemeters';

    for c = 1:size(BlockData.pos,1) %process all the movements in this block
        for d = 1:size(BlockData.pos,2)
            
            if isempty(BlockData.pos{c,d})
                BlockData.inds(c,d).Full = [NaN NaN];
                continue;
            end
            
            
            if ~isfield(BlockData,'inds') || (size(BlockData.inds,1) < c || size(BlockData.inds,2) < d) || isempty(BlockData.inds(c,d)) || ~isfield(BlockData.inds(c,d),'Full') || isempty(BlockData.inds(c,d).Full) || all(isnan(BlockData.inds(c,d).Full))
                inds{1} = [];
            else
                inds{1} = BlockData.inds(c,d).Full;
            end
            
            inds{2} = [];
            
            if isempty(inds{1}) || ClearMarks == 1
                inds = addmarks(BlockData.vel{c,d},'Nchan',3,'throwmid');
            end
            
            if isempty(inds{1})
                inds{1} = [];
            end
            
            inds = markdataGUI(BlockData.pos{c,d},'title',sprintf('Item %d: %s, VF=%d',c,BlockData.Items{c,d},BlockData.VisFdbk(c,d)),'ang',BlockData.ang{c,d},'mark',inds,'full',scaletype);
            
            if isempty(inds{1})
                BlockData.inds(c,d).Full = [NaN NaN];
            else
                BlockData.inds(c,d).Full = inds{1};
            end
            
        end  %end for d
        
%         loopcont = input('Continue to the next item? (1 = yes, 0 = no/quit): ');
%         if loopcont == 0
%             break;
%         end
        
    end %end for c

% %     if loopcont == 0
% %         break;
% %     end
% 
% end  %end for blk


%%


%resave the data to ensure that your marks have been recorded.
tosave = input('Save file (note this may overwrite existing marked data files)? (Y = 1, N = 0): ');
if (tosave)
    
    save(fullfile(fpath,strrep(fname,'_data','_markdata')),'BlockData');
    
end




