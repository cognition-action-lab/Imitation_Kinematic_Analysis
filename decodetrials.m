%This function decodes block and trial information. You can provide any of
%several inputs and obtain the remaining information as outputs.
% Specifically, this function operates with:
%  INPUT:                       OUTPUT:
%    file path, file name         BlockID, trialID, itemName, blockName, subjID, Group
%    blockID, trialID             itemName, blockName
%    blockID, itemName            trialID, blockName
%    blockName, itemName          blockID, trialID
%    blockName, trialID           blockID, trialName
%    blockName                    blockID
%
%

function [varargout] = decodetrials(varargin)

fpath = '';
fname = '';
blockID = -1;
blockname = '';
trialID = -1;
item = '';
subjID = -1;
group = '';
itemcode = [];

%parse the inputs and determine the desired outputs. we could make this
%more flexible using flags and a loop, but that makes it tricky to know
%what the desired output should be. so we will parse the inputs explicitly.
if nargin == 2 && (ischar(varargin{1}) && ischar(varargin{2})) && contains(varargin{1},filesep)
    %we have fpath and fname as input
    fpath = varargin{1};
    fname = varargin{2};
    
    %parse the file name
    isep = strfind(fname,'_');   %find all the separator indicies
    isep = [isep strfind(fname,'.')];
    isep = [isep strfind(fname,'-')];
    isep = sort(isep);
    
    %find the subject ID number
    isubjID = strfind(lower(fname),'sub');
    if ~isempty(isubjID)
        isepsubjID = isep(find(isep-isubjID>0,1,'first'));
        subjID = fname(isubjID+3:isepsubjID-1);
    else
        isubjID = strfind(fpath,'AB');
        if isempty(isubjID)
            isubjID = strfind(fpath,'CS');
        end
        
        if isempty(isubjID)
            subjID = '0';
        else
            islash = strfind(fpath,filesep);
            isepsubjID = islash(find(islash-isubjID>0,1,'first'));
            if isempty(isepsubjID)
                %subjID = [fpath(isubjID:isubjID+1) num2str(str2double(fpath(isubjID+2:end)))];
                subjID = fpath(isubjID:end);
            else
                %subjID = [fpath(isubjID:isubjID+1) num2str(str2double(fpath(isubjID+2:isepsubjID-1)))];
                subjID = fpath(isubjID:isepsubjID-1);
            end
        end
    end
    
    
    
    %find the trial number
    itrial = strfind(fname,'trial');
    irepsep = strfind(fname,'-');
    if ~isempty(itrial)
        iseptrial = isep(find(isep-itrial>0,1,'first'));
        trialID = str2double(fname(itrial+5:iseptrial-1));
    else
        itrial = strfind(lower(fname),'item');
        
        if ~isempty(itrial)
            iseptrial = irepsep(find(irepsep-itrial>0,1,'first'));
            trialID = str2double(fname(itrial+4:iseptrial-1));
        else
            item = lower(fname(isep(end-1)+1:isep(end)-1));
        end
    end
    
    %find the group
    if contains(fpath,'control','IgnoreCase',true) || contains(fpath,'CB')
        group = 'Control';
    elseif contains(fpath,'patient','IgnoreCase',true) || contains(fpath,'AB')
        group = 'Patient';
    elseif contains(fpath,'model','IgnoreCase',true)
        group = 'Model';
    end
    
    outputtype = 1;
    
elseif nargin == 2 && (isnumeric(varargin{1}) && isnumeric(varargin{2}))
    %we have blockID and trialID as input, and want to know the trial name
    blockID = varargin{1};
    trialID = varargin{2};
    
    outputtype = 2;
    
elseif nargin == 2 && (isnumeric(varargin{1}) && ischar(varargin{2}))
    %we have blockID and item name as input, and want to know the trialID
    blockID = varargin{1};
    item = varargin{2};
    
    outputtype = 3;
    
elseif nargin == 2 && (ischar(varargin{1}) && ischar(varargin{2}))
    %we have blockname and item name as input, and want to know the blockID and trialID 
    blockname = varargin{1};
    fname = blockname;
    item = varargin{2};
    
    outputtype = 4;
    
elseif nargin == 2 && (ischar(varargin{1}) && isnumeric(varargin{2}))
    %we have blockname and trialID as input, and want to know the blockID and trial name 
    blockname = varargin{1};
    fname = blockname;
    trialID = varargin{2};
    
    outputtype = 5;
    
elseif nargin == 1 && isnumeric(varargin{1})
    %we have blockInd and want to know the blockName
    blockID = varargin{1};
    
    outputtype = 6;
elseif nargin == 1 && ischar(varargin{1})
    %we have blockName and want to know the blockID
    blockname = varargin{1};
    
    outputtype = 7;
    
end
    
fpath = lower(fpath);
fname = lower(fname);

fullpath = fullfile(fpath,fname);

%condition on the available information to know what block we are in
if (contains(fullpath,'mful') && contains(fullpath,'unnamed')) || blockID == 2 || (contains(blockname,'ful','IgnoreCase',true) && contains(blockname,'unnamed','IgnoreCase',true))
    %meaningful unnamed
    blockID = 2;
    blockname = 'Mful-Unnamed';
   
elseif (contains(fullpath,'mful') && contains(fullpath,'named')) || (contains(fullpath,'mful') && ~contains(fullpath,'unnamed') && ~contains(fullpath,'awk') && ~contains(fullpath,'static')) || blockID == 5 || (contains(blockname,'ful','IgnoreCase',true) && contains(blockname,'named','IgnoreCase',true) && ~contains(blockname,'unnamed','IgnoreCase',true))
    %meaningful named (this is the default meaningful condition if it isn't labeled as something else)
    blockID = 5;
    blockname = 'Mful-Named';

elseif (contains(fullpath,'mful') && contains(fullpath,'awk')) || blockID == 9 || (contains(blockname,'ful','IgnoreCase',true) && contains(blockname,'awk','IgnoreCase',true))
    %meaningful awkward
    blockID = 9;
    blockname = 'Mful-Awk';
    
elseif (contains(fullpath,'mful') && contains(fullpath,'static')) || blockID == 11 || (contains(blockname,'ful','IgnoreCase',true) && contains(blockname,'static','IgnoreCase',true))
    %meaningful static
    blockID = 11;
    blockname = 'Mful-Static';
    
elseif (contains(fullpath,'mless') && ~contains(fullpath,'awk') && ~contains(fullpath,'static')) || blockID == 3 || (contains(blockname,'less','IgnoreCase',true) &&  ~contains(blockname,'awk','IgnoreCase',true) && ~contains(blockname,'static','IgnoreCase',true))
    %meaningless unnamed  (this is the default meaningful condition if it isn't labeled as something else)
    blockID = 3;
    blockname = 'Mless';
    
elseif (contains(fullpath,'mless') && contains(fullpath,'awk')) || blockID == 7 || (contains(blockname,'less','IgnoreCase',true) && contains(blockname,'awk','IgnoreCase',true))
    %meaningless awkward
    blockID = 7;
    blockname = 'Mless-Awk';
    
elseif (contains(fullpath,'mless') && contains(fullpath,'static')) || blockID == 10 || (contains(blockname,'less','IgnoreCase',true) && contains(blockname,'static','IgnoreCase',true))
    %meaningless static
    blockID = 10;
    blockname = 'Mless-Static';
    
elseif contains(fullpath,'rou') || blockID == 13 || (contains(blockname,'real','IgnoreCase',true) && contains(blockname,'object','IgnoreCase',true))
    %real object use
    blockID = 13;
    blockname = 'Real-Object';
    
    if isempty(trialID) || isnan(trialID)
        trialID = NaN;
    end
    
elseif contains(fullpath,'point') || blockID == 12  || contains(blockname,'point','IgnoreCase',true)
    %point to point
    blockID = 12;
    blockname = 'Point-to-Point';
    if isempty(trialID) || isnan(trialID)
        trialID = NaN;
    end
    item = 'ptp';

elseif contains(fpath,'sight','IgnoreCase',true) || (contains(fname,'tbl') && contains(fname,'VF')) || blockID == 14 || contains(blockname,'gesture-to-Sight','IgnoreCase',true)
    blockID = 14;
    blockname = 'Gesture-to-Sight';
    if contains(fname,'NVF','IgnoreCase',true)
        itemcode = 2;
    else
        itemcode = 1;
    end
end

%fill in the missing trial information
if trialID > 0
    %blockID
    %trialID
    item = trialcode(blockID,trialID);
elseif ~isnan(trialID)
    %blockID
    %item
    trialID = trialcode(blockID,item);
end


%set the function output
switch(outputtype)
    case 1
        %we have fpath and fname as input, output everything
        
        varargout{1} = blockID;
        varargout{2} = trialID;
        varargout{3} = item;
        varargout{4} = itemcode;
        varargout{5} = blockname;
        varargout{6} = subjID;
        varargout{7} = group;
    
    case 2
        %we have blockID and trialID as input, and want to know the trial name and blockname
        varargout{1} = item;
        varargout{2} = blockname;
        
    case 3
        %we have blockID and item name as input, and want to know the trialID
        varargout{1} = trialID;
        varargout{2} = blockname;
        
    case 4
        %we have blockname and item name as input, and want to know the blockID and trialID 
        varargout{1} = blockID;
        varargout{2} = trialID;
        
    case 5
        %we have blockname and trialID input, and want to know the trial name
        varargout{1} = blockID;
        varargout{2} = item;
        
    case 6
        varargout{1} = blockname;
        
    case 7
        varargout{1} = blockID;
        
end

%in case you ask for the wrong number of inputs, make the rest of them
%empty
for a = length(varargout)+1:nargout
    varargout{a} = [];
end



end





%*************************************************************************%

%the helper function below takes as input either a trialID or a trialName,
%and return the opposite.

function varargout = trialcode(blockid,input)

trials = cell(14,1);

%meaningless
trials{3} = {1, 'watch';
             2, 'toothbrush';
             3, 'screwdriver';
             4, 'scissors';
             5, 'saw';
             6, 'razor';
             7, 'nailclippers';
             8, 'lighter';
             9, 'hammer';
             10, 'fork';
             11, 'eraser';
             12, 'comb';
             13, 'bottleopener';
             14, 'drinking'
            };

%meaningful (named/unnamed)
trials{2} = {1, 'comb';
             2, 'eraser';
             3, 'fork';
             4, 'lighter';
             5, 'nailclippers';
             6, 'razor';
             7, 'scissors';
             8, 'toothbrush';
             9, 'watch';
             10, 'bottleopener';
             11, 'drinking';
             12, 'hammer';
             13, 'stapler';
             14, 'lipstick'
            };
trials{5} = trials{2};

%awkward (meaningful and meaningless are coded the same)
trials{7} = {1, 'comb';
             2, 'eraser';
             3, 'scissors';
             4, 'toothbrush';
             5, 'lipstick';
             6, 'stapler';
            };
trials{9} = trials{7};

%static (meaningful and meaningless are coded the same)
trials{10} = {1, 'book';
             2, 'camcorder';
             3, 'compactmirror';
             4, 'cupwstraw';
             5, 'hairdryer';
             6, 'magnifyingglass';
             7, 'microphone';
             8, 'phone';
             9, 'umbrella';
             10, 'whistle';
            };
trials{11} = trials{10};

%point-to-point
trials{12} = {1, 'horizontal';
              2, 'vertical'
             };

%real object use
trials{13} = {1, 'scissors';
             2, 'watch';
             3, 'toothbrush';
             4, 'comb';
             5, 'fork';
             6, 'bottleopener';
             7, 'lighter';
             8, 'razor';
             9, 'eraser';
             10, 'nailclippers';
            };

%gesture-to-sight
trials{14} = {1, 'beermug'
              2, 'pingpongpaddle'
              3, 'teapot'
              4, 'teabag'
              5, 'comb'
              6, 'eraser'
              7, 'paintroller'
              8, 'knife'
              9, 'axe'
              10, 'bottleopener'
              11, 'hammer'
              12, 'iron'
              13, 'saw'
              14, 'drill'
              15, 'razor'
              16, 'fork'
              17, 'woodenspoon'
              18, 'cigarette'
              19, 'deoderant'
              20, 'match'
              21, 'toothbrush'
              22, 'tambourine'
              23, 'wrench'
              24, 'tweezers'
              25, 'screwdriver'
              26, 'nailclippers'
              27, 'perfume'
              28, 'camera'
              29, 'remote'
              30, 'stapler'
              31, 'key'
              32, 'lipstick'
              33, 'scissors'
              34, 'corkscrew'
              35, 'lighter'
              36, 'syringe'
              37, 'soapdispenser'
              38, 'calculator'
              39, 'drum'
              40, 'keyboard'
              41, 'watch'
             };

        
% if blockid == 12
%     %if we are fed point-to-point, there's no "trial" information so we
%     %handle this exception
%     if isnumeric(input)
%         %input the trial number, output the trial name
%         varargout{1} = 'ptp';
%         
%     elseif ischar(input)
%         %input the trial name, output the trial number. Since this could be
%         %anything we just return a null value.
%         varargout{1} = NaN;
%     end
%else
if isempty(trials(blockid))
    %if this is a practice block we will have no information to return
    if isnumeric(input)
        %input the trial number, output the trial name
        varargout{1} = '';
        
    elseif ischar(input)
        %input the trial name, output the trial number. Since this could be
        %anything we just return a null value.
        varargout{1} = NaN;
    end
    
else
    %for all other trial types, figure out what the input was and return
    %the correct output
    
    
    if isnumeric(input)
        %input the trial number, output the trial name
        trialID = input;
        
        tmp = cell2mat(trials{blockid}(:,1));
        if trialID <= length(tmp)
            trialname = trials{blockid}{tmp == trialID,2};
        else
            trialname = 'undefined';
        end
        
        varargout{1} = trialname;
        
    elseif ischar(input)
        %input the trial name, output the trial number
        trialname = input; %don't be case sensitive
        trialname = strrep(trialname,'_',''); %remove all tokens
        trialname = strrep(trialname,' ',''); %remove all tokens
        trialname = strrep(trialname,'-',''); %remove all tokens
        
        tmp = contains(trials{blockid}(:,2),trialname);
        trialID = find(tmp == 1);
        if isempty(trialID)
            trialID = NaN;
        end
        
        varargout{1} = trialID;
    end
end


end


