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

%parse the inputs and determine the desired outputs. we could make this
%more flexible using flags and a loop, but that makes it tricky to know
%what the desired output should be. so we will parse the inputs explicitly.
if nargin == 2 && (ischar(varargin{1}) && ischar(varargin{2})) && (contains(varargin{1},'/') || contains(varargin{1},'\'))
    %we have fpath and fname as input
    fpath = lower(varargin{1});
    fname = lower(varargin{2});
    
    %parse the file name
    isep = strfind(fname,'_');   %find all the separator indicies
    isep = [isep strfind(fname,'.')];
    isep = [isep strfind(fname,'-')];
    isep = sort(isep);
    
    %find the subject ID number
    isubjID = strfind(fname,'sub');
    if ~isempty(isubjID)
        isepsubjID = isep(find(isep-isubjID>0,1,'first'));
        subjID = fname(isubjID+3:isepsubjID-1);
    else
        subjID = '0';
    end
    
    %find the trial number
    itrial = strfind(fname,'trial');
    if ~isempty(itrial)
        iseptrial = isep(find(isep-itrial>0,1,'first'));
        trialID = str2double(fname(itrial+5:iseptrial-1));
    else
        item = lower(fname(isep(end-1)+1:isep(end)-1));
    end
    
    %find the group
    if contains(fpath,'control')
        group = 'Control';
    elseif contains(fpath,'patient')
        group = 'Patient';
    elseif contains(fpath,'model')
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
    

%condition on the available information to know what block we are in
if (contains(fpath,'meaningful') && contains(fpath,'unnamed')) || (contains(fname,'meaningful') && contains(fname,'unnamed')) || blockID == 4 || (contains(blockname,'ful','IgnoreCase',true) && contains(blockname,'unnamed','IgnoreCase',true))
    %meaningful unnamed
    blockID = 4;
    blockname = 'Mful-Unnamed';
   
elseif (contains(fpath,'meaningful') && contains(fpath,'named')) || (contains(fname,'meaningful') && contains(fname,'named')) || (contains(fname,'meaningful') && ~contains(fname,'unnamed') && ~contains(fname,'awk') && ~contains(fname,'static')) || blockID == 7 || (contains(blockname,'ful','IgnoreCase',true) && contains(blockname,'named','IgnoreCase',true)) || (contains(blockname,'ful','IgnoreCase',true) && ~contains(blockname,'awk','IgnoreCase',true) && ~contains(blockname,'static','IgnoreCase',true))
    %meaningful named (this is the default meaningful condition if it isn't labeled as something else)
    blockID = 7;
    blockname = 'Mful-Named';

elseif (contains(fpath,'meaningful') && contains(fpath,'awkward')) || (contains(fname,'meaningful') && contains(fname,'awkward')) || blockID == 5 || (contains(blockname,'ful','IgnoreCase',true) && contains(blockname,'awk','IgnoreCase',true))
    %meaningful awkward
    blockID = 5;
    blockname = 'Mful-Awk';
    
elseif (contains(fpath,'meaningful') && contains(fpath,'static')) || (contains(fname,'meaningful') && contains(fname,'static')) || blockID == 6 || (contains(blockname,'ful','IgnoreCase',true) && contains(blockname,'static','IgnoreCase',true))
    %meaningful static
    blockID = 6;
    blockname = 'Mful-Static';
    
elseif (contains(fpath,'meaningless') && ~contains(fpath,'awk') && ~contains(fpath,'static')) || (contains(fname,'meaningless')  && ~contains(fname,'awk') && ~contains(fname,'static')) || blockID == 1 || (contains(blockname,'less','IgnoreCase',true) &&  ~contains(blockname,'awk','IgnoreCase',true) && ~contains(blockname,'static','IgnoreCase',true))
    %meaningless unnamed  (this is the default meaningful condition if it isn't labeled as something else)
    blockID = 1;
    blockname = 'Mless';
    
elseif (contains(fpath,'meaningless') && contains(fpath,'awkward')) || (contains(fname,'meaningless') && contains(fname,'awkward')) || blockID == 2 || (contains(blockname,'less','IgnoreCase',true) && contains(blockname,'awk','IgnoreCase',true))
    %meaningless awkward
    blockID = 2;
    blockname = 'Mless-Awk';
    
elseif (contains(fpath,'meaningless') && contains(fpath,'static')) || (contains(fname,'meaningless') && contains(fname,'static')) || blockID == 3 || (contains(blockname,'less','IgnoreCase',true) && contains(blockname,'static','IgnoreCase',true))
    %meaningless static
    blockID = 3;
    blockname = 'Mless-Static';
    
elseif (contains(fpath,'real') && contains(fpath,'object')) || (contains(fname,'real') && contains(fname,'object')) || blockID == 8 || (contains(blockname,'real','IgnoreCase',true) && contains(blockname,'object','IgnoreCase',true))
    %real object use
    blockID = 8;
    blockname = 'Real-Object';
    
elseif contains(fpath,'point') || contains(fname,'point') || blockID == 9  || contains(blockname,'point','IgnoreCase',true)
    %point to point
    blockID = 9;
    blockname = 'Point-to-Point';
    trialID = 0;
    item = 'ptp';
    
end

%fill in the missing trial information
if trialID > 0
    %blockID
    %trialID
    item = trialcode(blockID,trialID);
else
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
        varargout{4} = blockname;
        varargout{5} = subjID;
        varargout{6} = group;
    
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

trials = cell(9,1);

%meaningless
trials{1} = {1, 'watch';
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
trials{4} = {1, 'comb';
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
trials{7} = trials{4};

%awkward (meaningful and meaningless are coded the same)
trials{2} = {1, 'comb';
             2, 'eraser';
             3, 'scissors';
             4, 'toothbrush';
             5, 'lipstick';
             6, 'stapler';
            };
trials{5} = trials{2};

%static (meaningful and meaningless are coded the same)
trials{3} = {1, 'book';
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
trials{6} = trials{3};

%real object use
trials{8} = {1, 'scissors';
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


        
if blockid == 9
    %if we are fed point-to-point, there's no "trial" information so we
    %handle this exception
    if isnumeric(input)
        %input the trial number, output the trial name
        varargout{1} = 'ptp';
        
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


