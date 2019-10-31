%pull in data from all trials for one subject in one block and compile

%clear all;
%close all;

function [pathname] = ABblockcompile(varargin)

savepath = [];

if nargin > 0
    pathname = varargin{1};
    
    if nargin > 1
        savepath = varargin{2};
    end
else
    % %allow the user to select files, order them, and return the data to the
    % %  workspace as a char array pathname and a cell filename with the list of
    % %  selected files.
    % global exten
    % exten = '*.xls';
    % inputGUI;
    % uiwait(guih.figure1);
    % clear exten guih;
    
    fprintf('\nSelect data folder to analyze.\n');
    [pathname] = uigetdir('*.*','Select data folder to analyze');
end

%%

curdir = pwd;

cd(pathname);
%filename = ls;
filename = dir;

%ensure the files are sorted alphabetically
xlsfiles={filename.name};
[~,idx]=sort(xlsfiles);
filename=filename(idx);

cd(curdir);

% if ~strcmp(pathname(end),'/')
%     pathname = [pathname '/'];
% end

Data = [];
Items = {};
Fdbk = [];
BlockInd = [];
BlockName = '';
NTrials = zeros(14,2);

iptp = 0;
irou = 0;

%loadtype = input('Select transmitter configuration - (1) to the left, (2) upside_down, (3) wide-range :');
%loadtype = 3;

%for every trial in the folder, decode it and load it
for a = 1:size(filename,1)
    
    if ~contains(filename(a).name,'.xls') && ~contains(filename(a).name,'.txt') && ~contains(filename(a).name,'.dat')
        continue;
    end
    
    %don't process the practice trials
    if contains(filename(a).name,'practice','IgnoreCase',true) || contains(filename(a).name,'_p_','IgnoreCase',true)
        continue;
    end
    
    %decode the data and load it
    [blockind,trialID,item,itemcode,BlockName,SubjID,Group] = decodetrials(pathname,filename(a).name);
    
    %assume blockname, SubjID, and Group are the same for every trial in this folder.
%     
%     if isempty(SubjID)
%         fid = fopen(fullfile(pathname,filename(a).name),'r');
%         
%         temp = fgetl(fid);
%         cols = textscan(temp,'%s');
%         cols = cols{1};
%         
%         if length(cols) <= 2
%             hvals = textscan(temp,'%s %s');
%             if ~isempty(str2num(char(hvals{2}))) && contains(hvals{1},'sub','IgnoreCase',true)
%                 SubjID = hvals{2};
%             else
%                 doloop = 1;
%                 while (doloop)
%                     temp = fgetl(fid);
%                     if strcmp(temp,'--')
%                         %detect the end of the header
%                         doloop = 0;
%                     else
%                         hvals = textscan(temp,'%s %s');
%                         if ~isempty(str2num(char(hvals{2}))) && contains(hvals{1},'sub','IgnoreCase',true)
%                             SubjID = hvals{2};
%                             doloop = 0;
%                         end
%                     end
%                 end
%             end
%         end
%         if iscell(SubjID)
%             SubjID = SubjID{1};
%         end
%         
%         fclose(fid);
%     end

    if contains(filename(a).name,'multiple') || contains(filename(a).name,'sub41') || contains(filename(a).name,'sub116') || contains(filename(a).name,'sub345')  %we make the strong assumption that if the filename contains "multiple" it was an old recording with the upsidedown transmitter configuration
        trackerconfigtype = 2;
    else %we assume that otherwise it was recorded using the wide-range transmitter setup
        trackerconfigtype = 1;
    end
    
    
    if isempty(Data)
        %Data = ABLoad([pathname filename(a).name]);
        if trackerconfigtype == 2
            Data = ABLoad_upsidedowntransmitter(fullfile(pathname,filename(a).name));
        else
            Data = ABLoad(fullfile(pathname,filename(a).name));
        end
        if isempty(Data)
            continue;
        end
        if isempty(itemcode)
            itemCode = 1;
        elseif ~isempty(itemcode) || itemcode < 1
            itemCode = itemcode;
        else
            itemCode = 1;
        end
        if isnan(trialID) && blockind == 12
            iptp = iptp+1;
            trialID = iptp;
            item = decodetrials(blockind,trialID);
        elseif isnan(trialID) && blockind == 13
            irou = irou+1;
            trialID = irou;
            item = decodetrials(blockind,trialID);
        end
        
        if trialID > NTrials(blockind)
            NTrials(blockind,itemCode) = trialID;
        end
        
        Items{blockind,trialID,itemCode} = item;
        %BlockInd(blockind,trialID) = blockind;
        
        if (itemCode) == 2
            Fdbk(blockind,trialID,itemCode) = 0;
        else
            Fdbk(blockind,trialID,itemCode) = 1;
        end
        
        if trialID ~= 1 || itemCode ~= 1 || blockind ~= 1
            Data(blockind,trialID,itemCode) = Data;
        end
        
    else
        if ~isempty(itemcode)
            itemCode = itemcode;
        else
            itemCode = 1;
        end
        %BlockInd(trialID) = blockind;
        
        %Data(trialID) = ABLoad([pathname filename(a).name]);
        if trackerconfigtype == 2
            dat = ABLoad_upsidedowntransmitter(fullfile(pathname,filename(a).name));
        else
            dat = ABLoad(fullfile(pathname,filename(a).name));
        end
        
        if isempty(dat)
            continue;
        end
        
        if isnan(trialID) && blockind == 12
            iptp = iptp+1;
            trialID = iptp;
            item = decodetrials(blockind,trialID);
        elseif isnan(trialID) && blockind == 13
            irou = irou+1;
            trialID = irou;
            item = decodetrials(blockind,trialID);
        end
        
        if trialID > NTrials(blockind,itemCode)
            NTrials(blockind,itemCode) = trialID;
        end
        
        Items{blockind,trialID,itemCode} = item;
        
        if (itemCode) == 2
            Fdbk(blockind,trialID,itemCode) = 0;
        else
            Fdbk(blockind,trialID,itemCode) = 1;
        end
        
        Data(blockind,trialID,itemCode) = dat;
        
%         if loadtype == 3
%             Data(trialID) = ABLoad([pathname filename(a).name]);
%         elseif loadtype == 2
%             Data(trialID) = ABLoad_upsidedowntransmitter([pathname filename(a).name]);
%         else
%             disp('Not able to process requested transmitter configuration.');
%             return;
%         end
    end
    
%     Data.SubjID = SubjID;
%     Data.Group = Group;
    
end



%%

%rotate data into a fixed reference frame. to do that, we need to compute
%  two reference planes. The first should be along the frontal plane of the
%  body, which will contain the two shoulder markers. The orthogonal plane
%  should be the transverse plane (i.e., parallel to the surface of the
%  table). We can either assume that this plane is in line with the
%  transmitter's x-y plane, or we can estimate it by assuming that at the
%  time of trial initiation, the markers on the hand (thumb, index, middle,
%  and wrist) are all in the same plane, i.e. that they are all resting
%  relatively flat on the table surface. if that assumption may not
%  necesarily be true, we could end up with a very skewed transverse plane.
%  thus for now (and to keep the analysis simple) we will assume that the
%  transmitter is aligned with the table surface, which should be
%  approximately true.

%If we assume the transverse plane doesn't need to be rotated, we can just
%  focus on the frontal plane. To deal with this rotation, we will just
%  rotate the vector containing the two shoulders to align it with the x
%  axis. Because the shoulders continuously move in space, we will need to
%  compute a separate rotation for each time point and rotate all the
%  markers at that time point by that rotation angle.

fprintf('Data loaded. Processing...');

for a = 1:size(Data,1)
    for b = 1:size(Data,2)
        for z = 1:size(Data,3)
            
            if isempty(Items(a,b,z)) || isempty(Data(a,b,z).m)
                continue;
            end
            
            %calculate and remove the origin, which we will set at the left shoulder
            origin = [Data(a,b,z).m(7).x, Data(a,b,z).m(7).y, Data(a,b,z).m(7).z];
            
            for c = 1:length(Data(b).m)
                Data(a,b,z).m(c).x = Data(a,b,z).m(c).x-origin(:,1);
                Data(a,b,z).m(c).y = Data(a,b,z).m(c).y-origin(:,2);
                Data(a,b,z).m(c).z = Data(a,b,z).m(c).z-origin(:,3);
            end
            
            for d = 1:size(Data(a,b,z).m(1).x,1)
                
                %calculate the angle between the shoulders and the x axis; this
                %will tell us how far to rotate the shoulders to align the frontal
                %plane with the x-z plane.
                v1 = [Data(a,b,z).m(8).x(d)-Data(a,b,z).m(7).x(d), Data(a,b,z).m(8).y(d)-Data(a,b,z).m(7).y(d), Data(a,b,z).m(8).z(d)-Data(a,b,z).m(7).z(d)];
                theta = -acos(v1(1) / sqrt(sum(v1(1:2).^2,2)));
                thetas(d) = theta;
                
                %calculate the rotation matrix
                RotMat = [cos(theta) -sin(theta) 0;
                    sin(theta)  cos(theta) 0;
                    0           0          1];
                
                %rotate all the markers and the rotation matrix for this sample
                for c = 1:length(Data(a,b,z).m)
                    vec = [Data(a,b,z).m(c).x(d) Data(a,b,z).m(c).y(d) Data(a,b,z).m(c).z(d)];
                    rotvec = RotMat * vec';
                    Data(a,b,z).m(c).rotx(d,1) = rotvec(1);
                    Data(a,b,z).m(c).roty(d,1) = rotvec(2);
                    Data(a,b,z).m(c).rotz(d,1) = rotvec(3);
                    
                    %rotate the rotation matrix
                    Data(a,b,z).m(c).rotrotang(:,:,d) = Data(a,b,z).m(c).rotang(:,:,d)*RotMat;
                    
                    %recalculate the Euler angles from the rotated rotation matrix
                    Data(a,b,z).m(c).azim(d) = atan2d(Data(a,b,z).m(c).rotang(2,1,d),Data(a,b,z).m(c).rotang(1,1,d));
                    Data(a,b,z).m(c).elev(d) = -asind(Data(a,b,z).m(c).rotang(3,1,d));
                    Data(a,b,z).m(c).roll(d) = atan2d(Data(a,b,z).m(c).rotang(3,2,d),Data(a,b,z).m(c).rotang(3,3,d));
                    
                end
                
            end
        end
    end
end


%%

%compute joint angles from position data
for a = 1:size(Data,1)
    for b = 1:size(Data,2)
        for z = 1:size(Data,3)
            
            if isempty(Items(a,b,z)) || isempty(Data(a,b,z).m)
                continue;
            end
            
            tmp = Data(a,b,z).m;
            
            for c = 1:length(tmp)
                
                %we already set the origin above to be the left shoulder, so we'll
                %keep that, and create the position and angle matrix.
                Data(a,b,z).m(c).pos = [Data(a,b,z).m(c).x Data(a,b,z).m(c).y Data(a,b,z).m(c).z];
                Data(a,b,z).m(c).angs = [Data(a,b,z).m(c).azim Data(a,b,z).m(c).elev Data(a,b,z).m(c).roll];
                tmp(c).pos = Data(a,b,z).m(c).pos;
            end
            
            
            %position
            %1 = thumb
            %2 = index finger
            %3 = middle finger
            %4 = back of hand
            %5 = wrist
            %6 = elbow
            %7 = left shoulder
            %8 = right shoulder
            
            %calculate joint angles
            Data(a,b,z).ja(1,:) = acosd( sum( (tmp(5).pos-tmp(4).pos).*(tmp(1).pos-tmp(4).pos),2)./ (sqrt(sum((tmp(5).pos-tmp(4).pos).^2,2)).* sqrt(sum((tmp(1).pos-tmp(4).pos).^2,2)))); %thumb angle
            Data(a,b,z).ja(2,:) = acosd( sum( (tmp(5).pos-tmp(4).pos).*(tmp(2).pos-tmp(4).pos),2)./ (sqrt(sum((tmp(5).pos-tmp(4).pos).^2,2)).* sqrt(sum((tmp(2).pos-tmp(4).pos).^2,2)))); %index angle
            Data(a,b,z).ja(3,:) = acosd( sum( (tmp(5).pos-tmp(4).pos).*(tmp(3).pos-tmp(4).pos),2)./ (sqrt(sum((tmp(5).pos-tmp(4).pos).^2,2)).* sqrt(sum((tmp(3).pos-tmp(4).pos).^2,2)))); %middle angle
            Data(a,b,z).ja(4,:) = acosd( sum( (tmp(6).pos-tmp(5).pos).*(tmp(4).pos-tmp(5).pos),2)./ (sqrt(sum((tmp(6).pos-tmp(5).pos).^2,2)).* sqrt(sum((tmp(4).pos-tmp(5).pos).^2,2)))); %wrist angle
            Data(a,b,z).ja(5,:) = acosd( sum( (tmp(7).pos-tmp(6).pos).*(tmp(5).pos-tmp(6).pos),2)./ (sqrt(sum((tmp(7).pos-tmp(6).pos).^2,2)).* sqrt(sum((tmp(5).pos-tmp(6).pos).^2,2)))); %elbow angle
            Data(a,b,z).ja(6,:) = acosd( sum( (tmp(8).pos-tmp(7).pos).*(tmp(6).pos-tmp(7).pos),2)./ (sqrt(sum((tmp(8).pos-tmp(7).pos).^2,2)).* sqrt(sum((tmp(6).pos-tmp(7).pos).^2,2)))); %shoulder angle
            %joint angles:
            %1 = thumb angle
            %2 = index finger angle
            %3 = middle finger angle
            %4 = wrist angle
            %5 = elbow angle
            %6 = left shoulder angle w.r.t. a line between the two shoulders
            
            if all(isnan(Data(a,b,z).ja(1,:)))
                Data(a,b,z).ja(1,:)= acosd( sum( (tmp(6).pos-tmp(5).pos).*(tmp(1).pos-tmp(5).pos),2)./ (sqrt(sum((tmp(6).pos-tmp(5).pos).^2,2)).* sqrt(sum((tmp(1).pos-tmp(5).pos).^2,2)))); %approximated thumb angle
            end
            if all(isnan(Data(a,b,z).ja(2,:)))
                Data(a,b,z).ja(2,:)= acosd( sum( (tmp(6).pos-tmp(5).pos).*(tmp(2).pos-tmp(5).pos),2)./ (sqrt(sum((tmp(6).pos-tmp(5).pos).^2,2)).* sqrt(sum((tmp(2).pos-tmp(5).pos).^2,2)))); %approximated finger angle
            end
            
        end
    end
end

%eval(sprintf('%s = Data;',BlockName));  %give the data array the block name!
%clear Data;


%%

% 
% if length(blocktypes) > 1
% 	fprintf('\n\nFolder has multiple block types. They will be parsed.\n\n')
% end

for a = 1:size(Data,1)
    if sum(NTrials(a,:)) == 0
        continue;
    end
    
    iT1 = 0;
    iT2 = [];
    
    if NTrials(a,1) > 0
        iT1 = max(iT1,NTrials(a,1));
        iT2 = [iT2 1];
    end
    if NTrials(a,2) > 0
        iT1 = max(iT1,NTrials(a,2));
        iT2 = [iT2 2];
    end
    
    BlockData{a}.Data = Data(a,1:iT1,iT2);
    BlockData{a}.BlockID = a;
    BlockData{a}.BlockName = decodetrials(a);
    BlockData{a}.Items = Items(a,1:iT1,iT2);
    BlockData{a}.VisFdbk = Fdbk(a,1:iT1,iT2);
    BlockData{a}.SubjID = SubjID;
    BlockData{a}.Group = Group;
end

BD = BlockData;

%clear Data BlockName Items BlockInd;
    



%%


%save compiled data
clear a b c d curdir i1 i2 ind tmp* varargin

% tmpfname = strfind(pathname,'/');
% tmpfname = pathname(tmpfname(end-1)+1:tmpfname(end)-1);
% tmpfname = regexp(tmpfname,'\d*','Match');
% tmpfname = str2double(tmpfname{1});


% if exist([pathname 'S' tmpfname '_abdata.mat'],'file')
%     fappend = 1;
%     while exist([pathname 'S' tmpfname '_abdata' num2str(fappend) '.mat'],'file')
%         fappend = fappend+1;
%     end
%     
%     save([pathname 'S' tmpfname '_abdata' num2str(fappend) '.mat']);
%     fprintf('\n\nFile saved: %s\n\n',['S' tmpfname '_abdata' num2str(fappend) '.mat']);
% else
    if isempty(savepath)
        [savepath] = uigetdir(pathname,'Choose directory to save data into');
        if isempty(savepath)
            savepath = pathname;
        end
    elseif isnumeric(savepath) && savepath == 0
        savepath = fileparts(pathname); %we will go up one level, which should put us in the subject folder instead of the "Raw" folder
    end

    for a = 1:length(BD)
        if isempty(BD{a})
            continue;
        end
        
        BlockData = BD{a};
    
        save(fullfile(savepath, [SubjID '_' BlockData.BlockName '_data.mat']),'BlockData','filename','pathname','Group');  %BlockName 
        
        fprintf('\n\nFile saved: %s\n\n',[SubjID '_' BlockData.BlockName '_data.mat']);  %BlockName
    end
% end

