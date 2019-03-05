%pull in data from all trials for one subject and compile

%clear all;
%close all;

function [pathname] = ABblockcompile(varargin)

if nargin > 0
    pathname = varargin{1};
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
cd(curdir);

% if ~strcmp(pathname(end),'/')
%     pathname = [pathname '/'];
% end

Data = [];
Items = [];
BlockInd = [];
BlockName = '';

%loadtype = input('Select transmitter configuration - (1) to the left, (2) upside_down, (3) wide-range :');
%loadtype = 3;

%for every trial in the folder, decode it and load it
for a = 1:size(filename,1)
    
    if ~contains(filename(a).name,'.xls') && ~contains(filename(a).name,'.txt')
        continue;
    end
    
    %decode the data and load it
    [blockind,trialID,item,BlockName,SubjID,Group] = decodetrials(pathname,filename(a).name);
    %assume blockname, SubjID, and Group are the same for every trial in this folder.

    if contains(filename(a).name,'multiple')  %we make the strong assumption that if the filename contains "multiple" it was an old recording with the upsidedown transmitter configuration
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
        Items{trialID} = item;
        BlockInd(trialID) = blockind;
        
        if trialID ~= 1
            Data(trialID) = Data;
        end
        
    else
        BlockInd(trialID) = blockind;
        Items{trialID} = item;
        %Data(trialID) = ABLoad([pathname filename(a).name]);
        if trackerconfigtype == 2
            Data(trialID) = ABLoad_upsidedowntransmitter(fullfile(pathname,filename(a).name));
        else
            Data(trialID) = ABLoad(fullfile(pathname,filename(a).name));
        end
        
%         if loadtype == 3
%             Data(trialID) = ABLoad([pathname filename(a).name]);
%         elseif loadtype == 2
%             Data(trialID) = ABLoad_upsidedowntransmitter([pathname filename(a).name]);
%         else
%             disp('Not able to process requested transmitter configuration.');
%             return;
%         end
    end
    
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

for b = 1:length(Data)
    
    if isempty(Items(b))
        continue;
    end
    
    %calculate and remove the origin, which we will set at the left shoulder
    origin = [Data(b).m(7).x, Data(b).m(7).y, Data(b).m(7).z];
    
    for c = 1:length(Data(b).m)
        Data(b).m(c).x = Data(b).m(c).x-origin(:,1);
        Data(b).m(c).y = Data(b).m(c).y-origin(:,2);
        Data(b).m(c).z = Data(b).m(c).z-origin(:,3);
    end
    
    for d = 1:size(Data(b).m(1).x,1)
        
        %calculate the angle between the shoulders and the x axis; this
        %will tell us how far to rotate the shoulders to align the frontal
        %plane with the x-z plane.
        v1 = [Data(b).m(8).x(d)-Data(b).m(7).x(d), Data(b).m(8).y(d)-Data(b).m(7).y(d), Data(b).m(8).z(d)-Data(b).m(7).z(d)];
        theta = -acos(v1(1) / sqrt(sum(v1(1:2).^2,2)));
        thetas(d) = theta;
        
        %calculate the rotation matrix
        RotMat = [cos(theta) -sin(theta) 0;
                  sin(theta)  cos(theta) 0;
                  0           0          1];
        
        %rotate all the markers and the rotation matrix for this sample
        for c = 1:length(Data(b).m)
            vec = [Data(b).m(c).x(d) Data(b).m(c).y(d) Data(b).m(c).z(d)];
            rotvec = RotMat * vec';
            Data(b).m(c).rotx(d,1) = rotvec(1);
            Data(b).m(c).roty(d,1) = rotvec(2);
            Data(b).m(c).rotz(d,1) = rotvec(3);
            
            %rotate the rotation matrix
            Data(b).m(c).rotrotang(:,:,d) = Data(b).m(c).rotang(:,:,d)*RotMat;
            
            %recalculate the Euler angles from the rotated rotation matrix
            Data(b).m(c).azim(d) = atan2d(Data(b).m(c).rotang(2,1,d),Data(b).m(c).rotang(1,1,d));
            Data(b).m(c).elev(d) = -asind(Data(b).m(c).rotang(3,1,d));
            Data(b).m(c).roll(d) = atan2d(Data(b).m(c).rotang(3,2,d),Data(b).m(c).rotang(3,3,d));
            
        end
        
    end
end



%%

%compute joint angles from position data
for b = 1:length(Data)
    
    if isempty(Items(b))
        continue;
    end
    
    tmp = Data(b).m;
    
    for c = 1:length(tmp)
        
        %we already set the origin above to be the left shoulder, so we'll
        %keep that, and create the position and angle matrix.
        Data(b).m(c).pos = [Data(b).m(c).x Data(b).m(c).y Data(b).m(c).z];
        Data(b).m(c).angs = [Data(b).m(c).azim Data(b).m(c).elev Data(b).m(c).roll];
        tmp(c).pos = Data(b).m(c).pos;
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
    Data(b).ja(1,:) = acosd( sum( (tmp(5).pos-tmp(4).pos).*(tmp(1).pos-tmp(4).pos),2)./ (sqrt(sum((tmp(5).pos-tmp(4).pos).^2,2)).* sqrt(sum((tmp(1).pos-tmp(4).pos).^2,2)))); %thumb angle
    Data(b).ja(2,:) = acosd( sum( (tmp(5).pos-tmp(4).pos).*(tmp(2).pos-tmp(4).pos),2)./ (sqrt(sum((tmp(5).pos-tmp(4).pos).^2,2)).* sqrt(sum((tmp(2).pos-tmp(4).pos).^2,2)))); %index angle
    Data(b).ja(3,:) = acosd( sum( (tmp(5).pos-tmp(4).pos).*(tmp(3).pos-tmp(4).pos),2)./ (sqrt(sum((tmp(5).pos-tmp(4).pos).^2,2)).* sqrt(sum((tmp(3).pos-tmp(4).pos).^2,2)))); %middle angle
    Data(b).ja(4,:) = acosd( sum( (tmp(6).pos-tmp(5).pos).*(tmp(4).pos-tmp(5).pos),2)./ (sqrt(sum((tmp(6).pos-tmp(5).pos).^2,2)).* sqrt(sum((tmp(4).pos-tmp(5).pos).^2,2)))); %wrist angle
    Data(b).ja(5,:) = acosd( sum( (tmp(7).pos-tmp(6).pos).*(tmp(5).pos-tmp(6).pos),2)./ (sqrt(sum((tmp(7).pos-tmp(6).pos).^2,2)).* sqrt(sum((tmp(5).pos-tmp(6).pos).^2,2)))); %elbow angle
    Data(b).ja(6,:) = acosd( sum( (tmp(8).pos-tmp(7).pos).*(tmp(6).pos-tmp(7).pos),2)./ (sqrt(sum((tmp(8).pos-tmp(7).pos).^2,2)).* sqrt(sum((tmp(6).pos-tmp(7).pos).^2,2)))); %shoulder angle
    %joint angles:
    %1 = thumb angle
    %2 = index finger angle
    %3 = middle finger angle
    %4 = wrist angle
    %5 = elbow angle
    %6 = left shoulder angle w.r.t. a line between the two shoulders
    
    if all(isnan(Data(b).ja(1,:)))
        Data(b).ja(1,:)= acosd( sum( (tmp(6).pos-tmp(5).pos).*(tmp(1).pos-tmp(5).pos),2)./ (sqrt(sum((tmp(6).pos-tmp(5).pos).^2,2)).* sqrt(sum((tmp(1).pos-tmp(5).pos).^2,2)))); %approximated thumb angle
    end
    if all(isnan(Data(b).ja(2,:)))
        Data(b).ja(2,:)= acosd( sum( (tmp(6).pos-tmp(5).pos).*(tmp(2).pos-tmp(5).pos),2)./ (sqrt(sum((tmp(6).pos-tmp(5).pos).^2,2)).* sqrt(sum((tmp(2).pos-tmp(5).pos).^2,2)))); %approximated finger angle
    end
    
end

%eval(sprintf('%s = Data;',BlockName));  %give the data array the block name!
%clear Data;


%%

blocktypes = unique(BlockInd);

if length(blocktypes) > 1
	fprintf('\n\nFolder has multiple block types. They will be parsed.\n\n')
end
    
for a = 1:length(blocktypes)
    BlockData{blocktypes(a)}.Data = Data(BlockInd == blocktypes(a));
    BlockData{blocktypes(a)}.BlockID = blocktypes(a);
    BlockData{blocktypes(a)}.BlockName = decodetrials(blocktypes(a));
    BlockData{blocktypes(a)}.Items = Items(BlockInd == blocktypes(a));
    BlockData{blocktypes(a)}.SubjID = SubjID;
    BlockData{blocktypes(a)}.Group = Group;
end

clear blocktypes Data BlockName Items BlockInd;
    



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
    save([pathname 'S' SubjID '_data.mat'],'BlockData','filename','pathname','Group');  %BlockName 
    fprintf('\n\nFile saved: %s\n\n',['S' SubjID '_data.mat']);  %BlockName
% end
