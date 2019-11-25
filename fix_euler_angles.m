%this code fixes the problem that arises because some of the
%data files only recorded euler angles, not the rotation matrix.
%to fix this problem we will transform the data into rotation matrices and
%then visualize and fix the rotation data on a sphere

clear all;
%close all;

[fname,fpath] = uigetfile('*.*','Select data file to analyze');

if contains(fullfile(fpath,fname),'.xls','IgnoreCase',true)
    
    
    fid = fopen([fpath fname],'r');
    
    temp = fgetl(fid);
    cols = textscan(temp,'%s');
    cols = cols{1};
    
    %didtime = 0;
    
    readstr = '';
    for a = 1:length(cols)
        
        %     if didtime == 1 && strcmp(lower(char(cols{a})),'time') == 1
        %         continue;
        %     end
        
        switch(lower(char(cols{a})))
            case {'sub','session','block_type','trial','sample'}
                readstr = [readstr '%d '];
            otherwise
                readstr = [readstr '%f '];
        end
    end
    
    dvals = textscan(fid,readstr,inf);
    
    fclose(fid);
    
    %dvals quits as soon as it runs into an error. That means that if there are
    %an unequal number of columns, the last line is problably where the error
    %arose. So we will trim the remaining columns to match as a quick
    %workaround rather than having to diagnose the problem in the data file.
    minlength = inf;
    for a = 1:size(dvals,2)
        if length(dvals{a}) < minlength
            minlength = length(dvals{a});
        end
    end
    for a = 1:size(dvals,2)
        dvals{a} = dvals{a}(1:minlength);
    end
    
    
    didtime = 0;
    
    for a = 1:length(cols)
        if didtime == 1 && strcmp(cols{a},'time') == 1
            continue;
        elseif didtime == 1 && strcmp(cols{a},'time') == 1
            didtime = 1;
        end
        
        if contains(cols{a},'m') && (contains(cols{a},'_x') || contains(cols{a},'_y') || contains(cols{a},'_z') || contains(cols{a},'_qual') || contains(cols{a},'_azim') || contains(cols{a},'_elev') || contains(cols{a},'_roll') || contains(cols{a},'_rotang'))
            tmp1 = strfind(cols{a},'m');
            tmp2 = strfind(cols{a},'_');
            eval(sprintf('data.m(%s).%s = dvals{%d};',char(cols{a}(tmp1+1:tmp2-1)),char(cols{a}(tmp2+1:end)),a));
        elseif strcmp(cols{a},'sub') || strcmp(cols{a},'session') || strcmp(cols{a},'block_type') || strcmp(cols{a},'trial')
            eval(sprintf('data.%s = unique(dvals{%d});',char(cols{a}),a));
        else
            eval(sprintf('data.%s = dvals{%d};',char(cols{a}),a));
        end
    end
    
    
    %get the rotation matrix and Euler angles
    for a = 1:8
        if ~isfield(data.m(a),'azim') && ~isfield(data.m(a),'rotang11')
            %no rotation data is available, fill in NaNs
            data.m(a).azim = NaN*ones(size(data.m(a).x));
            data.m(a).elev = NaN*ones(size(data.m(a).x));
            data.m(a).roll = NaN*ones(size(data.m(a).x));
            
            data.m(a).rotang = NaN*ones(3,3,size(data.m(a).x));
            
        elseif ~isfield(data.m(a),'rotang11') || all(isnan(data.m(a).rotang11))
            %the rotation matrix is not available (but presumably we have Euler Angles)
            %
            %compute the rotation matrix from the azim, elev, and roll (i.e.,
            %yaw/pitch/roll)
            yaw = data.m(a).azim;
            pitch = data.m(a).elev;
            roll = data.m(a).roll;
            
            data.m(a).rotang(1,1,:) = cosd(yaw).*cosd(pitch);
            data.m(a).rotang(1,2,:) = cosd(yaw).*sind(pitch).*sind(roll)-sind(yaw).*cosd(roll);
            data.m(a).rotang(1,3,:) = cosd(yaw).*sind(pitch).*cosd(roll)+sind(yaw).*sind(roll);
            data.m(a).rotang(2,1,:) = sind(yaw).*cosd(pitch);
            data.m(a).rotang(2,2,:) = sind(yaw).*sind(pitch).*sind(roll)+cosd(yaw).*cosd(roll);
            data.m(a).rotang(2,3,:) = sind(yaw).*sind(pitch).*cosd(roll)-cosd(yaw).*sind(roll);
            data.m(a).rotang(3,1,:) = -sind(pitch);
            data.m(a).rotang(3,2,:) = cosd(pitch).*sind(roll);
            data.m(a).rotang(3,3,:) = cosd(pitch).*cosd(roll);
            
        else
            %the rotation matrix is available!
            data.m(a).rotang = [];
            data.m(a).rotang(1,1,:) = data.m(a).rotang11;
            data.m(a).rotang(1,2,:) = data.m(a).rotang12;
            data.m(a).rotang(1,3,:) = data.m(a).rotang13;
            data.m(a).rotang(2,1,:) = data.m(a).rotang21;
            data.m(a).rotang(2,2,:) = data.m(a).rotang22;
            data.m(a).rotang(2,3,:) = data.m(a).rotang23;
            data.m(a).rotang(3,1,:) = data.m(a).rotang31;
            data.m(a).rotang(3,2,:) = data.m(a).rotang32;
            data.m(a).rotang(3,3,:) = data.m(a).rotang33;
            
            if ~isfield(data.m(a),'azim') || all(isnan(data.m(a).azim))
                %if the Euler angles are not available, compute them
                data.m(a).azim = atan2d(data.m(a).rotang21,data.m(a).rotang11);
                data.m(a).elev = -asind(data.m(a).rotang31);
                data.m(a).roll = atan2d(data.m(a).rotang32,data.m(a).rotang33);
            end
            
        end
    end
    
elseif contains(fullfile(fpath,fname),'.dat','IgnoreCase',true)
    
    
    fid = fopen([fpath fname],'r');
    
    temp = fgetl(fid);
    cols = textscan(temp,'%s');
    cols = cols{1};
    
    if length(cols) <= 2
        
        filetype = 'aaron';
        
        %this file has a header, we have to parse the header
        hvals = textscan(temp,'%s %s');
        if ~isempty(str2num(char(hvals{2}))) && ~contains(hvals{1},'sub','IgnoreCase',true)
            hvals = textscan(temp,'%s %d');
            tmp = strrep(char(hvals{1}),':','');
            eval(sprintf('header.%s = %d;',tmp,hvals{2}));
        else
            tmp = strrep(char(hvals{1}),':','');
            tmp1 = strrep(char(hvals{2}),'.txt','');
            eval(sprintf('header.%s = ''%s'';',tmp,char(tmp1)));
        end
        
        doloop = 1;
        while (doloop)
            temp = fgetl(fid);
            if strcmp(temp,'--')
                %detect the end of the header
                temp = fgetl(fid);
                header.cols = textscan(temp,'%s');
                header.cols = header.cols{1};
                fgetl(fid);
                doloop = 0;
            else
                hvals = textscan(temp,'%s %s');
                if ~isempty(str2num(char(hvals{2}))) && ~contains(hvals{1},'sub','IgnoreCase',true)
                    hvals = textscan(temp,'%s %d');
                    tmp = strrep(char(hvals{1}),':','');
                    eval(sprintf('header.%s = %d;',tmp,hvals{2}));
                else
                    tmp = strrep(char(hvals{1}),':','');
                    tmp1 = strrep(char(hvals{2}),'.txt','');
                    eval(sprintf('header.%s = ''%s'';',tmp,char(tmp1)));
                end
                
            end
        end
        
        % %unknown columns, use this code for complete flexiblity
        readstr = '';
        for a = 1:length(header.cols)
            switch(char(header.cols{a}))
                case {'Device_Num','FakeTime','Time'}
                    readstr = [readstr '%d '];
                case {'HandX','HandY','HandZ','HandAzim','HandElev','HandRoll','Theta'}
                    readstr = [readstr '%f '];
                case {'StartX','StartY','TargetX','TargetY'}
                    readstr = [readstr '%f '];
                case {'Trial', 'Redo','Stimulus'}
                    readstr = [readstr '%d '];
                case {'Keymap' 'Stim'}
                    readstr = [readstr '%s '];
                otherwise
                    readstr = [readstr '%f '];
            end
        end
        
        dvals = textscan(fid,readstr,inf);
        
        %compute number of birds
        numBirds = length(unique(dvals{1}));
        
        for a = 1:numBirds
            
            %divide data structure into appropriate fields, based on column names
            for b = 1:length(header.cols)
                
                if contains(header.cols{b},'hand','IgnoreCase',true) || contains(header.cols{b},'velocity','IgnoreCase',true) || contains(header.cols{b},'latency','IgnoreCase',true) || contains(header.cols{b},'duration','IgnoreCase',true)
                    fieldname = lower(header.cols{b});
                    fieldname = strrep(fieldname,'hand','');
                    fieldname = strrep(fieldname,'rotmat','rotang');
                    eval(sprintf('data.m(%d).%s = dvals{%d}(%d:%d:end);',a,fieldname,b,a,numBirds));
                elseif strcmpi(header.cols{b},'trial') || strcmpi(header.cols{b},'redo') || strcmpi(header.cols{b},'session') || strcmpi(header.cols{b},'sub') || strcmpi(header.cols{b},'device_num')
                    eval(sprintf('data.%s = unique(dvals{%d}(%d:%d:end));',char(header.cols{b}),b,a,numBirds));
                else
                    eval(sprintf('data.%s = dvals{%d}(%d:%d:end);',char(header.cols{b}),b,a,numBirds));
                end
                
                
                %tmp1 = strfind(cols{a},'m');
                %tmp2 = strfind(cols{a},'_');
                %eval(sprintf('data.m(%s).%s = dvals{%d};',char(cols{a}(tmp1+1:tmp2-1)),char(cols{a}(tmp2+1:end)),a));
                
            end
            
        end %end for loop
        
        if ~exist('data','var') || isempty(data) || isempty(data.m) || length(data.m) < 8 || length(data.m(1).x) <= header.Sampling_Rate
            data = [];
            
            if nargout > 2
                varargout{1} = fname;
            end
            if nargout > 3
                varargout{2} = fpath;
            end
            if nargout > 4
                varargout{3} = [];
            end
            
            return;
            
        end
        
        if isfield(header,'Sampling_Rate') && isfield(data.m(1),'x')
            data.time = [0:1:length(data.m(1).x)-1]'/header.Sampling_Rate;  %create time stamp in msec
        elseif isfield(data.m(1),'x')
            data.time = [0:1:length(data.m(1).x)-1]';
        else  %unrecognized trial table type
            data.time = [];
        end
        
        data.sub = header.SubjectID;
        
        fclose(fid);
        
        
    else
        
        filetype = 'steve';
        
        %didtime = 0;
        
        readstr = '';
        for a = 1:length(cols)
            
            %     if didtime == 1 && strcmp(lower(char(cols{a})),'time') == 1
            %         continue;
            %     end
            
            switch(lower(char(cols{a})))
                case {'sub','session','block_type','trial','sample'}
                    readstr = [readstr '%d '];
                otherwise
                    readstr = [readstr '%f '];
            end
        end
        
        dvals = textscan(fid,readstr,inf);
        
        fclose(fid);
        
        %dvals quits as soon as it runs into an error. That means that if there are
        %an unequal number of columns, the last line is problably where the error
        %arose. So we will trim the remaining columns to match as a quick
        %workaround rather than having to diagnose the problem in the data file.
        minlength = inf;
        for a = 1:size(dvals,2)
            if length(dvals{a}) < minlength
                minlength = length(dvals{a});
            end
        end
        for a = 1:size(dvals,2)
            dvals{a} = dvals{a}(1:minlength);
        end
        
        
        didtime = 0;
        
        for a = 1:length(cols)
            if didtime == 1 && strcmp(cols{a},'time') == 1
                continue;
            elseif didtime == 1 && strcmp(cols{a},'time') == 1
                didtime = 1;
            end
            
            if contains(cols{a},'m') && (contains(cols{a},'_x') || contains(cols{a},'_y') || contains(cols{a},'_z') || contains(cols{a},'_qual') || contains(cols{a},'_azim') || contains(cols{a},'_elev') || contains(cols{a},'_roll') || contains(cols{a},'_rotang'))
                tmp1 = strfind(cols{a},'m');
                tmp2 = strfind(cols{a},'_');
                eval(sprintf('data.m(%s).%s = dvals{%d};',char(cols{a}(tmp1+1:tmp2-1)),char(cols{a}(tmp2+1:end)),a));
            elseif strcmp(cols{a},'sub') || strcmp(cols{a},'session') || strcmp(cols{a},'block_type') || strcmp(cols{a},'trial')
                eval(sprintf('data.%s = unique(dvals{%d});',char(cols{a}),a));
            else
                eval(sprintf('data.%s = dvals{%d};',char(cols{a}),a));
            end
        end
        
        if ~exist('data','var') || isempty(data) || isempty(data.m) || length(data.m) < 8 || length(data.m(1).x) <= 140
            data = [];
            
            if nargout > 2
                varargout{1} = fname;
            end
            if nargout > 3
                varargout{2} = fpath;
            end
            if nargout > 4
                varargout{3} = [];
            end
            
            return;
            
        end
        
    end
    
    
    %get the rotation matrix and Euler angles
    for a = 1:8
        if ~isfield(data.m(a),'azim') && ~isfield(data.m(a),'rotang11')
            %no rotation data is available, fill in NaNs
            data.m(a).azim = NaN*ones(size(data.m(a).x));
            data.m(a).elev = NaN*ones(size(data.m(a).x));
            data.m(a).roll = NaN*ones(size(data.m(a).x));
            
            data.m(a).rotang = NaN*ones(3,3,size(data.m(a).x));
            
        elseif ~isfield(data.m(a),'rotang11') || all(isnan(data.m(a).rotang11))
            %the rotation matrix is not available (but presumably we have Euler Angles)
            %
            %compute the rotation matrix from the azim, elev, and roll (i.e.,
            %yaw/pitch/roll)
            yaw = data.m(a).azim;
            pitch = data.m(a).elev;
            roll = data.m(a).roll;
            
            %***We need to fix the rotation angles to address the gimbal lock
            %problem. Once we do this, the resulting rotation matrices should
            %be good going forward and we will live in the rotation matrix
            %space.
            
            
            data.m(a).rotang(1,1,:) = cosd(yaw).*cosd(pitch);
            data.m(a).rotang(1,2,:) = cosd(yaw).*sind(pitch).*sind(roll)-sind(yaw).*cosd(roll);
            data.m(a).rotang(1,3,:) = cosd(yaw).*sind(pitch).*cosd(roll)+sind(yaw).*sind(roll);
            data.m(a).rotang(2,1,:) = sind(yaw).*cosd(pitch);
            data.m(a).rotang(2,2,:) = sind(yaw).*sind(pitch).*sind(roll)+cosd(yaw).*cosd(roll);
            data.m(a).rotang(2,3,:) = sind(yaw).*sind(pitch).*cosd(roll)-cosd(yaw).*sind(roll);
            data.m(a).rotang(3,1,:) = -sind(pitch);
            data.m(a).rotang(3,2,:) = cosd(pitch).*sind(roll);
            data.m(a).rotang(3,3,:) = cosd(pitch).*cosd(roll);
            
        else
            %the rotation matrix is available!
            data.m(a).rotang = [];
            data.m(a).rotang(1,1,:) = data.m(a).rotang11;
            data.m(a).rotang(1,2,:) = data.m(a).rotang12;
            data.m(a).rotang(1,3,:) = data.m(a).rotang13;
            data.m(a).rotang(2,1,:) = data.m(a).rotang21;
            data.m(a).rotang(2,2,:) = data.m(a).rotang22;
            data.m(a).rotang(2,3,:) = data.m(a).rotang23;
            data.m(a).rotang(3,1,:) = data.m(a).rotang31;
            data.m(a).rotang(3,2,:) = data.m(a).rotang32;
            data.m(a).rotang(3,3,:) = data.m(a).rotang33;
            
            if ~isfield(data.m(a),'azim') || all(isnan(data.m(a).azim))
                %if the Euler angles are not available, compute them
                data.m(a).azim = atan2d(data.m(a).rotang21,data.m(a).rotang11);
                data.m(a).elev = -asind(data.m(a).rotang31);
                data.m(a).roll = atan2d(data.m(a).rotang32,data.m(a).rotang33);
            end
            
        end
    end
    
end

EulerAng_vis(data.m(2).rotang)



%we want to reorient the axes to a standard configuration, where x points
%to the participant's right, y is straight ahead to the participant, and z
%points up.
%this is equivalent to applying a rotation matrix:
% RotMat = [ 0  1  0;
%            1  0  0;
%            0  0 -1;
%          ]; 
RotMat = [ 0 1 0;
          -1 0 0;
           0 0 1
         ];


%we will manually apply this below:

for a = 1:length(data.m)
    tmp1 = data.m(a).x; %x axis originally points to the patient's backward
    tmp2 = data.m(a).y; %y axis originally points to the patient's right
    tmp3 = data.m(a).z; %z axis originally points to the patient's up
    
    data.m(a).x = tmp2;  %x axis now points to the patient's right
    data.m(a).y = -tmp1;  %y axis now points to the patient's forward
    data.m(a).z = tmp3; %z axis now points to the patient's upward
    
    %we also need to apply this transformation to the rotation matrices
    for b = 1:size(data.m(a).rotang,3)
        data.m(a).rotang(:,:,b) = data.m(a).rotang(:,:,b)*RotMat;
    end
    
    
end




%we will assume data are in inches, so we need to convert units to m
for a = 1:length(data.m)
    data.m(a).x = (data.m(a).x)*0.0254;
    data.m(a).y = (data.m(a).y)*0.0254;
    data.m(a).z = (data.m(a).z)*0.0254;
end


%we will assume data are in inches, so we need to convert units to m
for a = 1:length(data.m)
    data.m(a).x = (data.m(a).x)*10;
    data.m(a).y = (data.m(a).y)*10;
    data.m(a).z = (data.m(a).z)*10;
end




