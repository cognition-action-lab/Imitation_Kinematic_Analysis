%this code does a quick comparison between two data sets.


clear all;
close all;

fprintf('\n\nSelect first data file to compare.\n');
[fname{1},fpath{1}] = uigetfile('*.*','Select first data file.');

fprintf('\nSelect second data file to compare.\n');
[fname{2},fpath{2}] = uigetfile('*.*','Select second data file.');

%%

Data = [];

%load in the raw data and transform to standard room coordinates
for ifile = 1:2
    
    clear data;
    
    fid = fopen(fullfile(fpath{ifile},fname{ifile}),'r');
    
    %the file could be the model, which used the upside-down transmitter
    %configuration (except in a couple special cases). So we have to account for that.
    if (contains(fullfile(fpath{ifile},fname{ifile}),'model','IgnoreCase',true) || contains(fullfile(fpath{ifile},fname{ifile}),'.xls')) && (~contains(fullfile(fpath{ifile},fname{ifile}),'hammer') || ~contains(fullfile(fpath{ifile},fname{ifile}),'drink'))
        
        
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
                data.converteuler = 1;
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
                data.converteuler = 0;
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
        
        %we want to reorient the axes to a standard configuration, where x points
        %to the participant's right, y is straight ahead to the participant, and z
        %points up.
        %this is equivalent to applying a rotation matrix:
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
                data.m(a).rotang(:,:,b) = RotMat*data.m(a).rotang(:,:,b);
            end
            
        end
        
        %we will assume data are in inches, so we need to convert units to m
        for a = 1:length(data.m)
            data.m(a).x = (data.m(a).x)*0.0254;
            data.m(a).y = (data.m(a).y)*0.0254;
            data.m(a).z = (data.m(a).z)*0.0254;
        end
        
        %there is a potential problem of data crossing out into the wrong
        %hemisphere. to avoid this, we will detect sudden discontinuities in two of
        %the three axes of the data, and when we find one we will systematically flip
        %the sign until the next discontinuity. these discontinuities should be
        %paired.
        for a = 1:length(data.m)
            velx = diff(data.m(a).x);
            vely = diff(data.m(a).y);
            velz = diff(data.m(a).z);
            
            accx = diff(velx);
            accy = diff(vely);
            accz = diff(velz);
            
            indx = find(abs(velx(1:end)) > 10);
            indy = find(abs(vely(1:end)) > 10);
            indz = find(abs(velz(1:end)) > 10);
            
            %find and throw out all the single-sample spikes
            indax = find(abs(accx(1:end)) > 10);
            inday = find(abs(accy(1:end)) > 10);
            indaz = find(abs(accz(1:end)) > 10);
            for b = 1:length(indax)-1
                if indax(b+1)-indax(b) < 2
                    indx(indx == indax(b+1)) = [];
                end
            end
            for b = 1:length(inday)-1
                if inday(b+1)-inday(b) < 2
                    indy(indy == inday(b+1)) = [];
                end
            end
            for b = 1:length(indaz)-1
                if indaz(b+1)-indaz(b) < 2
                    indz(indz == indaz(b+1)) = [];
                end
            end
            
            if isempty(indx)
                %y and z
                inds = intersect(indy,indz);
            elseif isempty(indy)
                %x and z
                inds = intersect(indx,indz);
            elseif isempty(indz)
                %x and y
                inds = intersect(indx,indy);
            else
                inds = [];
            end
            
            if isempty(inds)
                continue;
            end
            
            if size(inds,1) > size(inds,2) %make sure this is a row vector
                inds = inds';
            end
            
            %comment this line out to avoid automatically fixing odd-mark errors by
            %assuming the last point is missing - but sometimes it is the first point!
            %     if mod(length(inds),2) ~= 0
            %         inds(end+1) = length(data.m(a).x);
            %     end
            
            %     %comment these lines out to automatically process data without manual
            %     %verification
            
            %fprintf('   Subj: %s\tTrial: %d',data.sub,
            fpn = [fpath fname];
            tmpinds = strfind(fpn,filesep);
            fprintf('   %s\n',fpn(tmpinds(end)+1:end));
            
            figure(2)
            clf;
            torso = [(data.m(7).x(1)+data.m(8).x(1))/2 (data.m(7).y(1)+data.m(8).y(1))/2 (data.m(7).z(1)+data.m(8).z(1))/2];
            joint.refshoulder = [data.m(8).x(1) data.m(8).y(1) data.m(8).z(1)];
            joint.shoulder = [data.m(7).x(1) data.m(7).y(1) data.m(7).z(1)];
            joint.elbow = [data.m(6).x(1) data.m(6).y(1) data.m(6).z(1)];
            joint.wrist = [data.m(5).x(1) data.m(5).y(1) data.m(5).z(1)];
            joint.hand = [data.m(4).x(1) data.m(4).y(1) data.m(4).z(1)];
            joint.middlefinger = [data.m(3).x(1) data.m(3).y(1) data.m(3).z(1)];
            joint.indexfinger = [data.m(2).x(1) data.m(2).y(1) data.m(2).z(1)];
            joint.thumb = [data.m(1).x(1) data.m(1).y(1) data.m(1).z(1)];
            if contains(lower(fpath),'6markervers')
                %do simplified arm without wrist/hand
                [xSE,ySE,zSE] = cylinder2P(.45,8,joint.shoulder(1,:),joint.elbow(1,:));
                [xEW,yEW,zEW] = cylinder2P(.375,8,joint.elbow(1,:),joint.wrist(1,:));
                [xWF,yWF,zWF] = cylinder2P([.25 .25 .25 .1 .1 .1],8,joint.wrist(1,:),joint.indexfinger(1,:));
                [xWT,yWT,zWT] = cylinder2P([.25 .25 .25 .1 .1 .1],8,joint.wrist(1,:),joint.thumb(1,:));
            else%if contains(lower(fpath),'8markervers')
                [xSE,ySE,zSE] = cylinder2P(.45,8,joint.shoulder(1,:),joint.elbow(1,:));
                [xEW,yEW,zEW] = cylinder2P(.375,8,joint.elbow(1,:),joint.wrist(1,:));
                [xWH,yWH,zWH] = cylinder2P(.375,8,joint.wrist(1,:),joint.hand(1,:));
                [xHM,yHM,zHM] = cylinder2P([.25 .25 .25 .1 .1 .1],8,joint.hand(1,:),joint.middlefinger(1,:));
                [xWF,yWF,zWF] = cylinder2P([.25 .25 .25 .1 .1 .1],8,joint.hand(1,:),joint.indexfinger(1,:));
                [xWT,yWT,zWT] = cylinder2P([.25 .25 .25 .1 .1 .1],8,joint.hand(1,:),joint.thumb(1,:));
            end
            %rectangle representing the body, assuming a vertical plane
            torsolength = abs(joint.shoulder(1,1)-joint.refshoulder(1,1))*1.5;
            xtorso = [joint.shoulder(1,1); joint.refshoulder(1,1); joint.refshoulder(1,1); joint.shoulder(1,1)];
            ytorso = [joint.shoulder(1,2); joint.refshoulder(1,2); joint.refshoulder(1,2); joint.shoulder(1,2)];
            ztorso = [joint.shoulder(1,3); joint.refshoulder(1,3); joint.refshoulder(1,3)-torsolength; joint.shoulder(1,3)-torsolength];
            headrad = abs(joint.refshoulder(1,1)-joint.shoulder(1,1))/4;
            headc = torso + [0 0 headrad];
            theta = [0:pi/10:2*pi];
            headx = headrad*cos(theta)+headc(1);
            heady = zeros(size(headx))+headc(2);
            headz = headrad*sin(theta)+headc(3);
            
            [spherex,spherey,spherez] = sphere(10);
            
            h = patch('XData',xtorso,'YData',ytorso,'ZData',ztorso);
            set(h,'FaceColor',[.5 .5 .5]);
            hold on;
            h = patch('XData',headx,'YData',heady,'ZData',headz);
            set(h,'FaceColor',[.5 .5 .5]);
            h = surf(xSE, ySE, zSE);
            set(h,'FaceColor',[1 0 0]);
            h = surf(xEW, yEW, zEW);
            set(h,'FaceColor',[1 1 0]);
            if contains(lower(fpath),'8markervers')
                h = surf(xWH, yWH, zWH);
                set(h,'FaceColor',[1 1 0]);
                h = surf(xHM, yHM, zHM);
                set(h,'FaceColor',[0 0 1]);
            end
            h = surf(xWF, yWF, zWF);
            set(h,'FaceColor',[0 0 1]);
            h = surf(xWT, yWT, zWT);
            set(h,'FaceColor',[0 0 1]);
            h = surf(spherex+data.m(a).x(1),spherey+data.m(a).y(1),spherez+data.m(a).z(1));
            set(h,'FaceColor',[.8 0 0]);
            hold off;
            grid on;
            view(195,20);
            xlabel('x');
            ylabel('y');
            zlabel('z');
            
            figure(1)
            inds = markdata3d(data.m(a).x,data.m(a).y,data.m(a).z,length(data.m(a).x),inds,0,1,'Name',sprintf('Marker %d',a));
            
            inds = sort(inds);
            if inds(end) > length(data.m(a).x)
                inds(end) = length(data.m(a).x);
            end
            if inds(1) < 1
                inds(1) = 1;
            end
            
            while mod(length(inds),2) ~= 0  %prompt until we have an even number of marks
                disp('Odd number of marks.');
                inds = markdata3d(data.m(a).x,data.m(a).y,data.m(a).z,length(data.m(a).x),inds,0,1,'Name',sprintf('Marker %d',a));
                
                if inds(end) > length(data.m(a).x)
                    inds(end) = length(data.m(a).x);
                end
                if inds(1) < 1
                    inds(1) = 1;
                end
            end
            
            inds = sort(inds);
            
            for b = 1:2:length(inds)
                data.m(a).x(inds(b)+1:inds(b+1)) = -data.m(a).x(inds(b)+1:inds(b+1));
                data.m(a).y(inds(b)+1:inds(b+1)) = -data.m(a).y(inds(b)+1:inds(b+1));
                data.m(a).z(inds(b)+1:inds(b+1)) = -data.m(a).z(inds(b)+1:inds(b+1));
            end
            
        end
        
        Data(ifile).m = data.m;
        Data(ifile).time = data.time;
        Data(ifile).converteuler = data.converteuler;
        
    else %we assume that otherwise it was recorded using the wide-range transmitter setup
        
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
            
            if isfield(data,'TrackerTime')
                data.time = data.TrackerTime;
            elseif isfield(header,'Sampling_Rate') && isfield(data.m(1),'x')
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
                data.converteuler = 1;
                
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
                data.converteuler = 0;
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
        
        
        %we want to reorient the axes to a standard configuration, where x points
        %to the participant's right, y is straight ahead to the participant, and z
        %points up.
        %this is equivalent to applying a rotation matrix:
        RotMat = [ 0  1  0;
                   1  0  0;
                   0  0 -1;
                 ];
        
        %we will manually apply this below:
        
        for a = 1:length(data.m)
            tmp1 = data.m(a).x; %x axis originally points to the patient's forward
            tmp2 = data.m(a).y; %y axis originally points to the patient's right
            tmp3 = data.m(a).z; %z axis originally points to the patient's down
            
            data.m(a).x = tmp2;  %x axis now points to the patient's right
            data.m(a).y = tmp1;  %y axis now points to the patient's forward
            data.m(a).z = -tmp3; %z axis now points to the patient's upward
            
            %we also need to apply this transformation to the rotation matrices
            for b = 1:size(data.m(a).rotang,3)
                data.m(a).rotang(:,:,b) = RotMat*data.m(a).rotang(:,:,b);
            end
            
            
        end
        
        %we will assume data are in inches, so we need to convert units to cm. This
        % is true for the imitation data, but the gesture-to-sight data is already
        % collected in meters.
        if strcmpi(filetype,'aaron')
            %aaron's code writes the imitation data in inches/1000 and the
            %gesture-to-sight data in meters. so we have to convert only the
            %imitation data to meters.
            
            if ~contains(fullfile(fpath{ifile},fname{ifile}),'gesture','IgnoreCase',true) && ~contains(fullfile(fpath{ifile},fname{ifile}),'gts','IgnoreCase',true) && ~contains(fullfile(fpath{ifile},fname{ifile}),'VF')
                for a = 1:length(data.m)
                    data.m(a).x = (data.m(a).x)*0.0254*1000;  %multiply by 1000 since we are accidentally dividing by 1000 in the data file.
                    data.m(a).y = (data.m(a).y)*0.0254*1000;
                    data.m(a).z = (data.m(a).z)*0.0254*1000;
                end
            end
            
        elseif strcmpi(filetype,'steve')
            %steve wrote out all his data in inches so we just convert to meters
            
            for a = 1:length(data.m)
                data.m(a).x = (data.m(a).x)*0.0254;
                data.m(a).y = (data.m(a).y)*0.0254;
                data.m(a).z = (data.m(a).z)*0.0254;
            end
        end
        
        %there is a potential problem of data crossing out into the wrong
        %hemisphere. to avoid this, we will detect sudden discontinuities in two of
        %the three axes of the data, and when we find one we will systematically flip
        %the sign until the next discontinuity. these discontinuities should be
        %paired.
        for a = 1:length(data.m)
            velx = diff(data.m(a).x);
            vely = diff(data.m(a).y);
            velz = diff(data.m(a).z);
            
            accx = diff(velx);
            accy = diff(vely);
            accz = diff(velz);
            
            indx = find(abs(velx(1:end)) > 10);
            indy = find(abs(vely(1:end)) > 10);
            indz = find(abs(velz(1:end)) > 10);
            
            %find and throw out all the single-sample spikes
            indax = find(abs(accx(1:end)) > 10);
            inday = find(abs(accy(1:end)) > 10);
            indaz = find(abs(accz(1:end)) > 10);
            for b = 1:length(indax)-1
                if indax(b+1)-indax(b) < 2
                    indx(indx == indax(b+1)) = [];
                end
            end
            for b = 1:length(inday)-1
                if inday(b+1)-inday(b) < 2
                    indy(indy == inday(b+1)) = [];
                end
            end
            for b = 1:length(indaz)-1
                if indaz(b+1)-indaz(b) < 2
                    indz(indz == indaz(b+1)) = [];
                end
            end
            
            if isempty(indx)
                %y and z
                inds = intersect(indy,indz);
            elseif isempty(indy)
                %x and z
                inds = intersect(indx,indz);
            elseif isempty(indz)
                %x and y
                inds = intersect(indx,indy);
            else
                inds = [];
            end
            
            if isempty(inds)
                continue;
            end
            
            if size(inds,1) > size(inds,2) %make sure this is a row vector
                inds = inds';
            end
            
            %comment this line out to avoid automatically fixing odd-mark errors by
            %assuming the last point is missing - but sometimes it is the first point!
            %     if mod(length(inds),2) ~= 0
            %         inds(end+1) = length(data.m(a).x);
            %     end
            
            %     %comment these lines out to automatically process data without manual
            %     %verification
            
            %fprintf('   Subj: %s\tTrial: %d',data.sub,
            fpn = [fpath fname];
            tmpinds = strfind(fpn,filesep);
            fprintf('   %s\n',fpn(tmpinds(end)+1:end));
            
            figure(2)
            clf;
            torso = [(data.m(7).x(1)+data.m(8).x(1))/2 (data.m(7).y(1)+data.m(8).y(1))/2 (data.m(7).z(1)+data.m(8).z(1))/2];
            joint.refshoulder = [data.m(8).x(1) data.m(8).y(1) data.m(8).z(1)];
            joint.shoulder = [data.m(7).x(1) data.m(7).y(1) data.m(7).z(1)];
            joint.elbow = [data.m(6).x(1) data.m(6).y(1) data.m(6).z(1)];
            joint.wrist = [data.m(5).x(1) data.m(5).y(1) data.m(5).z(1)];
            joint.hand = [data.m(4).x(1) data.m(4).y(1) data.m(4).z(1)];
            joint.middlefinger = [data.m(3).x(1) data.m(3).y(1) data.m(3).z(1)];
            joint.indexfinger = [data.m(2).x(1) data.m(2).y(1) data.m(2).z(1)];
            joint.thumb = [data.m(1).x(1) data.m(1).y(1) data.m(1).z(1)];
            if contains(lower(fpath),'6markervers')
                %do simplified arm without wrist/hand
                [xSE,ySE,zSE] = cylinder2P(.45,8,joint.shoulder(1,:),joint.elbow(1,:));
                [xEW,yEW,zEW] = cylinder2P(.375,8,joint.elbow(1,:),joint.wrist(1,:));
                [xWF,yWF,zWF] = cylinder2P([.25 .25 .25 .1 .1 .1],8,joint.wrist(1,:),joint.indexfinger(1,:));
                [xWT,yWT,zWT] = cylinder2P([.25 .25 .25 .1 .1 .1],8,joint.wrist(1,:),joint.thumb(1,:));
            else%if contains(lower(fpath),'8markervers')
                [xSE,ySE,zSE] = cylinder2P(.45,8,joint.shoulder(1,:),joint.elbow(1,:));
                [xEW,yEW,zEW] = cylinder2P(.375,8,joint.elbow(1,:),joint.wrist(1,:));
                [xWH,yWH,zWH] = cylinder2P(.375,8,joint.wrist(1,:),joint.hand(1,:));
                [xHM,yHM,zHM] = cylinder2P([.25 .25 .25 .1 .1 .1],8,joint.hand(1,:),joint.middlefinger(1,:));
                [xWF,yWF,zWF] = cylinder2P([.25 .25 .25 .1 .1 .1],8,joint.hand(1,:),joint.indexfinger(1,:));
                [xWT,yWT,zWT] = cylinder2P([.25 .25 .25 .1 .1 .1],8,joint.hand(1,:),joint.thumb(1,:));
            end
            %rectangle representing the body, assuming a vertical plane
            torsolength = abs(joint.shoulder(1,1)-joint.refshoulder(1,1))*1.5;
            xtorso = [joint.shoulder(1,1); joint.refshoulder(1,1); joint.refshoulder(1,1); joint.shoulder(1,1)];
            ytorso = [joint.shoulder(1,2); joint.refshoulder(1,2); joint.refshoulder(1,2); joint.shoulder(1,2)];
            ztorso = [joint.shoulder(1,3); joint.refshoulder(1,3); joint.refshoulder(1,3)-torsolength; joint.shoulder(1,3)-torsolength];
            headrad = abs(joint.refshoulder(1,1)-joint.shoulder(1,1))/4;
            headc = torso + [0 0 headrad];
            theta = [0:pi/10:2*pi];
            headx = headrad*cos(theta)+headc(1);
            heady = zeros(size(headx))+headc(2);
            headz = headrad*sin(theta)+headc(3);
            
            [spherex,spherey,spherez] = sphere(10);
            
            h = patch('XData',xtorso,'YData',ytorso,'ZData',ztorso);
            set(h,'FaceColor',[.5 .5 .5]);
            hold on;
            h = patch('XData',headx,'YData',heady,'ZData',headz);
            set(h,'FaceColor',[.5 .5 .5]);
            h = surf(xSE, ySE, zSE);
            set(h,'FaceColor',[1 0 0]);
            h = surf(xEW, yEW, zEW);
            set(h,'FaceColor',[1 1 0]);
            if contains(lower(fpath),'8markervers')
                h = surf(xWH, yWH, zWH);
                set(h,'FaceColor',[1 1 0]);
                h = surf(xHM, yHM, zHM);
                set(h,'FaceColor',[0 0 1]);
            end
            h = surf(xWF, yWF, zWF);
            set(h,'FaceColor',[0 0 1]);
            h = surf(xWT, yWT, zWT);
            set(h,'FaceColor',[0 0 1]);
            h = surf(spherex+data.m(a).x(1),spherey+data.m(a).y(1),spherez+data.m(a).z(1));
            set(h,'FaceColor',[.8 0 0]);
            hold off;
            grid on;
            view(195,20);
            xlabel('x');
            ylabel('y');
            zlabel('z');
            
            figure(1)
            inds = markdata3d(data.m(a).x,data.m(a).y,data.m(a).z,length(data.m(a).x),inds,0,1,'Name',sprintf('Marker %d',a));
            
            inds = sort(inds);
            if inds(end) > length(data.m(a).x)
                inds(end) = length(data.m(a).x);
            end
            if inds(1) < 1
                inds(1) = 1;
            end
            
            while mod(length(inds),2) ~= 0  %prompt until we have an even number of marks
                disp('Odd number of marks.');
                inds = markdata3d(data.m(a).x,data.m(a).y,data.m(a).z,length(data.m(a).x),inds,0,1,'Name',sprintf('Marker %d',a));
                
                if inds(end) > length(data.m(a).x)
                    inds(end) = length(data.m(a).x);
                end
                if inds(1) < 1
                    inds(1) = 1;
                end
            end
            
            inds = sort(inds);
            
            for b = 1:2:length(inds)
                data.m(a).x(inds(b)+1:inds(b+1)) = -data.m(a).x(inds(b)+1:inds(b+1));
                data.m(a).y(inds(b)+1:inds(b+1)) = -data.m(a).y(inds(b)+1:inds(b+1));
                data.m(a).z(inds(b)+1:inds(b+1)) = -data.m(a).z(inds(b)+1:inds(b+1));
            end
            
        end
        
        Data(ifile).m = data.m;
        Data(ifile).time = data.time;
        Data(ifile).converteuler = data.converteuler;

        
    end %end upside-down or rightside-up transmitter transformation
    
end

%% transform into body-centered coordinates


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


for a = 1:length(Data)
    
    %calculate and remove the origin, which we will set at the left shoulder
    origin = [Data(a).m(7).x, Data(a).m(7).y, Data(a).m(7).z];
    
    for c = 1:length(Data(a).m)
        Data(a).m(c).x = Data(a).m(c).x-origin(:,1);
        Data(a).m(c).y = Data(a).m(c).y-origin(:,2);
        Data(a).m(c).z = Data(a).m(c).z-origin(:,3);
    end
    
    for d = 1:size(Data(a).m(1).x,1)
        
        
%         %calculate the angle between the shoulders and the x axis; this
%         %will tell us how far to rotate the shoulders to align the frontal
%         %plane with the x-z plane.
%         v1 = [Data(a).m(8).x(d)-Data(a).m(7).x(d), Data(a).m(8).y(d)-Data(a).m(7).y(d), Data(a).m(8).z(d)-Data(a).m(7).z(d)];
%         theta = -acos(v1(1) / sqrt(sum(v1(1:2).^2,2)));
%         thetas(d) = theta;
%         
%         %calculate the rotation matrix
%         RotMat = [cos(theta) -sin(theta) 0;
%                   sin(theta)  cos(theta) 0;
%                   0           0          1];

        %calculate the body-centered coordinate frame based on the
        %shoulder-shoulder vector and the gravity vector. first, compute
        %the normal vector (which points in the y_hat direction).
        v1 = [Data(a).m(8).x(d)-Data(a).m(7).x(d), Data(a).m(8).y(d)-Data(a).m(7).y(d), Data(a).m(8).z(d)-Data(a).m(7).z(d)];
        v2 = [0 0 -1]; %gravity vector starts at the shoulder and points straight down
        n = cross(v1,v2); %cross product of v1 and v2 is the y_hat vector
        z_hat = cross(v1,n); %cross product of v1 and n is the z_hat vector
        
        %now we assmble the rotation matrix that will take us from world
        %coordinates to body coordinates. I think this assumes we are pre-multiplying
        v1 = v1./sqrt(sum(v1.^2));
        n = n./sqrt(sum(n.^2));
        z_hat = z_hat./sqrt(sum(z_hat.^2));
        RotMat = [v1; 
                  n; 
                  z_hat];

        %rotate all the markers and the rotation matrix for this sample
        for c = 1:length(Data(a).m)
            
            vec = [Data(a).m(c).x(d) Data(a).m(c).y(d) Data(a).m(c).z(d)];
            rotvec = RotMat * vec';
            Data(a).m(c).x(d,1) = rotvec(1);
            Data(a).m(c).y(d,1) = rotvec(2);
            Data(a).m(c).z(d,1) = rotvec(3);
            
            %rotate the rotation matrix (U_hat = RUR')
            Data(a).m(c).rotang(:,:,d) = RotMat*Data(a).m(c).rotang(:,:,d);
            
            %recalculate the Euler angles from the rotated rotation matrix
            Data(a).m(c).azim(d) = atan2d(Data(a).m(c).rotang(2,1,d),Data(a).m(c).rotang(1,1,d));
            Data(a).m(c).elev(d) = -asind(Data(a).m(c).rotang(3,1,d));
            Data(a).m(c).roll(d) = atan2d(Data(a).m(c).rotang(3,2,d),Data(a).m(c).rotang(3,3,d));
            
        end
        
        
    end
    
    Thetas{a} = thetas;

end

%% cut the data

for a = 1:length(Data)
    
    %calculate the velocity if possible
    dt = mean(diff(Data(a).time));
    
    for c = 1:length(Data(a).m)
        Data(a).pos(:,1,c) = Data(a).m(c).x;
        Data(a).pos(:,2,c) = Data(a).m(c).y;
        Data(a).pos(:,3,c) = Data(a).m(c).z;
        
        Data(a).vel(:,1,c) = gradient(sgolayfilt(Data(a).m(c).x,2,min([19,2*floor((size(Data(a).m(c).x,1)-2)/2)+1])),dt);
        Data(a).vel(:,2,c) = gradient(sgolayfilt(Data(a).m(c).y,2,min([19,2*floor((size(Data(a).m(c).y,1)-2)/2)+1])),dt);
        Data(a).vel(:,3,c) = gradient(sgolayfilt(Data(a).m(c).z,2,min([19,2*floor((size(Data(a).m(c).z,1)-2)/2)+1])),dt);
    end
    
    inds = addmarks(Data(a).vel,'Nchan',3,'throwmid','vthresh',0.3,'vthreshMin',0.2);
    Data(a).inds = inds{1};
    
    %manually verify indices
    Data(a).inds = markdataNd(permute(Data(a).pos,[1 3 2]),size(Data(a).pos,1),0,Data(a).inds,0);

end


%% visually inspect the data

figure(1);
clf;
ylabels = {'thumb','index','middle','hand','wrist','elbow'};

for a = 1:6
    subplot(3,2,a)
    set(gca,'ColorOrder',[0 0 .7; .7 0 0; 0 .7 0; .4 .4 1; 1 .4 .4; .4 1 .4],'NextPlot','ReplaceChildren')
    inds = Data(1).inds;
    t = Data(1).time(inds(1):inds(2));
    t = t-t(1);
    h = plot(t,squeeze(Data(1).pos(inds(1):inds(2),:,a)),'-');
    set(h(1),'LineWidth',2);
    set(h(2),'LineWidth',2);
    set(h(3),'LineWidth',2);
    hold on;
    inds = Data(2).inds;
    t = Data(2).time(inds(1):inds(2));
    t = t-t(1);
    h = plot(t,squeeze(Data(2).pos(inds(1):inds(2),:,a)),'-');
    set(h(1),'LineWidth',1);
    set(h(2),'LineWidth',1);
    set(h(3),'LineWidth',1);
    ylabel([ylabels{a} 'pos']);
    set(gca,'ColorOrder','factory','NextPlot','replace');
    
    if a == 1
        legend('Data1-x','Data1-y','Data1-z','Data2-x','Data2-y','Data2-z');
    end
    
end


pd = [];
N = 200;

%compute Pdist
for a = 1:6
   
    inds = Data(1).inds;
    t = Data(1).time(inds(1):inds(2));
    t = t-t(1);
    dur = t(end);
    dN = dur/N;
    tinterp = [0:dN:dur];
    x = squeeze(Data(1).pos(inds(1):inds(2),:,a));
    x = csapi(t',x',tinterp);
    x = x';
    
    inds = Data(2).inds;
    t = Data(2).time(inds(1):inds(2));
    t = t-t(1);
    dur = t(end);
    dN = dur/N;
    tinterp = [0:dN:dur];
    y = squeeze(Data(2).pos(inds(1):inds(2),:,a));
    y = csapi(t',y',tinterp);
    y = y';
    
    pd(a) = procrustesMultiScale(x,y);
end

fprintf('\n\nProcrustes (sensors):\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n\n',pd(1),pd(2),pd(3),pd(4),pd(5),pd(6));


% 
% %plot the FT 
% figure(2); clf;
% figure(3); clf;
% figure(4); clf;
% figure(5); clf;
% figure(6); clf;
% figure(7); clf;
% %ylabels = {'thumb','index','middle','hand','wrist','elbow'};
% 
% iplot = 0;
% for a = 1:6
%     indsA = Data(1).inds;
%     tA = Data(1).time(indsA(1):indsA(2));
%     tA = tA-tA(1);
%     LA = length(tA);
%     TA = mean(diff(tA));
%     FsA = 1/TA;
%     YA = fft(squeeze(Data(1).pos(indsA(1):indsA(2),:,a)));
%     %P2A = abs(YA/LA);
%     %P1A = P2A(1:LA/2+1,:);
%     %P1A(2:end-1,:) = 2*P1A(2:end-1,:);
%     %fA = FsA*(0:(LA/2))/LA;
%     P1A = abs(YA);
%     P1A = P1A(1:ceil(LA/2),:);
%     fA = FsA*(0:(ceil(LA/2)-1))/LA;
%     
%     indsB = Data(2).inds;
%     tB = Data(2).time(indsB(1):indsB(2));
%     tB = tB-tB(1);
%     LB = length(tB);
%     TB = mean(diff(tB));
%     FsB = 1/TB;
%     YB = fft(squeeze(Data(2).pos(indsB(1):indsB(2),:,a)));
%     %P2B = abs(YB/LB);
%     %P1B = P2B(1:LB/2+1,:);
%     %P1B(2:end-1,:) = 2*P1B(2:end-1,:);
%     P1B = abs(YB);
%     P1B = P1B(1:ceil(LB/2),:);
%     fB = FsB*(0:(ceil(LB/2)-1))/LB;
%     
%     figure(a+1);
% 
%     for b = 1:3
%         %iplot = iplot+1;
%         plotlength = max([ceil(length(fA)/2)+1,ceil(length(fB)/2)+1]);
%         subplot(1,3,b)
%         h = plot(fA(1:plotlength),10*log10(P1A(1:plotlength,b)),'b-');
%         hold on;
%         h = plot(fB(1:plotlength),10*log10(P1B(1:plotlength,b)),'g-');
%         if b == 1
%             ylabel([ylabels{a} '|P1(f)|']);
%         end
%         
%     end
%     
%     if a == 1
%         legend('Data1','Data2');
%     end
%     if a >= 5
%         xlabel('f (Hz)');
%     end
%     
% end
