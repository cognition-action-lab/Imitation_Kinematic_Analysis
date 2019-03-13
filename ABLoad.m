%Generic function to load kinematic data from a data file from the apraxia battery

function [data,varargout] = ABLoad(varargin)

if nargin == 0
    [fname,fpath] = uigetfile('*.*','Select data file to analyze');
else
    fpath = varargin{1};
    fname = '';
end

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


%if this is an old data set, we want to re-align the marker numbers
if contains(lower(fpath),'6markervers')
    data.m(8) = data.m(6);  %right shoulder
    data.m(7) = data.m(5);  %left shoulder
    data.m(6) = data.m(4);  %elbow
    data.m(5) = data.m(3);  %wrist
    data.m(4).x = NaN*ones(size(data.m(1).x));         %back of hand, not recorded
    data.m(4).y = NaN*ones(size(data.m(1).x));         %back of hand, not recorded
    data.m(4).z = NaN*ones(size(data.m(1).x));         %back of hand, not recorded
    data.m(4).qual = NaN*ones(size(data.m(1).x));         %back of hand, not recorded
    data.m(4).azim = NaN*ones(size(data.m(1).x));         %back of hand, not recorded
    data.m(4).elev = NaN*ones(size(data.m(1).x));         %back of hand, not recorded
    data.m(4).roll = NaN*ones(size(data.m(1).x));         %back of hand, not recorded
    data.m(3).x = NaN*ones(size(data.m(1).x));         %middle finger, not recorded
    data.m(3).y = NaN*ones(size(data.m(1).x));         %middle finger, not recorded
    data.m(3).z = NaN*ones(size(data.m(1).x));         %middle finger, not recorded
    data.m(3).qual = NaN*ones(size(data.m(1).x));         %middle finger, not recorded
    data.m(3).azim = NaN*ones(size(data.m(1).x));         %middle finger, not recorded
    data.m(3).elev = NaN*ones(size(data.m(1).x));         %middle finger, not recorded
    data.m(3).roll = NaN*ones(size(data.m(1).x));         %middle finger, not recorded
end


%uncomment these lines to modify the sensor order...

%For patient 2221, we need to swap sensors 7 and 8
if contains(fpath,'2221') && contains(lower(fpath),'patient')
    tmp = data.m(8);
    data.m(8) = data.m(7);
    data.m(7) = tmp;
    clear tmp;
end

%for patient 2777, we have to swap sensor 4 and 8
if contains(fpath,'2777') && contains(lower(fpath),'patient')
    tmp = data.m(8);
    data.m(8) = data.m(4);
    data.m(4) = tmp;
    clear tmp;
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
        
        data.m(a).rotang(1,1,:) = cos(yaw).*cos(pitch);
        data.m(a).rotang(1,2,:) = cos(yaw).*sin(pitch).*sin(roll)-sin(yaw).*cos(roll);
        data.m(a).rotang(1,3,:) = cos(yaw).*sin(pitch).*cos(roll)+sin(yaw).*sin(roll);
        data.m(a).rotang(2,1,:) = sin(yaw).*cos(pitch);
        data.m(a).rotang(2,2,:) = sin(yaw).*sin(pitch).*sin(roll)+cos(yaw).*cos(roll);
        data.m(a).rotang(2,3,:) = sin(yaw).*sin(pitch).*cos(roll)-cos(yaw).*sin(roll);
        data.m(a).rotang(3,1,:) = -sin(pitch);
        data.m(a).rotang(3,2,:) = cos(pitch).*sin(roll);
        data.m(a).rotang(3,3,:) = cos(pitch).*cos(roll);
                            
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
        data.m(a).rotang(:,:,b) = data.m(a).rotang(:,:,b)*RotMat;
    end
    
    
end

%we will assume data are in inches, so we need to convert units to cm
for a = 1:length(data.m)
    data.m(a).x = (data.m(a).x)*2.54;
    data.m(a).y = (data.m(a).y)*2.54;
    data.m(a).z = (data.m(a).z)*2.54;
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

    indx = find(abs(velx(1:end)) > 10);
    indy = find(abs(vely(1:end)) > 10);
    indz = find(abs(velz(1:end)) > 10);
    
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
    tmpinds = strfind(fpn,'/');
    fprintf('   %s\n',fpn(tmpinds(end):end));
    
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
    elseif contains(lower(fpath),'8markervers')
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
    
    


%there is an analogous problem in the angle data where we hit 180 deg and
%reset to -180. we can fix this problem in the same way, by adding 360 to
%the angle within these discontinuous regions.
%we don't have to worry about this for the rotation matrix (which we will
%preferentially use going forward...)
for a = 1:length(data.m)
    velaz = diff(data.m(a).azim);
    velel = diff(data.m(a).elev);
    velro = diff(data.m(a).roll);

    indaz = find(abs(velaz(1:end)) > 100);
    indel = find(abs(velel(1:end)) > 100);
    indro = find(abs(velro(1:end)) > 100);
    
    %comment these lines to avoid automatically fixing odd-numbers of marks
    %by assuming the last mark is missing -- but sometimes it is the first
    %mark that is missing!
    if mod(length(indaz),2) ~= 0
        indaz(end+1) = length(data.m(a).azim);
    end
    if mod(length(indel),2) ~= 0
        indel(end+1) = length(data.m(a).elev);
    end
    if mod(length(indro),2) ~= 0
        indro(end+1) = length(data.m(a).roll);
    end

    
    if ~isempty(indaz)
%         %comment these lines out to automatically process data without manual
%         %verification 
%         figure(1)
%         indaz = markdata3d(data.m(a).azim,data.m(a).elev,data.m(a).roll,length(data.m(a).x),indaz,0,1,'Name',sprintf('Marker %d azim',a));

        for b = 1:2:length(indaz)
            if sign(data.m(a).azim(indaz(b))) == 1  %positive to negative discontinuity
                data.m(a).azim(indaz(b)+1:indaz(b+1)) = data.m(a).azim(indaz(b)+1:indaz(b+1))+360;
            else %negative to positive discontinuity
                data.m(a).azim(indaz(b)+1:indaz(b+1)) = data.m(a).azim(indaz(b)+1:indaz(b+1))-360;
            end
        end
        
    end
    
    if ~isempty(indel)
%         %comment these lines out to automatically process data without manual
%         %verification 
%         figure(1)
%         indel = markdata3d(data.m(a).elev,data.m(a).elev,data.m(a).roll,length(data.m(a).x),indel,0,1,'Name',sprintf('Marker %d elev',a));
        for b = 1:2:length(indel)
            if sign(data.m(a).elev(indel(b))) == 1  %positive to negative discontinuity
                data.m(a).elev(indel(b)+1:indel(b+1)) = data.m(a).elev(indel(b)+1:indel(b+1))+360;
            else %negative to positive discontinuity
                data.m(a).elev(indel(b)+1:indel(b+1)) = data.m(a).elev(indel(b)+1:indel(b+1))-360;
            end
        end
    end
    
    if ~isempty(indro)
%         %comment these lines out to automatically process data without manual
%         %verification 
%         figure(1)
%         indro = markdata3d(data.m(a).roll,data.m(a).roll,data.m(a).roll,length(data.m(a).x),indro,0,1,'Name',sprintf('Marker %d roll',a));
        for b = 1:2:length(indro)
            if sign(data.m(a).roll(indro(b))) == 1  %positive to negative discontinuity
                data.m(a).roll(indro(b)+1:indro(b+1)) = data.m(a).roll(indro(b)+1:indro(b+1))+360;
            else %negative to positive discontinuity
                data.m(a).roll(indro(b)+1:indro(b+1)) = data.m(a).roll(indro(b)+1:indro(b+1))-360;
            end
        end
    end
    
end


% %compute the velocity
% for a = 1:length(data.m)
%     data.m(a).velx = gradient(sgolayfilt(data.m(a).x,2,min([19,2*floor((length(data.m(a).x)-2)/2)+1])));
%     data.m(a).vely = gradient(sgolayfilt(data.m(a).x,2,min([19,2*floor((length(data.m(a).x)-2)/2)+1])));
%     data.m(a).velz = gradient(sgolayfilt(data.m(a).x,2,min([19,2*floor((length(data.m(a).x)-2)/2)+1])));
% end


if nargout > 2
    varargout{1} = fname;
end
if nargout > 3
    varargout{2} = fpath;
end
if nargout > 4
    varargout{3} = cols;
end

