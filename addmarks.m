%code to automatically identify some initial marks for a gesture based on a
%velocity threshold. 

function ind = addmarks(vel,varargin)


vthresh = 0.3;      %m/s threshold  %40*0.0254
vthreshMin = 0.2;    %m/s threshold  %10*0.0254
NchanAgree = 4;     %number of channels that have to meet the threshold
tcombine = 8;      %duration within which sections will be combined
throwfirst = 0;     %flag to throw out the first set of marks
throwlast = 0;      %flag to throw out the last set of marks
throwmid = 0;
%NchanWeight = [0.75 0.75 0.75 1.25 1.5 1.75 1 1]; %weighting on the channels, to weight agreement between the larger joints more than the fingers
NchanWeight = ones(1,8); %weighting on the channels, to weight agreement between the larger joints more than the fingers

a = 1;
while a <= length(varargin)
    switch(varargin{a})
        case 'vthresh'
            vthresh = varargin{a+1};
            a = a+2;
        case 'vthreshMin'
            vthreshMin = varargin{a+1};
            a = a+2;
        case 'Nchan'
            NchanAgree = varargin{a+1};
            a = a+2;
        case 'NchanWeight'
            NchanWeight = varargin{a+1};
            a = a+2;
        case 'tcombine'
            tcombine = varargin{a+1};
            a = a+2;
        case 'throw'
            throwfirst = 1;
            throwlast = 1;
            a = a+1;
        case 'throwfirst'
            throwfirst = 1;
            a = a+1;
        case 'throwlast'
            throwlast = 1;
            a = a+1;
        case 'throwmid'
            throwmid = 1;
            a = a+1;
        otherwise
            disp('Unrecognized input')
            a = a+1;
    end
end
    

vel3d = squeeze(sqrt(sum(vel.^2,2)));

%agreement across at least N channels
vless = vel3d > vthresh;
%vless = vless.*NchanWeight;
vless = sum(vless,2) >= NchanAgree;
if ~any(vless)
    vless = vel3d > vthresh;
    vless = sum(vless,2) > 1;
    if ~any(vless)
        vless = vel3d > vthresh;
        vless = sum(vless,2) > 0;
    end
end
veliless = find(vless == 1);
velidless = find(diff([0; veliless])>1); %now i have the indices into veli of all the starts of all reaches exceeding vthreshMax
if isempty(velidless)
    if ~isempty(veliless)
        velidless = veliless(1);
    else
        velidless = [];
    end
end
velidless = sort([velidless; velidless-1]);
velidless = [velidless(2:end); length(veliless)];
if ~isempty(veliless)
    veliless = veliless(velidless); %now i have the start and end indices of all reaches exceeding vthreshMax
else
    veliless = [];
end
for a = 1:length(veliless)-1
    if abs(veliless(a+1)-veliless(a)) < tcombine %if the start/end marks are too close, throw them out! alternatively maybe merge them?
        veliless(a) = NaN;
        veliless(a+1) = NaN;
    end
end
veliless = veliless(~isnan(veliless));
if mod(length(veliless),2) ~= 0
    %we hopefully now have an even number of marks...  if not we have to
    %fix this! we could get an isolated mark if the mark has both high or
    %both low peak velocities on both sides. if we find this case, we will
    %duplicate it we'll look for the first case where this happens and fix
    %it.

    for a = 2:length(veliless)-1
        if (all(vel3d(veliless(a-1)+1:veliless(a)-1) > vthresh) && all(vel3d(veliless(a)+1:veliless(a+1)-1) > vthresh)) || (all(vel3d(veliless(a-1)+1:veliless(a)-1) < vthresh) && all(vel3d(veliless(a)+1:veliless(a+1)-1) < vthresh))
            veliless = sort([veliless; veliless(a)]);
            break;
        end
    end
end
if mod(length(veliless),2) ~= 0
    %if we still haven't fixed the problem, we'll just return what we can
    if length(veliless) > 3
        ind{1} = [veliless(1) veliless(end)];
        ind{2} = [veliless(2:end-1)];
    elseif length(veliless) > 2
        ind{1} = [veliless(1) veliless(end)];
        ind{2} = [];
    else
        ind = cell(2);
    end
    return;
end


veli = [1; veliless; size(vel3d,1)];

%now we have an even number of marks, representing regions where the
%velocity < vthresh.

%within each pair, we will search for the index of the minimum velocity
i = 1;
clear velimin;
for a = 1:2:length(veli)
    
    vsection = vel3d(veli(a):veli(a+1),:);
    [~,imin] = min(sum(vsection,2));
    if (a < length(veli)-1) %if we're not on the last pair, find a start
        %iminfirst = find(sum(vsection < vthreshMin,2)>=NchanAgree,1,'last');
        iminfirst = find(sum(vsection(imin:end,:) > vthreshMin,2)>=NchanAgree,1,'first');
        if isempty(iminfirst)
            [~,iminfirst] = min(vsection);
            iminfirst = min(iminfirst);
        end
        velimin(i,2) = veli(a)+ imin + iminfirst-1;
    else
        velimin(i,2) = NaN;
    end
    
    if (a > 1)
        %iminlast = find(sum(vsection < vthreshMin,2)>=NchanAgree,1,'first');
        iminlast = find(sum(vsection(1:imin,:) > vthreshMin,2)>=NchanAgree,1,'last');
        if isempty(iminlast)
            [~,iminlast] = min(vsection);
            iminlast = min(iminlast);
        end
        velimin(i,1) = veli(a)+iminlast-1;
    else
        velimin(i,1) = NaN;
    end
    i = i+1;
    
end
velimin = [velimin(1:end-1,2) velimin(2:end,1)]; %these are the indices of velocity starts and ends
for a = 2:size(velimin,1)
    if abs(velimin(a,1)-velimin(a-1,2)) < (tcombine/3) %if successive start/end marks are too close, throw them out! alternatively maybe merge them?
        velimin(a,1) = NaN;
        velimin(a-1,2) = NaN;
    end
end
tmp = velimin;
velimin = tmp(~isnan(tmp(:,1)),1);
velimin(:,2) = tmp(~isnan(tmp(:,2)),2);

if throwfirst && size(velimin,1) > 2
    velimin(1,:) = [];
end
if throwlast && size(velimin,1) > 2
    velimin(end,:) = [];
end


% if size(velimin,1) > 2
% %     %we will throw out the first pair and last pair, which represent the
% %     %initial and final phase of the action and not the "core" action
% %     ind{1} = [velimin(2,1); velimin(end-1,2)];
% %     ind{2} = unique(sort(reshape([velimin(3:end-1,1); velimin(2:end-2,2)],[],1)));
%     
%     %if we want to preserve the first pair
%     ind{1} = [velimin(1,1); velimin(end,2)];
%     ind{2} = unique(sort(reshape([velimin(2:end,1); velimin(1:end-1,2)],[],1)));
%     
% else
if size(velimin,1) > 1
    ind{1} = [velimin(1,1); velimin(end,2)];
    ind{2} = unique(sort(reshape([velimin(2:end,1); velimin(1:end-1,2)],[],1)));
    if length(ind{2}) > 2
        ind{2}(1) = [];
        ind{2}(end) = [];
    end
    
else
    ind{1} = velimin';  %even this might be empty!
    ind{2} = [];
end

if throwmid == 1
    ind{2} = [];
end

% %at the individual sensor level
% for b = 1:size(vel3d,2)
%     %find all the indices of when the data meet the threshold
%     veli = find(abs(vel3d(:,b)) > vthresh);
%     velid = find(diff([0; veli])>1); %now i have the indices into veli of all the starts of all reaches exceeding vthreshMax
%     if isempty(velid)
%         velid = veli(1);
%     end
%     velid = sort([velid; velid-1]);
%     velid = [velid(2:end); length(veli)];
%     veli = veli(velid); %now i have the start and end indices of all reaches exceeding vthreshMax
%     if veli(1) > 1
%         veli = [1; veli];
%     end
%     if veli(end) < size(vel3d,1)
%         veli = [veli; size(vel3d,1)];
%     end
%     if mod(length(veli),2) ~= 0
%         %odd number of marks....
%         if any(vel3d(veli(1):veli(2)) > vthresh)
%             veli(1) = [];
%         elseif any(vel3d(veli(end-1):veli(end)) > vthresh)
%             veli(end) = [];
%         elseif (veli(2)-veli(1)) < (veli(end)-veli(end-1)) %if we can't figure it out, we'll just have to guess...
%             veli(1) = [];
%         else
%             veli(end) = [];
%         end
%     end
% 	
%     %now we have an even number of marks, representing regions where the
%     %velocity < vthresh.
% 
%     %within each pair, we will search for the index of the minimum velocity
%     i = 1;
%     for a = 1:2:length(veli)
%         
%         vsection = vel3d(veli(a):veli(a+1));
%         if (a < length(veli)-1) %if we're not on the last pair, find a start
%             iminfirst = find(vsection < vthreshMin,1,'last');
%             if isempty(iminfirst)
%                 [~,iminfirst] = min(vsection);
%             end
%             velimin(i,2) = veli(a)+ iminfirst-1;
%         else
%             velimin(i,2) = NaN;
%         end
%         
%         if (a > 1)
%             iminlast = find(vsection < vthreshMin,1,'first');
%             if isempty(iminlast)
%                 [~,iminlast] = min(vsection);
%             end
%             velimin(i,1) = veli(a)+iminlast-1;
%         else
%             velimin(i,1) = NaN;
%         end
%         i = i+1;
% 
%     end
%     velimin = [velimin(1:end-1,2) velimin(2:end,1)]; %these are the indices of velocity starts and ends
%     inds{b} = velimin;
% end

