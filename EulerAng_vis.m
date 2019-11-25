function varargout = EulerAng_vis(varargin)
% EULERANG_VIS MATLAB code for EulerAng_vis.fig
%      EULERANG_VIS, by itself, creates a new EULERANG_VIS or raises the existing
%      singleton*.
%
%      H = EULERANG_VIS returns the handle to a new EULERANG_VIS or the handle to
%      the existing singleton*.
%
%      EULERANG_VIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EULERANG_VIS.M with the given input arguments.
%
%      EULERANG_VIS('Property','Value',...) creates a new EULERANG_VIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EulerAng_vis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EulerAng_vis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EulerAng_vis

% Last Modified by GUIDE v2.5 15-Nov-2019 18:31:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EulerAng_vis_OpeningFcn, ...
                   'gui_OutputFcn',  @EulerAng_vis_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before EulerAng_vis is made visible.
function EulerAng_vis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EulerAng_vis (see VARARGIN)

rotmat = varargin{1}; %first entry is the rotation matrix. Dims are 3x3xN
orig_rotmat = rotmat; %keep a record of the unmodified rotation matrix so we can revert sections back as requested
setappdata(handles.figure1,'rotmat',rotmat);
setappdata(handles.figure1,'orig_rotmat',orig_rotmat);


%set the axis view angles
views.az = -37.5;
views.el = 30;
views.allowupdate = 0;
setappdata(handles.figure1,'views',views);

%set up the slider
set(handles.slider1,'Min',1);
set(handles.slider1,'Max',size(rotmat,3));
set(handles.slider1,'SliderStep',[1/size(rotmat,3) , 10/size(rotmat,3)]);

c = 1;
setappdata(handles.figure1,'c',c);

indStarts = [];  %this will be a 2d array of start and end marks
indEnds = [];
setappdata(handles.figure1,'indStarts',indStarts);
setappdata(handles.figure1,'indEnds',indEnds);

%sets up the play timer. Change the Period to speed up/slow down playback
playtimer = timer('Name','DataPlayer','Period',0.5,'StartDelay',0.1,'ExecutionMode','fixedRate','TimerFcn',{@timerCallbackFcn,handles});
stop(playtimer);
setappdata(handles.figure1,'playtimer',playtimer);

updateplot(handles);

views.allowupdate = 1;
setappdata(handles.figure1,'views',views);


% Choose default command line output for EulerAng_vis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EulerAng_vis wait for user response (see UIRESUME)
uiwait(handles.figure1);



% --- Executes when user attempts to close figure1.
function varargout = EulerAng_vis_OutputFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%output the modified rotation matrix
rotmat = getappdata(handles.figure1,'rotmat');
varargout{1} = rotmat;

%output the indices of the modified sections
inds = getappdata(handles.figure1,'inds');
varargout{2} = inds;

varargout{3} = handles;

delete(handles.figure1);



% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'figure1')
    playtimer = getappdata(handles.figure1,'playtimer');
    if isvalid(playtimer)
        stop(playtimer);
    end
end

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end


% --- Executes on button press in pushbtnAddStart (Add Start).
function pushbtnAddStart_Callback(hObject, eventdata, handles)
% hObject    handle to pushbtnAddStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

indStarts = getappdata(handles.figure1,'indStarts');
indEnds = getappdata(handles.figure1,'indEnds');
c = getappdata(handles.figure1,'c');

if ~any(indEnds == c)  %we will take this as a start as long as it isn't an end!
    indStarts = [indStarts; c];
    indStarts = sort(unique(indStarts)); %make sure we don't double mark, and sort in order
end

setappdata(handles.figure1,'indStarts',indStarts);

updateplot(handles);



% --- Executes on button press in pushbtnRmvStart (Remove Start).
function pushbtnRmvStart_Callback(hObject, eventdata, handles)
% hObject    handle to pushbtnRmvStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

indStarts = getappdata(handles.figure1,'indStarts');
c = getappdata(handles.figure1,'c');

if ~isempty(indStarts)
    diffinds = indStarts-c;
    
    [~,idiff] = min(abs(diffinds));
    indStarts(idiff) = [];
end

setappdata(handles.figure1,'indStarts',indStarts);

updateplot(handles);


% --- Executes on button press in pushbtnAddEnd (Add End).
function pushbtnAddEnd_Callback(hObject, eventdata, handles)
% hObject    handle to pushbtnDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

indEnds = getappdata(handles.figure1,'indEnds');
indStarts = getappdata(handles.figure1,'indStarts');
c = getappdata(handles.figure1,'c');

if ~any(indStarts == c)  %we will take this as an end as long as it isn't a start! 
    indEnds = [indEnds; c];
    indEnds = sort(unique(indEnds)); %make sure we don't double mark, and sort in order
end

setappdata(handles.figure1,'indEnds',indEnds);

updateplot(handles);


% --- Executes on button press in pushbtnRmvEnd.
function pushbtnRmvEnd_Callback(hObject, eventdata, handles)
% hObject    handle to pushbtnRmvEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

indEnds = getappdata(handles.figure1,'indEnds');
c = getappdata(handles.figure1,'c');

if ~isempty(indEnds)
    diffinds = indEnds-c;
    
    [~,idiff] = min(abs(diffinds));
    indEnds(idiff) = [];
end

setappdata(handles.figure1,'indEnds',indEnds);

updateplot(handles);



% --- Executes on button press in pushbtnInterp (Interpolate Section).
function pushbtnInterp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbtnInterp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

c = getappdata(handles.figure1,'c');
indStarts = getappdata(handles.figure1,'indStarts');
indEnds = getappdata(handles.figure1,'indEnds');
rotmat = getappdata(handles.figure1,'rotmat');
orig_rotmat = getappdata(handles.figure1,'orig_rotmat');

%find the segment we are going to be modifying. We only do this if there is
%at least 1 marked segment, and we are IN that segment
if ~isempty(indStarts) && ~isempty(indEnds)

    diffinds = indStarts-c;
    [~,idiff] = min(abs(diffinds));
    ind1 = indStarts(idiff);
    
    diffinds = indEnds-c;
    [~,idiff] = min(abs(diffinds));
    ind2 = indEnds(idiff);
    
    if c >= ind1 && c <= ind2
        %do the interpolation
        
        
        
        
        
    end
    
end

setappdata(handles.figure1,'rotmat',rotmat);




% --- Executes on button press in pushbtnRevert (Revert Section).
function pushbtnRevert_Callback(hObject, eventdata, handles)
% hObject    handle to pushbtnRevert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

c = getappdata(handles.figure1,'c');
indStarts = getappdata(handles.figure1,'indStarts');
indEnds = getappdata(handles.figure1,'indEnds');
rotmat = getappdata(handles.figure1,'rotmat');
orig_rotmat = getappdata(handles.figure1,'orig_rotmat');

%find the segment we are going to be modifying. We only do this if there is
%at least 1 marked segment, and we are IN that segment
if ~isempty(indStarts) && ~isempty(indEnds)

    diffinds = indStarts-c;
    [~,idiff] = min(abs(diffinds));
    ind1 = indStarts(idiff);
    
    diffinds = indEnds-c;
    [~,idiff] = min(abs(diffinds));
    ind2 = indEnds(idiff);
    
    if c >= ind1 && c <= ind2
        %revert back this section
        rotmat(:,:,ind1:ind2) = orig_rotmat(:,:,ind1:ind2);
    end
    
end

setappdata(handles.figure1,'rotmat',rotmat);



% --- Executes on button press in pushbtnDone.
function pushbtnDone_Callback(hObject, eventdata, handles)
% hObject    handle to pushbtnDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on mouse press over axes background.
function axes2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
doClickCallback(handles,handles.axes2)


% --- Executes on mouse press over axes background.
function axes3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
doClickCallback(handles,handles.axes3)


% --- Executes on mouse press over axes background.
function axes4_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
doClickCallback(handles,handles.axes4)



% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject,'Min',1);
set(hObject,'Max',2000);
set(hObject,'Value',1);
set(hObject,'SliderStep',[1/2000 , 10/2000]);


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

if isfield(handles,'figure1')
    playtimer = getappdata(handles.figure1,'playtimer');
    stop(playtimer);
end

rotmat = getappdata(handles.figure1,'rotmat');

set(handles.pushbtnAddStart,'enable','on');
set(handles.pushbtnInterp,'enable','on');
set(handles.pushbtnRevert,'enable','on');
set(handles.pushbtnRmvStart,'enable','on');
set(handles.pushbtnRevert,'enable','on');
c = get(handles.slider1,'Value');
c = round(c);

if c < 1
    c = 1;
end
if c > size(rotmat,3)
    c = size(rotmat,3);
end

setappdata(handles.figure1,'c',c);

updateplot(handles);


% --- Executes on button press in pushbtnPlay (Play).
function pushbtnPlay_Callback(hObject, eventdata, handles)
% hObject    handle to pushbtnPlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pushbtnAddStart,'enable','off');
set(handles.pushbtnInterp,'enable','off');
set(handles.pushbtnRevert,'enable','off');
set(handles.pushbtnRmvStart,'enable','off');
set(handles.pushbtnRevert,'enable','off');

playtimer = getappdata(handles.figure1,'playtimer');
pl = get(playtimer);
if strcmpi(pl.Running,'off')  %only ask to start the timer if it was off before.
    start(playtimer);
end




% --- Executes on button press in pushbtnStop (Stop).
function pushbtnStop_Callback(hObject, eventdata, handles)
% hObject    handle to pushbtnStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

playtimer = getappdata(handles.figure1,'playtimer');
stop(playtimer);

set(handles.pushbtnAddStart,'enable','on');
set(handles.pushbtnInterp,'enable','on');
set(handles.pushbtnRevert,'enable','on');
set(handles.pushbtnRmvStart,'enable','on');
set(handles.pushbtnRevert,'enable','on');




function timerCallbackFcn(obj, event, handles)
c = getappdata(handles.figure1,'c');
rotmat = getappdata(handles.figure1,'rotmat');
c = c+10;  %number of frames to step forward on each tick

if c > size(rotmat,3)
    c = 1;
    playtimer = getappdata(handles.figure1,'playtimer');
    stop(playtimer);
    
    set(handles.pushbtnAddStart,'enable','on');
    set(handles.pushbtnInterp,'enable','on');
    set(handles.pushbtnRevert,'enable','on');
    set(handles.pushbtnRmvStart,'enable','on');
    set(handles.pushbtnRevert,'enable','on');

end

setappdata(handles.figure1,'c',c);

set(handles.slider1,'Value',c);

updateplot(handles);




% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)

rotmat = getappdata(handles.figure1,'rotmat');

c = getappdata(handles.figure1,'c');

c = c + eventdata.VerticalScrollCount;

if c < 1
    c = 1;
end
if c > size(rotmat,3)
    c = size(rotmat,3);
end

setappdata(handles.figure1,'c',c);
    
set(handles.slider1,'Value',c);
    
updateplot(handles);






% --- create new function to update plot
function updateplot(handles)

c = getappdata(handles.figure1,'c');
indStarts = getappdata(handles.figure1,'indStarts');
indEnds = getappdata(handles.figure1,'indEnds');
rotmat = getappdata(handles.figure1,'rotmat');
views = getappdata(handles.figure1,'views');

vecX = [1 0 0]'; %vector aligned along the x axis
vecY = [0 1/3 0]'; %vector aligned along the y axis
vecZ = [0 0 1/3]'; %vector aligned along the z axis

%we will plot the data in the current segment or, if no segments are
%marked, the data in a 100-sample window.

%find the segment we are going to be visualizing
if isempty(indStarts)
%     if c < size(rotmat,3)-99
%         ind1 = max(1,c-49);
%     else
%         ind1 = size(rotmat,3)-99;
%     end
ind1 = 1;
else
    diffinds = indStarts-c;
    [~,idiff] = min(abs(diffinds));
    ind1 = indStarts(idiff);
end
if isempty(indEnds)
%     if c > size(rotmat,3)-99
%         ind2 = size(rotmat,3);
%     else
%         ind2 = c+49;
%     end
ind2 = size(rotmat,3);
else
    diffinds = indEnds-c;
    [~,idiff] = min(abs(diffinds));
    ind2 = indEnds(idiff);
end

%if there is another end between the indStart and c, drop the start
tmp = c-indEnds;
tmp(tmp < 0) = [];
if any(tmp < (c-ind1))
    ind1 = max(1,c-49);
end
tmp = indStarts - c;
tmp(tmp < 0) = [];
if any(tmp < (ind2-c))
    ind2 = min(c+49,size(rotmat,3));
end

vecplotX = [];
vecplotY = [];
vecplotZ = [];
for d = 1:size(rotmat,3)
    vecplotX(:,d) = rotmat(:,:,d)*vecX;
    vecplotY(:,d) = rotmat(:,:,d)*vecY;
    vecplotZ(:,d) = rotmat(:,:,d)*vecZ;
end

axes(handles.axes1)
[az,el] = view;
if (az ~= views.az && views.allowupdate == 1)
    views.az = az;
    setappdata(handles.figure1,'views',views);
end
if (el ~= views.el && views.allowupdate == 1)
    views.el = el;
    setappdata(handles.figure1,'views',views);
end
cla(handles.axes1);
%set(gca,'xlim',[-1.2 1.2],'ylim',[-1.2 1.2],'zlim',[-1.2 1.2]);
sphere(handles.axes1);
hsphere = get(handles.axes1,'Children');
view(views.az,views.el);
axis('equal') %set the axes to be equal so that the sphere looks like a sphere
set(hsphere,'FaceAlpha',0.5,'EdgeAlpha',0)
pbaspect('manual');
hold on;
plot3(vecplotX(1,ind1:ind2),vecplotX(2,ind1:ind2),vecplotX(3,ind1:ind2),'k.-')
plot3(vecplotX(1,c),vecplotX(2,c),vecplotX(3,c),'ro')
arrow3([0 0 0],vecplotX(:,c)','r2')
arrow3([0 0 0],vecplotY(:,c)','y2')
arrow3([0 0 0],vecplotZ(:,c)','b2')

if c > ind2-9
    h = plot3(vecplotX(1,ind2:min(c+10,size(rotmat,3))),vecplotX(2,ind2:min(c+10,size(rotmat,3))),vecplotX(3,ind2:min(c+10,size(rotmat,3))),'k.-');
    set(h,'Color',[0.5 0.5 0.5]);
elseif c < ind1+9
    h = plot3(vecplotX(1,max(c-10,1):ind1),vecplotX(2,max(c-10,1):ind1),vecplotX(3,max(c-10,1):ind1),'k.-');
    set(h,'Color',[0.5 0.5 0.5]);
end
hold off;
set(gca,'xlimmode','auto','ylimmode','auto','zlimmode','auto')

xlabel('x');
ylabel('y');
zlabel('z');

azim = atan2d(squeeze(rotmat(2,1,:)),squeeze(rotmat(1,1,:)));
elev = -asind(squeeze(rotmat(3,1,:)));
roll = atan2d(squeeze(rotmat(3,2,:)),squeeze(rotmat(3,3,:)));

axes(handles.axes2)
cla(handles.axes2);
set(handles.axes2,'NextPlot','add');
h = plot(azim,'-');
set(h,'HitTest','off')
hold on
for a = 1:length(indStarts)
    h = plot(indStarts(a)*[1 1],get(gca,'ylim'),'r-');
    set(h,'HitTest','off')
end
for a = 1:length(indEnds)
    h = plot(indEnds(a)*[1 1],get(gca,'ylim'),'b-');
    set(h,'HitTest','off')
end
h = plot([c c],get(gca,'ylim'),'k:');
set(h,'HitTest','off','LineWidth',1.5)
hold off;
set(handles.axes2,'Ytick',[-360:90:360]);
ylabel('azim')

axes(handles.axes3)
cla(handles.axes3);
set(handles.axes3,'NextPlot','add');
h = plot(elev,'-');
set(h,'HitTest','off')
hold on
for a = 1:length(indStarts)
    h = plot(indStarts(a)*[1 1],get(gca,'ylim'),'r-');
    set(h,'HitTest','off')
end
for a = 1:length(indEnds)
    h = plot(indEnds(a)*[1 1],get(gca,'ylim'),'b-');
    set(h,'HitTest','off')
end
h = plot([c c],get(gca,'ylim'),'k:');
set(h,'HitTest','off','LineWidth',1.5)
hold off;
set(handles.axes3,'Ytick',[-360:90:360]);
ylabel('elev')

axes(handles.axes4)
cla(handles.axes4);
set(handles.axes4,'NextPlot','add');
h = plot(roll,'-');
set(h,'HitTest','off')
hold on
for a = 1:length(indStarts)
    h = plot(indStarts(a)*[1 1],get(gca,'ylim'),'r-');
    set(h,'HitTest','off')
end
for a = 1:length(indEnds)
    h = plot(indEnds(a)*[1 1],get(gca,'ylim'),'b-');
    set(h,'HitTest','off')
end
h = plot([c c],get(gca,'ylim'),'k:');
set(h,'HitTest','off','LineWidth',1.5)
hold off;
set(handles.axes4,'Ytick',[-360:90:360]);
ylabel('roll')

set(handles.text1,'String',num2str(c));



function doClickCallback(handles,haxes)

playtimer = getappdata(handles.figure1,'playtimer');
stop(playtimer);

set(handles.pushbtnAddStart,'enable','on');
set(handles.pushbtnInterp,'enable','on');
set(handles.pushbtnRevert,'enable','on');
set(handles.pushbtnRmvStart,'enable','on');
set(handles.pushbtnRevert,'enable','on');

coordinates = get(haxes,'CurrentPoint');
coords = coordinates(1,1:2);
c = round(coords(1));

set(handles.slider1,'Value',c);

setappdata(handles.figure1,'c',c);

% stype = get(handles.figure1,'SelectionType');
% if strcmpi(stype,'alt')  %right mouse click, delete the closest mark
%     
%     inds = getappdata(handles.figure1,'inds');
%     
%     if ~isempty(inds{MarkType})
%         diffinds = inds{MarkType}-c;
%         
%         [~,idiff] = min(abs(diffinds));
%         inds{MarkType}(idiff) = [];
%     end
%     
%     setappdata(handles.figure1,'inds',inds);
% 
% end

updateplot(handles)
