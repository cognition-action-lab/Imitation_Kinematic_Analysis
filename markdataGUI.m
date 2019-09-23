function varargout = markdataGUI(varargin)
% markdataGUI    visualization program for movement data
%      markdataGUI(dat) creates a user interface for visualizing
%      and/or marking the data in dat from a single trial (movement). dat
%      can have one of two formats; it may be an array of structures where
%      each structure contains the x, y, and z coordinates of the movement
%      in the format m(i).x, m(i).y, and m(i).z where i is the index of the
%      marker (up to 8 trackers can be provided). alternatively, dat may be
%      a 3D array with dimensions [samples, axes, markers]; this should
%      therefore be a [: by 3 by 8] array. If a structure or array with
%      fewer markers is provided, the remaining dimensions are assumed to
%      be NaN; any additional markers beyond 8 will be ignored.
%
%      markdataGUI assumes that the markers are as follows:
%          1: thumb*
%          2: index finger*
%          3: middle finger
%          4: back of the hand
%          5: wrist*
%          6: elbow*
%          7: shoulder*
%          8: opposite shoulder
%      where the markers indicated with an asterisk are actually displayed
%      in the GUI. Note, markers 7 and 8 are used to align/center the other
%      markers as a reference point, so using these indices for other data
%      may cause confusion in the visualization.
%
%      [inds] = markdataGUI(dat) provides as output the array of indices
%      that were marked during the call to markdataGUI. This output is
%      only populated when the GUI window is closed.
%
%      Clicking on the slider or on any of the axes on the right side of
%      the figure will position a dashed line on the right-hand axes
%      indicating the current index, and will update the 3D plot. When an
%      index of interest has been identified, it can be saved by using the
%      "Mark" button. All currently saved indices will be denoted by a
%      vertical red line on the right-hand axes. To remove a mark, choose
%      a position close to the the mark to be removed and click the
%      "Unmark" button, which will remove the closest marked index. 
%
%      The "play" and "stop" features will automatically step through the
%      data set like a movie. Note that in some cases (i.e., when the
%      processor or Matlab are very slow) this can cause unintended
%      blocking, preventing additional updates to be made. In this case, it
%      may become necessary to force-quit out of the GUI.



% MARKDATAGUI MATLAB code for markdataGUI.fig
%      MARKDATAGUI, by itself, creates a new MARKDATAGUI or raises the existing
%      singleton*.
%
%      H = MARKDATAGUI returns the handle to a new MARKDATAGUI or the handle to
%      the existing singleton*.
%
%      MARKDATAGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARKDATAGUI.M with the given input arguments.
%
%      MARKDATAGUI('Property','Value',...) creates a new MARKDATAGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before markdataGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to markdataGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help markdataGUI

% Last Modified by GUIDE v2.5 04-Apr-2019 17:02:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @markdataGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @markdataGUI_OutputFcn, ...
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


% --- Executes just before markdataGUI is made visible.
function markdataGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to markdataGUI (see VARARGIN)

%create the data plot
m = varargin{1};

inds = cell(1,2);
nmarkers = [];
handles.ctp = [];
ma = [];
subenable = 1;

i = 2;
while i <= length(varargin)
    if ~ischar(varargin{i}) && ~isscalar(varargin{i}) %catch cases where invalid input into switch
        continue;
    end
    switch(varargin{i})
        case 'rot'
            if isstruct(m)
                for a = 1:length(m)
                    m(a).x = m(a).rotx;
                    m(a).y = m(a).roty;
                    m(a).z = m(a).rotz;
                end
            end
            i = i+1;
        case 'title'
            ttlstr = varargin{i+1};
            set(gcf,'Name',sprintf('markdataGUI: %s',ttlstr));
            i = i+2;
        case 'mark'
            inds = varargin{i+1};
            if ~iscell(inds)
                inds = {inds};
                inds{2} = [];
            end
            i = i+2;
        case 'nmarkers'
            nmarkers = varargin{i+1};
            i = i+2;
        case 'ctp' %curvature,torsion,armplane data
            handles.ctp = varargin{i+1};
            i = i+2;
        case 'ang'
            angm = varargin{i+1};
            if isstruct(angm) && isstruct(m)
                for a = 1:length(m)
                    m(a).azim = angm(a).azim;
                    m(a).elev = angm(a).elev;
                    m(a).roll = angm(a).roll;
                end
            elseif isstruct(angm) && ~isstruct(m)
                for a = 1:length(m)
                    ma(a).azim = angm(a).azim;
                    ma(a).elev = angm(a).elev;
                    ma(a).roll = angm(a).roll;
                end
            elseif ~isstruct(angm) && isstruct(m)
                for a = 1:size(angm,3)
                    m(a).azim = angm(:,1,a);
                    m(a).elev = angm(:,2,a);
                    m(a).roll = angm(:,3,a);
                end
            else  %is matrix
                for a = 1:size(angm,3)
                    ma(a).azim = angm(:,1,a);
                    ma(a).elev = angm(:,2,a);
                    ma(a).roll = angm(:,3,a);
                end
            end
            i = i+2;
        case 'full' %disable the subaction option
            subenable = 0;
        otherwise
            i = i+1;
    end
end

if subenable == 0
    set(handles.radiobutton2,'enable','off');
else%if subenable == 1
    set(handles.radiobutton2,'enable','on');
end

if ~isstruct(m)
    for a = 1:size(m,3)
        n(a).x = m(:,1,a);
        n(a).y = m(:,2,a);
        n(a).z = m(:,3,a);
    end
    m = n;
    
    for a = 1:length(m)
        m(a).azim = ma(a).azim;
        m(a).elev = ma(a).elev;
        m(a).roll = ma(a).roll;
    end
end

nm = 0;
for a = 1:length(m)
    if sum(isnan(m(a).x)) < 10
        nm = nm+1;
    end
end
nmarkers = nm;
if nmarkers ~= 6 && nmarkers ~= 8
    nmarkers = 6;
end

if ~isfield(m,'azim')
    for a = 1:length(m)
        m(a).azim = NaN*ones(size(m(a).x));
        m(a).elev = NaN*ones(size(m(a).x));
        m(a).roll = NaN*ones(size(m(a).x));
    end
end

for a = length(m)+1:7
    m(a).x = zeros(size(m(1).x));
    m(a).y = zeros(size(m(1).y));
    m(a).z = zeros(size(m(1).z));
    m(a).azim = NaN*ones(size(m(a).x));
    m(a).elev = NaN*ones(size(m(a).x));
    m(a).roll = NaN*ones(size(m(a).x));
end

if length(m) < 8
    m(8).x = m(7).x-30;
    m(8).y = m(7).y;
    m(8).z = m(7).z;
    m(8).azim = NaN*ones(size(m(8).x));
    m(8).elev = NaN*ones(size(m(8).x));
    m(8).roll = NaN*ones(size(m(8).x));
end

handles.m = m;
handles.ma = ma;


% for a = 1:size(m,3)
%     m(a).y = -m(a).y;
% end

handles.torso = [(m(7).x+m(8).x)/2 (m(7).y+m(8).y)/2 (m(7).z+m(8).z)/2];
%handles.torso = [zeros(size(m(7).x)) zeros(size(m(7).y)) zeros(size(m(7).z))];

handles.joint.refshoulder = [m(8).x m(8).y m(8).z] - handles.torso;
handles.joint.shoulder = [m(7).x m(7).y m(7).z] - handles.torso;
handles.joint.elbow =    [m(6).x m(6).y m(6).z] - handles.torso;
handles.joint.wrist =    [m(5).x m(5).y m(5).z] - handles.torso;
handles.joint.hand =     [m(4).x m(4).y m(4).z] - handles.torso;
handles.joint.middlefinger = [m(3).x m(3).y m(3).z] - handles.torso;
handles.joint.indexfinger = [m(2).x m(2).y m(2).z] - handles.torso;
handles.joint.thumbfinger = [m(1).x m(1).y m(1).z] - handles.torso;
handles.nmarkers = nmarkers;

%also store the angle data if it is available
handles.angle.refshoulder = [m(8).azim m(8).elev m(8).roll];
handles.angle.shoulder = [m(7).azim m(7).elev m(7).roll];
handles.angle.elbow =    [m(6).azim m(6).elev m(6).roll];
handles.angle.wrist =    [m(5).azim m(5).elev m(5).roll];
handles.angle.hand =     [m(4).azim m(4).elev m(4).roll];
handles.angle.middlefinger = [m(3).azim m(3).elev m(3).roll];
handles.angle.indexfinger = [m(2).azim m(2).elev m(2).roll];
handles.angle.thumbfinger = [m(1).azim m(1).elev m(1).roll];

%estimate the lengths of the fingers when they are fully extended
handles.finglen.middlefinger = handles.joint.middlefinger - handles.joint.hand;
handles.finglen.middlefinger = sqrt( sum(handles.finglen.middlefinger.^2,2) );
handles.maxlen.middlefinger = max(handles.finglen.middlefinger);
handles.finglen.indexfinger = handles.joint.indexfinger - handles.joint.hand;
handles.finglen.indexfinger = sqrt( sum(handles.finglen.indexfinger.^2,2) );
handles.maxlen.indexfinger = max(handles.finglen.indexfinger);
handles.finglen.thumbfinger = handles.joint.thumbfinger - handles.joint.wrist;
handles.finglen.thumbfinger = sqrt( sum(handles.finglen.thumbfinger.^2,2) );
handles.maxlen.thumbfinger = max(handles.finglen.thumbfinger);

% handles.fingratios.middlefinger = [1 1.5 2.6];
% handles.fingratios.middlefinger = handles.fingratios.middlefinger/sum(handles.fingratios.middlefinger);
% handles.fingratios.indexfinger = [1 1.4 2.5];
% handles.fingratios.indexfinger = handles.fingratios.indexfinger/sum(handles.fingratios.indexfinger);
% handles.fingratios.thumbfinger = [1 1.5 2.1];
% handles.fingratios.thumbfinger = handles.fingratios.thumbfinger/sum(handles.fingratios.thumbfinger);
% %finger length proportions taken from Alexander B and Viktor K "Proportions
% % of Hand Segments" Int J Morphol 28(3) 755-758, 2010.

%estimate finger angles - for computational ease, assume the same angle at
%  every finger joint and that every finger joint is equal length

handles.fingangs.middlefinger = -asind(0.5 - 0.5*3*handles.finglen.middlefinger/handles.maxlen.middlefinger) + 90;
handles.fingangs.indexfinger = -asind(0.5 - 0.5*3*handles.finglen.indexfinger/handles.maxlen.indexfinger) + 90;
handles.fingangs.thumbfinger = -asind(0.5 - 0.5*3*handles.finglen.thumbfinger/handles.maxlen.thumbfinger) + 90;

xmin = min([handles.joint.thumbfinger(:,1); handles.joint.indexfinger(:,1); handles.joint.wrist(:,1); handles.joint.elbow(:,1); handles.joint.shoulder(:,1); handles.joint.refshoulder(:,1)]);
xmax = max([handles.joint.thumbfinger(:,1); handles.joint.indexfinger(:,1); handles.joint.wrist(:,1); handles.joint.elbow(:,1); handles.joint.shoulder(:,1); handles.joint.refshoulder(:,1)]);
ymin = min([handles.joint.thumbfinger(:,2); handles.joint.indexfinger(:,2); handles.joint.wrist(:,2); handles.joint.elbow(:,2); handles.joint.shoulder(:,2); handles.joint.refshoulder(:,2)]);
ymax = max([handles.joint.thumbfinger(:,2); handles.joint.indexfinger(:,2); handles.joint.wrist(:,2); handles.joint.elbow(:,2); handles.joint.shoulder(:,2); handles.joint.refshoulder(:,2)]);
zmin = min([handles.joint.thumbfinger(:,3); handles.joint.indexfinger(:,3); handles.joint.wrist(:,3); handles.joint.elbow(:,3); handles.joint.shoulder(:,3); handles.joint.refshoulder(:,3)]);
zmax = max([handles.joint.thumbfinger(:,3); handles.joint.indexfinger(:,3); handles.joint.wrist(:,3); handles.joint.elbow(:,3); handles.joint.shoulder(:,3); handles.joint.refshoulder(:,3)]);

handles.plot.xmin = floor(xmin/5)*5;
handles.plot.xmax = ceil(xmax/5)*5;
handles.plot.ymin = floor(ymin/5)*5;
handles.plot.ymax = ceil(ymax/5)*5;
handles.plot.zmin = floor(zmin/5)*5;
handles.plot.zmax = ceil(zmax/5)*5;
if( handles.plot.zmax < (max([handles.joint.shoulder(:,3); handles.joint.refshoulder(:,3)]) + abs(handles.joint.refshoulder(1,1)-handles.joint.shoulder(1,1))/4) )
    handles.plot.zmax = (max([handles.joint.shoulder(:,3); handles.joint.refshoulder(:,3)]) + abs(handles.joint.refshoulder(1,1)-handles.joint.shoulder(1,1))/2*1.5);
    %handles.plot.zmax = ceil(zmax/5)*5;
end

% handles.plot.xmin = -80;
% handles.plot.xmax = 80;
% handles.plot.ymin = -70;
% handles.plot.ymax = 70;
% handles.plot.zmin = -60;
% handles.plot.zmax = 60;

%handles.viewaz = -25;
%handles.viewel = 25;
views.az = 195;
views.el = 20;
views.allowupdate = 0;
setappdata(handles.figure1,'views',views);

set(handles.slider1,'Min',1);
set(handles.slider1,'Max',size(handles.joint.shoulder,1));
set(handles.slider1,'SliderStep',[1/size(handles.joint.shoulder,1) , 10/size(handles.joint.shoulder,1)]);

%handles.c = 1;
c = 1;
setappdata(handles.figure1,'c',c);

%handles.inds = [];
%inds = [];
setappdata(handles.figure1,'inds',inds);

%sets up the play timer. Change the Period to speed up/slow down playback
playtimer = timer('Name','DataPlayer','Period',0.5,'StartDelay',0.1,'ExecutionMode','fixedRate','TimerFcn',{@timerCallbackFcn,handles});
stop(playtimer);
setappdata(handles.figure1,'playtimer',playtimer);

updateplot(handles);

views.allowupdate = 1;
setappdata(handles.figure1,'views',views);

MarkType = 1;%Set to 1 for full action, set to 2 for best subaction
setappdata(handles.figure1,'MarkType',MarkType);



% Choose default command line output for markdataGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes markdataGUI wait for user response (see UIRESUME)
uiwait(handles.figure1);


% function [] = timerCallback(~,~,guiHandle)
% 
% if ~isempty(guiHandle)
%     
%     % get the handles for the GUI/figure
%     handles = guidata(guiHandle);
% 
% end


% --- Outputs from this function are returned to the command line.
function varargout = markdataGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;
inds = getappdata(handles.figure1,'inds');
varargout{1} = inds;

varargout{2} = handles;

delete(handles.figure1);

% --- Executes on button press in pushbutton1 (play).
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pushbutton4,'enable','off');
set(handles.pushbutton5,'enable','off');
set(handles.pushbutton6,'enable','off');

playtimer = getappdata(handles.figure1,'playtimer');
pl = get(playtimer);
if strcmpi(pl.Running,'off')  %only ask to start the timer if it was off before.
    start(playtimer);
end


% --- Executes on button press in pushbutton2 (stop).
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
playtimer = getappdata(handles.figure1,'playtimer');
stop(playtimer);

set(handles.pushbutton4,'enable','on');
set(handles.pushbutton5,'enable','on');
set(handles.pushbutton6,'enable','on');



% --- Executes on button press in pushbutton4 (mark).
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
inds = getappdata(handles.figure1,'inds');
c = getappdata(handles.figure1,'c');
MarkType = getappdata(handles.figure1,'MarkType');

%inds = [inds; c];
%inds = unique(inds); %get rid of duplicate marks
%inds = sort(inds);   %resort the indices in ascending order
inds{MarkType} = [inds{MarkType}; c];
inds{MarkType} = unique(inds{MarkType}); %get rid of duplicate marks
inds{MarkType} = sort(inds{MarkType});   %resort the indices in ascending order

setappdata(handles.figure1,'inds',inds);

updateplot(handles);


% --- Executes on button press in pushbutton5 (unmark).
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

inds = getappdata(handles.figure1,'inds');
c = getappdata(handles.figure1,'c');
MarkType = getappdata(handles.figure1,'MarkType');

if ~isempty(inds{MarkType})
    diffinds = inds{MarkType}-c;
    
    [~,idiff] = min(abs(diffinds));
    inds{MarkType}(idiff) = [];
end

setappdata(handles.figure1,'inds',inds);

updateplot(handles);



% --- Executes on button press in pushbutton6 (snap to mark).
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
c = getappdata(handles.figure1,'c');
inds = getappdata(handles.figure1,'inds');

playtimer = getappdata(handles.figure1,'playtimer');
stop(playtimer);

set(handles.pushbutton4,'enable','on');
set(handles.pushbutton5,'enable','on');


allinds = [inds{1}; inds{2}];
if ~isempty(allinds)
    %allinds = sort(allinds);
    [~,nearestind] = min(abs(allinds-c));
    c = allinds(nearestind);
    
    setappdata(handles.figure1,'c',c);
    
    set(handles.slider1,'Value',c);
    
    updateplot(handles);
    
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

inds = getappdata(handles.figure1,'inds');
MarkType = getappdata(handles.figure1,'MarkType');

inds{MarkType} = [];

setappdata(handles.figure1,'inds',inds);

updateplot(handles);





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

set(handles.pushbutton4,'enable','on');
set(handles.pushbutton5,'enable','on');
set(handles.pushbutton6,'enable','on');

c = get(handles.slider1,'Value');
c = round(c);
%cmin = get(handles.slider1,'Min');
%cmax = get(handles.slider1,'Max');
%c = round((c-cmin)/(cmax-cmin) * size(handles.joint.shoulder,1));

if c < 1
    c = 1;
end
if c > size(handles.joint.shoulder,1)
    c = size(handles.joint.shoulder,1);
end

setappdata(handles.figure1,'c',c);

updateplot(handles);


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


% --- Executes on mouse press over axes background.
function axes5_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
doClickCallback(handles,handles.axes5)


% --- Executes on mouse press over axes background.
function axes6_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
doClickCallback(handles,handles.axes6)




% --- create new function to update plot
function updateplot(handles)

c = getappdata(handles.figure1,'c');
inds = getappdata(handles.figure1,'inds');
views = getappdata(handles.figure1,'views');



%compute finger segment coordinates
ssc = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
RU = @(A,B) eye(3) + ssc(cross(A,B)) + ssc(cross(A,B))^2*(1-dot(A,B))/(norm(cross(A,B))^2); %function to align the "standard" vector A to the reference hand-finger vector b
RA = @(A,b,gamma) eye(3) + sind(gamma)*ssc(b) + 2*((sind(gamma/2))^2)*(ssc(b))^2; %funtion to rotate a vector A around a reference vector b by an angle gamma;
PlaneAng = @(A1,A2,B1,B2) acosd(dot(cross(A1,A2),cross(B1,B2))/(norm(cross(A1,A2))*norm(cross(B1,B2))));
REuler = @(e,a,r) [  cosd(e)*cosd(a)                              cosd(e)*sind(a)                             -sind(e); 
                   -(cosd(r)*sind(a))+(sind(r)*sind(e)*cosd(a))  (cosd(r)*cosd(a))+(sind(r)*sind(e)*sind(a)) sind(r)*cosd(e); 
                    (sind(r)*sind(a))+(cosd(r)*sind(e)*cosd(a)) -(sind(r)*cosd(a))+(cosd(r)*sind(e)*sind(a)) cosd(r)*cosd(e)];

Raxisswitch = REuler(90,180,0);
                
basev = [-1 0 0];

mfv = [handles.joint.middlefinger(c,:)-handles.joint.hand(c,:)];
mfv = mfv/norm(mfv); %normalize the hand-finger vector
RotMat = RU(basev,mfv); %rotation matrix between (-1 0 0) and hand-finger vectors
mfbeta = handles.fingangs.middlefinger(c);
mfd = handles.finglen.middlefinger(c);
mflen = handles.maxlen.middlefinger;
mf1 = [-(mflen/3)*sind(mfbeta-90) 0 (mflen/3)*cosd(mfbeta-90) ];
mf2 = [-(mflen/3)*(sind(mfbeta-90)+1) 0 (mflen/3)*cosd(mfbeta-90)];
mf3 = [-(mflen/3)*(2*sind(mfbeta-90)+1) 0 0];
mf1 = RotMat*mf1';
mf2 = RotMat*mf2';
mf3 = RotMat*mf3'; %sanity check: this should be equal to mfv...
mfall = [mf1'; mf2'; mf3'];
% %estimate the angle between the finger plane and the hand plane, and rotate if necessary
% mfh = [handles.joint.middlefinger(c,:)-handles.joint.hand(c,:)];
% mwh = [handles.joint.wrist(c,:)-handles.joint.hand(c,:)];
% finghandang = PlaneAng(mf1,mf3,mfh,mwh)-90
% %if (abs(finghandang) > 90)
%     mf1 = RA(mf1,mfv,finghandang)*mf1;
%     mf2 = RA(mf2,mfv,finghandang)*mf2;
%     mf3 = RA(mf3,mfv,finghandang)*mf3;
% %end
%rotate the plane of the generated finger to align with the XZ plane of the finger sensor.
%mfang = [handles.angle.middlefinger(c,1)-handles.angle.middlefinger(1,1),handles.angle.middlefinger(c,2)-handles.angle.middlefinger(1,2),handles.angle.middlefinger(c,3)-handles.angle.middlefinger(1,3)];
mfang = [handles.angle.middlefinger(c,1),handles.angle.middlefinger(c,2),handles.angle.middlefinger(c,3)];
mfang = mfang/norm(mfang);
RotMatE = REuler(mfang(1),mfang(2),mfang(3));
%sensorang = RotMatE*[0 1 0]';
%fingang = cross(mf1,mf2);
fingang = PlaneAng(mf1/norm(mf1),mf2/norm(mf2),RotMatE*[0 -1 0]',RotMatE*[0 0 1]');
mf1 = RA(mf1/norm(mf1),mfv/norm(mfv),fingang)*mf1;
mf2 = RA(mf2/norm(mf2),mfv/norm(mfv),fingang)*mf2;
mf1 = RotMatE*mf1;
mf2 = RotMatE*mf2;
mf1 = mf1' + handles.joint.hand(c,:); %translate the origin back to the hand
mf2 = mf2' + handles.joint.hand(c,:); 
mf3 = mf3' + handles.joint.hand(c,:); %sanity check: this should be the middle finger marker pos
mfalltrans = [mf1; mf2; mf3];

ifv = [handles.joint.indexfinger(c,:)-handles.joint.hand(c,:)];
ifv = ifv/norm(ifv); %normalize the hand-finger vector
RotMat = RU(basev,ifv); %rotation matrix between (-1 0 0) and hand-finger vectors
ifbeta = handles.fingangs.indexfinger(c);
ifd = handles.finglen.indexfinger(c);
iflen = handles.maxlen.indexfinger;
if1 = [-(iflen/3)*sind(ifbeta-90) 0 (iflen/3)*cosd(ifbeta-90)];
if2 = [-(iflen/3)*(sind(ifbeta-90)+1) 0 (iflen/3)*cosd(ifbeta-90)];
if3 = [-(iflen/3)*(2*sind(ifbeta-90)+1) 0 0];
if1 = RotMat*if1';
if2 = RotMat*if2';
if3 = RotMat*if3'; %sanity check: this should be equal to mfv...
ifall = [if1; if2; if3];
ifang = [handles.angle.indexfinger(c,1),handles.angle.indexfinger(c,2),handles.angle.indexfinger(c,3)];
ifang = ifang/norm(ifang);
RotMatE = REuler(ifang(1),ifang(2),ifang(3));
fingang = PlaneAng(if1/norm(if1),if2/norm(if2),RotMatE*[0 -1 0]',RotMatE*[0 0 1]');
if1 = RA(if1/norm(if1),ifv/norm(ifv),fingang)*if1;
if2 = RA(if2/norm(if2),ifv/norm(ifv),fingang)*if2;
if1 = RotMatE*if1;
if2 = RotMatE*if2;
if1 = if1' + handles.joint.hand(c,:); %translate the origin back to the hand
if2 = if2' + handles.joint.hand(c,:); 
if3 = if3' + handles.joint.hand(c,:); %sanity check: this should be the middle finger marker pos
ifalltrans = [if1; if2; if3];


% %estimate the angle between the fingers and re-align then if needed
% fingfingang = PlaneAng(mfall(1,:),mfall(3,:),ifall(1,:),ifall(3,:))
% if (fingfingang > 90)
%     %we will assume the middle finger is off...
%     

tfv = [handles.joint.thumbfinger(c,:)-handles.joint.wrist(c,:)];
tfv = tfv/norm(tfv); %normalize the hand-finger vector
RotMat = RU(basev,tfv); %rotation matrix between (-1 0 0) and hand-finger vectors
tfbeta = handles.fingangs.thumbfinger(c);
tfd = handles.finglen.thumbfinger(c);
tflen = handles.maxlen.thumbfinger;
tf1 = [-(tflen/3)*sind(tfbeta-90) 0 (tflen/3)*cosd(tfbeta-90)];
tf2 = [-(tflen/3)*(sind(tfbeta-90)+1) 0 (tflen/3)*cosd(tfbeta-90)];
tf3 = [-(tflen/3)*(2*sind(tfbeta-90)+1) 0 0];
tf1 = RotMat*tf1';
tf2 = RotMat*tf2';
tf3 = RotMat*tf3'; %sanity check: this should be equal to mfv...
tfang = [handles.angle.thumbfinger(c,1),handles.angle.thumbfinger(c,2),handles.angle.thumbfinger(c,3)];
tfang = tfang/norm(tfang);
RotMatE = REuler(tfang(1),tfang(2),tfang(3));
fingang = PlaneAng(tf1/norm(tf1),tf2/norm(tf2),RotMatE*[0 1 0]',RotMatE*[0 0 -1]');
tf1 = RA(tf1/norm(tf1),tfv/norm(tfv),fingang)*tf1;
tf2 = RA(tf2/norm(tf2),tfv/norm(tfv),fingang)*tf2;
tf1 = RotMatE*tf1;
tf2 = RotMatE*tf2;
tf1 = tf1' + handles.joint.wrist(c,:); %translate the origin back to the hand
tf2 = tf2' + handles.joint.wrist(c,:); 
tf3 = tf3' + handles.joint.wrist(c,:); %sanity check: this should be the middle finger marker pos
tfall = [tf1; tf2; tf3];


armrad = abs(handles.joint.refshoulder(1,1)-handles.joint.shoulder(1,1))/16;

if handles.nmarkers == 6
    %do simplified arm without wrist/hand
    [xSE,ySE,zSE] = cylinder2P(armrad,8,handles.joint.shoulder(c,:),handles.joint.elbow(c,:));
    [xEW,yEW,zEW] = cylinder2P(armrad*.75,8,handles.joint.elbow(c,:),handles.joint.wrist(c,:));
    [xWF,yWF,zWF] = cylinder2P(armrad*[.25 .25 .25 .1 .1 .1],8,handles.joint.wrist(c,:),handles.joint.indexfinger(c,:));
    [xWT,yWT,zWT] = cylinder2P(armrad*[.25 .25 .25 .1 .1 .1],8,handles.joint.wrist(c,:),handles.joint.thumbfinger(c,:));
elseif handles.nmarkers == 8
    [xSE,ySE,zSE] = cylinder2P(armrad,8,handles.joint.shoulder(c,:),handles.joint.elbow(c,:));
    [xEW,yEW,zEW] = cylinder2P(armrad*.75,8,handles.joint.elbow(c,:),handles.joint.wrist(c,:));
    [xWH,yWH,zWH] = cylinder2P(armrad*.75,8,handles.joint.wrist(c,:),handles.joint.hand(c,:));
    [xHM1,yHM1,zHM1] = cylinder2P(armrad*[.3 .2],6,handles.joint.hand(c,:),mf1);
    [xHM2,yHM2,zHM2] = cylinder2P(armrad*[.2 .2],6,mf1,mf2);
    [xHM3,yHM3,zHM3] = cylinder2P(armrad*[.2 .2],6,mf2,mf3);
    [xHI1,yHI1,zHI1] = cylinder2P(armrad*[.3 .2],6,handles.joint.hand(c,:),if1);
    [xHI2,yHI2,zHI2] = cylinder2P(armrad*[.2 .2],6,if1,if2);
    [xHI3,yHI3,zHI3] = cylinder2P(armrad*[.2 .2],6,if2,if3);
    [xWT1,yWT1,zWT1] = cylinder2P(armrad*[.3 .2],6,handles.joint.wrist(c,:),tf1);
    [xWT2,yWT2,zWT2] = cylinder2P(armrad*[.2 .2],6,tf1,tf2);
    [xWT3,yWT3,zWT3] = cylinder2P(armrad*[.2 .2],6,tf2,tf3);
    
    %[xHM,yHM,zHM] = cylinder2P([.25 .25 .25 .1 .1 .1],8,handles.joint.hand(c,:),handles.joint.middlefinger(c,:));
    %[xHM,yHM,zHM] = cylinder2P(.1,8,handles.joint.hand(c,:),handles.joint.middlefinger(c,:));
    %[xWF,yWF,zWF] = cylinder2P([.25 .25 .25 .1 .1 .1],8,handles.joint.hand(c,:),handles.joint.indexfinger(c,:));
    %[xWT,yWT,zWT] = cylinder2P([.25 .25 .25 .1 .1 .1],8,handles.joint.hand(c,:),handles.joint.thumbfinger(c,:));
end




%rectangle representing the body, assuming a vertical plane
torsolength = abs(handles.joint.shoulder(1,1)-handles.joint.refshoulder(1,1))*1.5;
xtorso = [handles.joint.shoulder(c,1); handles.joint.refshoulder(c,1); handles.joint.refshoulder(c,1); handles.joint.shoulder(c,1)];
ytorso = [handles.joint.shoulder(c,2); handles.joint.refshoulder(c,2); handles.joint.refshoulder(c,2); handles.joint.shoulder(c,2)];
ztorso = [handles.joint.shoulder(c,3); handles.joint.refshoulder(c,3); handles.joint.refshoulder(c,3)-torsolength; handles.joint.shoulder(c,3)-torsolength];

%headrad = 9;
headrad = abs(handles.joint.refshoulder(1,1)-handles.joint.shoulder(1,1))/4;
headc = handles.torso(c,:) + [0 0 headrad];
theta = [0:pi/10:2*pi];
headx = headrad*cos(theta);
heady = zeros(size(headx));
headz = headrad*sin(theta)+headrad;

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
set(gca,'xlim',[handles.plot.xmin handles.plot.xmax],'ylim',[handles.plot.ymin handles.plot.ymax],'zlim',[handles.plot.zmin handles.plot.zmax]);
view(views.az,views.el);
h = patch('XData',xtorso,'YData',ytorso,'ZData',ztorso);
set(h,'FaceColor',[.5 .5 .5]);
hold on;
h = patch('XData',headx,'YData',heady,'ZData',headz);
set(h,'FaceColor',[.5 .5 .5]);
h = surf(xSE, ySE, zSE);
set(h,'FaceColor',[1 0 0]);
h = surf(xEW, yEW, zEW);
set(h,'FaceColor',[1 1 0]);
if handles.nmarkers == 8
    h = surf(real(xWH), real(yWH), real(zWH));
    set(h,'FaceColor',[1 1 0]);
    h = surf(real(xHM1), real(yHM1), real(zHM1));
    set(h,'FaceColor',[0 0 1]);
    h = surf(real(xHM2), real(yHM2), real(zHM2));
    set(h,'FaceColor',[0 0 1]);
    h = surf(real(xHM3), real(yHM3), real(zHM3));
    set(h,'FaceColor',[0 0 1]);
    h = surf(real(xHI1), real(yHI1), real(zHI1));
    set(h,'FaceColor',[0 .5 1]);
    h = surf(real(xHI2), real(yHI2), real(zHI2));
    set(h,'FaceColor',[0 .5 1]);
    h = surf(real(xHI3), real(yHI3), real(zHI3));
    set(h,'FaceColor',[0 .5 1]);
    h = surf(real(xWT1), real(yWT1), real(zWT1));
    set(h,'FaceColor',[1 0 1]);
    h = surf(real(xWT2), real(yWT2), real(zWT2));
    set(h,'FaceColor',[1 0 1]);
    h = surf(real(xWT3), real(yWT3), real(zWT3));
    set(h,'FaceColor',[1 0 1]);
else
    h = surf(real(xWF), real(yWF), real(zWF));
    set(h,'FaceColor',[0 0 1]);
    h = surf(real(xWT), real(yWT), real(zWT));
    set(h,'FaceColor',[0 0 1]);
end
hold off;
grid on;
%set(gca,'xlim',[handles.plot.xmin handles.plot.xmax],'ylim',[handles.plot.ymin handles.plot.ymax],'zlim',[handles.plot.zmin handles.plot.zmax]);
%view(handles.viewaz,handles.viewel);
xlabel('x');
ylabel('y');
zlabel('z');

if isempty(handles.ctp)
    
    axes(handles.axes2)
    cla(handles.axes2);
    set(handles.axes2,'NextPlot','add');
    h = plot(handles.joint.elbow,'-');
    set(h,'HitTest','off')
    hold on
    for a = 1:length(inds{1})
        h = plot([inds{1}(a) inds{1}(a)],get(gca,'ylim'),'r-');
        set(h,'HitTest','off')
    end
    for a = 1:length(inds{2})
        if mod(a,2) ~= 0
            h = plot([inds{2}(a) inds{2}(a)],get(gca,'ylim'),'b-.');
        else
            h = plot([inds{2}(a) inds{2}(a)],get(gca,'ylim'),'b--');
            set(h,'Color',[0 0 .5]);
        end
        set(h,'HitTest','off')
    end
    h = plot([c c],get(gca,'ylim'),'k:');
    set(h,'HitTest','off','LineWidth',1.5)
    hold off;
    ylabel('elbow')
    
    axes(handles.axes3)
    cla(handles.axes3);
    set(handles.axes3,'NextPlot','add');
    h = plot(handles.joint.wrist,'-');
    set(h,'HitTest','off');
    hold on
    for a = 1:length(inds{1})
        h = plot([inds{1}(a) inds{1}(a)],get(gca,'ylim'),'r-');
        set(h,'HitTest','off')
    end
    for a = 1:length(inds{2})
        %h = plot([inds{2}(a) inds{2}(a)],get(gca,'ylim'),'b-.');
        if mod(a,2) ~= 0
            h = plot([inds{2}(a) inds{2}(a)],get(gca,'ylim'),'b-.');
        else
            h = plot([inds{2}(a) inds{2}(a)],get(gca,'ylim'),'b--');
            set(h,'Color',[0 0 .5]);
        end
        set(h,'HitTest','off')
    end
    h = plot([c c],get(gca,'ylim'),'k:');
    set(h,'HitTest','off','LineWidth',1.5);
    hold off;
    ylabel('wrist')
    
    axes(handles.axes4)
    cla(handles.axes4);
    set(handles.axes4,'NextPlot','add');
    h = plot(handles.joint.hand,'-');
    set(h,'HitTest','off');
    hold on
    for a = 1:length(inds{1})
        h = plot([inds{1}(a) inds{1}(a)],get(gca,'ylim'),'r-');
        set(h,'HitTest','off')
    end
    for a = 1:length(inds{2})
        if mod(a,2) ~= 0
            h = plot([inds{2}(a) inds{2}(a)],get(gca,'ylim'),'b-.');
        else
            h = plot([inds{2}(a) inds{2}(a)],get(gca,'ylim'),'b--');
            set(h,'Color',[0 0 .5]);
        end
        set(h,'HitTest','off')
    end
    h = plot([c c],get(gca,'ylim'),'k:');
    set(h,'HitTest','off','LineWidth',1.5);
    hold off;
    ylabel('hand')
    
    axes(handles.axes5)
    cla(handles.axes5);
    set(handles.axes5,'NextPlot','add');
    h = plot(handles.joint.indexfinger,'-');
    set(h,'HitTest','off');
    hold on
    for a = 1:length(inds{1})
        h = plot([inds{1}(a) inds{1}(a)],get(gca,'ylim'),'r-');
        set(h,'HitTest','off')
    end
    for a = 1:length(inds{2})
        if mod(a,2) ~= 0
            h = plot([inds{2}(a) inds{2}(a)],get(gca,'ylim'),'b-.');
        else
            h = plot([inds{2}(a) inds{2}(a)],get(gca,'ylim'),'b--');
            set(h,'Color',[0 0 .5]);
        end
        set(h,'HitTest','off')
    end
    h = plot([c c],get(gca,'ylim'),'k:');
    set(h,'HitTest','off','LineWidth',1.5);
    hold off;
    ylabel('index')
    
    axes(handles.axes6)
    cla(handles.axes6);
    set(handles.axes6,'NextPlot','add');
    h = plot(handles.joint.thumbfinger,'-');
    set(h,'HitTest','off');
    hold on
    for a = 1:length(inds{1})
        h = plot([inds{1}(a) inds{1}(a)],get(gca,'ylim'),'r-');
        set(h,'HitTest','off')
    end
    for a = 1:length(inds{2})
        if mod(a,2) ~= 0
            h = plot([inds{2}(a) inds{2}(a)],get(gca,'ylim'),'b-.');
        else
            h = plot([inds{2}(a) inds{2}(a)],get(gca,'ylim'),'b--');
            set(h,'Color',[0 0 .5]);
        end
        set(h,'HitTest','off')
    end
    h = plot([c c],get(gca,'ylim'),'k:');
    set(h,'HitTest','off','LineWidth',1.5);
    hold off;
    ylabel('thumb')
    
else %plot transformed arm data instead of raw data
    
    axes(handles.axes2)
    cla(handles.axes2);
    set(handles.axes2,'NextPlot','add');
    h = plot(handles.ctp(:,5:7),'-');
    set(h,'HitTest','off');
    hold on
    for a = 1:length(inds{1})
        h = plot([inds{1}(a) inds{1}(a)],get(gca,'ylim'),'r-');
        set(h,'HitTest','off')
    end
    for a = 1:length(inds{2})
        h = plot([inds{2}(a) inds{2}(a)],get(gca,'ylim'),'b-.');
        set(h,'HitTest','off')
    end
    h = plot([c c],get(gca,'ylim'),'k:');
    set(h,'HitTest','off','LineWidth',1.5);
    hold off;
    ylabel('arm plane (phi,theta,psi)')
    
    axes(handles.axes3)
    cla(handles.axes3);
    set(handles.axes3,'NextPlot','add');
    h = plot(handles.ctp(:,8),'-');
    set(h,'HitTest','off');
    hold on
    for a = 1:length(inds{1})
        h = plot([inds{1}(a) inds{1}(a)],get(gca,'ylim'),'r-');
        set(h,'HitTest','off')
    end
    for a = 1:length(inds{2})
        h = plot([inds{2}(a) inds{2}(a)],get(gca,'ylim'),'b-.');
        set(h,'HitTest','off')
    end
    h = plot([c c],get(gca,'ylim'),'k:');
    set(h,'HitTest','off','LineWidth',1.5);
    hold off;
    ylabel('elbow ang')
    
    axes(handles.axes4)
    cla(handles.axes4);
    set(handles.axes4,'NextPlot','add');
    h = plot(handles.ctp(:,9),'-');
    set(h,'HitTest','off');
    hold on
    for a = 1:length(inds{1})
        h = plot([inds{1}(a) inds{1}(a)],get(gca,'ylim'),'r-');
        set(h,'HitTest','off')
    end
    for a = 1:length(inds{2})
        h = plot([inds{2}(a) inds{2}(a)],get(gca,'ylim'),'b-.');
        set(h,'HitTest','off')
    end
    h = plot([c c],get(gca,'ylim'),'k:');
    set(h,'HitTest','off','LineWidth',1.5);
    hold off;
    ylabel('wrist ang')
    
    
    axes(handles.axes5)
    cla(handles.axes5);
    set(handles.axes5,'NextPlot','add');
    h = plot(handles.ctp(:,3:4),'-');
    set(h,'HitTest','off');
    hold on
    for a = 1:length(inds{1})
        h = plot([inds{1}(a) inds{1}(a)],get(gca,'ylim'),'r-');
        set(h,'HitTest','off')
    end
    for a = 1:length(inds{2})
        h = plot([inds{2}(a) inds{2}(a)],get(gca,'ylim'),'b-.');
        set(h,'HitTest','off')
    end
    h = plot([c c],get(gca,'ylim'),'k:');
    set(h,'HitTest','off','LineWidth',1.5);
    hold off;
    ylabel('index curv/tors')
    
    axes(handles.axes6)
    cla(handles.axes6);
    set(handles.axes6,'NextPlot','add');
    h = plot(handles.ctp(:,1:2),'-');
    set(h,'HitTest','off');
    hold on
    for a = 1:length(inds{1})
        h = plot([inds{1}(a) inds{1}(a)],get(gca,'ylim'),'r-');
        set(h,'HitTest','off')
    end
    for a = 1:length(inds{2})
        h = plot([inds{2}(a) inds{2}(a)],get(gca,'ylim'),'b-.');
        set(h,'HitTest','off')
    end
    h = plot([c c],get(gca,'ylim'),'k:');
    set(h,'HitTest','off','LineWidth',1.5);
    hold off;
    ylabel('thumb curv/tors')
    
end

set(handles.text1,'String',num2str(c));



function doClickCallback(handles,haxes)

playtimer = getappdata(handles.figure1,'playtimer');
stop(playtimer);

set(handles.pushbutton4,'enable','on');
set(handles.pushbutton5,'enable','on');
set(handles.pushbutton6,'enable','on');

coordinates = get(haxes,'CurrentPoint');
coords = coordinates(1,1:2);
c = round(coords(1));

%cmin = get(handles.slider1,'Min');
%cmax = get(handles.slider1,'Max');
%set(handles.slider1,'Value',(c/size(handles.joint.shoulder,1)*(cmax-cmin)+cmin));
set(handles.slider1,'Value',c);

setappdata(handles.figure1,'c',c);

updateplot(handles)


function timerCallbackFcn(obj, event, handles)
c = getappdata(handles.figure1,'c');
c = c+10;  %number of frames to step forward on each tick

if c > size(handles.joint.shoulder,1)
    c = 1;
    playtimer = getappdata(handles.figure1,'playtimer');
    stop(playtimer);
    
    set(handles.pushbutton4,'enable','on');
    set(handles.pushbutton5,'enable','on');
    set(handles.pushbutton6,'enable','on');
end

setappdata(handles.figure1,'c',c);

%cmin = get(handles.slider1,'Min');
%cmax = get(handles.slider1,'Max');
%set(handles.slider1,'Value',(c/size(handles.joint.shoulder,1)*(cmax-cmin)+cmin));
set(handles.slider1,'Value',c);

updateplot(handles);




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

% Hint: delete(hObject) closes the figure
%delete(hObject);


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MarkType = 1;%Set to 1 for full action
setappdata(handles.figure1,'MarkType',MarkType);
%updateplot(handles);

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MarkType = 2;%Set to 2 for best subaction
setappdata(handles.figure1,'MarkType',MarkType);

% Hint: get(hObject,'Value') returns toggle state of radiobutton2



% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)

c = getappdata(handles.figure1,'c');

c = c + eventdata.VerticalScrollCount;

if c < 1
    c = 1;
end
if c > size(handles.joint.shoulder,1)
    c = size(handles.joint.shoulder,1);
end


setappdata(handles.figure1,'c',c);
    
set(handles.slider1,'Value',c);
    
updateplot(handles);



