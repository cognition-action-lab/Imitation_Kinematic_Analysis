
playtimer = timer('Name','DataPlayer','Period',0.25,'StartDelay',0.1,'ExecutionMode','fixedRate','TimerFcn',{@timerCallbackFcn});
stop(playtimer);




function timerCallbackFcn(obj, event)
c = getappdata(handles.figure1,'c');
c = c+10;

if c > size(handles.joint.shoulder,1)
    c = 1;
    playtimer = getappdata(handles.figure1,'playtimer');
    stop(playtimer);
    
    set(handles.pushbutton4,'enable','on');
    set(handles.pushbutton5,'enable','on');
end

setappdata(handles.figure1,'c',c);

%cmin = get(handles.slider1,'Min');
%cmax = get(handles.slider1,'Max');
%set(handles.slider1,'Value',(c/size(handles.joint.shoulder,1)*(cmax-cmin)+cmin));
set(handles.slider1,'Value',c);

updateplot(handles);

end