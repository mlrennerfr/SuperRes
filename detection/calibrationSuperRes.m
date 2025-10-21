function varargout = calibrationSuperRes(varargin)
%CALIBRATIONSUPERRES MATLAB code file for calibrationSuperRes.fig
% Calibration GUI for SuperRes programs
%
% Marianne Renner feb 2025 version for SuperRes_v4_app
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last Modified by GUIDE v2.5 21-Feb-2025 18:09:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @calibrationSuperRes_OpeningFcn, ...
                   'gui_OutputFcn',  @calibrationSuperRes_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calibrationSuperRes_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
set(handles.output,'userdata',varargin{1}(1));
set(handles.quit,'value',1);
guidata(hObject, handles);

%-------------------------------------------------------------------------
function varargout = calibrationSuperRes_OutputFcn(hObject, eventdata, handles)

varargout{1} = handles.output;

detoptions=get(handles.output,'userdata');

set(handles.text2,'Userdata',detoptions.file); %thresold text (just a variable)
lastframe=detoptions.lastframe;
set(handles.text8,'value',lastframe); %nro frames

set(handles.thresh,'String',num2str(detoptions.seuil_detec_1vue));
set(handles.windsize,'String',num2str(detoptions.wn));
set(handles.alpha,'String',num2str(detoptions.seuil_alpha));
set(handles.gaussrad,'String',num2str(detoptions.r0));
set(handles.defloops,'String',num2str(detoptions.nb_defl));

%options(1)=detoptions.seuil_alpha  ;
%options(2)=detoptions.seuil_detec_1vue;         % threshold
%options(3)=detoptions.wn;                       % window size
%options(4)=detoptions.r0;                       % gaussian radius
%options(5)=detoptions.nb_defl;                  % number deflation loops


axes(handles.axes1);
set(handles.factor,'value',1); %frame #1

detect_Callback(hObject, eventdata, handles)

guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function thresh_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function windsize_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function width_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function edit4_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function skip_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function factor_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function gaussrad_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function defloops_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function alpha_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function thresh_Callback(hObject, eventdata, handles)
handles.threshold=get(hObject,'String');
guidata(hObject, handles);

function windsize_Callback(hObject, eventdata, handles)
handles.maxintensity=get(hObject,'String');
guidata(hObject, handles);

function skip_Callback(hObject, eventdata, handles)
handles.skipframes=get(hObject,'String');
guidata(hObject, handles);

function factor_Callback(hObject, eventdata, handles)
handles.valfactor=str2num(get(hObject,'String'));
guidata(hObject, handles);

function defloops_Callback(hObject, eventdata, handles)
function alpha_Callback(hObject, eventdata, handles)
function gaussrad_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in detect.
function detect_Callback(hObject, eventdata, handles)

valfactor=str2num(get(handles.factor,'string'));
first=get(handles.factor,'value'); %frame #1

actualframe=get(handles.quit,'value');
lastframe=get(handles.text8,'value');
if actualframe<1
    actualframe=1;
elseif actualframe>lastframe
    actualframe=lastframe;
end

moviename=get(handles.text2,'Userdata'); 

 %   k=strfind(moviename,'nd2');
   % if k>0
    %    data = bfopen(moviename);
    %    stack = data{1, 1};
     %   Nb_image_ds_stack = size(stack, 1);
   % else
        stack=[];
     %   info=imfinfo(name_stk);
     %   Nb_image_ds_stack =length(info);
   % end


seuil_alpha = str2num(get(handles.alpha,'String'));
seuil_detec_1vue = str2num(get(handles.thresh,'String'));
wn = str2num(get(handles.windsize, 'String'));
r0 = str2num(get(handles.gaussrad, 'String'));
nb_defl = str2num(get(handles.defloops,'String'));
activation_sig_fit= 1;


[matrice_results, radiusmean,sommevalid] = detectionMTToneimageSR(moviename,stack, actualframe,1, seuil_detec_1vue, wn, r0, nb_defl, seuil_alpha,activation_sig_fit);  

if first==1;
    %first
    showimage(handles, stack, 1);
    set(handles.factor,'value',0); 
else
    showimage(handles, stack);
end

axes(handles.axes1);
set(handles.text8,'string',['Frame = ',num2str(actualframe),' (of ',num2str(lastframe),')']);

if isempty(matrice_results)==0
    peak=[matrice_results(1,:)' matrice_results(3,:)' matrice_results(2,:)'  matrice_results(4,:)'];
    plot (peak(:,2)+1, peak(:,3)+1,'o','markeredgecolor',[1 0 0],'markersize',8);
end

        
set(handles.gaussrad,'string',num2str(radiusmean));
set(handles.result,'string',['Valid particles :',num2str(sommevalid)]);



disp(' ')
disp(['Frame # ',num2str(actualframe)]);
disp(['Threshold alpha: ',num2str(seuil_alpha)]);
disp(['Threshold detection: ',num2str(seuil_detec_1vue)]);
disp(['Window size: ',num2str(wn)]);
if activation_sig_fit==1  
    disp(['Fit of gaussian free, radius obtained: ',num2str(radiusmean)])
else
    disp(['Radius of gaussian free: ',num2str(r0)])
end
disp(['Deflation loops: ',num2str(nb_defl)]);
disp(' ')
disp(['Valid particles: ',num2str(sommevalid)]);
disp(' ')


%clear all

guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in accept.
function accept_Callback(hObject, eventdata, handles)

%detoptions=get(handles.output,'userdata');
detoptions.seuil_alpha = str2num(get(handles.alpha,'String'));
detoptions.seuil_detec_1vue = str2num(get(handles.thresh,'String'));
detoptions.wn = str2num(get(handles.windsize, 'String'));
detoptions.r0 = str2num(get(handles.gaussrad, 'String'));
detoptions.nb_defl = str2num(get(handles.defloops,'String'));

detoptions.image=[];
detoptions.imagepar=[];

pathdet=['detecoptions.mat'];
save(pathdet,'detoptions','-mat');
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function showimage(handles, stack, first)

valfactor=str2num(get(handles.factor,'string'));
datamatrix=[];
detoptions=get(handles.output,'userdata');
actualframe=get(handles.quit,'value');
lastframe=get(handles.text8,'value');
if actualframe<1
    actualframe=1;
elseif actualframe>lastframe
    actualframe=lastframe;
end
%Image=detoptions.image;
%ImagePar=detoptions.imagepar;
name_stk=get(handles.text2,'userdata');

if isempty(stack)==1
  %  datamatrix = imread(name_stk,actualframe);
    datamatrix = imread_big(name_stk,actualframe);
    datamatrix = double(datamatrix);
end

%figure;
axes(handles.axes1);
        ax = gca;
        ax.Box = 'on';
        ax.BoxStyle = 'full';
        ax.XTickLabel = {''};
        ax.YTickLabel = {''};
       % hold on

set(handles.text8,'string',['Frame = ',num2str(actualframe),' (of ',num2str(lastframe),')']);
[Xdim,Ydim]=size(datamatrix);

prop=get(handles.text11,'value') ;  %correcion contraste
if prop==0
    prop=0.5;
end
stackmin=(min(min(min(datamatrix))));
stackmax=(max(max(max(datamatrix))));

if nargin>2 %first
   conval=1;
   valfactor=1000/(stackmax*0.5);
   val=num2str(valfactor);
   set(handles.factor,'string',val);
else
   valfactor=str2num(get(handles.factor,'string'));
   conval=round(prop*(stackmax*valfactor))/1000;
   if conval==0
      conval=0.00001;
   end
end
datamatrix=datamatrix*conval;

imshow(datamatrix,[stackmin stackmax],'InitialMagnification','fit');
axis([0,Ydim,0,Xdim]);
hold on
datamatrix=[];
guidata(gcbo,handles) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)

frames=get(handles.text8,'value');
slidervalue=get(hObject,'Value');
minvalue=get(hObject,'Min');
maxvalue=get(hObject,'Max');
prop=slidervalue/(minvalue+maxvalue);
frame=round(prop*frames);
if frame<1
    frame=1;
end
set(handles.quit,'value',frame);
moviename=get(handles.text2,'Userdata'); 

  %  k=strfind(moviename,'nd2');
 %   if k>0
  %      data = bfopen(moviename);
  %      stack = data{1, 1};
  %  else
        stack=[];
  %  end

showimage(handles, stack);

guidata(hObject, handles);

%--------------------------------------------------------------------------
function slider1_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
slidervalue=get(hObject,'Value');
minvalue=get(hObject,'Min');
maxvalue=get(hObject,'Max');

prop=slidervalue/(minvalue+maxvalue);
set(handles.text11,'value',prop);
valfactor=str2num(get(handles.factor,'string'));
moviename=get(handles.text2,'Userdata'); 

  %  k=strfind(moviename,'nd2');
 %   if k>0
  %      data = bfopen(moviename);
  %      stack = data{1, 1};
     %   Nb_image_ds_stack = size(stack, 1);
   % else
        stack=[];
     %   info=imfinfo(name_stk);
     %   Nb_image_ds_stack =length(info);
  %  end

showimage(handles, stack);
guidata(hObject, handles);

%-------------------------------------------------------------------------
function slider2_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in quit.
function quit_Callback(hObject, eventdata, handles)

%detoptions=get(handles.output,'userdata');
detoptions.seuil_alpha = str2num(get(handles.alpha,'String'));
detoptions.seuil_detec_1vue = str2num(get(handles.thresh,'String'));
detoptions.wn = str2num(get(handles.windsize, 'String'));
detoptions.r0 = str2num(get(handles.gaussrad, 'String'));
detoptions.nb_defl = str2num(get(handles.defloops,'String'));

detoptions.image=[];
detoptions.imagepar=[];

pathdet=['detecoptions.mat'];
save(pathdet,'detoptions','-mat');

close
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%end of file

