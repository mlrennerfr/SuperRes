function varargout = ROIshiftcorrection(varargin)
%
% GUI for color shift correction
%
% Marianne Renner, SuperRes progrms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last Modified by GUIDE v2.5 13-Mar-2014 15:48:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ROIshiftcorrection_OpeningFcn, ...
                   'gui_OutputFcn',  @ROIshiftcorrection_OutputFcn, ...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ROIshiftcorrection_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

output=[0 0];
set(handles.acceptpushbutton,'userdata',output);

save('auxiliar.mat','output','-mat')

set(handles.panel,'userdata',zeros(4,4)); %shifts

varargout{1} = [];
   
handles.file=cell2mat(varargin{1}(1)); % varargin{1}=namefile; %
handles.Xdim=cell2mat(varargin{1}(2)); % varargin{2}=dimx; %inverted
handles.Ydim=cell2mat(varargin{1}(3));   % varargin{3}=dimy;
handles.traces1=cell2mat(varargin{1}(4)); % varargin{4}=newpoints;
handles.traces2=cell2mat(varargin{1}(5));       % varargin{5}=newpoints2; %  
handles.newtraces=handles.traces2;

showframeshift(handles)

guidata(hObject,handles) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = ROIshiftcorrection_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function selectobjectpushbutton_Callback(hObject, eventdata, handles)


% ROI
axes(handles.axes1);
[areaselect,xi,yi]=roipolyold;    %seleccion ROI   

aux = find(inpolygon( handles.traces1(:,1), handles.traces1(:,2),xi,yi) );
newpoints(:,1)= handles.traces1(aux,1);
newpoints(:,2)= handles.traces1(aux,2);
aux = find(inpolygon( handles.newtraces(:,1), handles.newtraces(:,2),xi,yi)) ;
newpoints2(:,1)= handles.newtraces(aux,1);
newpoints2(:,2)= handles.newtraces(aux,2); 

if isempty(newpoints)==0 && isempty(newpoints2)==0 
    % center of points cloud  
    centerx1=mean(newpoints(:,1));
    centery1=mean(newpoints(:,2));
    centerx2=mean(newpoints2(:,1));
    centery2=mean(newpoints2(:,2));       
    difx=centerx1-centerx2;
    dify=centery1-centery2;
    handles.newtraces(:,1)=handles.newtraces(:,1)+difx;
    handles.newtraces(:,2)=handles.newtraces(:,2)+dify;   
    
    % plot
    showframeshift(handles,handles.newtraces);      
    
    ktev=str2num(get(handles.trajverticalshift,'String'));
    kteh=str2num(get(handles.trajhorizontalshift,'String'));
    output(1)=difx+kteh;
    output(2)=dify+ktev;
    set(handles.trajhorizontalshift,'String',num2str(output(1)));
    set(handles.trajverticalshift,'String',num2str(output(2)));
    set(handles.acceptpushbutton,'userdata',output)   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function acceptpushbutton_Callback(hObject, eventdata, handles)

output=get(handles.acceptpushbutton,'userdata');

save('auxiliar.mat','output','-mat')
close


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% correct shift
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trajUbutton_Callback(hObject, eventdata, handles)

status=get(handles.trajverticalshift,'enable');
output=get(handles.acceptpushbutton,'userdata');

if size(status,2)==2 %on
    
%UP
handles.thold=get(handles.trajverticalshift,'String');
ktev=str2num(handles.thold);
ktev=ktev-0.05;
set(handles.trajverticalshift,'string',num2str(ktev));
datashift=get(handles.panel,'userdata');
datashift(1,4)=ktev;
set(handles.panel,'userdata',datashift);

output(2)=ktev;
set(handles.acceptpushbutton,'userdata',output)
handles.newtraces=shiftpoints(handles,ktev,datashift(2,4),handles.traces2);

indexzero=find(handles.newtraces(:,2)>0);
handles.newtraces=handles.newtraces(indexzero,:);
%set(handles.finishedpushbutton,'value',1);

points2=handles.newtraces;
showframeshift(handles,points2)

end %status

guidata(hObject, handles);

%-----------------------------------------------------------------
function trajDbutton_Callback(hObject, eventdata, handles)

status=get(handles.trajverticalshift,'enable');
output=get(handles.acceptpushbutton,'userdata');

if size(status,2)==2 %on

%DOWN
handles.thold=get(handles.trajverticalshift,'String');
ktev=str2num(handles.thold);

ktev=ktev+0.05;
set(handles.trajverticalshift,'string',num2str(ktev));
datashift=get(handles.panel,'userdata');
datashift(1,4)=ktev;
set(handles.panel,'userdata',datashift);

output(2)=ktev;
set(handles.acceptpushbutton,'userdata',output)
handles.newtraces=shiftpoints(handles,ktev,datashift(2,4),handles.traces2);

indexzero=find(handles.newtraces(:,2)<handles.Ydim);
handles.newtraces=handles.newtraces(indexzero,:);
%set(handles.finishedpushbutton,'value',1);

points2=handles.newtraces;
showframeshift(handles,points2)

end

guidata(hObject, handles);

%-----------------------------------------------------------------
function trajLbutton_Callback(hObject, eventdata, handles)

status=get(handles.trajhorizontalshift,'enable');
output=get(handles.acceptpushbutton,'userdata');

if size(status,2)==2 %on

%LEFT
handles.thold=get(handles.trajhorizontalshift,'String');
kteh=str2num(handles.thold);
kteh=kteh-0.05;
set(handles.trajhorizontalshift,'string',num2str(kteh));
datashift=get(handles.panel,'userdata');
datashift(2,4)=kteh;
set(handles.panel,'userdata',datashift);

output(1)=kteh;
set(handles.acceptpushbutton,'userdata',output)

handles.newtraces=shiftpoints(handles,datashift(1,4),kteh,handles.traces2);
indexzero=find(handles.newtraces(:,2)>0);
handles.newtraces=handles.newtraces(indexzero,:);
%set(handles.finishedpushbutton,'value',1);

points2=handles.newtraces;
showframeshift(handles,points2)

end

guidata(hObject, handles);

%-----------------------------------------------------------------
function trajRbutton_Callback(hObject, eventdata, handles)

status=get(handles.trajhorizontalshift,'enable');
output=get(handles.acceptpushbutton,'userdata');

if size(status,2)==2 %on

%RIGHT
handles.thold=get(handles.trajhorizontalshift,'String');
kteh=str2num(handles.thold);
kteh=kteh+0.05;
set(handles.trajhorizontalshift,'string',num2str(kteh));
datashift=get(handles.panel,'userdata');
datashift(2,4)=kteh;
set(handles.panel,'userdata',datashift);

output(1)=kteh;
set(handles.acceptpushbutton,'userdata',output)

handles.newtraces=shiftpoints(handles,datashift(1,4),kteh,handles.traces2);

indexzero=find(handles.newtraces(:,2)<handles.Xdim);
handles.newtraces=handles.newtraces(indexzero,:);
%set(handles.finishedpushbutton,'value',1);

points2=handles.newtraces;
showframeshift(handles,points2)

end

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function showframeshift(handles,newpoints)

if nargin<2
    newpoints=handles.traces2;
end

axes(handles.axes1);
imshow(ones(handles.Xdim,handles.Ydim),'InitialMagnification','fit');
axis off

hold on

% plot positions
plot(handles.traces1(:,1),handles.traces1(:,2),'.','MarkerSize',3,'Color','r');

if isempty(newpoints)==0
    hold on
    plot(newpoints(:,1),newpoints(:,2),'.','MarkerSize',3,'Color','b');
end


clear rectangle yaelegidos actualimagen datamatrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function trajverticalshift_Callback(hObject, eventdata, handles)
valor=get(handles.trajverticalshift,'String');
set(handles.trajverticalshift,'string',valor);
guidata(hObject, handles);

function trajhorizontalshift_Callback(hObject, eventdata, handles)
valor=get(handles.trajhorizontalshift,'String');
set(handles.trajhorizontalshift,'string',valor);
guidata(hObject, handles);

function trajverticalshift_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function trajhorizontalshift_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



