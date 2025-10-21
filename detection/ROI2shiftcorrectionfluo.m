function varargout = ROI2shiftcorrectionfluo(varargin)
%
% GUI for color correction
%
% Marianne Renner, SuperRes_v3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last Modified by GUIDE v2.5 03-Mar-2020 16:01:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ROI2shiftcorrectionfluo_OpeningFcn, ...
                   'gui_OutputFcn',  @ROI2shiftcorrectionfluo_OutputFcn, ...
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


% --- Executes just before ROI2shiftcorrectionfluo is made visible.
function ROI2shiftcorrectionfluo_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;


output=[0 0];
set(handles.acceptpushbutton,'userdata',output);

save('auxiliar.mat','output','-mat')

set(handles.panel,'userdata',zeros(4,4)); %shifts

varargout{1} = [];

%        varargin{1}=namefile; %
%        varargin{2}=dimy; %inverted!
 %       varargin{3}=dimx;
 %       varargin{4}=newpoints;
 %       varargin{5}=newpoints2; 
 %       varargin{6}=data2zoom; 

   
handles.file=cell2mat(varargin{1}(1)); % varargin{1}=namefile; %
handles.Xdim=cell2mat(varargin{1}(2)); % varargin{2}=dimx;
handles.Ydim=cell2mat(varargin{1}(3));   % varargin{3}=dimy;
handles.traces1=cell2mat(varargin{1}(4)); % varargin{4}=newpoints;
handles.traces2=cell2mat(varargin{1}(5)); % varargin{4}=newpoints;
handles.imagefluo=cell2mat(varargin{1}(6));       % varargin{5}=data2; %  
handles.newtraces1=handles.traces1;
handles.newtraces2=handles.traces2;

showframeshiftfluo2(handles)
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = ROI2shiftcorrectionfluo_OutputFcn(hObject, eventdata, handles)

varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
ktev=ktev-0.5;
set(handles.trajverticalshift,'string',num2str(ktev));
datashift=get(handles.panel,'userdata');
datashift(1,4)=ktev;
set(handles.panel,'userdata',datashift);

output(2)=ktev;
set(handles.acceptpushbutton,'userdata',output)
handles.newtraces=shiftpoints(handles,ktev,datashift(2,4),handles.traces1);
indexzero=find(handles.newtraces(:,2)>0);
handles.newtraces=handles.newtraces(indexzero,:);
points2=handles.newtraces;

handles.newtraces2=shiftpoints(handles,ktev,datashift(2,4),handles.traces2);
indexzero=find(handles.newtraces2(:,2)>0);
handles.newtraces2=handles.newtraces2(indexzero,:);
points22=handles.newtraces2;

showframeshiftfluo2(handles,points2,points22)

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

ktev=ktev+0.5;
set(handles.trajverticalshift,'string',num2str(ktev));
datashift=get(handles.panel,'userdata');
datashift(1,4)=ktev;
set(handles.panel,'userdata',datashift);

output(2)=ktev;
set(handles.acceptpushbutton,'userdata',output)
handles.newtraces=shiftpoints(handles,ktev,datashift(2,4),handles.traces1);
indexzero=find(handles.newtraces(:,2)<handles.Ydim);
handles.newtraces=handles.newtraces(indexzero,:);
points2=handles.newtraces;

handles.newtraces2=shiftpoints(handles,ktev,datashift(2,4),handles.traces2);
indexzero=find(handles.newtraces2(:,2)>0);
handles.newtraces2=handles.newtraces2(indexzero,:);
points22=handles.newtraces2;

showframeshiftfluo2(handles,points2,points22)

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
kteh=kteh-0.5;
set(handles.trajhorizontalshift,'string',num2str(kteh));
datashift=get(handles.panel,'userdata');
datashift(2,4)=kteh;
set(handles.panel,'userdata',datashift);

output(1)=kteh;
set(handles.acceptpushbutton,'userdata',output)

handles.newtraces=shiftpoints(handles,datashift(1,4),kteh,handles.traces1);
indexzero=find(handles.newtraces(:,2)>0);
handles.newtraces=handles.newtraces(indexzero,:);
points2=handles.newtraces;

handles.newtraces2=shiftpoints(handles,datashift(1,4),kteh,handles.traces2);
indexzero=find(handles.newtraces2(:,2)>0);
handles.newtraces2=handles.newtraces2(indexzero,:);
points22=handles.newtraces2;

showframeshiftfluo2(handles,points2,points22)

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
kteh=kteh+0.5;
set(handles.trajhorizontalshift,'string',num2str(kteh));
datashift=get(handles.panel,'userdata');
datashift(2,4)=kteh;
set(handles.panel,'userdata',datashift);

output(1)=kteh;
set(handles.acceptpushbutton,'userdata',output)

handles.newtraces=shiftpoints(handles,datashift(1,4),kteh,handles.traces1);
indexzero=find(handles.newtraces(:,2)<handles.Xdim);
handles.newtraces=handles.newtraces(indexzero,:);
points2=handles.newtraces;

handles.newtraces2=shiftpoints(handles,datashift(1,4),kteh,handles.traces2);
indexzero=find(handles.newtraces2(:,2)>0);
handles.newtraces2=handles.newtraces2(indexzero,:);
points22=handles.newtraces2;

showframeshiftfluo2(handles,points2,points22)

end

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function showframeshiftfluo2(handles,newpoints,newpoints2)

if nargin<2
    newpoints=handles.traces1;
    newpoints2=handles.traces2;
end


axes(handles.axes1);
datamatrix=handles.imagefluo;
stackmin=min(min(datamatrix));
stackmax=max(max(datamatrix));
imshow(datamatrix,[stackmin stackmax],'InitialMagnification','fit');
axis off

hold on

% plot positions
%if isempty(newpoints)==0
   plot(newpoints(:,1),newpoints(:,2),'.','MarkerSize',3,'Color','r');
%end
%if isempty(newpoints2)==0
    plot(newpoints2(:,1),newpoints2(:,2),'.','MarkerSize',3,'Color','b');
%end


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
