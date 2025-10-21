function varargout = testspeakSM(varargin)
%TESTSPEAKSM M-file for testspeakSM.fig
% Edit the above text to modify the response to help testspeakSM

% Last Modified by GUIDE v2.5 20-Jul-2021 09:36:23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @testspeakSM_OpeningFcn, ...
                   'gui_OutputFcn',  @testspeakSM_OutputFcn, ...
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
function testspeakSM_OpeningFcn(hObject, eventdata, handles, varargin)

% default output
handles.output = hObject;
handles.type= 1; %intensity
handles.low= '0';
handles.up='0';
handles.bin='0';
handles.mean=0;
handles.nro=0;
handles.file='';

handles.output = hObject;

guidata(hObject, handles);

function varargout = testspeakSM_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------
function lowlimit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function lowlimit_Callback(hObject, eventdata, handles)
handles.low=get(hObject,'String');
guidata(hObject, handles);
%-------------------------------------------------------------------
function uplimit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function uplimit_Callback(hObject, eventdata, handles)
handles.up=get(hObject,'String');
guidata(hObject, handles);
%-------------------------------------------------------------------
function binsize_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function binsize_Callback(hObject, eventdata, handles)
handles.bin=get(hObject,'String');
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function loadpk_Callback(hObject, eventdata, handles)

% input data
[file,path] = uigetfile('*.pk','Load peak data file');
filename = [path,file];
if file==0
    return
end

handles.file=filename;
set (handles.filepk,'string',file);
set (handles.pushbutton5,'enable','on');
set(handles.uplimit,'string','0');
handles.up='0';
set(handles.lowlimit,'string','0');
handles.low='0';
set(handles.binsize,'string','');
handles.bin='1';

guidata(hObject, handles);

%-------------------------------------------------------------------
function pushbutton5_Callback(hObject, eventdata, handles)

switch handles.type
    case 1
        column=5;
    case 2
        column=4;
end
resp=1;
firstentry=1;
filename=handles.file;
peak =load(filename);
[nropeak,colu]=size(peak);
[filas,col]=size(peak);
minval=ceil(min(peak(:,column)))-1
if minval<0; minval=0; end
maxval=ceil(max(peak(:,column)))

if maxval>minval    %only if good values
    %initial values
    limsup=str2num(handles.up)
    if limsup==0
        limsup=maxval
    end
    liminf=str2num(handles.low)
    if liminf<limsup
    else
        liminf=minval
    end
    
    indexpeak=find(peak(:,column)>liminf & peak(:,column)<limsup);
    for counter=1:size(indexpeak,1)
        y(counter)=peak(indexpeak(counter),column);
    end
    bin=str2num(handles.bin);
    if bin<2
        bin=ceil((limsup-liminf)/50);
    end
    nbins=ceil((limsup-liminf)/bin);
    bin=ceil(bin);
    set(handles.uplimit,'string',num2str(limsup));
    set(handles.lowlimit,'string',num2str(liminf));
    set(handles.binsize,'string',num2str(bin));
    
    %disp('Calculating the histogram...')
    [n,xout]=hist(y(:),nbins);
    [fil,col]=size(n);
    count=1;
    if liminf>min(y(:))
        for i=1:col
            if xout(i)>liminf
                xnew(count)=xout(i);
                nnew(count)=n(i);
                count=count+1;
            end
        end
        xout=xnew;
        n=nnew;
    end
    xnew=[];
    nnew=[];
    [fil,col]=size(n);
    count=1;
    if limsup<(xout(1))
        limsup=xout(1)+1;
    end
    if limsup<max(y(:))
        for i=1:col
            if xout(i)<limsup
                xnew(count)=xout(i);
                nnew(count)=n(i);
                count=count+1;
            end
        end
        xout=xnew;
        n=nnew;
    end
    sizen=max(n(:))+max(n(:))/10;
    linf=liminf-liminf/20;
    lsup=limsup+limsup/20;
    
    axes(handles.axes1)
    axis([linf lsup 0 sizen]);
    bar(xout,n)
    m=mean(peak(:,column));
    switch handles.type
        case 1
            xlabel('Intensity'); ylabel('Counts');
        case 2
            xlabel('Width'); ylabel('Counts');
    end
    set (handles.meanvalue,'string',num2str(m));
    set (handles.number,'string',num2str(nropeak));
else
    disp('Cannot make the histogram, something is wrong with the data')
end % good values

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function popupmenu1_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function popupmenu1_Callback(hObject, eventdata, handles)

handles.type=get(hObject,'Value');
set(handles.uplimit,'string','0');
handles.up='0';
set(handles.lowlimit,'string','0');
handles.low='0';
set(handles.binsize,'string','');
handles.bin='0';

guidata(hObject, handles);

%-------------------------------------------------------------------
% --- Executes on button press in quit.
function quit_Callback(hObject, eventdata, handles)

close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
