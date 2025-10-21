function varargout = Segment(varargin)
% function varargout = Segment(varargin)
% Diinamic package
% GUI to create ROIs and segmented images (masks)
%
% Marianne Renner may 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last Modified by GUIDE v2.5 09-Nov-2022 10:23:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Segment_OpeningFcn, ...
                   'gui_OutputFcn',  @Segment_OutputFcn, ...
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes just before Segment is made visible.
function Segment_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

set(handles.displaypushbutton,'Enable','off');
set(handles.nextpushbutton,'Enable','off');
set(handles.createroipushbutton,'Enable','off');
set(handles.createmaskpushbutton,'Enable','off');

% to fit screen definition/changes in size
h3 = findobj('Type','figure');
txtHand = findall(h3, '-property', 'FontUnits');
%set(txtHand, 'FontUnits', 'normalized');
set(txtHand, 'FontUnits', 'centimeters');
    
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = Segment_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadpushbutton_Callback(hObject, eventdata, handles)

% load data simplified version
px_mu = str2num(get(handles.szpx,'String'));
mu=get(handles.muradiobutton,'Value');

dialog_title=['Select data folder'];
directory_name = uigetdir(cd,dialog_title);
if directory_name==0
    return
end
cd(directory_name);

d = dir('*.mat'); % movie
st={d.name};
cuenta=1;
if isempty(st)==1
   msgbox('No files!!','Select files','error');
   return
else
    for jj=1:size(st,2)
        k1 = strfind(st{jj},'roidata');
        k2 = strfind(st{jj},'Irend');
        k3 = strfind(st{jj},'auxdens');
        k4 = strfind(st{jj},'mask');
        if isempty(k1)==1 && isempty(k2)==1 && isempty(k3)==1 && isempty(k4)==1
            stok{cuenta}=st{jj};
            cuenta=cuenta+1;
        end
    end
end
[files,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',stok);
if v==0
   return
end

for i=1:size(files,2)
      listafiles{i}=stok{files(i)};
end

filename=listafiles{1}; %first file
set (handles.filename, 'Userdata',filename);

if size(files,2)>1 %batch
  set (handles.filename, 'string',['Batch: File ',filename,' (1/',num2str(size(files,2)),')']) ;
  set(handles.nextpushbutton,'Enable','on');
else
  set (handles.filename, 'string',filename) ;
end

set(handles.displaypushbutton,'Enable','on');
set(handles.createroipushbutton,'Enable','on');
set(handles.createmaskpushbutton,'Enable','on');

[namefile,~]=strtok(filename,'.'); %!!!!!!!!!!!!!!!!!!!!!!!!!!!!

handles.ffname = fullfile(directory_name,filename);
set(handles.savename,'String',namefile);
set(handles.textlistafiles,'userdata',listafiles);
set(handles.nextpushbutton,'userdata',directory_name);
set(handles.textnrofiles,'userdata',1);

if mu==0
    handles=loadselectPALMdata(filename,px_mu,handles);
else
    handles=loadselectPALMdata(filename,1,handles); %data already in microns
end

set(handles.alpha1,'string','0')
set(handles.alpha2,'string','100')

% pointillistic
showimagesauto = get(handles.autoimradiobutton,'value');
if showimagesauto==1
    pointillistic(handles)
end

% Update handles structure
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function createroipushbutton_Callback(hObject, eventdata, handles)

listafiles=get(handles.textlistafiles,'userdata');

if isempty(listafiles)==0
    varargin{1}=handles; 
    disp('')
    disp('Draw ROIs')
    disp('Reading previous data, please wait...')
    
    % window to create ROIs
    varargout=CreateROI(varargin); % 
else
    disp('Choose the files first')
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function createmaskpushbutton_Callback(hObject, eventdata, handles)

docreateMask(handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function displaypushbutton_Callback(hObject, eventdata, handles)

pointillistic(handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function previouspushbutton_Callback(hObject, eventdata, handles)

px_mu= str2num(get(handles.szpx,'String'));
mu=get(handles.muradiobutton,'Value');
listafiles=get(handles.textlistafiles,'userdata');
nrofile=get(handles.textnrofiles,'userdata');
directory_name=get(handles.nextpushbutton,'userdata');
next=nrofile-1;

% go for previous file
if next>0
    set(handles.textnrofiles,'userdata',next);
    filename=listafiles{next}
    handles.ffname = fullfile(directory_name,filename);
    [namefile,~]=strtok(filename,'.'); %sin extension
    set(handles.savename,'String',namefile);

    set(handles.filename,'string',[filename,' (',num2str(next),' of ',num2str(size(listafiles,2)),')']);
    
    set(handles.textlistafiles,'userdata',listafiles);
    set(handles.textnrofiles,'userdata',next);
    set(handles.filename,'userdata',filename);
    
    if mu==0
        handles=loadselectPALMdata(filename,px_mu,handles);
    else
        handles=loadselectPALMdata(filename,1,handles); %data already in microns
    end
    
    % pointillistic
    showimagesauto = get(handles.autoimradiobutton,'value');
    if showimagesauto==1
        pointillistic(handles)
    end
    set(handles.nextpushbutton,'Enable','on')
    if next==1
        set(handles.previouspushbutton,'Enable','off')
    end
else
    set(handles.previouspushbutton,'Enable','off')
end %there is previous
 
clear file listafiles 

% Update handles structure
guidata(gcbo, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nextpushbutton_Callback(hObject, eventdata, handles)

px_mu= str2num(get(handles.szpx,'String'));
mu=get(handles.muradiobutton,'Value');
listafiles=get(handles.textlistafiles,'userdata');
nrofile=get(handles.textnrofiles,'userdata');
directory_name=get(handles.nextpushbutton,'userdata');
next=nrofile+1;

% go for next file
if size(listafiles,2)<next
    set(handles.nextpushbutton,'Enable','off')
else
    set(handles.textnrofiles,'userdata',next);
    filename=listafiles{next}
    handles.ffname = fullfile(directory_name,filename);
    [namefile,~]=strtok(filename,'.'); %sin extension
    set(handles.savename,'String',namefile);

    set(handles.filename,'string',[filename,' (',num2str(next),' of ',num2str(size(listafiles,2)),')']);
    
    set(handles.textlistafiles,'userdata',listafiles);
    set(handles.textnrofiles,'userdata',next);
    set(handles.filename,'userdata',filename);
    
    if mu==0
        handles=loadselectPALMdata(filename,px_mu,handles);
    else
        handles=loadselectPALMdata(filename,1,handles); %data already in microns
    end
    
    % pointillistic
    showimagesauto = get(handles.autoimradiobutton,'value');
    if showimagesauto==1
        pointillistic(handles)
    end
    
    set(handles.previouspushbutton,'Enable','on')
    if next==size(listafiles,2)
        set(handles.nextpushbutton,'Enable','off')
    end

end %there is next
 
clear file listafiles 

% Update handles structure
guidata(gcbo, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function helptag_Callback(hObject, eventdata, handles)

%-------------------------------------------------------------------------

function openhelp_Callback(hObject, eventdata, handles)

helpsegment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function extras_Callback(hObject, eventdata, handles)

%-------------------------------------------------------------------------

function batchfull_Callback(hObject, eventdata, handles)

disp('Creating ROIs of the full image')

fullimageroi(handles)

disp('Done')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pointillistic(handles)
% plots pointillistic images
 
X_mu=handles.x;
Y_mu=handles.y;
alpha=handles.alpha;
alpha(1:10);
fr=handles.fr;
gausssmooth = 1;
dx = str2num(get(handles.PALMszpx,'String'));
sigmaloc=dx*2;
    
vect2remove=selectremovealpha2(handles,X_mu,Y_mu,alpha); %selection by absolute intensities
if ~isempty(vect2remove)
    [X_mu,Y_mu,alpha,fr]=removealpha(vect2remove,X_mu,Y_mu,alpha,fr);
end
 
% clean data
xmax = max(X_mu);
ymax = max(Y_mu);
 
% figure
fig_points = figure('Name','Pointillist image','Toolbar','figure');
axis equal;
pointsize = 2;
set(gca,'YDir','reverse')

xlim([0 xmax]); ylim([0 ymax]);
color = 'k';    
hold on;
title([' Selected positions =',num2str(length(X_mu))]);
handles.pointsplotted =  plot(X_mu,Y_mu,'.','MarkerSize',pointsize,'Color',color);
hold on

figure(fig_points);
handles.fig_points = fig_points;

%rendered
xdim=ceil(max(X_mu)/dx);
ydim=ceil(max(Y_mu)/dx);

%% Create rendered image
[I,xxi,yyi] = PALM_rendering3( X_mu,Y_mu,alpha,dx*2,dx,0,xdim, ydim, gausssmooth);
figselected = figure('Name','Rendered image','Toolbar','figure');
hold on;

% rendered in false colors, pixel 0,0 top left
imshow(I,'InitialMagnification','fit')
colormap(gca,'hot')
hold on
 
rendx=X_mu/dx;
rendy=Y_mu/dx;
handles.rendx=rendx;
handles.rendy=rendy;

set(handles.displaypushbutton,'userdata',I); %also to recover size of rendered
clear x y X_mu Y_mu alpha fr I

% Update handles structure
guidata(gcbo, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = flagmatrixelements(M,vectorindices)

siz = size(M);
aux = M(:);
aux(vectorindices) = 1;
M = reshape(aux,siz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function figurapoints=plotpoints(x,y,fr,handles)
 
figurapoints = figure('Name','Pointillist image','Toolbar','figure');
 
% set boundaries
if isfield(handles,'xlim')
    [xmin,xmax] = deal(handles.xlim(1),handles.xlim(2));
    [ymin,ymax] = deal(handles.ylim(1),handles.ylim(2));
else
    xmin = min(x); xmax = max(x);
    ymin = min(y); ymax = max(y);
end
axis equal;
pointsize = str2num(get(handles.dotsize,'String'));
xlim([xmin xmax]); ylim([ymin ymax]);
 
color = 'k';    
hold on;
 
xlabel('X (mu)');
ylabel('Y (mu)');   
handles.pointsplotted =  plot(x,y,'.','MarkerSize',pointsize,'Color',color);
title([' Total nb of positions =',num2str(length(x))]);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I = image_ch(varargin)
xx1 = varargin{1};
yy1 = varargin{2};
I = varargin{3};
pct_sat = 0;
if ndims(I)==2
    if nargin>3
        clims = varargin{4};
    else
        Imax = prctile(I(:),100-pct_sat);
        Imin = min(I(:));
        clims = [Imin Imax];
    end
    if any(isnan(clims)) || clims(2)<=clims(1)
        warning(['clims (=',num2str(clims),')! Using [0,1] instead !']);
        clims = [0 1];
    end
    I = imagesc(xx1,yy1,I,clims);
else
    image(xx1,yy1,I);
end
colorbar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function szpx_Callback(hObject, eventdata, handles)
function szpx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function savename_Callback(hObject, eventdata, handles)
function savename_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function alpha1_Callback(hObject, eventdata, handles)
function alpha1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function alpha2_Callback(hObject, eventdata, handles)
function alpha2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function PALMszpx_Callback(hObject, eventdata, handles)
function PALMszpx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function minsizeclumask_Callback(hObject, eventdata, handles)
function minsizeclumask_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function intensthresh_Callback(hObject, eventdata, handles)
function intensthresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function muradiobutton_Callback(hObject, eventdata, handles)
function lengthradiobutton_Callback(hObject, eventdata, handles)
function autoimradiobutton_Callback(hObject, eventdata, handles)
function savematradiobutton_Callback(hObject, eventdata, handles)
function smallradiobutton_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
