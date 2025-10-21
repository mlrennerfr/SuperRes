function varargout = Cluster(varargin)
% CLUSTER MATLAB code for Cluster.fig
% Diinamic package for clustering analysis
% GUI for selecting files 
% Marianne Renner oct22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last Modified by GUIDE v2.5 29-Sep-2022 16:11:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Cluster_OpeningFcn, ...
                   'gui_OutputFcn',  @Cluster_OutputFcn, ...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cluster_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
set(handles.selectfilepushbutton,'userdata',[]);
set(handles.previouspushbutton,'Enable','off')
set(handles.nextpushbutton,'Enable','off')

h3 = findobj('Type','figure');
txtHand = findall(h3, '-property', 'FontUnits'); 
set(txtHand, 'FontUnits', 'normalized');

guidata(hObject, handles);

%--------------------------------------------------------------------------
function varargout = Cluster_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selectfilepushbutton_Callback(hObject, eventdata, handles)

% load data simplified version
szpx = str2num(get(handles.szpx,'String'));
mu=get(handles.muradiobutton,'Value');

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
  set(handles.nextpushbutton,'Enable','on')
else
  set (handles.filename, 'string',filename) ;
end

[namefile,~]=strtok(filename,'.'); %!!!!!!!!!!!!!!!!!!!!!!!!!!!!

set(handles.namefile,'String',namefile);
set(handles.analyzepushbutton,'Enable','on');
set(handles.analyzepresentpushbutton,'Enable','on');
set(handles.refreshpushbutton,'Enable','on');
set(handles.textlistafiles,'userdata',listafiles);
set(handles.textnrofiles,'userdata',1);

%load data, in µm
disp(' ')
disp('Reading data, please wait...')
if mu==0
    handles=loadselectPALMdata(filename,szpx,handles);
else
    handles=loadselectPALMdata(filename,1,handles); %data already in microns
end

set(handles.selectfilepushbutton,'enable','on')
set(handles.alpha1,'string','0')
set(handles.alpha2,'string','100')

% pointillistic
showimagesauto = get(handles.autoimradiobutton,'value');
if showimagesauto==1
    pointillistic(handles)
end

disp('Done')
disp(' ')

% Update handles structure
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refreshpushbutton_Callback(hObject, eventdata, handles)

pointillistic(handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function analyzepresentpushbutton_Callback(hObject, eventdata, handles)

% Values from window
dx = str2num(get(handles.PALMszpx,'String'));
szpx = str2num(get(handles.szpx,'String'));
alphamin = str2num(get(handles.alpha1,'String'));
alphamax = str2num(get(handles.alpha2,'String'));
mu=get(handles.muradiobutton,'Value');

nrofile=get(handles.textnrofiles,'userdata');
listafiles=get(handles.textlistafiles,'userdata');

%check for previous parameters
if length(dir('auxdens.mat'))>0 
    load('auxdens.mat','-mat')
else
    Densparam.maxdiamcluster=0;
end

%load data
filename=listafiles{nrofile};
handles.file=filename;

[namefile,~]=strtok(filename,'.')  ;   %sin extension
set(handles.namefile,'String', namefile);

% control existance of ROIs
if exist([namefile,'.rgn'],'file')
    
        names{1}=filename;  % for report
        
        if mu==0
            handles=loadselectPALMdata(filename,szpx,handles);
        else
            handles=loadselectPALMdata(filename,1,handles); %data already in microns
        end
        
        X_mu=handles.x-min(handles.x);
        Y_mu=handles.y-min(handles.y);
        alpha = handles.alpha(:);
        alpha(1:10);   
        xdim=ceil(max(X_mu)/dx);
        ydim=ceil(max(Y_mu)/dx);
        
        % remove by intensities
        vect2remove=selectremovealpha2(handles,X_mu,Y_mu,alpha);
        if ~isempty(vect2remove)
            [X_mu,Y_mu,alpha,fr]=removealpha(vect2remove,X_mu,Y_mu,alpha,handles.fr);
        end
       
        % Create rendered image
        [handles.I,xxi,yyi] = PALM_rendering3(X_mu,Y_mu,alpha,dx*2,dx,0,xdim, ydim, 1);
        
      %  figure
      %  imshow(handles.I,'InitialMagnification','fit')
        
        % clustering analysis
        typeanalysis=1; % only one file
        
        varargin{1} = handles.I;
        varargin{2} = []; %roifile;
        varargin{3} = handles.fr;
        varargin{4} = handles.alpha;
        varargin{5} = filename;  
        varargin{6} = dx;
        varargin{7} = szpx;
        varargin{8} = Densparam;
        varargin{9} = typeanalysis;
        
        if length(dir('auxiliar.mat'))>0  
            delete auxiliar.mat        %clears previous results
        end

        varargout=DetectClusters(varargin);
        uiwait;   
              
        %read results (in auxiliar.mat) for each file (all rois)
        control=0;
        aux='auxiliar.mat'; 
       
        if length(dir(aux))>0  
            load(aux,'-mat')
            load('auxdens.mat','-mat')
            
            %abort signal?
            if abortsignal==1
                return                
            else
                savenameout=[namefile,'-allclust.txt'];
                savenamenano=[namefile,'-nanoclust.txt'];
                savenamedist=[namefile,'-clustperdist.txt'];
                resultsclu=allresult.clu;
                resultsnano=allresult.nano;
                resdistance=allresult.dist;
                
                if isempty(resultsclu)==0
                    save(savenameout,'resultsclu','-ascii')
                end
                if isempty(resdistance)==0 % distance data
                    save(savenamedist,'resdistance','-ascii')
                end
                if isempty(resultsnano)==0
                    save(savenamenano,'resultsnano','-ascii')
                end
                delete 'auxiliar.mat'
                
                %save report
                reportclustering(Densparam,names);

            end
        else
            return %something wrong happened
        end % empty results (aux)
else % no ROI for this file
    disp('No ROI found for this file : no analysis done')
    disp(' ')
end % control existance ROI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function analyzepushbutton_Callback(hObject, eventdata, handles)

% Values from window
dx = str2num(get(handles.PALMszpx,'String'));
szpx = str2num(get(handles.szpx,'String'));
alphamin = str2num(get(handles.alpha1,'String'));
alphamax = str2num(get(handles.alpha2,'String'));
mu=get(handles.muradiobutton,'Value');
batch=get(handles.batchradiobutton,'Value');

currentfolder=cd;

listafiles=get(handles.textlistafiles,'userdata');
for i=1:size(listafiles,2)
    newlistafiles(i).name=listafiles{i};
end

if batch==1 %batch mode
    %check for previous parameters
    if length(dir('auxdens.mat'))>0 
        load('auxdens.mat','-mat')
    else
        Densparam.maxdiamcluster=0;
    end
    
    varargin{1} = newlistafiles;
    varargin{2} = dx;
    varargin{3} = szpx;
    varargin{4} = mu;
    varargin{5} = alphamin;
    varargin{6} = alphamax;
    varargin{7} = Densparam;

    %DetectClustersBatchOptim(varargin); % detects clusters and saves results in batch mode (parameters set at the begining)
    DetectClustersBatch(varargin); % detects clusters and saves results in batch mode (parameters set at the begining)
    uiwait;   
   
else

  set(handles.textnrofiles,'userdata',0); %starts from file 1
  names={};
  count=1;
  
  logical=1;
  
  while logical  % loop files
      
    %check for previous parameters
    if length(dir('auxdens.mat'))>0 
        load('auxdens.mat','-mat')
    else
        Densparam.maxdiamcluster=0;
    end
    
    %load data
    cd(currentfolder)
    nrofile=get(handles.textnrofiles,'userdata');
    next=nrofile+1;
    if next>size(listafiles,2)
        logical=0;
        break
    end
    
    set(handles.textnrofiles,'userdata',next);
    filename=listafiles{next};
    handles.file=filename;
    
    [namefile,~]=strtok(filename,'.')  ;   %sin extension
    set(handles.namefile,'String', namefile);
    
    % control existance of ROIs
    if exist([namefile,'.rgn'],'file')
        
        names{count}=filename;  % for report
        set(handles.filename,'string',[filename,' (',num2str(next),' of ',num2str(size(listafiles,2)),')']);
        
        if mu==0
            handles=loadselectPALMdata(filename,szpx,handles);
        else
            handles=loadselectPALMdata(filename,1,handles); %data already in microns
        end
        
        X_mu=handles.x-min(handles.x);
        Y_mu=handles.y-min(handles.y);
        alpha = handles.alpha(:);
        alpha(1:10);   
        xdim=ceil(max(X_mu)/dx);
        ydim=ceil(max(Y_mu)/dx);
        
        % remove by intensities
        vect2remove=selectremovealpha2(handles,X_mu,Y_mu,alpha);
        if ~isempty(vect2remove)
            [X_mu,Y_mu,alpha,fr]=removealpha(vect2remove,X_mu,Y_mu,alpha,handles.fr);
        end
       
        % Create rendered image
        [handles.I,xxi,yyi] = PALM_rendering3(X_mu,Y_mu,alpha,dx*2,dx,0,xdim, ydim, 1);
        
        % clustering analysis
        typeanalysis=2; % loop aver all files
        
        varargin{1} = handles.I;
        varargin{2} = []; %roifile;
        varargin{3} = handles.fr;
        varargin{4} = handles.alpha;
        varargin{5} = filename;  
        varargin{6} = dx;
        varargin{7} = szpx;
        varargin{8} = Densparam;
        varargin{9} = typeanalysis;
        
        if length(dir('auxiliar.mat'))>0  
            delete auxiliar.mat        %clears previous results
        end

        varargout=DetectClusters(varargin);
        uiwait;   
              
        %read results (in auxiliar.mat) for each file (all rois)
        control=0;
        aux='auxiliar.mat'; 
       
        if length(dir(aux))>0  
            load(aux,'-mat')
            load('auxdens.mat','-mat')
            
            if abortsignal==1
                logical=0;
                break
            else
                savenameout=[namefile,'-allclust.txt'];
                savenamenano=[namefile,'-nanoclust.txt'];
                savenamedist=[namefile,'-clustperdist.txt'];
                resultsclu=allresult.clu;
                resultsnano=allresult.nano;
                resdistance=allresult.dist;
                
                if isempty(resultsclu)==0
                    save(savenameout,'resultsclu','-ascii')
                end
                if isempty(resdistance)==0 % distance data
                    save(savenamedist,'resdistance','-ascii')
                end
                if isempty(resultsnano)==0
                    save(savenamenano,'resultsnano','-ascii')
                end
                delete 'auxiliar.mat'
                
                %save report
                reportclustering(Densparam,names);
                
            end
        else
            break %something wrong happened
        end % empty results (aux)
    else % no ROI for this file
        disp('No ROI found for this file : no analysis done')
        disp(' ')
    end % control existance ROI
    count=count+1; %number file

  end %while loop files
  
end %batch mode

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function previouspushbutton_Callback(hObject, eventdata, handles)

szpx= str2num(get(handles.szpx,'String'));
mu=get(handles.muradiobutton,'Value');
PALMszpx= str2num(get(handles.PALMszpx,'String'));
listafiles=get(handles.textlistafiles,'userdata');
nrofile=get(handles.textnrofiles,'userdata');

next=nrofile-1;

% go for previous file
if next>0
    set(handles.textnrofiles,'userdata',next);
    filename=listafiles{next}
    [namefile,~]=strtok(filename,'.'); %sin extension
    set(handles.namefile,'String',namefile);
    set(handles.filename,'string',[filename,' (',num2str(next),' of ',num2str(size(listafiles,2)),')']);
    set(handles.textlistafiles,'userdata',listafiles);
    set(handles.textnrofiles,'userdata',next);
    set(handles.filename,'userdata',filename);
    
    if mu==0
        handles=loadselectPALMdata(filename,szpx,handles);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nextpushbutton_Callback(hObject, eventdata, handles)

szpx= str2num(get(handles.szpx,'String'));
mu=get(handles.muradiobutton,'Value');
PALMszpx= str2num(get(handles.PALMszpx,'String'));
listafiles=get(handles.textlistafiles,'userdata');
nrofile=get(handles.textnrofiles,'userdata');

next=nrofile+1;

% go for next file
if size(listafiles,2)<next
    set(handles.nextpushbutton,'Enable','off')
else
    set(handles.textnrofiles,'userdata',next);
    filename=listafiles{next}
    [namefile,~]=strtok(filename,'.'); %sin extension
    set(handles.namefile,'String',namefile);
    set(handles.filename,'string',[filename,' (',num2str(next),' of ',num2str(size(listafiles,2)),')']);
    set(handles.textlistafiles,'userdata',listafiles);
    set(handles.textnrofiles,'userdata',next);
    set(handles.filename,'userdata',filename);
    
    if mu==0
        handles=loadselectPALMdata(filename,szpx,handles);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function optim_Callback(hObject, eventdata, handles)

function menuoptim_Callback(hObject, eventdata, handles)

% Values from window
dx = str2num(get(handles.PALMszpx,'String'));
szpx = str2num(get(handles.szpx,'String'));
alphamin = str2num(get(handles.alpha1,'String'));
alphamax = str2num(get(handles.alpha2,'String'));
mu=get(handles.muradiobutton,'Value');
batch=get(handles.batchradiobutton,'Value');

currentfolder=cd;

listafiles=get(handles.textlistafiles,'userdata');
for i=1:size(listafiles,2)
    newlistafiles(i).name=listafiles{i};
end

%check for previous parameters
    if length(dir('auxdens.mat'))>0 
        load('auxdens.mat','-mat')
    else
        Densparam.maxdiamcluster=0;
    end
    
%    varargin{1} = newlistafiles; %not needed
    varargin{2} = dx;
    varargin{3} = szpx;
    varargin{4} = mu;
    varargin{5} = alphamin;
    varargin{6} = alphamax;
    varargin{7} = Densparam;

    OptimDetectClusters(varargin); % detects clusters and saves results in batch mode (parameters set at the begining)
    uiwait;   
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function namefile_Callback(hObject, eventdata, handles)
function namefile_CreateFcn(hObject, eventdata, handles)
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

function szpx_Callback(hObject, eventdata, handles)
function szpx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function PALMszpx_Callback(hObject, eventdata, handles)
function PALMszpx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function muradiobutton_Callback(hObject, eventdata, handles)
function autoimradiobutton_Callback(hObject, eventdata, handles)
function batchradiobutton_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



