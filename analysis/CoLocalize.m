function varargout = CoLocalize(varargin)
% Last Modified by GUIDE v2.5 27-Jun-2022 21:23:44
%
% GUI for colocalization of SMLM data
%
% Marianne Renner, SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CoLocalize_OpeningFcn, ...
                   'gui_OutputFcn',  @CoLocalize_OutputFcn, ...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes just before CoLocalize is made visible.
function CoLocalize_OpeningFcn(hObject, eventdata, handles, varargin)

set(handles.localizepushbutton,'enable','off')
set(handles.locfilename,'String','')
set(handles.displaypushbutton,'enable','off')
%set(handles.nextpushbutton,'enable','off')

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

function varargout = CoLocalize_OutputFcn(hObject, eventdata, handles)

varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selectfilepushbutton_Callback(hObject, eventdata, handles)

% load data simplified version
px_mu = str2num(get(handles.szpx,'String'));
mu=get(handles.muradiobutton,'Value');

dialog_title='Select data folder';
directory_name = uigetdir(cd,dialog_title);
if directory_name==0
    return
end
cd(directory_name);

d = dir('*.mat'); % file 1: always .mat
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
[files,v] = listdlg('PromptString','Select one file:','SelectionMode','multiple','ListString',stok);
if v==0
   return
end

%for i=1:size(files,2)
%listafiles{1}=stok{files(1)};
%end

filename=stok{files(1)}; %first file
set (handles.filename, 'Userdata',filename);

if size(files,2)>1 %batch
  set (handles.filename, 'string',['Batch: File ',filename,' (1/',num2str(size(files,2)),')']) ;
else
  set (handles.filename, 'string',filename) ;
end

%[namefile,rem]=strtok(filename,'.'); 

handles.ffname = fullfile(directory_name,filename);

set(handles.filename,'userdata',filename);

set(handles.displaypushbutton,'userdata',directory_name);

if mu==0
    handles=loadselectPALMdata(filename,px_mu,handles);
else
    handles=loadselectPALMdata(filename,1,handles); %data already in microns
end

set(handles.displaypushbutton,'enable','on')

set(handles.alpha1,'string','0')
set(handles.alpha2,'string','100')

% pointillistic
%onlypointillistic(handles)

% Update handles structure
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localizepushbutton_Callback(hObject, eventdata, handles)

clusteringdata=get(handles.clusteringradiobutton,'Value');
dx = str2num(get(handles.PALMszpx,'String')); %in Âµm
erode=str2num(get(handles.erodevalue,'String'));
dilate=str2num(get(handles.dilatevalue,'String'));
distcent=str2num(get(handles.maxdistcentroids,'String')); % in Âµm
savemat=get(handles.savematradiobutton,'value');
mu=get(handles.muradiobutton,'Value');

detecincluster=get(handles.detecincluradiobutton,'value');
maxdistancecentroids=str2num(get(handles.maxdistcentroids,'String')); %in µm

detecfluo=get(handles.fluoradiobutton,'Value');
detecclu=get(handles.clusteringradiobutton,'Value');
%cludetec=get(handles.cludetradiobutton,'Value');
%clufluo=get(handles.cluvsfluoradiobutton,'Value');
cluclu=get(handles.clucluradiobutton,'Value');

%listafiles=get(handles.filename,'userdata'); % only one!!!
filename=get(handles.filename,'userdata'); % only one!!!
enterlocfilename=get(handles.locfilename,'String');

disp('Reading data, please wait...')


%for i=1:size(listafiles,2)
    
    %filename=listafiles{i}
    handles.file=filename;
    [namefile,~]=strtok(filename,'.')  ; 
    
    %ROIs
    %    ROI.coord{count}=coord;
    %    ROI.in{count}=in; %index of selected points
    %    ROI.x{count}=x*dx;
    %    ROI.y{count}=y*dx;
    %    ROI.xselec{count}=roiselectedx; %index of selected points
    %    ROI.yselec{count}=roiselectedy; %index of selected points
    %    ROI.imagerend{count}=imageroirend; %crop of rendered image
    %    ROI.dimx=[1 max(X_mu)];
    %    ROI.dimy=[1 max(Y_mu)];
    %    ROI.dim=[newxdim newydim];
    %    ROI.dx=dx;
    
    % ROIs ------------------------------------------------------------
    if exist([namefile,'.rgn'])==2
        dataroi=importdata([namefile,'.rgn']);
        if isstruct(dataroi)
            %roifile=dataroi.coord;
          %  indata=dataroi.in; %index points
            numberroi=size(dataroi.coord,2);
            dx=dataroi.dx;
           
            %if isfield(dataroi,'dist')
                % distance
           %     distance=dataroi.dist;
           % end
        else
            disp('ROI file not found')
            roifile=[];
        end %if sstruct
        disp(['Loading regions from ',namefile,'.rgn'])
        disp(' ')
    else
        % no ROI
       msgbox('Please define ROIs before','Error','error')
    end % of exist 
    
    %---------------------------------------------------------------------
    % comparison to a localization image
    %---------------------------------------------------------------------
   if detecfluo==1 %|| clufluo==1 
       
        % loc file
        % check if already loaded
        locfilename=get(handles.locfilename,'String')
        if isempty(locfilename)==1      
           % locfilename=[namefile,ident,maskid,'.tif']
            control=1;
            if exist(locfilename)==2
                set(handles.locfilename,'String',locfilename);
            else
                locfilename=enterlocfilename;
                if exist(locfilename)==2
                    % else
                    disp(['Localization file ',locfilename,' not found. Load the file manually'])
                    control=0;
                end
            end
            if control==0
                return
            end
        end
        
        % check that locfilename is .tif!!!
        if ~contains(locfilename,'.tif')
            msgbox('The file #2 has to be an image file (.tif)','Error','error')
            return
        end

        % modification of size
        synimage=double(imread(locfilename));
        
        %binarize
        level = graythresh(synimage);
        synimage = im2bw(synimage,level);
   
        %dilate/erode
        if dilate>0
                se = strel('disk',dilate,0);
                I = imdilate(synimage,se);
        else
                I=synimage;
        end
        if erode>0
                se = strel('disk',erode,0);
                I = imerode(I,se);
        end
        synimage=I;  
        
        %create a rendered image with data to localize to adjust image size
        synimage=checkfluoimagesize(handles, synimage,dx,1,mu);

        if detecfluo==1 % detections with respect to pixels
            
            detecvsfluo(handles, namefile, synimage,dx, dataroi,numberroi)
            
        else  % clusters with respect to pixels
            
            clustersvsfluo(synimage, namefile,dx,dataroi, numberroi)
           
        end %types coloc
        
   end
    
   %----------------------------------------------------------------------
    % detections vs clusters
    %---------------------------------------------------------------------
    
    if detecclu==1 % || cludetec==1

        % mask name
        locfilename=get(handles.locfilename,'String');
        
        % check that locfilename is .mat!!!
        if ~contains(locfilename,'.mat')
            msgbox('The file #2 has to be of .mat type','Error','error')
            return
        end

        [namefilemask,~]=strtok(locfilename,'.')  ; 
       % if detecclu==1 
            detecvsclusters(handles, namefilemask, namefile, dx,dataroi,numberroi,detecincluster,distcent)
       % else
       %     clustersvsdetec(handles, namefilemask, namefile, dx,dataroi,numberroi,detecincluster,distcent)
       % end
    
    end %type analysis det vs clusters
            
     %----------------------------------------------------------------------
    % clusters vs clusters
    %----------------------------------------------------------------------
    
    if cluclu==1    
        
         % file 2 name
         locfilename=get(handles.locfilename,'String');
         [namefilemask,~]=strtok(locfilename,'.'); 

        clustersvsclusters(locfilename, namefile,dx, dataroi, numberroi,maxdistancecentroids)
        
    end %type analysis clusters vs clusters
 
    
%end % loop files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function displaypushbutton_Callback(hObject, eventdata, handles)

detecfluo=get(handles.fluoradiobutton,'Value');
%detecclu=get(handles.clusteringradiobutton,'Value');
fluorend=get(handles.fluorendradiobutton,'Value');
%cludetec=get(handles.cludetradiobutton,'Value');
%clufluo=get(handles.cluvsfluoradiobutton,'Value');
cluclu=get(handles.clucluradiobutton,'Value');
mu=get(handles.muradiobutton,'Value');
erode=str2num(get(handles.erodevalue,'String'));
dilate=str2num(get(handles.dilatevalue,'String'));

dx = str2num(get(handles.PALMszpx,'String')); %in µm
szpx=str2num(get(handles.szpx,'String'));


%listafiles=get(handles.filename,'userdata'); % only one!!!
filename=get(handles.filename,'userdata'); % only one!!!
enterlocfilename=get(handles.locfilename,'String');

%filename=listafiles{1};
[namefile,~]=strtok(filename,'.')  ; 
locfilename=get(handles.locfilename,'String');

disp('Reading data, please wait...')

% ROIs ------------------------------------------------------------
if exist([namefile,'.rgn'])==2
        dataroi=importdata([namefile,'.rgn']);
        
        if isstruct(dataroi)
            roifile=dataroi.coord;
           % indata=dataroi.in; %index points
            numberroi=size(dataroi.coord,2);
            dx=dataroi.dx;
           
            %if isfield(dataroi,'dist')
                % distance
           %     distance=dataroi.dist;
           % end
        else
            disp('ROI file not found')
            roifile=[];
        end %if sstruct
        disp(['Loading regions from ',namefile,'.rgn'])
        disp(' ')
else
        % no ROI
       msgbox('Please define ROIs before','Error','error')
end % of exist 

%------------------------------------------------------------------------
% type colocalization

figure

if detecfluo==1 %|| clufluo==1           % with fluo
    
    % check that locfilename is .tif!!!
    if ~contains(locfilename,'.tif')
        msgbox('Check the colocalization options and the type of files','Error','error')
        return
    end

    synimage=double(imread(locfilename));
    %synimage=checkfluoimagesize(handles, synimage,dx,szpx,mu);
    synimage=checkfluoimagesize(handles, synimage,dx);
    
    % erode and dilate
    if dilate>0
        se = strel('disk',dilate,0);
        I2 = imdilate(synimage,se);
        I2 = imfill(I2,'holes');
    else
        I2=synimage;
    end
    
    if erode>0
        se = strel('disk',erode,0);
        I2 = imerode(I2,se);
    end
    synimage=I2;
    
    imshow(synimage,'InitialMagnification','fit')
    hold on
    
    % print detections
    X_mu=handles.x;
    if fluorend==1
        Y_mu=handles.y;
    else
        Y_mu=handles.y-szpx;
    end
    
    plot(X_mu/dx,Y_mu/dx,'.','MarkerSize',2,'Color',[1 0 0]);

    % plot ROIs
    if isempty(roifile)==0
        % for all rois
        for roicounter=1:numberroi
            coordroi=roifile{roicounter};
            plot(coordroi(:,1)/dx,coordroi(:,2)/dx,'-g')
        end
    end

else
    
     % print first set detections
    X_mu=handles.x;
    Y_mu=handles.y-szpx;
    plot(X_mu,Y_mu,'.','MarkerSize',3,'Color',[1 0 0]);
    hold on
    axis equal
    axis ij

   
    % second detections set
    if mu==0
        handles2=loadselectPALMdata(locfilename,szpx,handles);
    else
        handles2=loadselectPALMdata(locfilename,1,handles); %data already in microns
    end
    
    X_mu2=handles2.x;  %%%%%%%%%%%%%%%%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Y_mu2=handles2.y-szpx;
    plot(X_mu2,Y_mu2,'.','MarkerSize',3,'Color',[0 1 0]);
    xlim([0 max([max(X_mu), max(X_mu2)])]); ylim([0 max([max(Y_mu), max(Y_mu2)])]);

    % plot ROIs
    if isempty(roifile)==0
        % for all rois
        for roicounter=1:numberroi
            coordroi=roifile{roicounter};
            plot(coordroi(:,1),coordroi(:,2),'-b')
        end
        
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadlocfilepushbutton_Callback(hObject, eventdata, handles)

% options
detecfluo=get(handles.fluoradiobutton,'Value');
detecclu=get(handles.clusteringradiobutton,'Value');
%clufluo=get(handles.cluvsfluoradiobutton,'Value');
cluclu=get(handles.clucluradiobutton,'Value');


% loc file
%ATT mask .tif or .mat
if detecfluo==1 %|| clufluo==1
    % tif
    [locfilename,PathName,FilterIndex] = uigetfile('.tif','Choose the localization file','');
else
   % [locfilename,PathName,FilterIndex] = uigetfile('.mat','Choose the localization file','');
   
   d = dir('*.mat'); % file 1: always .mat
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
   [files,v] = listdlg('PromptString','Select one file:','SelectionMode','multiple','ListString',stok);
   if v==0
       return
   end
   
   locfilename=stok{files(1)}; %first file
    
end

set(handles.locfilename,'String',locfilename);
set(handles.localizepushbutton,'enable','on')


% Update handles structure
guidata(gcbo, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function quitpushbutton_Callback(hObject, eventdata, handles)

clear all
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function muradiobutton_Callback(hObject, eventdata, handles)
function epifluoradiobutton_Callback(hObject, eventdata, handles)

%options type of colocalization
function fluoradiobutton_Callback(hObject, eventdata, handles)
code=get(hObject, 'value');
if code==1
    set(handles.clusteringradiobutton,'Value',0);
 %   set(handles.cluvsfluoradiobutton,'Value',0);
    set(handles.clucluradiobutton,'Value',0);
  %  set(handles.cludetradiobutton,'Value',0);
end
guidata(hObject, handles);

function clusteringradiobutton_Callback(hObject, eventdata, handles)
code=get(hObject, 'value');
if code==1
    set(handles.fluoradiobutton,'Value',0);
%    set(handles.cluvsfluoradiobutton,'Value',0);
    set(handles.clucluradiobutton,'Value',0);
  %  set(handles.cludetradiobutton,'Value',0);
end
guidata(hObject, handles);

%function cluvsfluoradiobutton_Callback(hObject, eventdata, handles)
%val=get(hObject,'Value');
%if val==1
%    set(handles.clucluradiobutton,'Value',0);
%    set(handles.fluoradiobutton,'Value',0);
%    set(handles.clusteringradiobutton,'Value',0);
  %  set(handles.cludetradiobutton,'Value',0);
%end

function clucluradiobutton_Callback(hObject, eventdata, handles)
val=get(hObject,'Value');
if val==1
  %  set(handles.cluvsfluoradiobutton,'Value',0);
    set(handles.fluoradiobutton,'Value',0);
    set(handles.clusteringradiobutton,'Value',0);
  %  set(handles.cludetradiobutton,'Value',0);
end

%only detections inside clusters
function detecincluradiobutton_Callback(hObject, eventdata, handles)


%save .mat
function savematradiobutton_Callback(hObject, eventdata, handles)

%fluo from rendered
function fluorendradiobutton_Callback(hObject, eventdata, handles)

function distcentroid_Callback(hObject, eventdata, handles)
function distcentroid_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxdistcentroids_Callback(hObject, eventdata, handles)
function maxdistcentroids_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function distborder_Callback(hObject, eventdata, handles)
function distborder_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function szpx_Callback(hObject, eventdata, handles)
function szpx_CreateFcn(hObject, eventdata, handles)
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

function erodevalue_Callback(hObject, eventdata, handles)
function erodevalue_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function minoverlap_Callback(hObject, eventdata, handles)
function minoverlap_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dilatevalue_Callback(hObject, eventdata, handles)
function dilatevalue_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
