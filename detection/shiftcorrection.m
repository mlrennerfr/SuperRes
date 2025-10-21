function varargout = shiftcorrection(varargin)
% Last Modified by GUIDE v2.5 13-Mar-2014 17:23:28
%
% GUI shift correction
%
% Marianne Renner, SuperRes programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @shiftcorrection_OpeningFcn, ...
                   'gui_OutputFcn',  @shiftcorrection_OutputFcn, ...
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes just before shiftcorrection is made visible.
function shiftcorrection_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
set(handles.visupushbutton,'enable','off');
set(handles.donebutton,'enable','off');
set(handles.trypushbutton,'enable','off');
set(handles.tagtwocolor,'userdata',' '); % name second file
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = shiftcorrection_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadpushbutton_Callback(hObject, eventdata, handles)

option=get(handles.popupmenu1,'value');
st=[]; 
ident=get(handles.tagtwocolor,'string');
px_mu = str2num(get(handles.szpx,'String'));
dx = str2num(get(handles.accu,'String'));
set(handles.trypushbutton,'string','Try');

% carga files y loop general 
dialog_title=['Select data folder'];
directory_name = uigetdir(cd,dialog_title);
if directory_name==0
    return
end
cd(directory_name);

d = dir('*mat*'); % movie
st={d.name};
if isempty(st)==1
   msgbox(['No files!!'],'Select files','error');
   return
end
[files,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',st);
if v==0
   return
end

for i=1:size(files,2)
      listafiles{i}=st{files(i)};
end

file=listafiles{1}; %first file
firststr=file(1:2) ; %numbering!!!
%filename=file; 
if size(files,2)>1 %batch
  set (handles.filename, 'string',['Batch: File ',file,' (1/',num2str(size(files,2)),')']) ;
else
  set (handles.filename, 'string',file) ;
end
set(handles.file2,'string',[]);


handles=loadPALMdata(file,px_mu,0,handles);

%disp('first')
%disp(handles.x(1))

set(handles.loadpushbutton,'userdata',listafiles);
set(handles.trypushbutton,'userdata',1);

if option==3
    control=0;
    % name second super res image
    d=dir(['*',ident,'*']); %  files
    st={d.name} ;
    for j=1:size(st,2);
        k = findstr(st{j},firststr);
        if k>0      
            file2=st{j};
            if isempty(file2)==0
                control=1;  
                set(handles.tagtwocolor,'userdata',file2);
                set(handles.file2,'string',file2);
            end
        end
    end
    if control==0
       % disp('Second image not found')
        %not found: enter second image manually
        [file2,path] = uigetfile('*.*','Choose second image file');
        if file2==0
            return
        end
        set(handles.tagtwocolor,'userdata',file2);
        set(handles.file2,'string',file2);
       % return
    end
elseif option==4
    % fluo image
    
    
end % option color


%pushbuttons & radiobuttons
set (handles.donebutton, 'Enable','on');
set (handles.visupushbutton, 'Enable','on');
set(handles.trypushbutton,'enable','on');

guidata(gcbo,handles) ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function secondpushbutton_Callback(hObject, eventdata, handles)

[file2,path] = uigetfile('*.*','Choose second image file');
if file2==0
    return
end
set(handles.tagtwocolor,'userdata',file2);
set(handles.file2,'string',file2);
guidata(gcbo,handles) ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function donebutton_Callback(hObject, eventdata, handles)

% save data
listafiles=get(handles.loadpushbutton,'userdata');
nrofile=get(handles.trypushbutton,'userdata');
next=nrofile+1;
option=get(handles.popupmenu1,'value');
file=listafiles{nrofile};
[namefile,rem]=strtok(file,'.'); %sin extension
px_mu = str2num(get(handles.szpx,'String'));
dx = str2num(get(handles.accu,'String'));
ident=get(handles.tagtwocolor,'string');


currentdir=cd;
k = strfind(currentdir, 'driftcorrected');
if isempty(k)==1
    if isdir('driftcorrected'); else mkdir('driftcorrected'); end
    cd('driftcorrected')  
end

if option==3
    % second file
    file2=get(handles.tagtwocolor,'userdata');
    [namefile2,rem]=strtok(file2,'.'); %sin extension
end

matrice_results=get(handles.donebutton,'userdata');

% save new data
if isempty(matrice_results)==0
    auxy= matrice_results(2,:);
    auxx= matrice_results(3,:);  
    matrice_results(3,:)=auxy/px_mu; %reconversion
    matrice_results(2,:)=auxx/px_mu;
    
    xdim=ceil(max(matrice_results(2,:)));
    ydim=ceil(max(matrice_results(3,:)));
    
    %correction coordinates!!
    % att: µm!
    minx1=min(handles.x)/px_mu;
    miny1=min(handles.y)/px_mu;
    maxx1=max(handles.x)/px_mu;
    maxy1=max(handles.y)/px_mu;
    
    indexneg=find(matrice_results(2,:)<minx1);
    if isempty(indexneg)==0
        matrice_results(:,indexneg)=NaN;
    end
    indexneg=find(matrice_results(3,:)<miny1);
    if isempty(indexneg)==0
        matrice_results(:,indexneg)=NaN;
    end
    indexextra=find(matrice_results(3,:)>maxx1);
    if isempty(indexextra)==0
        matrice_results(:,indexextra)=NaN;
    end
    indexextra=find(matrice_results(2,:)>maxy1);
    if isempty(indexextra)==0
        matrice_results(:,indexextra)=NaN;
    end
    
    % false points for size
    last=size(matrice_results,2);
    aux= [matrice_results(1,last) miny1 minx1  matrice_results(4,last)  matrice_results(5,last)];
    matrice_results=[matrice_results, aux'];
    aux=[matrice_results(1,last) maxy1 maxx1  matrice_results(4,last)  matrice_results(5,last)];
    matrice_results=[matrice_results, aux'];
    
    if option==2 || option==4
        savename=[namefile, '.mat'];
        save(savename, 'matrice_results'); % one image, corrected for stage drift
        
    elseif option==3
        savename=[namefile,'.mat'];
        savename2=[namefile2,'.mat'];
        save(savename2, 'matrice_results'); % second image, corrected for color drift
        clear x y alpha fr radius matrice_results
        
    end 
    
end

cd(currentdir)

% habilita pasar a la siguiente
if size(listafiles,2)<next
    set(handles.trypushbutton,'enable','off');
else
        set(handles.trypushbutton,'userdata',next);
        file=listafiles{next}; %next file
        set(handles.filename,'string',[file,' (',num2str(next),' of ',num2str(size(listafiles,2)),')']);
        set(handles.file2,'string',[]);

        if option==3   
            firststr=file(1:2) ; %numbering!!!

            control=0;
            % name second image
            d=dir(['*',ident,'*']); %  files
            st={d.name} ;
            for j=1:size(st,2);
                k = findstr(st{j},firststr);
                if k>0      
                    file2=st{j};
                    if isempty(file2)==0
                        control=1;  
                        set(handles.tagtwocolor,'userdata',file2);
                        set(handles.file2,'string',file2);
                    end
                end
            end
            if control==0
                %not found: enter second image manually
                [file2,path] = uigetfile('*.*','Choose second image file');
                if file2==0
                    return
                end
                set(handles.tagtwocolor,'userdata',file2);
                set(handles.file2,'string',file2);
            end
        end % option
end %there is next

clear file listafiles  matrice_results
guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function visupushbutton_Callback(hObject, eventdata, handles)

listafiles=get(handles.loadpushbutton,'userdata');
nrofile=get(handles.trypushbutton,'userdata');
file=listafiles{nrofile};
px_mu = str2num(get(handles.szpx,'String'));
dx = str2num(get(handles.accu,'String'));
file2=get(handles.tagtwocolor,'userdata'); % name second file
option=get(handles.popupmenu1,'value');

handles=loadPALMdata(file,px_mu,0,handles);  
        
xdim=ceil(max(handles.x));
ydim=ceil(max(handles.y));
   
if option==3
    if isempty(file2)==0
        handles2=loadPALMdata(file2,px_mu,0,handles);  
       % showframe(handles,xdim,ydim,[],handles2);
        showframe(handles,ydim,xdim,[],handles2);
    end
else
   % showframe(handles,xdim,ydim);
    showframe(handles,ydim,xdim);
end

clear handles2

guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trypushbutton_Callback(hObject, eventdata, handles)

option=get(handles.popupmenu1,'value');
warning off all
px_mu = str2num(get(handles.szpx,'String'));
dx = str2num(get(handles.accu,'String'));
matrice_results=[];

listafiles=get(handles.loadpushbutton,'userdata');
nrofile=get(handles.trypushbutton,'userdata');

file=listafiles{nrofile};

set(handles.filename,'string',[file,' (',num2str(nrofile),' of ',num2str(size(listafiles,2)),')']);

if option==2 %stage drift
    matrice_results=stagedrift(file,px_mu,handles);
        
elseif option==3
    matrice_results=colorshift(file,px_mu,dx, handles);
    
else
    msgbox('Choose the analysis type')
end
set(handles.donebutton,'userdata',matrice_results);

 guidata(gcbo,handles) ;
   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [matrice_results]=stagedrift(file,px_mu,handles)

handles=loadPALMdata(file,px_mu,0,handles); 
handles.smoothframes=str2num(get(handles.smooth,'string'));

%disp('2nd')
%disp(handles.x(1))

matrice_results=[];
disp(' ')
disp(['Correcting stage drift in file ', file])

matrice_results =correctshiftstage2(handles);  

set(handles.donebutton,'userdata',matrice_results);

if isempty(matrice_results)==0
    data.x=matrice_results(2,:);
    data.y=matrice_results(3,:);
    data.alpha=matrice_results(4,:);
    data.fr=matrice_results(1,:);
    xdim=ceil(max(data.x));
    ydim=ceil(max(data.y));
   % showframe(handles,xdim,ydim,data)
    showframe(handles,ydim,xdim,data)
   % disp('3rd')
%disp(data.x(1))

end
clear data


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [matrice_results]=colorshift(file,px_mu,dx, handles)

firststr=file(1:2) ; %numbering!!!
nroregion=0;  
ident=get(handles.tagtwocolor,'string');
[namefile,rem]=strtok(file,'.');
matrice_results=[];

% name second image
file2=get(handles.tagtwocolor,'userdata'); % name second file
disp(' ')
disp(['Correcting color shift between files ',file, ' and ',file2])

data2=loadPALMdata(file2,px_mu,0,handles);
handles=loadPALMdata(file,px_mu,0,handles);

xdim=ceil(max(handles.x));
ydim=ceil(max(handles.y));

%showframe(handles,xdim,ydim,[],data2);  
showframe(handles,ydim,xdim,[],data2);  
set(handles.textzoom,'string','Please select the region to zoom')

% ROI
[areaselect,xi,yi]=roipolyold;    %seleccion ROI
control=0; 

%data en ROI
if isempty(data2)==0        
        aux = find(inpolygon( handles.x, handles.y,xi,yi) );
        newpoints(:,1)= handles.x(aux);
        newpoints(:,2)= handles.y(aux);
        aux = find(inpolygon( data2.x, data2.y,xi,yi)) ;
        newpoints2(:,1)= data2.x(aux);
        newpoints2(:,2)= data2.y(aux);
        
        if isempty(newpoints)==0
            control=1;
        end
end

set(handles.textzoom,'string',[' '])

if control>0
        % define extremos con las pos min y max de las traj
        minposx=ceil(min(newpoints(:,1)));
        minposy=ceil(min(newpoints(:,2)));
        maxposx=floor(max(newpoints(:,1)));
        maxposy=floor(max(newpoints(:,2)));
        dimx=maxposx-minposx+1;
        dimy=maxposy-minposy+1;  
        
        % correccion coordenadas
        newpoints(:,1)= newpoints(:,1)-minposx+1;
        newpoints(:,2)= newpoints(:,2)-minposy+1;  
        newpoints2(:,1)= newpoints2(:,1)-minposx+1;
        newpoints2(:,2)= newpoints2(:,2)-minposy+1;  

        varargin{1}=namefile; %
        varargin{2}=dimy; %inverted!
        varargin{3}=dimx;
        varargin{4}=newpoints;
        varargin{5}=newpoints2; 
        
        % ventana ROI
        varargout=ROIshiftcorrection(varargin);
        uiwait;    
        
        aux=['auxiliar.mat'];    
        
        if length(dir(aux))>0
            det=load(aux);
            detopt = struct2cell(det);
            % result: shift    
            shift=detopt{1};
            
            if isempty(shift)==0
                data2.x=data2.x+shift(1); 
                data2.y=data2.y+shift(2);
            end
            
            matrice_results(1,:)=data2.fr;
            matrice_results(2,:)= data2.x; 
            matrice_results(3,:)= data2.y;
            matrice_results(4,:)=data2.alpha;
            matrice_results(5,:)=data2.radius;  

        end % dir aux  
else  % control
    msgbox('Error','No data files','error')
end  % control

%showframe(handles,xdim,ydim,[],data2)
showframe(handles,ydim,xdim,[],data2)

clear data2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function showframe(handles,xdim,ydim,data,data2)

if nargin<4
    data=handles;
end
if nargin<5
    data2=[];
end

figure
imshow(zeros(xdim,ydim),'InitialMagnification','fit');
hold on

if isempty(data)==1 && isempty(data2)==1
    return
elseif isempty(data)==1 && isempty(data2)==0;
    data=handles;
end

pctmin = str2num(get(handles.alpha1,'String'));
pctmax = str2num(get(handles.alpha2,'String'));

[X_mu1,Y_mu1]=selectremovealphaquartiles(pctmin,pctmax,data.x,data.y,data.alpha);

% plot positions
if isempty(data2)==0
    [X_mu,Y_mu]=selectremovealphaquartiles(pctmin,pctmax,data2.x,data2.y,data2.alpha);
    plot(X_mu,Y_mu,'.','MarkerSize',2,'Color','r'); 
    hold on
    plot(X_mu1,Y_mu1,'.','MarkerSize',2,'Color','y');
else
    plot(X_mu1,Y_mu1,'.','MarkerSize',2,'Color','w');
end

clear X_mu1 Y_mu1 X_mu Y_mu data data2 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles=loadPALMdata(file,px_mu,dx,handles)

points=[];
S = load(file);

if isfield(S,'Xmatrix')
           % disp(' File with tracking information');
           % Xmatrix = S.Xmatrix * px_mu;
           % Ymatrix = S.Ymatrix * px_mu;
            Xmatrix = S.Xmatrix;
            Ymatrix = S.Ymatrix;
            alphamatrix = S.alphamatrix;
            frmatrix = S.frmatrix;
            if isfield(S,'radiusmatrix')
                radiusmatrix = S.radiusmatrix;
            else
                radiusmatrix=zeros(size(Xmatrix,1),1);
            end
            x = Xmatrix(:);
            y = Ymatrix(:);
            alpha = alphamatrix(:);
            fr = frmatrix(:);
            radius = radiusmatrix(:); %????
            clear S;
            
elseif isfield(S,'matrice_results')
            aux = S.matrice_results;
            clear S;
            fr = aux(1,:); fr = fr(:);
           % x = aux(3,:) * px_mu;
            x = aux(3,:);
            x = x(:);
           % y = aux(2,:) * px_mu; 
            y = aux(2,:); 
            y = y(:);
            alpha = aux(4,:); alpha = alpha(:);
            radius=(aux(5,:));
            
elseif isfield(S,'alpha') && isfield(S,'fr') && isfield(S,'x') && isfield(S,'y')
           % x = S.x * px_mu;
           % y = S.y * px_mu;
            x = S.x ;
            x = x(:);
            y = S.y ;
            y = y(:);
            alpha = S.alpha; alpha = alpha(:);
            fr = S.fr; fr = fr(:);
            if isfield(S,'radius')
                radius=S.radius;
            else
                radius=zeros(size(x,1),1);
            end
            radius=radius(:);
else
            %validfile = 0;
            disp('Invalid file!');
            return
end


if dx>0
    x=x/dx;
    y=y/dx;
end

clear S

handles.x = x;
handles.y = y;
handles.alpha = alpha;
handles.fr = fr;
handles.radius=radius;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function popupmenu1_Callback(hObject, eventdata, handles)

handles.option=get(handles.popupmenu1,'value');

switch handles.option
    case 1 %mock
        %Select item
    case 2 %stage drift
        set(handles.secondpushbutton,'enable','off')
        set(handles.tagtwocolor,'enable','off')
        set(handles.filedefin,'enable','off')
        set(handles.textmerged,'enable','off')
        set(handles.file2,'enable','off')
        set(handles.textsmooth,'enable','on')
        set(handles.smooth,'enable','on')
    case 3 % color drift
        set(handles.secondpushbutton,'enable','on')
        set(handles.tagtwocolor,'enable','on')
        set(handles.filedefin,'enable','on')
        set(handles.textmerged,'enable','on')
        set(handles.file2,'enable','on')
        set(handles.textsmooth,'enable','off')
        set(handles.smooth,'enable','off')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tagtwocolor_Callback(hObject, eventdata, handles)
function szpx_Callback(hObject, eventdata, handles)
function accu_Callback(hObject, eventdata, handles)
function alpha1_Callback(hObject, eventdata, handles)
function alpha2_Callback(~, eventdata, handles)
function smooth_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tagtwocolor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function szpx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function accu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function alpha1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function alpha2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function smooth_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function figure1_ResizeFcn(hObject, eventdata, handles)
%mock
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
