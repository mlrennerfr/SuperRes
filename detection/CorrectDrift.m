function varargout = CorrectDrift(varargin)
% Last Modified by GUIDE v2.5 28-Oct-2015 11:50:56
%
% GUI for stage drift correction
%
% Marianne Renner 2015, for SuperRes programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CorrectDrift_OpeningFcn, ...
                   'gui_OutputFcn',  @CorrectDrift_OutputFcn, ...
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

% --- Executes just before CorrectDrift is made visible.
function CorrectDrift_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
set(handles.visupushbutton,'enable','off');
set(handles.donebutton,'enable','off');
set(handles.trypushbutton,'enable','off');
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = CorrectDrift_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadpushbutton_Callback(hObject, eventdata, handles)

st=[]; 
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
if size(files,2)>1 %batch
  set (handles.filename, 'string',['Batch: File ',file,' (1/',num2str(size(files,2)),')']) ;
else
  set (handles.filename, 'string',file) ;
end

handles=loadPALMdata(file,0,handles);

set(handles.loadpushbutton,'userdata',listafiles);
set(handles.trypushbutton,'userdata',1);

%pushbuttons & radiobuttons
set (handles.donebutton, 'Enable','on');
set (handles.visupushbutton, 'Enable','on');
set(handles.trypushbutton,'enable','on');

guidata(gcbo,handles) ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function donebutton_Callback(hObject, eventdata, handles)

listafiles=get(handles.loadpushbutton,'userdata');
nrofile=get(handles.trypushbutton,'userdata');
next=nrofile+1;
file=listafiles{nrofile};

% habilita pasar a la siguiente
if size(listafiles,2)<next
    set(handles.trypushbutton,'enable','off');
else
    set(handles.trypushbutton,'userdata',next);
    file=listafiles{next}; %next file
    set(handles.filename,'string',[file,' (',num2str(next),' of ',num2str(size(listafiles,2)),')']);
end %there is next

clear file listafiles  x y alpha fr radius matrice_results
guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function visupushbutton_Callback(hObject, eventdata, handles)

listafiles=get(handles.loadpushbutton,'userdata');
nrofile=get(handles.trypushbutton,'userdata');
file=listafiles{nrofile};
handles=loadPALMdata(file,0,handles);  
        
xdim=ceil(max(handles.x));
ydim=ceil(max(handles.y));
   
showframe(handles,ydim,xdim);

clear handles2
guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trypushbutton_Callback(hObject, eventdata, handles)

warning off all
matrice_results=[];

listafiles=get(handles.loadpushbutton,'userdata');
nrofile=get(handles.trypushbutton,'userdata');

file=listafiles{nrofile};

set(handles.filename,'string',[file,' (',num2str(nrofile),' of ',num2str(size(listafiles,2)),')']);

% correction
stagedrift(file,handles);
        
 guidata(gcbo,handles) ;
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stagedrift(file,handles)

handles=loadPALMdata(file,0,handles); 
handles.smoothframes=str2num(get(handles.smooth,'string'));
listafiles=get(handles.loadpushbutton,'userdata');
nrofile=get(handles.trypushbutton,'userdata');
file=listafiles{nrofile};
[namefile,rem]=strtok(file,'.'); %sin extension

matrice_results=[];
disp(' ')
disp(['Correcting stage drift in file ', file])

[matrice_results, controlsave] =correctshiftstage2(handles);  

currentdir=cd;
k = strfind(currentdir, 'driftcorrected');
if isempty(k)==1
    if isdir('driftcorrected'); else mkdir('driftcorrected'); end
end
cd('driftcorrected')  

if isempty(matrice_results)==0
end

if controlsave==1
    
    % save new data
    if isempty(matrice_results)==0
        
        %data.x=matrice_results(2,:);
       % data.y=matrice_results(3,:);
        data.x=matrice_results(3,:);
        data.y=matrice_results(2,:);
        
        data.alpha=matrice_results(4,:);
        data.fr=matrice_results(1,:);
        xdim=ceil(max(data.x));
        ydim=ceil(max(data.y));
        showframe(handles,ydim,xdim,data)
        
       % auxy= matrice_results(2,:);
      %  auxx= matrice_results(3,:);  
      %  matrice_results(3,:)=auxy; 
      %  matrice_results(2,:)=auxx;  
      
      % plot transposed!
        xdim=ceil(max(matrice_results(2,:)));
        ydim=ceil(max(matrice_results(3,:)))  ;
        
        %correction coordinates!! VOIR
        minx1=min(handles.x);
        miny1=min(handles.y);
        maxx1=max(handles.x);
        maxy1=max(handles.y);  
        
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
        
        % false points for size ?????????????????????????????????
         last=size(matrice_results,2);
         if size(matrice_results,1)>7
             aux= [matrice_results(1,last) miny1 minx1  matrice_results(4,last)  matrice_results(5,last) matrice_results(6,last) matrice_results(7,last)...
                 matrice_results(8,last) matrice_results(9,last) matrice_results(10,last) matrice_results(11,last)];
             matrice_results=[matrice_results, aux'];
             aux=[matrice_results(1,last) maxy1 maxx1  matrice_results(4,last)  matrice_results(5,last) matrice_results(6,last) matrice_results(7,last)...
                 matrice_results(8,last) matrice_results(9,last) matrice_results(10,last) matrice_results(11,last)];
             matrice_results=[matrice_results, aux'];  
         else
             aux= [matrice_results(1,last) miny1 minx1  matrice_results(4,last)  matrice_results(5,last) matrice_results(6,last) matrice_results(7,last)];
             matrice_results=[matrice_results, aux'];
             aux=[matrice_results(1,last) maxy1 maxx1  matrice_results(4,last)  matrice_results(5,last) matrice_results(6,last) matrice_results(7,last)];
             matrice_results=[matrice_results, aux'];  
         end
        
        clear pk
        pk(:,1)=matrice_results(1,:); 
     %   pk(:,2)=matrice_results(2,:);  
     %   pk(:,3)=matrice_results(3,:);    
        pk(:,2)=matrice_results(3,:);  
        pk(:,3)=matrice_results(2,:);  
        
        pk(:,5)= matrice_results(4,:); 
        pk(:,4)=matrice_results(5,:)*2; 
        pk(:,7)=matrice_results(6,:);  
        pk(:,6)=matrice_results(7,:);
        
        if size(matrice_results,1)>7
            pk(:,8)=matrice_results(8,:); %ratio
            pk(:,9)=matrice_results(9,:); %z
            pk(:,10)=matrice_results(10,:); %test
            pk(:,11)=matrice_results(11,:); %test

            if max(pk(:,9))==0 && min(pk(:,9))==0 % no values for z
                if isdir('pk'); else mkdir('pk');end
                pkpath=['pk'];
                savename=[namefile, '.pk'];
            else
                if isdir('pk3'); else mkdir('pk3');end
                pkpath=['pk3'];
                savename=[namefile, '.pk3'];
                % prepare file for visp
                indexoui=find(pk(:,8)<1000);
                disp('Pixel size for Visp : 160 nm')
                szpxum=160/1000;
                if isempty(indexoui)==0
                    pkvisp=[pk(indexoui,2)*szpxum pk(indexoui,3)*szpxum pk(indexoui,9)*szpxum pk(indexoui,5) pk(indexoui,1)];    %[x y z alpha fr];
                    save([namefile, '.3d'], 'pkvisp','-ascii');
                    clear pkvisp indexoui
                end
            end
        else
            if isdir('pk'); else mkdir('pk');end
            pkpath=['pk'];
            savename=[namefile, '.pk'];
        end

        save([namefile, '.mat'], 'matrice_results'); % one image, corrected for stage drift
        
        %currentdir=cd;
        
        cd(pkpath);
        save(savename, 'pk','-ascii');
       % cd(currentdir)

    end %empty matrice_results
end % controlsave


cd(currentdir)

clear pk data matrice_results

 
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

%function handles=loadPALMdata(file,px_mu,dx,handles)
function handles=loadPALMdata(file,dx,handles)

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
            sigma=(aux(6,:));
            blink=(aux(7,:));
            if size(aux,1)>7
                ratio=aux(8,:);
                z=aux(9,:);
                test1=aux(10,:);
                test2=aux(11,:);
            else
                ratio=[];
                z=[];
                test1=[];
                test2=[];
            end

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
handles.sigma=sigma;
handles.blink=blink;
handles.ratio=ratio;
handles.z=z;
handles.test1=test1;
handles.test2=test2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function quitpushbutton_Callback(hObject, eventdata, handles)

close
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function szpx_Callback(hObject, eventdata, handles)
function accu_Callback(hObject, eventdata, handles)
function alpha1_Callback(hObject, eventdata, handles)
function alpha2_Callback(~, eventdata, handles)
function smooth_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
function smooth_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function figure1_ResizeFcn(hObject, eventdata, handles)
%mock
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

