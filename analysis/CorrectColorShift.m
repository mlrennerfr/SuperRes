function varargout = CorrectColorShift(varargin)
%CORRECTCOLORSHIFT M-file for CorrectColorShift.fig
% 
% GUI for color shift correction, on SMLM data
% Marianne Renner, for SuperRes programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last Modified by GUIDE v2.5 01-Jul-2019 10:19:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CorrectColorShift_OpeningFcn, ...
                   'gui_OutputFcn',  @CorrectColorShift_OutputFcn, ...
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
function CorrectColorShift_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
set(handles.visupushbutton,'enable','off');
set(handles.donebutton,'enable','off');
set(handles.trypushbutton,'enable','off');
set(handles.tagcolor1,'userdata',' '); % name second file
set(handles.tagcolor2,'userdata',' '); % name second file
set(handles.tagcolor1,'value',0);
set(handles.alreadycorrradiobutton,'value',0);
set(handles.alreadycorrradiobutton,'enable','off');


guidata(hObject, handles);

function varargout = CorrectColorShift_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function figure1_SizeChangedFcn(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tagcolor2_Callback(hObject, eventdata, handles)
function tagcolor2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function szpx_Callback(hObject, eventdata, handles)
function szpx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function accu_Callback(hObject, eventdata, handles)
function accu_CreateFcn(hObject, eventdata, handles)
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

function popupmenu1_Callback(hObject, eventdata, handles)
function popupmenu1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tagcolor1_Callback(hObject, eventdata, handles)
function tagcolor1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadpushbutton_Callback(hObject, eventdata, handles)

option=get(handles.popupmenu1,'value');
st=[]; 
ident1=get(handles.tagcolor1,'string');
ident2=get(handles.tagcolor2,'string');
px_mu = str2num(get(handles.szpx,'String'));
dx = str2num(get(handles.accu,'String'));
set(handles.trypushbutton,'string','Try');
set(handles.alreadycorrradiobutton,'value',0);


% carga files y loop general 
dialog_title=['Select data folder'];
directory_name = uigetdir(cd,dialog_title);
if directory_name==0
    return
end
cd(directory_name);

%d = dir('*mat*'); % movie
d = dir(['*',ident1,'.mat']); % movie
st={d.name};

if isempty(st)==1
   %msgbox(['No files!!'],'Select files','error');
   %return
   d = dir('*mat*'); % movie
   st={d.name};

end

[files,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',st);
if v==0
   return
end

for i=1:size(files,2)
      listafiles{i}=st{files(i)};
end

file=listafiles{1}; %first file
if size(files,2)>1 %batch
  set (handles.filename, 'string',['Batch: File ',file,' (1/',num2str(size(files,2)),')']) ;
else
  set (handles.filename, 'string',file) ;
end
set(handles.file2,'string',[]);

data1=loadPALMdatanew(file,0); %first file
set(handles.tagcolor1,'userdata',data1);
set(handles.loadpushbutton,'userdata',listafiles);
set(handles.trypushbutton,'userdata',1);

%name second file
[namefile,rem]=strtok(file,'.');
k = findstr(namefile, ident1);
rootname=namefile(1:k-1);

if option==1
    control=0;
    % name second image
    file2=[rootname,ident2,'.mat']
    if isempty(dir(file2))==0
        control=1;  
    end
    if control==0
        %not found: enter second image manually
        [file2,path] = uigetfile('*.mat','Choose second image file');
        if file2==0
            return
        end
       % return
    end
    %read second image
    data2=loadPALMdatanew(file2,0); %first file
    set(handles.tagcolor2,'userdata',data2);
    
elseif option==2    
    
    control=0;
    % name second image
    file2=[rootname,ident2,'.tif']
    if isempty(dir(file2))==0
        control=1;  
    end

    if control==0
        %not found: enter second image manually
        [file2,path] = uigetfile('*.tif','Choose second image file');
        if file2==0
            return
        end
       % return
    end
    
   % xdim=ceil(max(data1.x)/dx)
  %  ydim=ceil(max(data1.y)/dx)

    [~,datamatrix] = tifdataread(file2);


    if isfield(datamatrix,'data')
        dat=datamatrix.data;
        datamatrix=double(dat);
    end

    data2 = datamatrix;
    
    %  disp('originals')
    xdim=ceil(max(data1.x));
    ydim=ceil(max(data1.y));
   % disp('fluo')
   % disp(size(datamatrix))
    
    if size(datamatrix,1)>max(data1.y) %data detections in µm
        data1.x=data1.x/px_mu;
        data1.y=data1.y/px_mu;
        set(handles.tagcolor1,'userdata',data1);
        set(handles.tagcolor1,'value',1);
        data2 = imresize(datamatrix,[ceil(max(data1.y)),ceil(max(data1.x))]);
    end
    
    %  disp('converted')
    %    disp(max(data1.x))
    %    disp(max(data1.y))
    %    disp('fluo')
    %    disp(size(data2))

end % option color

set(handles.tagcolor2,'userdata',data2);
set(handles.file2,'String',file2);

%pushbuttons & radiobuttons
if size(listafiles,2)>1
    set (handles.donebutton, 'Enable','on');
end
set (handles.visupushbutton, 'Enable','on');
set(handles.trypushbutton,'enable','on');

clear data1 data2 datamatrix BWmatrix newBWmatrix

guidata(gcbo,handles) ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function visupushbutton_Callback(hObject, eventdata, handles)

listafiles=get(handles.loadpushbutton,'userdata');
nrofile=get(handles.trypushbutton,'userdata');
file=listafiles{nrofile};
px_mu = str2num(get(handles.szpx,'String'));
dx = str2num(get(handles.accu,'String'));
data1=get(handles.tagcolor1,'userdata'); % data first file
data2=get(handles.tagcolor2,'userdata'); % data second file
option=get(handles.popupmenu1,'value');

previouscorr=get(handles.alreadycorrradiobutton,'value');
if previouscorr==1
    if option==1 %two sets of detections
        data2=get(handles.alreadycorrradiobutton,'userdata'); % data second file
    elseif option==2 %one set of detections and one fluo file
        data1=get(handles.alreadycorrradiobutton,'userdata'); % data first file
    end
end


xdim=ceil(max(data1.x));
ydim=ceil(max(data1.y));
   
if isempty(data2)==0
   showframe(handles,ydim,xdim,option,data1,data2);
end

clear data1 data2

guidata(gcbo,handles) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trypushbutton_Callback(hObject, eventdata, handles)

option=get(handles.popupmenu1,'value');
warning off all
px_mu = str2num(get(handles.szpx,'String'));
dx = str2num(get(handles.accu,'String'));
matrice_results=[];
data1=get(handles.tagcolor1,'userdata'); % data first file
data2=get(handles.tagcolor2,'userdata'); % data second file

file2=get(handles.file2,'String');
[namefile2,~]=strtok(file2,'.'); 

listafiles=get(handles.loadpushbutton,'userdata');
nrofile=get(handles.trypushbutton,'userdata');

file=listafiles{nrofile};
[namefile,~]=strtok(file,'.'); 

set(handles.filename,'string',[file,' (',num2str(nrofile),' of ',num2str(size(listafiles,2)),')']);
clear matrice_results

% routines correction (below) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if option==1 %two sets of detections
    
    matrice_results=colorshift(file,px_mu,dx, handles);
    
elseif option==2 %one set of detections and one fluo file
    
    matrice_results=colorshiftfluo(file,px_mu,dx, handles);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %   disp('converted corrected')
   %     disp(max(matrice_results(2,:)))
   %     disp(max(matrice_results(3,:)))
   %     disp('fluo')
   %     disp(size(data2))

currentdir=cd;

% save new data
if isempty(matrice_results)==0
    
    k = strfind(currentdir, 'colorcorrected');
    if isempty(k)==1
        if isdir('colorcorrected'); else mkdir('colorcorrected'); end
        cd('colorcorrected')  
    end

    % image that moved
    auxy= matrice_results(2,:);
    auxx= matrice_results(3,:);  
    matrice_results(3,:)=auxy; %reconversion
    matrice_results(2,:)=auxx;
    
    %correction coordinates!!
    minx1=min(data1.x);
    miny1=min(data1.y);
    maxx1=max(data1.x);
    maxy1=max(data1.y);
    
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
    
    if option==2    %correct for convertion to µm
        matrice_results(2,:)= matrice_results(2,:).* px_mu; %data1.x; %%%%%%%%%%%%
        matrice_results(3,:)= matrice_results(3,:).* px_mu; %data1.y;
        
       % disp('converted re-corrected')
       % disp(max(matrice_results(2,:)))
       % disp(max(matrice_results(3,:)))
       % disp('fluo')
       % disp(size(data2))
    end
    
    
    % false points for size
    last=size(matrice_results,2);
    aux= [matrice_results(1,last) miny1*px_mu minx1*px_mu  matrice_results(4,last)  matrice_results(5,last)...
        matrice_results(6,last) matrice_results(7,last) matrice_results(8,last) matrice_results(9,last) matrice_results(10,last) matrice_results(11,last)];
    matrice_results=[matrice_results, aux'];
    aux=[matrice_results(1,last) maxy1*px_mu maxx1*px_mu  matrice_results(4,last)  matrice_results(5,last)...
                matrice_results(6,last) matrice_results(7,last) matrice_results(8,last) matrice_results(9,last) matrice_results(10,last) matrice_results(11,last)];
    matrice_results=[matrice_results, aux'];

    clear pk
    pk(:,1)=matrice_results(1,:); 
    pk(:,2)=matrice_results(3,:);  %!!!!!!!!!!!
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
            pkpath='pk';
            extpksavename='.pk';
        else
            if isdir('pk3'); else mkdir('pk3');end
            pkpath='pk3';
            extpksavename='.pk3';
        end
    else
        if isdir('pk'); else mkdir('pk');end
        pkpath='pk';
        extpksavename='.pk';
    end
    
    if option==1 %second image moved
        savename2=[namefile2,'.mat'];
         pksavename=[namefile2, extpksavename];
    elseif option==2 %first image move
        savename2=[namefile,'.mat'];
        pksavename=[namefile, extpksavename];
    end
    
    save(savename2, 'matrice_results'); %  image corrected for color drift
  %  set(handles.donebutton,'userdata',matrice_results);

    cd(pkpath);
    save(pksavename, 'pk','-ascii');
    clear x y alpha fr radius matrice_results pk
    
    %image 1 not moved
    if option==1
        cd(currentdir)
        savename=[namefile,'.mat'];
        
        cd('colorcorrected')  
        matrice_results(1,:)=data1.fr;
        matrice_results(2,:)= data1.y; 
        matrice_results(3,:)= data1.x;
        matrice_results(4,:)=data1.alpha;
        matrice_results(6,:)=data1.radius;  
        matrice_results(7,:)=data1.sigma;  
        matrice_results(8,:)=data1.blink;  
        matrice_results(9,:)=data1.z;  
        matrice_results(10,:)=data1.test1;  
        matrice_results(11,:)=data1.test2;  
        save(savename, 'matrice_results'); % second image, corrected for color drift
        
        %pk
        pk(:,1)=matrice_results(1,:); 
        pk(:,2)=matrice_results(3,:);  %%%
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
                pkpath='pk';
                pksavename=[namefile, '.pk'];
            else
                if isdir('pk3'); else mkdir('pk3');end
                pkpath='pk3';
                pksavename=[namefile, '.pk3'];
            end
        else
            if isdir('pk'); else mkdir('pk');end
            pkpath='pk';
            pksavename=[namefile, '.pk'];
        end
        cd(pkpath);
        save(pksavename, 'pk','-ascii');
        
    else
        % image fluo
        cd(currentdir)
        %cd('colorcorrected')  
        copyfile(file2,[currentdir,'\colorcorrected'])
        
    end %if option
        
end

cd(currentdir)

clear x y alpha fr radius matrice_results pk

 guidata(gcbo,handles) ;
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function donebutton_Callback(hObject, eventdata, handles)
% next file

% save data
ident1=get(handles.tagcolor1,'string');
ident2=get(handles.tagcolor2,'string');
listafiles=get(handles.loadpushbutton,'userdata');
nrofile=get(handles.trypushbutton,'userdata');
next=nrofile+1;
option=get(handles.popupmenu1,'value');
px_mu = str2num(get(handles.szpx,'String'));
dx = str2num(get(handles.accu,'String'));
set(handles.alreadycorrradiobutton,'enable','off');
set(handles.alreadycorrradiobutton,'value',0);


% habilita pasar a la siguiente
if size(listafiles,2)<next
    set(handles.trypushbutton,'enable','off');
else
    set(handles.trypushbutton,'userdata',next);
    set(handles.tagcolor1,'userdata',[]); % data first file
    set(handles.tagcolor2,'userdata',[]); % data second file

    %first file
    file1=listafiles{next} %
    data1=loadPALMdatanew(file1,0); 
    set(handles.tagcolor1,'userdata',data1);

    set(handles.filename,'string',[file1,' (',num2str(next),' of ',num2str(size(listafiles,2)),')']);
    set(handles.file2,'string',[]);  
    
    %name second file
    [namefile,rem]=strtok(file1,'.');
    k = findstr(namefile, ident1);
    rootname=namefile(1:k-1);
    
    if option==1
        control=0;
        % name second image
        file2=[rootname,ident2,'.mat']
        if isempty(dir(file2))==0
            control=1;  
        end
        if control==0
            %not found: enter second image manually
            [file2,path] = uigetfile('*.mat','Choose second image file');
            if file2==0
                return
            end
            % return
        end
        %read second image
        data2=loadPALMdatanew(file2,0); 
        set(handles.tagcolor2,'userdata',data2);
        
    elseif option==2    
        
        control=0;
        % name second image
        file2=[rootname,ident2,'.tif']
        if isempty(dir(file2))==0
            control=1;  
        end
        
        if control==0
            %not found: enter second image manually
            [file2,path] = uigetfile('*.tif','Choose second image file');
            if file2==0
                return
            end
            % return
        end
        
      %  xdim=ceil(max(data1.x)/dx);
       % ydim=ceil(max(data1.y)/dx);
        
      [~,datamatrix] = tifdataread(file2);
      datamatrix=double(datamatrix.data);
      data2 = datamatrix;
      
      xdim=ceil(max(data1.x));
      ydim=ceil(max(data1.y));

      if size(datamatrix,1)>max(data1.y) %data detections in µm
          data1.x=data1.x/px_mu;
          data1.y=data1.y/px_mu;
          set(handles.tagcolor1,'userdata',data1);
          set(handles.tagcolor1,'value',1);
          data2 = imresize(datamatrix,[ceil(max(data1.y)),ceil(max(data1.x))]);
      end
        
    end % option color
    
    set(handles.file2,'string',file2);
    set(handles.tagcolor2,'userdata',data2);

end %there is next

disp(' ')

clear file listafiles  matrice_results data1 data2 datamatrix BWmatrix newBWmatrix
guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function secondpushbutton_Callback(hObject, eventdata, handles)

option=get(handles.popupmenu1,'value');

if option==1
   [file2,path] = uigetfile('*.mat','Choose second image file');
   if file2==0
       return
   end
elseif option==2
    [file2,path] = uigetfile('*.tif','Choose second image file');
    if file2==0
        return
    end
    [stack_info,datamatrix] = tifdataread(file2);
end

set(handles.file2,'string',file2);


guidata(gcbo,handles) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [matrice_results]=colorshift(file,px_mu,dx, handles)

data1=get(handles.tagcolor1,'userdata'); % data first file
previouscorr=get(handles.alreadycorrradiobutton,'value');
if previouscorr==0
    data2=get(handles.tagcolor2,'userdata'); % data second file
else
    data2=get(handles.alreadycorrradiobutton,'userdata');%already corrected
end

[namefile,~]=strtok(file,'.');
matrice_results=[];

xdim=ceil(max(data1.x));
ydim=ceil(max(data2.y));

showframe(handles,ydim,xdim,1,data1,data2);  
title('Please select the region to zoom');

% ROI
%[~,xi,yi]=roipolyold;    %seleccion ROI
[~,xi,yi]=roipoly;    %seleccion ROI
control=0; 

%data en ROI
if isempty(data2)==0        
        aux = find(inpolygon( data1.x, data1.y,xi,yi) );
        newpoints(:,1)= data1.x(aux);
        newpoints(:,2)= data1.y(aux);
        aux = find(inpolygon( data2.x, data2.y,xi,yi)) ;
        newpoints2(:,1)= data2.x(aux);
        newpoints2(:,2)= data2.y(aux);
        if isempty(newpoints)==0
            control=1;
        end
end

if control>0
        % define extremos con las pos min y max de las traj
      %  minposx=ceil(min(newpoints(:,1)));
      %  minposy=ceil(min(newpoints(:,2)));
      %  maxposx=floor(max(newpoints(:,1)));
      %  maxposy=floor(max(newpoints(:,2)));
        minposx=floor(min(newpoints(:,1)));
        minposy=floor(min(newpoints(:,2)));
        maxposx=ceil(max(newpoints(:,1)));
        maxposy=ceil(max(newpoints(:,2)));
        
       
        dimx=maxposx-minposx+1;
        dimy=maxposy-minposy+1;  
        
        % correccion coordenadas
        newpoints(:,1)= newpoints(:,1)-minposx+1;
        newpoints(:,2)= newpoints(:,2)-minposy+1;  
        newpoints2(:,1)= newpoints2(:,1)-minposx+1;
        newpoints2(:,2)= newpoints2(:,2)-minposy+1;  

        varargin{1}=namefile; %
      %  varargin{2}=dimy; %inverted!
       % varargin{3}=dimx;
        varargin{2}=dimy; %inverted!
        varargin{3}=dimx;
        varargin{4}=newpoints;
        varargin{5}=newpoints2; 
        
        % ventana ROI
        varargout=ROIshiftcorrection(varargin);
        uiwait;    
        
        aux='auxiliar.mat';
        clear newdata2
        
        if length(dir(aux))>0
            det=load(aux);
            detopt = struct2cell(det);
            % result: shift    
            shift=detopt{1};
            newdata2=data2;
            
            if isempty(shift)==0
                newdata2.x=data2.x+shift(1); 
                newdata2.y=data2.y+shift(2);
            end
            
            showframe(handles,ydim,xdim,1,data1,newdata2)
            
            %Accept?
            button = questdlg('Accept?','','Yes','No','Yes');
            resp = strcmp(button,'No');
            if resp==0 
                
                limitxy = inpolygon(newdata2.x,newdata2.y,[0 max(data1.x)],[0 max(data1.y)]); % eliminates detections out of data1 range
                clear data2
                
                data2.fr=newdata2.fr(limitxy);
                data2.x=newdata2.x(limitxy); 
                data2.y=newdata2.y(limitxy);
                data2.alpha=newdata2.alpha(limitxy);
                data2.radius=newdata2.radius(limitxy); 
                data2.sigma=newdata2.sigma(limitxy);
                data2.blink=newdata2.blink(limitxy);  
                data2.z=newdata2.z(limitxy);  
                data2.test1=newdata2.test1(limitxy); 
                data2.test2=newdata2.test2(limitxy);  
                
            end
            
            matrice_results(1,:)=data2.fr;
          %  matrice_results(2,:)= data2.y; %%%%%%%%
          %  matrice_results(3,:)= data2.x;
            matrice_results(2,:)= data2.x; %%%%%%%%
            matrice_results(3,:)= data2.y;
            matrice_results(4,:)=data2.alpha;
            matrice_results(6,:)=data2.radius;  
            matrice_results(7,:)=data2.sigma;  
            matrice_results(8,:)=data2.blink;  
            matrice_results(9,:)=data2.z;  
            matrice_results(10,:)=data2.test1;  
            matrice_results(11,:)=data2.test2;  
            
            set(handles.alreadycorrradiobutton,'userdata',data2);
            set(handles.alreadycorrradiobutton,'enable','on');


        end % dir aux  
else  % control
    msgbox('Error','No data files','error')
end  % control

clear data2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [matrice_results]=colorshiftfluo(file,px_mu,dx, handles)

data2=get(handles.tagcolor2,'userdata'); % data second file (fluo)
previouscorr=get(handles.alreadycorrradiobutton,'value');
if previouscorr==0
    data1=get(handles.tagcolor1,'userdata'); % data first file
else
    data1=get(handles.alreadycorrradiobutton,'userdata'); %already corrected
end



[namefile,rem]=strtok(file,'.');
matrice_results=[];

xdim=ceil(max(data1.x));
ydim=ceil(max(data1.y));

showframe(handles,ydim,xdim,2,data1,data2);  
title('Please select the region to zoom');

% ROI
%[areaselect,xi,yi]=roipolyold;    %seleccion ROI
[areaselect,xi,yi]=roipoly;    %seleccion ROI
control=0; 
minposx=max(floor(min(xi)),1);
minposy=max(floor(min(yi)),1);
maxposx=min(ceil(max(xi)), max(data1.x));
maxposy=min(ceil(max(yi)), max(data1.y));
dimx=maxposx-minposx+1;
dimy=maxposy-minposy+1;  

%data en ROI
if isempty(data1)==0        
        aux = find(inpolygon( data1.x, data1.y,xi,yi));
        newpoints(:,1)= data1.x(aux);
        newpoints(:,2)= data1.y(aux);
        if isempty(newpoints)==0
            control=1;
        end
end
% image fluo
data2zoom=data2(minposy:maxposy,minposx:maxposx);

if control>0
        
        % correccion coordenadas
        newpoints(:,1)= newpoints(:,1)-minposx+1;
        newpoints(:,2)= newpoints(:,2)-minposy+1;  

        varargin{1}=namefile; %
        varargin{2}=dimy; %inverted!
        varargin{3}=dimx;
        varargin{4}=newpoints;
        varargin{5}=data2zoom; 
        
        % ventana ROI
        varargout=ROIshiftcorrectionfluo(varargin);
        uiwait;    
        
        aux=['auxiliar.mat'];    
        
        if length(dir(aux))>0
            det=load(aux);
            detopt = struct2cell(det);
            % result: shift    
            shift=detopt{1};
            newdata1=data1;
            
            if isempty(shift)==0
                newdata1.x=data1.x+shift(1); 
                newdata1.y=data1.y+shift(2);
            end
            
            showframe(handles,ydim,xdim,2,newdata1,data2)
            
            %Accept?
            button = questdlg('Accept?','','Yes','No','Yes');
            resp = strcmp(button,'No');
            
            if resp==0  
                data1=newdata1; 
            end
            
            matrice_results(1,:)=data1.fr;
         %   matrice_results(2,:)= data1.y; %%%%%%%%%%%%
         %   matrice_results(3,:)= data1.x;
            matrice_results(2,:)= data1.x; %%%%%%%%%%%%
            matrice_results(3,:)= data1.y;
            matrice_results(4,:)=data1.alpha;
            matrice_results(6,:)=data1.radius;  
            matrice_results(7,:)=data1.sigma;  
            matrice_results(8,:)=data1.blink;  
            matrice_results(9,:)=data1.z;  
            matrice_results(10,:)=data1.test1;  
            matrice_results(11,:)=data1.test2; 
            
            set(handles.alreadycorrradiobutton,'userdata',data1);
            set(handles.alreadycorrradiobutton,'enable','on');


        end % dir aux  
else  % control
    msgbox('Error','No data files','error')
end  % control

clear data2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function showframe(handles,xdim,ydim,option,data,data2)

ident1=get(handles.tagcolor1,'string');
ident2=get(handles.tagcolor2,'string');
ident1=ident1(2:size(ident1,2));
ident2=ident2(2:size(ident2,2));

if nargin<4
    data=handles;
end
if nargin<5
    data2=[];
end

if isempty(data)==1 && isempty(data2)==1
    return
elseif isempty(data)==1 && isempty(data2)==0;
    data=handles;
end
pctmin = str2num(get(handles.alpha1,'String'));
pctmax = str2num(get(handles.alpha2,'String'));

[X_mu1,Y_mu1]=selectremovealphaquartiles(pctmin,pctmax,data.x,data.y,data.alpha);

figure
hold on

if option==1
    imshow(zeros(xdim,ydim),'InitialMagnification','fit');
    % plot positions
    if isempty(data2)==0
        [X_mu,Y_mu]=selectremovealphaquartiles(pctmin,pctmax,data2.x,data2.y,data2.alpha);
        plot(X_mu,Y_mu,'.','MarkerSize',2,'Color','y'); %second image
      %  hold on
    end
elseif option==2
    stackmin=double(min(min(data2)));
    stackmax=double(max(max(data2)));
    imshow(data2,[stackmin stackmax],'InitialMagnification','fit');
end
plot(X_mu1,Y_mu1,'.','MarkerSize',2,'Color','r'); %data first image

%if option==1
   % legend([ident1,' in red'],[ident2,' in yellow'],'Location','southoutside','Orientation','horizontal');legend('boxoff')
     %   xlabel(['x (mu) (',ident1,' in red; ',ident2,' in yellow)']);    

%else
    xlabel('x (mu)');

%end

ylabel('y (mu)');    axis off ;
xlim([0 max(X_mu1)]);    ylim([0 max(Y_mu1)]);    
hold on;

clear X_mu1 Y_mu1 X_mu Y_mu data data2 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dataread=loadPALMdatanew(file,dx)

points=[];
S = load(file);
dataread.x=[];
dataread.y=[];
           
if isfield(S,'matrice_results')
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
    
else
    disp('Invalid file!');
    return
end


if dx>0
    x=x/dx;
    y=y/dx;
end

dataread.x = x;
dataread.y = y;
dataread.alpha = alpha;
dataread.fr = fr;
dataread.radius=radius;
dataread.sigma=sigma;
dataread.blink=blink;
dataread.ratio=ratio;
dataread.z=z;
dataread.test1=test1;
dataread.test2=test2;

clear aux


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function alreadycorrradiobutton_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%