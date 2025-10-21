function varargout = ConvertData(varargin)
% function ConvertData
%converts .csv file with SMLM detections into .mat format usable by Diinamic
% 
%input: 
%.csv file 
%
% output:
%.mat file with the structure (all sizes in µm)
%
%    matrice_results(1,:)=fr;
%    matrice_results(3,:)= x; 
%    matrice_results(2,:)= y;
%    matrice_results(4,:)=alpha;
%    matrice_results(5,:)=radius;  
%    matrice_results(6,:)=sigma;  
%    matrice_results(7,:)=blink;  
%    matrice_results(8,:)=ratio;  
%    matrice_results(9,:)=z;  
%    matrice_results(10,:)=test1;  
%    matrice_results(11,:)=test2;
%
% Note:
% for example ThunderSTORM-type data :
%    col1=id
%    col2=fr
%    col3=x
%    col4=y
%    col5=z    
%    col6=sigma1
%    col7=sigma2
%    col8:intensity
%    col9:offset
%    col10
%    col11:chi
%    col12: uncertainty

% Marianne Renner 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ConvertData_OpeningFcn, ...
                   'gui_OutputFcn',  @ConvertData_OutputFcn, ...
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ConvertData_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

% to fit screen definition/changes in size
h3 = findobj('Type','figure');
txtHand = findall(h3, '-property', 'FontUnits');
%set(txtHand, 'FontUnits', 'normalized');
set(txtHand, 'FontUnits', 'centimeters');
    
guidata(hObject, handles);

function varargout = ConvertData_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadpushbutton_Callback(hObject, eventdata, handles)

disp(' ')
disp('Convertion of SMLM data from .csv to .mat format')
disp(' ')
disp('Indicate the data columns in your file. Coordinates (x,y) are compulsory, the other values are optional. If a given type of data is not present in your file, put zero as number of column')

d = dir('*.csv'); %
st={d.name};
if isempty(st)==1
   msgbox('No files!!','Select files','error');
   return
end
[listafiles,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',st);
if v==0
   return
end

set(handles.loadpushbutton,'userdata',listafiles);
set(handles.matfilepushbutton,'userdata',st);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function matfilepushbutton_Callback(hObject, eventdata, handles)


files=get(handles.loadpushbutton,'userdata');
st=get(handles.matfilepushbutton,'userdata');
framecol=str2num(get(handles.framecol,'string'));
xcol=str2num(get(handles.xcol,'string'));
ycol=str2num(get(handles.ycol,'string'));
intcol=str2num(get(handles.intcol,'string'));
offsetcol=str2num(get(handles.offsetcol,'string'));
radcol=str2num(get(handles.radcol,'string'));
sigmacol=str2num(get(handles.sigmacol,'string'));
testcol=str2num(get(handles.testcol,'string'));
factor=str2num(get(handles.factor,'string'));

for i=1:size(files,2)   
    
    file=st{files(i)};
    disp(' ')
    disp(['Reading data in ',file,', please wait...'])

    [namefile,rem]=strtok(file,'.');
    savename=[namefile,'.mat'];
    
    data = csvread(file,1,0);
    
    
    x=(data(:,xcol)+abs(min(data(:,xcol))))*factor; % in µm with correction for negative coordinates
    y=(data(:,ycol)+abs(min(data(:,ycol))))*factor; % in µm
    
    if intcol>0
        alpha=data(:,intcol);
    else
        alpha=1000*ones(size(data,1),1); 
    end
    if framecol>0
        fr=data(:,framecol);
    else
        fr=1:size(data,1);
    end
    if radcol>0
        radius=data(:,radcol)/2;
    else
        radius=zeros(size(data,1),1); 
    end
    if sigmacol>0
        sigma=data(:,sigmacol);
    else
        sigma=zeros(size(data,1),1); 
    end
    blink=ones(size(data,1),1); 
    ratio=ones(size(data,1),1);  %data(:,6)/data(:,7);
    %z=data(:,5)*factor; % in µm
    z=ones(size(data,1),1); 
    if testcol>0
        test1=data(:,testcol);
    else
        test1=zeros(size(data,1),1);
    end
    test2=zeros(size(data,1),1);

    matrice_results=zeros(11,size(data,1));
    matrice_results(1,:)=fr;
    matrice_results(3,:)= x; 
    matrice_results(2,:)= y;
    matrice_results(4,:)=alpha;
    matrice_results(5,:)=radius;  
    matrice_results(6,:)=sigma;  
    matrice_results(7,:)=blink;  
    matrice_results(8,:)=ratio;  
    matrice_results(9,:)=z;  
    matrice_results(10,:)=test1;  
    matrice_results(11,:)=test2;


    save(savename, 'matrice_results');
    disp(['Data saved as ',savename])
    disp(' ')

    clear data
    clear matrice_results

end

disp('Done')

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function framecol_Callback(hObject, eventdata, handles)
function framecol_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function xcol_Callback(hObject, eventdata, handles)
function xcol_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ycol_Callback(hObject, eventdata, handles)
function ycol_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function intcol_Callback(hObject, eventdata, handles)
function intcol_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function offsetcol_Callback(hObject, eventdata, handles)
function offsetcol_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function radcol_Callback(hObject, eventdata, handles)
function radcol_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function sigmacol_Callback(hObject, eventdata, handles)
function sigmacol_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function testcol_Callback(hObject, eventdata, handles)
function testcol_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function factor_Callback(hObject, eventdata, handles)
function factor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



