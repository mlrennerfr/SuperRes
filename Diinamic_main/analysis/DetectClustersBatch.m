function varargout = DetectClustersBatch(varargin)
%function varargout = DetectClustersBatch(varargin)
% GUI for setting threshold and select clusters of detections in batch mode
%
% called by Cluster.m
%
% Marianne Renner oct22
% Last Modified by GUIDE v2.5 29-Sep-2022 12:03:42
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DetectClustersBatch_OpeningFcn, ...
                   'gui_OutputFcn',  @DetectClustersBatch_OutputFcn, ...
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

function DetectClustersBatch_OpeningFcn(hObject, eventdata, handles, varargin)
%
handles.output = hObject;

% to fit screen definition/changes in size
h3 = findobj('Type','figure');
txtHand = findall(h3, '-property', 'FontUnits');
%set(txtHand, 'FontUnits', 'normalized');
set(txtHand, 'FontUnits', 'centimeters');

disp(' ')
disp('Reading data, please wait...')
disp(' ')

handles.listafiles=cell2mat(varargin{1}(1));
handles.dx=cell2mat(varargin{1}(2));  
handles.szpx=cell2mat(varargin{1}(3));  
handles.mu=cell2mat(varargin{1}(4));  
handles.alpha1=cell2mat(varargin{1}(5));  
handles.alpha2=cell2mat(varargin{1}(6));  
Densparam=cell2mat(varargin{1}(7));  

filename=handles.listafiles(1).name;

%set window values using last used values (Densparam variable)
if length(dir('auxdens.mat'))>0  
    load('auxdens.mat','-mat');
    handles=setdefaults(Densparam, handles);
end

nano=get(handles.nanoradiobutton,'Value');
if nano==1
    set(handles.epsilon,'Enable','on');
    set(handles.mindetecnano,'Enable','on');
    set(handles.autoepsilonradiobutton,'Enable','on');
else
    set(handles.epsilon,'Enable','off');
    set(handles.mindetecnano,'Enable','off');
    set(handles.autoepsilonradiobutton,'Enable','off');
end
voro=get(handles.vororadiobutton,'Value');
if voro==1
    set(handles.polsize,'Enable','on');
    set(handles.mindenspx2,'Enable','on');
    set(handles.intenthresh,'Enable','off');
    set(handles.mindenspx,'Enable','off');
else
    set(handles.polsize,'Enable','off');
    set(handles.mindenspx2,'Enable','off');
    set(handles.intenthresh,'Enable','on');
    set(handles.mindenspx,'Enable','on');
end

autoep=get(handles.autoepsilonradiobutton,'Value');
if autoep==1
    set(handles.epsilon,'String',' ');
end

set(handles.filename,'String',[filename,' (1 of ',num2str(size(handles.listafiles,2)),')']);
        set(handles.textresults,'String',' #1'); 

guidata(hObject, handles);

%--------------------------------------------------------------------------
function varargout = DetectClustersBatch_OutputFcn(hObject, eventdata, handles)

varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function startpushbutton_Callback(hObject, eventdata, handles)

disp(' ')
disp('Batch analysis')

pxdilate=str2num(get(handles.pxdilate,'String'));
pxerode=str2num(get(handles.pxerode,'String'));
limitdens=str2num(get(handles.mindens,'String'));
mindenspx=str2num(get(handles.mindenspx,'String'));
mindenspx2=str2num(get(handles.mindenspx2,'String'));
minnrodetect=str2num(get(handles.mindetect,'String'));
mindiamcluster=str2num(get(handles.minsize,'String'));
maxdiamcluster=str2num(get(handles.maxsize,'String'));
inthresh=str2num(get(handles.intenthresh,'String'));
voro=get(handles.vororadiobutton,'Value');
if voro==0
    set(handles.polsize,'Enable','off')
    mindenspx=str2num(get(handles.mindenspx,'String'));
    vorosize=0;
else
    set(handles.polsize,'Enable','on')
    mindenspx=str2num(get(handles.mindenspx2,'String'));
    vorosize=str2num(get(handles.polsize,'String'));
end

donano=get(handles.nanoradiobutton,'Value');
if donano==0
    minpointsnano=0;
    set(handles.mindetecnano,'Enable','off')
else
    set(handles.polsize,'Enable','on')
    disp('With nanodomain analysis')
    minpointsnano=str2num(get(handles.mindetecnano,'String'));
end
autoep=get(handles.autoepsilonradiobutton,'Value');
if autoep==1
    epsilon=0;
else
    epsilon=str2num(get(handles.epsilon,'String'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Loop files
%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(handles.filename,'userdata',0); %starts from file 1
names={};
count=1;
currentdir=cd;
next=0;

 logical=1;
 while logical

     % initialize results files
     allresultsclu=[];
     allresultsnano=[];
     allresdistance=[];

     %load data
     cd(currentdir)
     
     next=next+1;
     if next>size(handles.listafiles,2)
        logical=0;
        break
     end
     
     filename=handles.listafiles(next).name
     [namefile,~]=strtok(filename,'.'); %sin extension

     set(handles.filename,'String',[filename,' (',num2str(next),' of ',num2str(size(handles.listafiles,2)),')']); 
     names{count}=filename;
     count=count+1;
   
     if handles.mu==0
        handles=loadselectPALMdata(filename,handles.szpx,handles);    
     else
        handles=loadselectPALMdata(filename,1,handles); %data already in microns
     end
     
     X_mu=handles.x-min(handles.x);
     Y_mu=handles.y-min(handles.y);
     alpha = handles.alpha(:);
     alpha(1:10);   
     xdim=ceil(max(X_mu)/handles.dx);
     ydim=ceil(max(Y_mu)/handles.dx);
     
     % remove by intensities
     vect2remove=selectremovealpha3(handles,X_mu,Y_mu,alpha);
     if ~isempty(vect2remove)
         [X_mu,Y_mu,alpha,fr]=removealpha(vect2remove,X_mu,Y_mu,alpha,handles.fr);
     end
     
     % Create rendered image
     gausssmooth = 1;
     [I,xxi,yyi] = PALM_rendering3(X_mu,Y_mu,alpha,handles.dx*2,handles.dx,0,xdim, ydim, 1);
     BW=zeros(size(I,1),size(I,2)); % for mask
    
     %ROIs
     if exist([namefile,'.rgn'])==2
        handles.data=importdata([namefile,'.rgn']);
        disp(['Loading regions from ',namefile,'.rgn'])
     else
         disp('ROI file not found')
         return
     end % of exist

     %read data
     data=handles.data;
     if isstruct(data)
        roifile=data.coord;
        handles.dx=data.dx;
        if isfield(data,'dist')
            handles.controldistanceroi=1;
        else
            handles.controldistanceroi=0;
        end
        handles.totalroi=size(roifile,2);
     else
        disp('Wrong ROI file')
     end %if sstruct
    
     % all ROIs
     for roij=1:size(roifile,2)
        
        set(handles.textresults,'String',[' #',num2str(roij)]); 
        disp(' ')
        disp(['ROI #',num2str(roij)])
        
        filenameIrend=[namefile,'-',num2str(roij),'-Irend.mat'];
        load(filenameIrend,'-mat') ;

        
        roiselectedx=[];
        roiselectedy=[];
        frselected=[];
        alphaselected =[];
        resultsclu=[];
        resultsnano=[];

        
        % points ROI
        roiselectedx=data.xselec{roij}; %index of selected points
        roiselectedy=data.yselec{roij}; %index of selected points
        handles.datax=roiselectedx/handles.dx;
        handles.datay=roiselectedy/handles.dx;
        
        if isempty(roiselectedx)==1
           disp('No data')
           return
        else
           
           if voro==0
               valintthresh=inthresh/100*max(max(Irend));
               BW=zeros(size(Irend,1),size(Irend,2));
               
               for i=1:size(Irend,1)
                   indexj=find(Irend(i,:)>valintthresh);
                   if isempty(indexj)==0
                       BW(i,indexj)=1;
                   end
               end
               
               if pxdilate>0
                se = strel('disk',pxdilate,0);
                I2 = imdilate(BW,se);
                I2 = imfill(I2,'holes');
               else
                I2=BW;
               end
            
               BW=I2;
               if pxerode>0
                se = strel('disk',pxerode,0);
                I2 = imerode(BW,se);
               end
               
           else % voronoi
               handles.vorosize= str2num(get(handles.polsize,'String'));
               fondo=zeros(size(Irend,1),size(Irend,2));
               [~, I2]=MaskVoronoi2(fondo,roiselectedx/handles.dx,roiselectedy/handles.dx,vorosize); %ATT disable figures!!!
           end
       end
       
       % Cluster analysis
       % keep clusters with more than mindens density of points
       [labeledall,~] = bwlabel(I2,4);
       
       maxdiamcluster=(((maxdiamcluster/2)^2)*pi)/(handles.dx*1000)^2; %!!!!!!!
       mindiamcluster=(((mindiamcluster/2)^2)*pi)/(handles.dx*1000)^2; %!!!!!!!
       
       [output,epsilon,~] =selectclusters(labeledall,handles.datax ,handles.datay,handles.fr,handles.alpha,limitdens,mindenspx,minnrodetect,mindiamcluster,maxdiamcluster,minpointsnano,epsilon);

       if autoep==1
           set(handles.epsilon,'String', num2str(epsilon)); % auto epsilon
       end

       rendxmask=output.rendxmask;
       countclusters=output.countclusters;
       disp([num2str(countclusters),' clusters found'])
       disp(' ')
       
       if isempty(output)==1
           disp('Nothing to save!')
           return
       end

       if isempty(rendxmask)==0
           rendymask=output.rendymask;
           
           sigmaloc=handles.dx*2;
           frmask=output.frmask;
           alphamask=output.alphamask;
           resultsclu=output.results ;  % nro cluster - # detect - #detect roi - area cluster - density
           listeboundary=output.listeboundary;
           countclusters=output.countclusters;
           auxoutx=output.auxoutx;
           auxouty=output.auxouty;
           rendxmasknano=output.rendxmasknano;
           rendymasknano=output.rendymasknano;
           resultsnano=output.resultsnano;
           countclustersnano=output.countclustersnano;
           if handles.controldistanceroi==1
               distance=data.dist{nroroi};
               allresdistance=[allresdistance; roij countclusters distance countclusters/distance]; % ATT units distance??
               %output.resdistance=resdistance;
           else
               %output.resdistance=[];
           end
 
           %convertion area from pixels to µm2
           if isempty(resultsclu)==0
               % 1 nro cluster - 2 # detect - 3 #detect roi - 4 area cluster - 5 density
               resultsclu(:,4)=resultsclu(:,4).*handles.dx^2; %µm2
               resultsclu(:,5)=resultsclu(:,2)./resultsclu(:,4);
           end
           if isempty(resultsnano)==0
               % nro cluster - nro nanocluster- #detect - #detect cluster - area nanocluster - density
               resultsnano(:,5)=resultsnano(:,4).*handles.dx^2; %µm2
               resultsnano(:,6)=resultsnano(:,2)./resultsnano(:,4);
           end
           
           % save data detections
           xdimmask=0; %size(imagemask,1);
           ydimmask=0; %size(imagemask,2);
           imagemask=[];
           dx=handles.dx;
           
           save([namefile,'-roidata_',num2str(roij),'.mat'],'auxoutx','auxouty','imagemask','limitdens','mindenspx','xdimmask', 'ydimmask',...
           'rendxmask','rendymask','rendxmasknano','rendymasknano','alphamask','frmask','minnrodetect','mindiamcluster',...
           'maxdiamcluster','dx','sigmaloc','voro','vorosize','-mat');
     
            if isempty(resultsclu)==0
                allresultsclu=[allresultsclu; roij*ones(size(resultsclu,1),1) resultsclu];
            end
            if isempty(resultsnano)==0
                allresultsnano=[allresultsnano; roij*ones(size(resultsnano,1),1)  resultsnano];
            end
           % if isempty(resdistance)==0
           %     allresdistance=[allresdistance; roij*ones(size(resdistance,1),1)  resdistance];
           % end
       
       end %empty data
     
    end %loop ROIs
   
    %save
    savename=[namefile,'-allclust.txt'];
    savenamenano=[namefile,'-allclustnano.txt'];
    savenamedist=[namefile,'-clustperdist.txt'];
    
    if isempty(allresultsclu)==0
        save(savename,'allresultsclu','-ascii')
    end
    if isempty(allresultsnano)==0
        save(savenamenano,'allresultsnano','-ascii')
    end
    if isempty(allresdistance)==0 % distance data
        save(savenamedist,'allresdistance','-ascii')
    end
    
    
 end %loop files


%go back to original values
mindiamcluster=str2num(get(handles.minsize,'String'));
maxdiamcluster=str2num(get(handles.maxsize,'String'));
   
Densparam.inthresh=inthresh;
Densparam.limitdens=limitdens;
Densparam.pxdilate=pxdilate;
Densparam.pxerode=pxerode;
Densparam.minnrodetect=minnrodetect;
Densparam.mindiamcluster=mindiamcluster;
Densparam.maxdiamcluster=maxdiamcluster;
Densparam.vorosize=vorosize;
Densparam.voro=voro;
if voro==0
    Densparam.mindenspx=mindenspx; 
else
    Densparam.mindenspx=mindenspx2; 
end
Densparam.epsilon=epsilon;
Densparam.autoepsilon=autoep;
Densparam.minpointsnano=minpointsnano;
Densparam.nano=donano;

save('auxdens.mat','Densparam','-mat') ;

%save report
reportclustering(Densparam,names);

disp('Done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nanoradiobutton_Callback(hObject, eventdata, handles)
val=get(hObject,'Value');
if val==1
    set(handles.epsilon,'Enable','on');
    set(handles.mindetecnano,'Enable','on');
    set(handles.autoepsilonradiobutton,'Enable','on');
else
    set(handles.epsilon,'Enable','off');
    set(handles.mindetecnano,'Enable','off');
    set(handles.autoepsilonradiobutton,'Enable','off');
end

function vororadiobutton_Callback(hObject, eventdata, handles)
val=get(hObject,'Value');
if val==1
    set(handles.polsize,'Enable','on');
    set(handles.mindenspx2,'Enable','on');
    set(handles.intenthresh,'Enable','off');
    set(handles.mindenspx,'Enable','off');
else
    set(handles.polsize,'Enable','off');
    set(handles.mindenspx2,'Enable','off');
    set(handles.intenthresh,'Enable','on');
    set(handles.mindenspx,'Enable','on');
end

function autoepsilonradiobutton_Callback(hObject, eventdata, handles)
val=get(hObject,'Value');
if val==1
    set(handles.epsilon,'Enable','off');
else
    set(handles.epsilon,'Enable','on');
end

function polsize_Callback(hObject, eventdata, handles)
function polsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mindens_Callback(hObject, eventdata, handles)
function mindens_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function mindenspx_Callback(hObject, eventdata, handles)
function mindenspx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function mindenspx2_Callback(hObject, eventdata, handles)
function mindenspx2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function minsize_Callback(hObject, eventdata, handles)
function minsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mindetect_Callback(hObject, eventdata, handles)
function mindetect_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxsize_Callback(hObject, eventdata, handles)
function maxsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function intenthresh_Callback(hObject, eventdata, handles)
function intenthresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pxdilate_Callback(hObject, eventdata, handles)
function pxdilate_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pxerode_Callback(hObject, eventdata, handles)
function pxerode_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mindetecnano_Callback(hObject, eventdata, handles)
function mindetecnano_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function epsilon_Callback(hObject, eventdata, handles)
function epsilon_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
