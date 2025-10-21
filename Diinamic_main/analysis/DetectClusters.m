function varargout = DetectClusters(varargin)
%function varargout = DetectClusters(varargin)
% GUI for setting thresholds and select clusters of detections
%
% called by Cluster.m
%
% Marianne Renner oct22
% Last Modified by GUIDE v2.5 29-Sep-2022 11:24:57
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DetectClusters_OpeningFcn, ...
                   'gui_OutputFcn',  @DetectClusters_OutputFcn, ...
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

function DetectClusters_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

% to fit screen definition/changes in size
h3 = findobj('Type','figure');
txtHand = findall(h3, '-property', 'FontUnits');
%set(txtHand, 'FontUnits', 'normalized');
set(txtHand, 'FontUnits', 'centimeters');

disp(' ')
disp('Reading data, please wait...')
disp(' ')
disp('Present ROI: #1') 

handles.I=cell2mat(varargin{1}(1));
handles.fr=cell2mat(varargin{1}(3));
handles.alpha=cell2mat(varargin{1}(4));
handles.file=cell2mat(varargin{1}(5));   
handles.PALMpx=cell2mat(varargin{1}(6));  
handles.szpx=cell2mat(varargin{1}(7));  
Densparam=cell2mat(varargin{1}(8));  
typeanalysis=cell2mat(varargin{1}(9));  

if typeanalysis==2
    set(handles.quitpushbutton,'String','Next file')
end
set(handles.quitpushbutton,'Enable','off')

[namefile,~]=strtok(handles.file,'.');
set(handles.filename,'String',handles.file);

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

%if ROIs, plot ROIs
if exist([namefile,'.rgn'])==2
    handles.data=importdata([namefile,'.rgn']);
    disp(['Loading regions from ',namefile,'.rgn'])
    data=handles.data;
    if isstruct(data)
        %coordroi=data.coord;
        handles.totalroi=size(data.coord,2);
    else
        disp('Error reading ROI data')
        return
    end
    clear data 
else
    % roifile= whole image
    % not implemented    
end % of exist    

%roi number=1
nroroi=1;
set(handles.nextroipushbutton,'Userdata',nroroi);
set(handles.textresults,'String',['#',num2str(nroroi),' of ',num2str(handles.totalroi)]);
if handles.totalroi==1
    set(handles.nextroipushbutton,'Enable','off');
    set(handles.previousroipushbutton,'Enable','off');
end

allresult.clu=[];
allresult.nano=[];
allresult.dist=[];
set(handles.acceptpushbutton,'Userdata',allresult);

guidata(hObject, handles);

%--------------------------------------------------------------------------
function varargout = DetectClusters_OutputFcn(hObject, eventdata, handles)

varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trypushbutton_Callback(hObject, eventdata, handles)
% launches analysis

disp(' ')
disp('Clustering analysis')

[namefile,rem]=strtok(handles.file,'.');
nroroi=get(handles.nextroipushbutton,'Userdata');

pxdilate=str2num(get(handles.pxdilate,'String'));
pxerode=str2num(get(handles.pxerode,'String'));
mindens=str2num(get(handles.mindens,'String'));
mindetect=str2num(get(handles.mindetect,'String'));
minsize=str2num(get(handles.minsize,'String'));
maxsize=str2num(get(handles.maxsize,'String'));
autoep=get(handles.autoepsilonradiobutton,'Value');
if autoep==1
    epsilon=0;
else
    epsilon=str2num(get(handles.epsilon,'String'));
end
inthresh=str2num(get(handles.intenthresh,'String'));

voro=get(handles.vororadiobutton,'Value');
if voro==0
    set(handles.polsize,'Enable','off')
    mindenspx=str2num(get(handles.mindenspx,'String'));

else
    set(handles.polsize,'Enable','on')
    vorosize=str2num(get(handles.polsize,'String'));
    mindenspx=str2num(get(handles.mindenspx2,'String'));
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

%read data
disp(' ');
disp('Analysing clustering by density');
disp(' ')
disp(['Mask in ROI # ', num2str(nroroi)])  
filenameIrend=[namefile,'-',num2str(nroroi),'-Irend.mat'];
load(filenameIrend,'-mat') ;
handles.dim=[size(Irend,1) size(Irend,2)];

data=handles.data;
if isstruct(data)
    handles.dx=data.dx;
    if isfield(data,'dist')
        handles.controldistanceroi=1;
    else
        handles.controldistanceroi=0;
    end
else
    disp('Wrong ROI file')
end %if sstruct

roiselectedx=[];
roiselectedy=[];
frselected=[];
alphaselected =[];
 
% points ROI
roiselectedx=data.xselec{nroroi}; %index of selected points
roiselectedy=data.yselec{nroroi}; %index of selected points
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
        [~, I2]=MaskVoronoi2(fondo,roiselectedx/handles.dx,roiselectedy/handles.dx,vorosize); %add nanoclusters!!!!!
    end
end

figurerend= figure('Name','Rendered image of ROI','Toolbar','figure');
imshow(Irend,'InitialMagnification','fit')
title([namefile,', ROI #',num2str(nroroi)],'Interpreter','none')

       
[labeledall,~] = bwlabel(I2,4);

maxsize=(((maxsize/2)^2)*pi)/(handles.PALMpx*1000)^2; 
minsize=(((minsize/2)^2)*pi)/(handles.PALMpx*1000)^2; 

% keep clusters with more than mindens density of points
[output,epsilon,labeled2] =selectclusters(labeledall,handles.datax ,handles.datay,handles.fr,handles.alpha,mindens,mindenspx,mindetect,minsize,maxsize,minpointsnano,epsilon);

if autoep==1
    set(handles.epsilon,'String', num2str(epsilon)); % auto epsilon
end
countclusters=output.countclusters;
disp([num2str(countclusters),' clusters found'])
disp(' ')

%figurerend= figure('Name','Segmented rendered image','Toolbar','figure');
%imshow(labeled2,'InitialMagnification','fit')
%title([namefile,', ROI #',num2str(nroroi)],'Interpreter','none')

rendxmask=output.rendxmask;
rendymask=output.rendymask;
countclusters=output.countclusters;
rendxmasknano=output.rendxmasknano;
rendymasknano=output.rendymasknano;
countclustersnano=output.countclustersnano;
if handles.controldistanceroi==1
    distance=data.dist{nroroi};
    resdistance=[nroroi countclusters distance countclusters/distance];
    output.resdistance=resdistance;
else
    output.resdistance=[];
end

imagemask=get(handles.textresults,'userdata');

if isempty(rendxmask)==0
    
  % visualization
  imageclu=zeros(size(I2,1),size(I2,2));
  fond=ones(size(I2,1),size(I2,2));
  
  for z=1:max(rendxmask(:,3)) %all selected clusters
    index=find(rendxmask(:,3)==z);
    if isempty(index)==0
        for t=1:size(index,1)
            posx=floor(rendxmask(index(t),1));
            if posx==0; posx=1; end
            posy=floor(rendymask(index(t),1));
            if posy==0; posy=1; end
            imageclu(posy,posx)=1; %new image with pixels that have enough dens
        end
    end
  end
  
  se = strel('ball',5,5);
  dilatedI = imdilate(imageclu,se);
  imageclu=imerode(dilatedI,se);
  for z=1:size(imageclu,1)
    index=find(imageclu(z,:)>0);
    if isempty(index)==0
        imageclu(z,index)=1;
    end
  end
  imageclu = imfill(imageclu,'holes');
  
  % figure clusters
  
  fondo=ones(ceil(max(handles.datay)),ceil(max(handles.datax)));

  figuremask=figure;
  hold on
  imshow(fondo,'InitialMagnification','fit')
  title([namefile,', ROI #',num2str(nroroi)],'Interpreter','none')
  axis off
  
  plot(handles.datax,handles.datay,'.','MarkerSize',6,'Color',[0.5 0.55 0.55]);
  if isempty(rendxmask)==0
        
        %possibility to plot with different colors!!
        colorcode=(jet(max(rendxmask(:,3))+1));
        for ll=1:size(colorcode,1)
            order=ceil(rand*size(colorcode,1));
            newcolorcode(ll,:)=colorcode(order,:);
        end
        for jj=1:max(rendxmask(:,3))
            index=find(rendxmask(:,3)==jj);
            if isempty(index)==0
               plot(rendxmask(index,1),rendymask(index,1),'.','MarkerSize',6,'Color',newcolorcode(jj,:)); % all clustered points or clusterd points out loc
               hold on
               polyin=boundary(rendxmask(index,1),rendymask(index,1));
               plot(rendxmask(index(polyin),1),rendymask(index(polyin),1),'LineWidth',0.75,'Color',newcolorcode(jj,:)); % all clustered points or clusterd points out loc
            end
        end
        set(handles.textresults,'string',[num2str(nroroi),': ',num2str(countclusters),' clusters']);
  end

  hold off
  
  if isempty(rendxmasknano)==0   % there are nanodomains
      
      % figure clusters
      figurenano=figure;
      hold on
      imshow(fond,'InitialMagnification','fit')
      title(['Nanodomains in ROI #',num2str(nroroi)],'Interpreter','none')
      axis off
      
      plot(handles.datax,handles.datay,'.','MarkerSize',5,'Color',[0.5 0.55 0.55]);
      colorcode=jet(max(rendxmasknano(:,3))+1);
      newcolorcode=colorcode;
      for ll=1:size(colorcode)
          order=ceil(rand*size(colorcode,1));
          newcolorcode(ll,:)=colorcode(order,:);
      end
      
      for jj=1:max(rendxmask(:,3))
          index=find(rendxmask(:,3)==jj);
          if isempty(index)==0
              polyin=boundary(rendxmask(index,1),rendymask(index,1));
              plot(rendxmask(index(polyin),1),rendymask(index(polyin),1),'Color',[0.2 0.2 0.2]); % all clustered points or clusterd points out loc
          end
      end
      
      for jj=1:max(rendxmasknano(:,3))
            index=find(rendxmasknano(:,3)==jj);
            if isempty(index)==0
               plot(rendxmasknano(index,1),rendymasknano(index,1),'.','MarkerSize',5,'Color',newcolorcode(jj,:)); % all clustered points or clusterd points out loc
               hold on
               polyinnano=boundary(rendxmasknano(index,1),rendymasknano(index,1));
               plot(rendxmasknano(index(polyinnano),1),rendymasknano(index(polyinnano),1),'Color',newcolorcode(jj,:)); % all clustered points or clusterd points out loc
            end
      end
      set(handles.textresults,'string',[num2str(nroroi),': ',num2str(countclusters-1),' clusters and ',num2str(countclustersnano-1),' nanodomains']);
  end
  set(handles.textresults,'userdata',I2);
  set(handles.trypushbutton,'Userdata',output)
else
    disp('No clusters left')
end

disp('Done')

set(handles.acceptpushbutton,'Enable','on')
set(handles.quitpushbutton,'Enable','on')

clear I mask rendxmask rendymask frmask alphamask results listeboundary countclusters maximage

 guidata(gcbo, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nextroipushbutton_Callback(hObject, eventdata, handles)

nroroi=get(handles.nextroipushbutton,'Userdata');
nroroi=nroroi+1;

%check if new ROI exists
if nroroi>handles.totalroi
    set(handles.nextroipushbutton,'Enable','off');
else
    set(handles.acceptpushbutton,'Enable','off');
    set(handles.nextroipushbutton,'Userdata',nroroi);
    set(handles.textresults,'String',['#',num2str(nroroi),' of ',num2str(handles.totalroi)]);
end
if nroroi>1
    set(handles.previousroipushbutton,'Enable','on');
end
if nroroi==handles.totalroi
    set(handles.nextroipushbutton,'Enable','off');
end

disp(['Present ROI: #',num2str(nroroi)]) 

 guidata(gcbo, handles);
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function previousroipushbutton_Callback(hObject, eventdata, handles)

nroroi=get(handles.nextroipushbutton,'Userdata');
nroroi=nroroi-1;

%check if new ROI exists
if nroroi<1
    set(handles.previousroipushbutton,'Enable','off');
else
    set(handles.acceptpushbutton,'Enable','off');
    set(handles.nextroipushbutton,'Userdata',nroroi);
    set(handles.textresults,'String',['#',num2str(nroroi),' of ',num2str(handles.totalroi)]);
end
if nroroi==1
    set(handles.previousroipushbutton,'Enable','off');
end
if nroroi<handles.totalroi
    set(handles.nextroipushbutton,'Enable','on');
end

disp(['Present ROI: #',num2str(nroroi)]) 

 guidata(gcbo, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function acceptpushbutton_Callback(hObject, eventdata, handles)

imagemask=get(handles.textresults,'userdata');
[namefile,~]=strtok(handles.file,'.');
nroroi=get(handles.nextroipushbutton,'Userdata');
output=get(handles.trypushbutton,'Userdata');

allresult=get(handles.acceptpushbutton,'Userdata');
allresultsclu=allresult.clu;
allresultsnano=allresult.nano;
allresultsdist=allresult.dist;

pxdilate=str2num(get(handles.pxdilate,'String'));
pxerode=str2num(get(handles.pxerode,'String'));
limitdens=str2num(get(handles.mindens,'String'));
mindenspx=str2num(get(handles.mindenspx,'String'));
mindenspx2=str2num(get(handles.mindenspx2,'String'));
minnrodetect=str2num(get(handles.mindetect,'String'));
mindiamcluster=str2num(get(handles.minsize,'String'));
maxdiamcluster=str2num(get(handles.maxsize,'String'));
inthresh=str2num(get(handles.intenthresh,'String'));
autoep=get(handles.autoepsilonradiobutton,'Value');
if autoep==0
    epsilon=0;
else
    epsilon=str2num(get(handles.epsilon,'String'));
end
inthresh=str2num(get(handles.intenthresh,'String'));

voro=get(handles.vororadiobutton,'Value');
donano=get(handles.nanoradiobutton,'Value');
minpointsnano=str2num(get(handles.mindetecnano,'String'));

dx=handles.dx;
sigmaloc=handles.dx*2;

if isempty(output)==1
    disp('Nothing to save!')
    return
end

rendxmask=output.rendxmask;
rendymask=output.rendymask;
frmask=output.frmask;
alphamask=output.alphamask;
resultsclu=output.results;   
resultsnano=output.resultsnano;
resdistance=output.resdistance;
auxoutx=output.auxoutx;
auxouty=output.auxouty;
rendxmasknano=output.rendxmasknano;
rendymasknano=output.rendymasknano;

if voro==1 % voronoi
    vorosize=str2num(get(handles.polsize,'String'));
else
    vorosize=[];
end

%convertion area from pixels to �m2
if isempty(resultsclu)==0
    % 1 nro cluster - 2 # detect - 3 #detect roi - 4 area cluster - 5 density
    resultsclu(:,4)=resultsclu(:,4).*handles.dx^2; %µm2
    resultsclu(:,5)=resultsclu(:,2)./resultsclu(:,4);
end
if isempty(resultsnano)==0  
    % nro cluster - nro nanocluster- #detect - #detect cluster - area nanocluster - density
    resultsnano(:,5)=resultsnano(:,5).*handles.dx^2; %µm2
    resultsnano(:,6)=resultsnano(:,3)./resultsnano(:,5);
end

%fond=ones(handles.dim(1),handles.dim(2));
fondo=ones(ceil(max(handles.datay)),ceil(max(handles.datax)));

finalfig=figure;
hold on
imshow(fondo,'InitialMagnification','fit')
title([namefile,', ROI #',num2str(nroroi)],'Interpreter','none')
axis off
    
plot(handles.datax,handles.datay,'.','MarkerSize',5,'Color',[0.5 0.55 0.55]);

if isempty(rendxmask)==0
    %possibility to plot with different colors!!
    colorcode=(jet(max(rendxmask(:,3))+1));
    for ll=1:size(colorcode)
        order=ceil(rand*size(colorcode,1));
        newcolorcode(ll,:)=colorcode(order,:);
    end
    
    for jj=1:max(rendxmask(:,3))
            index=find(rendxmask(:,3)==jj);
            if isempty(index)==0
               plot(rendxmask(index,1),rendymask(index,1),'.','MarkerSize',5,'Color',newcolorcode(jj,:)); % all clustered points or clusterd points out loc
               hold on
               polyin=boundary(rendxmask(index,1),rendymask(index,1));
               plot(rendxmask(index(polyin),1),rendymask(index(polyin),1),'Color',newcolorcode(jj,:)); 
            end
    end
    
    % att plot nanodomains?
     
end
hold off

presentimage2=getframe(gca);  % gets the mask image 
Savename=[namefile,'-roi',num2str(nroroi),'-selection.tif'];
[imagemask2,Map] = frame2im(presentimage2);  
imwrite(imagemask2,Savename,'tif','Compression','none');   

% nanodomains lacking!!!!

% save data detections
xdimmask=size(imagemask,1);
ydimmask=size(imagemask,2);

save([namefile,'-roidata_',num2str(nroroi),'.mat'],'auxoutx','auxouty','imagemask','limitdens','mindenspx','xdimmask', 'ydimmask',...
    'rendxmask','rendymask','rendxmasknano','rendymasknano','alphamask','frmask','minnrodetect','mindiamcluster',...
    'maxdiamcluster','dx','sigmaloc','voro','vorosize','-mat');
close(finalfig)

%save results
if isempty(resultsclu)==0
    allresultsclu=[allresultsclu; nroroi*ones(size(resultsclu,1),1) resultsclu];
end
if isempty(resultsnano)==0
    allresultsnano=[allresultsnano; nroroi*ones(size(resultsnano,1),1) resultsnano];
end
if isempty(resdistance)==0
    allresultsdist=[allresultsdist; resdistance];
end

allresult.clu=allresultsclu;
allresult.nano=allresultsnano;
allresult.dist=allresultsdist;

set(handles.acceptpushbutton,'Userdata',allresult);

% parameters
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

%detection results
abortsignal=0;
save('auxiliar.mat','imagemask','allresult','abortsignal','-mat') 

%close % DetectClusters window
if nroroi+1>handles.totalroi
    set(handles.nextroipushbutton,'Enable','off')
else
    set(handles.nextroipushbutton,'Enable','on')
end
disp('Results saved for the ROI')
disp(' ')

%clear  'imagemask' 'mindens'  'xdimmask'  'ydimmask'  'rendxmask'  'rendymask'  'alphamask'  'frmask' 'minnrodetect'  'minsize' 
%clear 'maxsize'  'dx'  'sigmaloc' 'presentimage2' 'presentimage' 'listeclucercle'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function quitpushbutton_Callback(hObject, eventdata, handles)

%just close the window
clear handles
close % DetectClusters window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function abortpushbutton_Callback(hObject, eventdata, handles)

% interrupt analysis
abortsignal=1;
save('auxiliar.mat','abortsignal','-mat') 

clear handles
close % DetectClusters window

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
