function varargout = OptimDetectClusters(varargin)
% Last Modified by GUIDE v2.5 14-Apr-2023 21:50:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OptimDetectClusters_OpeningFcn, ...
                   'gui_OutputFcn',  @OptimDetectClusters_OutputFcn, ...
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OptimDetectClusters_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

disp(' ');
disp('Optimization of parameters');
disp(' ')

% to fit screen definition/changes in size
h3 = findobj('Type','figure');
txtHand = findall(h3, '-property', 'FontUnits');
%set(txtHand, 'FontUnits', 'normalized');
set(txtHand, 'FontUnits', 'centimeters');

%handles.listafiles=cell2mat(varargin{1}(1));
handles.dx=cell2mat(varargin{1}(2));  
handles.szpx=cell2mat(varargin{1}(3));  
handles.mu=cell2mat(varargin{1}(4));  
handles.alpha1=cell2mat(varargin{1}(5));  
handles.alpha2=cell2mat(varargin{1}(6));  
Densparam=cell2mat(varargin{1}(7));  

%set window values using last used values (Densparam variable)
if length(dir('auxdens.mat'))>0  
    load('auxdens.mat','-mat');
    handles=setdefaults(Densparam, handles);
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

guidata(hObject, handles);

%-------------------------------------------------------------------------

function varargout = OptimDetectClusters_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function loadpushbutton_Callback(hObject, eventdata, handles)

%load data
[filename,path] = uigetfile('*.mat','Load detections file (.mat)');
if path>0
      handles.presentfolder=path;
else
    return
end
cd(path)

set(handles.filename,'String',filename); 
[namefile,~]=strtok(filename,'.'); %sin extension

if handles.mu==0
    handles=loadselectPALMdata(filename,handles.szpx,handles);    
else
    handles=loadselectPALMdata(filename,1,handles); %data already in microns
end

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
    handles.totalroi=size(roifile,2);
else
    disp('Wrong ROI file')
    return
end %if sstruct

nroroi=1;
set(handles.nextroipushbutton,'Userdata',nroroi);
set(handles.loadpushbutton,'Userdata',data);
set(handles.textresults,'String',['#',num2str(nroroi),' of ',num2str(handles.totalroi)]);
if handles.totalroi==1
    set(handles.nextroipushbutton,'Enable','off');
    set(handles.previousroipushbutton,'Enable','off');
else
    set(handles.nextroipushbutton,'Enable','on');
end

filenameIrend=[namefile,'-1-Irend.mat'];
load(filenameIrend,'-mat') ;
handles.dim=[size(Irend,1) size(Irend,2)];
%handles.Irend=Irend;
set(handles.text33,'Userdata',Irend);
set(handles.visuintpushbutton,'Enable','on')
set(handles.visudenspushbutton,'Enable','on')
set(handles.visudensvoropushbutton,'Enable','on')
set(handles.startpushbutton,'Enable','on')

disp(['Mask in ROI # ', num2str(nroroi)])  
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nextroipushbutton_Callback(hObject, eventdata, handles)

filename=get(handles.filename,'String'); 
[namefile,~]=strtok(filename,'.'); %sin extension
nroroi=get(handles.nextroipushbutton,'Userdata');
nroroi=nroroi+1;

%check if new ROI exists
if nroroi>handles.totalroi
    set(handles.nextroipushbutton,'Enable','off');
else
    set(handles.nextroipushbutton,'Userdata',nroroi);
end
if nroroi>1
    set(handles.previousroipushbutton,'Enable','on');
end
if nroroi==handles.totalroi
    set(handles.nextroipushbutton,'Enable','off');
end

%ROIs
if exist([namefile,'.rgn'])==2
    handles.data=importdata([namefile,'.rgn']);
    disp(['Loading regions from ',namefile,'.rgn'])
else
    disp('ROI file not found')
    return
end % of exist

filenameIrend=[namefile,'-',num2str(nroroi),'-Irend.mat'];
load(filenameIrend,'-mat') ;
handles.dim=[size(Irend,1) size(Irend,2)];
set(handles.text33,'Userdata',Irend);

%read data
data=handles.data;
if isstruct(data)
    roifile=data.coord;
else
    disp('Wrong ROI file')
    return
end %if sstruct

set(handles.nextroipushbutton,'Userdata',nroroi);
set(handles.loadpushbutton,'Userdata',data);

disp(['Present ROI: #',num2str(nroroi)]) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function previousroipushbutton_Callback(hObject, eventdata, handles)

filename=get(handles.filename,'String'); 
[namefile,~]=strtok(filename,'.'); %sin extension
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

filenameIrend=[namefile,'-',num2str(nroroi),'-Irend.mat'];
load(filenameIrend,'-mat') ;
handles.dim=[size(Irend,1) size(Irend,2)];
set(handles.text33,'Userdata',Irend);

%read data
data=handles.data;
if isstruct(data)
    roifile=data.coord;
else
    disp('Wrong ROI file')
    return
end %if sstruct

set(handles.nextroipushbutton,'Userdata',nroroi);
set(handles.loadpushbutton,'Userdata',data);

disp(['Present ROI: #',num2str(nroroi)]) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function visuintpushbutton_Callback(hObject, eventdata, handles)

filename=get(handles.filename,'String'); 
[namefile,~]=strtok(filename,'.'); %sin extension
nroroi=get(handles.nextroipushbutton,'Userdata');

inthresh=str2num(get(handles.intenthresh,'String'));
pxdilate=str2num(get(handles.pxdilate,'String'));
pxerode=str2num(get(handles.pxerode,'String'));
mindenspx=str2num(get(handles.mindenspx,'String'));

Irend=get(handles.text33,'Userdata');

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

% plot
figurerend= figure('Name','Rendered image segmented by intensity','Toolbar','figure');
axis equal
imshow(I2,'InitialMagnification','fit')
title([namefile,', ROI #',num2str(nroroi)],'Interpreter','none')
%savefig([namefile,', ROI #',num2str(nroroi),'.fig'])


%-------------------------------------------------------------------------
function visudenspushbutton_Callback(hObject, eventdata, handles)

filename=get(handles.filename,'String'); 
[namefile,~]=strtok(filename,'.'); %sin extension
nroroi=get(handles.nextroipushbutton,'Userdata');
data=get(handles.loadpushbutton,'Userdata');

inthresh=str2num(get(handles.intenthresh,'String'));
pxdilate=str2num(get(handles.pxdilate,'String'));
pxerode=str2num(get(handles.pxerode,'String'));
mindenspx=str2num(get(handles.mindenspx,'String'));

% points ROI

roiselectedx=data.xselec{nroroi}; %index of selected points
roiselectedy=data.yselec{nroroi}; %index of selected points
datax=roiselectedx/handles.dx;
datay=roiselectedy/handles.dx;

Irend=get(handles.text33,'Userdata');

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
[labeled,~] = bwlabel(I2,4);

stats = regionprops(labeled,'PixelList','Area') ; 
nb_clust=size(stats,1);
listpoints=[];

disp('')
disp('Visualization of density segmentation. Processing... please wait')
disp(' ')

BW=zeros(size(Irend,1),size(Irend,2));

if mindenspx>0 % only intensity threshold
    for n=1:nb_clust % chaque cluster
    
       clust=stats(n).PixelList;
        
       for cc=1:size(clust,1)
           index2 = intersect(find(round(datax(:))==clust(cc,1)),find(round(datay(:))==clust(cc,2)));
           if isempty(index2)==0
               if size(index2,1)> mindenspx % minimum density per pixel
                   BW(clust(cc,2), clust(cc,1))=1; %new image with pixels that have enough dens
                   listpoints=[listpoints index2'];
               end
           end
       end
    end
else
    BW=labeled;
end

se = strel('disk',1,0);
BW = imdilate(BW,se);
BW = imfill(BW,'holes');
se = strel('disk',1,0);
BW = imerode(BW,se);
[labeled2,nb_clust] = bwlabel(BW,4);

disp('Done')

% plot
figurerend= figure('Name','Rendered image segmented by intensity and density','Toolbar','figure');
axis equal
imshow(labeled2,'InitialMagnification','fit')
title([namefile,', ROI #',num2str(nroroi)],'Interpreter','none')
%savefig([namefile,', ROI #',num2str(nroroi),'.fig'])

%----------------------------------------------------------------------
function visudensvoropushbutton_Callback(hObject, eventdata, handles)

filename=get(handles.filename,'String'); 
[namefile,~]=strtok(filename,'.'); %sin extension
nroroi=get(handles.nextroipushbutton,'Userdata');
data=get(handles.loadpushbutton,'Userdata');

inthresh=str2num(get(handles.intenthresh,'String'));
pxdilate=str2num(get(handles.pxdilate,'String'));
pxerode=str2num(get(handles.pxerode,'String'));
mindenspx=str2num(get(handles.mindenspx2,'String'));

% points ROI

roiselectedx=data.xselec{nroroi}; %index of selected points
roiselectedy=data.yselec{nroroi}; %index of selected points
datax=roiselectedx/handles.dx;
datay=roiselectedy/handles.dx;

Irend=get(handles.text33,'Userdata');

% Voronoi tesselation
vorosize= str2num(get(handles.polsize,'String'));
fondo=zeros(size(Irend,1),size(Irend,2));
[~, I2]=MaskVoronoi2(fondo,roiselectedx/handles.dx,roiselectedy/handles.dx,vorosize); %add nanoclusters!!!!!

[labeled,~] = bwlabel(I2,4);

stats = regionprops(labeled,'PixelList','Area') ; 
nb_clust=size(stats,1);
listpoints=[];

disp('')
disp('Visualization of density segmentation after voronoi tesselation. Processing... please wait')
disp(' ')

BW=zeros(size(Irend,1),size(Irend,2));

if mindenspx>0 % only intensity threshold
    for n=1:nb_clust % chaque cluster
    
       clust=stats(n).PixelList;
        
       for cc=1:size(clust,1)
           index2 = intersect(find(round(datax(:))==clust(cc,1)),find(round(datay(:))==clust(cc,2)));
           if isempty(index2)==0
               if size(index2,1)> mindenspx % minimum density per pixel
                   BW(clust(cc,2), clust(cc,1))=1; %new image with pixels that have enough dens
                   listpoints=[listpoints index2'];
               end
           end
       end
    end
else
    BW=labeled;
end

se = strel('disk',1,0);
BW = imdilate(BW,se);
BW = imfill(BW,'holes');
se = strel('disk',1,0);
BW = imerode(BW,se);
[labeled2,nb_clust] = bwlabel(BW,4);

disp('Done')

% plot
figurerend= figure('Name','Rendered image segmented by intensity and density','Toolbar','figure');
axis equal
imshow(labeled2,'InitialMagnification','fit')
title([namefile,', ROI #',num2str(nroroi)],'Interpreter','none')
%savefig([namefile,', ROI #',num2str(nroroi),'.fig'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function startpushbutton_Callback(hObject, eventdata, handles)

disp(' ')
disp('Batch optimization of parameters')
nroroi=get(handles.nextroipushbutton,'Userdata');
data=get(handles.loadpushbutton,'Userdata');
Irend=get(handles.text33,'Userdata');
plotok=get(handles.plotradiobutton,'Value');
mindenspx=str2num(get(handles.mindenspx,'String'));
mindenspx2=str2num(get(handles.mindenspx2,'String'));

pxdilate=str2num(get(handles.pxdilate,'String'));
pxerode=str2num(get(handles.pxerode,'String'));

limitdens=str2num(get(handles.mindens,'String'));
stopdens=str2num(get(handles.mindensstop,'String'));
intervdens=str2num(get(handles.mindensinterv,'String'));

minnrodetect=str2num(get(handles.mindetect,'String'));
stopdet=str2num(get(handles.mindetstop,'String'));
intervdet=str2num(get(handles.mindetinterv,'String'));

mindiamcluster=str2num(get(handles.minsize,'String'));
stopsize=str2num(get(handles.minsizestop,'String'));
intervsize=str2num(get(handles.minsizeinterv,'String'));

%error check
errorinterval=0;
if stopdens>0 && intervdens>stopdens; errorinterval=1; end
if stopdet>0 && intervdet>stopdet; errorinterval=1; end
if stopsize>0 && intervsize>stopsize; errorinterval=1; end
if errorinterval==1
    msgbox('Check interval values','Error','error')
    return
end

if stopdens>0 && intervdens>0
    Vmindens=limitdens:intervdens:stopdens; %
else
    Vmindens=limitdens;
end

if stopdet>0 && intervdet>0
    Vmindet=minnrodetect:intervdet:stopdet; %
else
    Vmindet=minnrodetect;
end

if stopsize>0 && intervsize>0
    Vminsize=mindiamcluster:intervsize:stopsize;
else
    Vminsize=mindiamcluster;
end

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

filename=get(handles.filename,'String'); 
names{1}=filename;
[namefile,~]=strtok(filename,'.'); %sin extension

results=[];

% points ROI
roiselectedx=data.xselec{nroroi}; %index of selected points
roiselectedy=data.yselec{nroroi}; %index of selected points
handles.datax=roiselectedx/handles.dx;
handles.datay=roiselectedy/handles.dx;
handles.datafr=data.fr{nroroi};
handles.dataalpha=data.alpha{nroroi};

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
[labeledall,~] = bwlabel(I2,4);

for ii=1:size(Vmindens,2)
    limitdens=Vmindens(ii);
    %set(handles.mindens,'String',num2str(limitdens))
        
    for mm=1:size(Vmindet,2)
        minnrodetect=Vmindet(mm);    
        %set(handles.mindetect,'String',num2str(minnrodetect))

        for jj=1:size(Vminsize,2)
            mindiamcluster=Vminsize(jj);
            
            disp(['Min density: ',num2str(limitdens),' - Min # detections: ',num2str(minnrodetect),' - Min size : ',num2str(mindiamcluster)])
            % keep clusters with more than mindens density of points
            [labeledall,~] = bwlabel(I2,4);
            maxdiamcluster=(((maxdiamcluster/2)^2)*pi)/(handles.dx*1000)^2; %!!!!!!!
            mindiamcluster=(((mindiamcluster/2)^2)*pi)/(handles.dx*1000)^2; %!!!!!!!
            minpointsnano=0;
            epsilon=1;
            
            [output,epsilon,~] =selectclusters(labeledall,handles.datax ,handles.datay,handles.datafr,handles.dataalpha,limitdens,mindenspx,minnrodetect,mindiamcluster,maxdiamcluster,minpointsnano,epsilon);
  
            rendxmask=output.rendxmask;
            countclusters=output.countclusters;
            disp([num2str(countclusters),' clusters found'])
            disp(' ')

            if isempty(rendxmask)==0
                   rendymask=output.rendymask;
                   
                   resultsclu=output.results ;  % nro cluster - # detect - #detect roi - area cluster - density
                   aux=resultsclu;
                   count=1;
                   for nro=1:max(resultsclu(:,1))
                       index=find(resultsclu(:,1)==nro);
                       if isempty(index)==0
                           aux(index,1)=count;
                           count=count+1;
                       end
                   end
                   resultsclu=aux;
                   
                   results=[results; limitdens minnrodetect mindiamcluster max(resultsclu(:,1)) ]; %results optim
            else
                results=[results; limitdens minnrodetect mindiamcluster 0 ];
               
            end %empty data
               
            if plotok   % att only for ROI all image
                   
                   fondo=ones(ceil(max(handles.datay)),ceil(max(handles.datax)));

                   finalfig=figure;
                   hold on
                   imshow(fondo,'InitialMagnification','fit')
                   title([namefile,', ROI #',num2str(nroroi)],'Interpreter','none')
                   axis off
                   axis equal
                   
                   plot(handles.datax,handles.datay,'.','MarkerSize',5,'Color',[0.5 0.55 0.55]);
           
                   if isempty(rendxmask)==0
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
                   end
                   title([namefile,'-mindens',num2str(limitdens),'-mindet',num2str(minnrodetect)],'Interpreter','none')
                   hold off

                   savefig([namefile,'-mindens',num2str(limitdens),'-mindet',num2str(minnrodetect),'.fig'])
                   %close(finalfig)
                   
            end % plotok
        end %loop minsize
    end % loop mindet
end % loop mindens

%save results optimization
save([namefile,'-optimDiinamic.txt'],'results','-ascii')
    
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
Densparam.autoepsilon=0;
Densparam.minpointsnano=minpointsnano;
Densparam.nano=0;

save('auxdens.mat','Densparam','-mat') ;

%save report
%reportclustering(Densparam,names);

disp('Done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

function mindensstop_Callback(hObject, eventdata, handles)
function mindensstop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mindensradiobutton_Callback(hObject, eventdata, handles)
function mindetradiobutton_Callback(hObject, eventdata, handles)
function minsizeradiobutton_Callback(hObject, eventdata, handles)
function plotradiobutton_Callback(hObject, eventdata, handles)

function mindetinterv_Callback(hObject, eventdata, handles)
function mindetinterv_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mindetstop_Callback(hObject, eventdata, handles)
function mindetstop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function minsizestop_Callback(hObject, eventdata, handles)
function minsizestop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function minsizeinterv_Callback(hObject, eventdata, handles)
function minsizeinterv_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mindensinterv_Callback(hObject, eventdata, handles)
function mindensinterv_CreateFcn(hObject, eventdata, handles)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
