function varargout = MaskRendered(varargin)
% function varargout = MaskRendered(varargin)
% 
% GUI for setting threshold and select clusters of detections manually
%
% called by Clustering_v2.m and Diinamic.m
%
% Marianne Renner oct19
% modif 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Last Modified by GUIDE v2.5 10-Mar-2016 12:39:47
 
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MaskRendered_OpeningFcn, ...
                   'gui_OutputFcn',  @MaskRendered_OutputFcn, ...
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
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MaskRendered_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
 
set(handles.output,'userdata',varargin{1}(1));  %roimovie

handles.I=cell2mat(varargin{1}(1));
handles.fr=cell2mat(varargin{1}(3));
handles.alpha=cell2mat(varargin{1}(4));
handles.filename=cell2mat(varargin{1}(5));   
handles.PALMpx=cell2mat(varargin{1}(6));  
handles.szpx=cell2mat(varargin{1}(7));  

[namefile,rem]=strtok(handles.filename,'.');

%if ROIs, plot ROIs
if exist([namefile,'.rgn'])==2
    data=importdata([namefile,'.rgn']);
    if isstruct(data)
        handles.roifile=data.coord;
        handles.indata=data.in; %indexes points
        handles.Irend=data.imagerend; %rendered of the ROI
        handles.dx=data.dx;
        if isfield(data,'dist')
            % distance
            handles.distance=data.dist;
        end
    else
        handles.roifile=data;
    end %if sstruct
    disp(['Loading regions from ',namefile,'.rgn'])
else
    % roifile= whole image
    % not implemented    
end % of exist    

%filename


%roi number=1



guidata(hObject, handles);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = MaskRendered_OutputFcn(hObject, eventdata, handles) 
 
varargout{1} = handles.output;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


set(handles.mindetect,'String',num2str(handles.minnrodetect));
set(handles.mindens,'String',num2str(handles.limitdens));
set(handles.minsize,'String',num2str(handles.mindiamcluster));
set(handles.maxsize,'String',num2str(handles.maxdiamcluster));
set(handles.pxdilate,'String',num2str(handles.dilate));
set(handles.pxerode,'String',num2str(handles.erode));
set(handles.intenthresh,'String',num2str(handles.percinthresh));

if handles.voronoi==1 % voronoi
    set(handles.polsize,'String',num2str(handles.vorosize));
    set(handles.text17,'Enable','on');
    set(handles.polsize,'Enable','on');
else
    set(handles.text17,'Enable','off');
    set(handles.polsize,'Enable','off');
end

set(handles.roinro,'string',['ROI # ',num2str(handles.nroroi)]);

handles.maxdiamcluster=(((handles.maxdiamcluster/2)^2)*pi)/(handles.PALMpx*1000)^2; %!!!!!!!
handles.mindiamcluster=(((handles.mindiamcluster/2)^2)*pi)/(handles.PALMpx*1000)^2; %!!!!!!!

%---------------------------------------------------------------------------
%intensity threshold
%--------------------------------------------------------------------------
I=handles.data;

if handles.voronoi==0
    set(handles.polsize,'Enable','off')

    valintthresh=handles.percinthresh/100*max(max(I));
    BW=zeros(size(I,1),size(I,2));
    
    for i=1:size(I,1)
        indexj=find(I(i,:)>valintthresh);
        if isempty(indexj)==0
            BW(i,indexj)=1;
        end
    end
    
    if handles.dilate>0
        se = strel('disk',handles.dilate,0);
        I2 = imdilate(BW,se);
        I2 = imfill(I2,'holes');
    else
        I2=BW;
    end
    
    BW=I2;
    if handles.erode>0
        se = strel('disk',handles.erode,0);
        I2 = imerode(BW,se);
    end
    
    figurerend= figure('Name','Segmented rendered image','Toolbar','figure');
    imshow(I2,'InitialMagnification','fit')
    title([handles.namefile,', ROI #',num2str(handles.nroroi)])

else % voronoi
    
    set(handles.polsize,'Enable','on')
    
    fondo=zeros(size(I,1),size(I,2));
    [~, I2]=MaskVoronoi2(fondo,handles.datax,handles.datay,handles.vorosize);  
    
end

% keep clusters with more than mindens density of points
[labeledall,nb_clust] = bwlabel(I2,4);

rendxmask=[];
rendymask=[];
frmask=[];
alphamask=[];
results=[];
listeboundary=[];
countclusters=1;
rendxmask2=[];
rendymask2=[];
frmask2=[];
alphamask2=[];
results2=[];
listeboundary2=[];
countclusters2=1;

loc=handles.synroi;

[rendxmask, rendymask,frmask,alphamask, results,listeboundary,countclusters,auxoutx,auxouty] =makemask5(labeledall, handles.datax,handles.datay,handles.fr,...
    handles.alpha,handles.limitdens, handles.minnrodetect,handles.mindiamcluster,handles.maxdiamcluster,handles.maxnumberdet,handles.synroi,handles.voronoi);

imageclu=zeros(size(I2,1),size(I2,2));
fond=ones(size(I2,1),size(I2,2));

if isempty(rendxmask)==0 %there are clusters

 for z=1:max(rendxmask(:,4)) %all selected clusters
    index=find(rendxmask(:,4)==z);
    if isempty(index)==0
        for t=1:size(index)
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

    figuremask=figure;
    hold on
    imshow(fond,'InitialMagnification','fit')

    title([handles.namefile,', ROI #',num2str(handles.nroroi)])
    axis off

    plot(handles.datax,handles.datay,'.','MarkerSize',5,'Color',[0.5 0.55 0.55]);
    if isempty(rendxmask)==0
        
        %possibility to plot with different colors!!
        colorcode=jet(max(rendxmask(:,4))+1);
        for jj=1:max(rendxmask(:,4))
            index=find(rendxmask(:,4)==jj);
            if isempty(index)==0
               plot(rendxmask(index,1),rendymask(index,1),'.','MarkerSize',5,'Color',colorcode(jj,:)); % all clustered points or clusterd points out loc
               hold on
               polyin=boundary(rendxmask(index,1),rendymask(index,1));
               plot(rendxmask(index(polyin),1),rendymask(index(polyin),1),'Color',colorcode(jj,:)); % all clustered points or clusterd points out loc
            end
        end
        
        %plot(rendxmask(:,1),rendymask(:,1),'.','MarkerSize',5,'Color','g'); % all clustered points or clusterd points out loc
    end

   % imcontour(imageclu,1) % NON: boundary

hold off

set(handles.textresults,'string',['# of clusters: ',num2str(countclusters-1)]);
set(handles.acceptpushbutton,'userdata',results);
set(handles.quitpushbutton,'userdata',figuremask);
set(handles.textresults,'userdata',I2);
set(handles.trypushbutton,'userdata',loc);
set(handles.text11,'userdata',rendxmask);
set(handles.text2,'userdata',rendymask);
set(handles.text9,'userdata',frmask);
set(handles.text12,'userdata',alphamask);
set(handles.text7,'userdata',listeboundary);
set(handles.text8,'userdata',countclusters);

set(handles.text15,'userdata',auxoutx); % detections out of clusters
set(handles.text17,'userdata',auxouty);

 
else % empty result
    disp('No clusters left after thresholds')
end % 

clear mask rendxmask rendymask frmask alphamask results listeboundary countclusters maximage
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trypushbutton_Callback(hObject, eventdata, handles)
 
pxdilate=str2num(get(handles.pxdilate,'String'));
pxerode=str2num(get(handles.pxerode,'String'));

roifile=handles.roifile;

mindens=str2num(get(handles.mindens,'String'));
mindetect=str2num(get(handles.mindetect,'String'));
minsize=str2num(get(handles.minsize,'String'));
maxsize=str2num(get(handles.maxsize,'String'));

inthresh=str2num(get(handles.intenthresh,'String'));

voro=str2num(get(handles.vororadiobutton,'Value'));
vorosize=str2num(get(handles.polsize,'String'));

donano=str2num(get(handles.nanoradiobutton,'Value'));
if donano==0
    minpointsnano=0;
else
    minpointsnano=str2num(get(handles.epsilon,'String'));
end

%ROI!!!
%get the number
%roij

%read data

disp(' ')
disp(['Mask in ROI # ', num2str(roij)])  

roiselectedx=[];
roiselectedy=[];
frselected=[];
alphaselected =[];
            
% crop image by ROI
BW=zeros(size(handles.I,1),size(handles.I,2));

xi=roifile{roij}(:,1)/dx; % coordinates of the ROI in PALM pixels size
yi=roifile{roij}(:,2)/dx; 
%rendered ROI
imageroirend=Irend{roij}; 

in= indata{roij};

minxi=floor(min(xi));
maxxi=ceil(max(xi));
minyi=floor(min(yi));
maxyi=ceil(max(yi));  
minyi=max(1,floor(min(yi)));
maxyi=min(size(handles.I,1),ceil(max(yi)));

if minxi<1; minxi=1; end % to avoid points out of the image size
if maxxi>size(handles.I,2); maxxi=size(handles.I,2);  end
if minyi<1; minyi=1; end
if maxyi>size(handles.I,1); maxyi=size(handles.I,1) ; end  
            
roi=roipolyold(BW,xi,yi);   % mask of the ROI (the same size than the rendered image)
imageroi=immultiply(roi,handles.I); % picks up pixels of the rendered image that corresponds to the ROI
imageroi=imageroi(minyi:maxyi,minxi:maxxi); % crops to the size of the ROI

%pick points ROI
roiselectedx=data.xselec{roij}; %index of selected points
roiselectedy=data.yselec{roij}; %index of selected points
            
% frames and intensities
frselected=[frselected; handles.fr(in)];   
alphaselected=[alphaselected;handles.alpha(in)];  

if isempty(roiselectedx)==1
    disp('No data')
    return
else
    I=handles.I;
    
    if handles.voronoi==0
        
        set(handles.polsize,'Enable','off')
        valintthresh=inthresh/100*max(max(I));
        BW=zeros(size(I,1),size(I,2));
        
        for i=1:size(I,1)
            indexj=find(I(i,:)>valintthresh);
            if isempty(indexj)==0
                BW(i,indexj)=1;
            end
        end
        
        if handles.dilate>0
            se = strel('disk',pxdilate,0);
            I2 = imdilate(BW,se);
            I2 = imfill(I2,'holes');
        else
            I2=BW;
        end
        
        BW=I2;
        if handles.erode>0
            se = strel('disk',pxerode,0);
            I2 = imerode(BW,se);
        end
        
        figurerend= figure('Name','Segmented rendered image','Toolbar','figure');
        imshow(I2,'InitialMagnification','fit')
        title([handles.namefile,', ROI #',num2str(handles.nroroi)])
        
        
    else % voronoi
        
        handles.vorosize= str2num(get(handles.polsize,'String'));
        
        fondo=zeros(size(I,1),size(I,2));
        [~, I2]=MaskVoronoi2(fondo,handles.datax,handles.datay,handles.vorosize); %add nanoclusters!!!!!
    end
end

       
% keep clusters with more than mindens density of points
[labeledall,nb_clust] = bwlabel(I2,4);
maxsize=(((maxsize/2)^2)*pi)/(handles.PALMpx*1000)^2; %!!!!!!!
minsize=(((minsize/2)^2)*pi)/(handles.PALMpx*1000)^2; %!!!!!!!
labeledinloc=[];

% improve variables!
%add rendxmasknano, rendymasknano, results nano, countclustersnano
[rendxmask, rendymask,frmask,alphamask, results,listeboundary,countclusters,auxoutx,auxouty] =makemask5(labeledall, handles.datax,handles.datay,handles.fr,...
    handles.alpha,mindens, mindetect,minsize,maxsize,handles.maxnumberdet,handles.synroi,handles.voronoi);

imagemask=get(handles.textresults,'userdata');

if isempty(rendxmask)==0
  % visualization
  imageclu=zeros(size(I2,1),size(I2,2));
  fond=ones(size(I2,1),size(I2,2));

  for z=1:max(rendxmask(:,4)) %all selected clusters
    index=find(rendxmask(:,4)==z);
    if isempty(index)==0
        for t=1:size(index)
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

    figuremask=figure;
    hold on
    imshow(fond,'InitialMagnification','fit')

    title([handles.namefile,', ROI #',num2str(handles.nroroi)])
    axis off

    plot(handles.datax,handles.datay,'.','MarkerSize',5,'Color',[0.5 0.55 0.55]);
    plot(rendxmask(:,1),rendymask(:,1),'.','MarkerSize',5,'Color','g'); % all clustered points or clusterd points out loc
    imcontour(imageclu,1)

  hold off

  set(handles.textresults,'string',['# of clusters: ',num2str(countclusters-1)]);
  set(handles.acceptpushbutton,'userdata',results);
  set(handles.quitpushbutton,'userdata',figuremask);

  set(handles.textresults,'userdata',I2);

  set(handles.text11,'userdata',rendxmask);
  set(handles.text2,'userdata',rendymask);
  set(handles.text9,'userdata',frmask);
  set(handles.text12,'userdata',alphamask);
  set(handles.text7,'userdata',listeboundary);
  set(handles.text8,'userdata',countclusters);
  
  set(handles.text15,'userdata',auxoutx); % detections out of clusters
  set(handles.text17,'userdata',auxouty);

else
    disp('No clusters left')
end


clear I mask rendxmask rendymask frmask alphamask results listeboundary countclusters maximage

 guidata(gcbo, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function acceptpushbutton_Callback(hObject, eventdata, handles)
 
imagemask=get(handles.textresults,'userdata');
resultsclu=get(handles.acceptpushbutton,'userdata');
rendxmask=get(handles.text11,'userdata');
rendymask= get(handles.text2,'userdata');
frmask= get(handles.text9,'userdata');
alphamask= get(handles.text12,'userdata');
listeboundary= get(handles.text7,'userdata');
countclusters= set(handles.text8,'userdata');

auxoutx=get(handles.text15,'userdata'); % detections out of clusters
auxouty=get(handles.text17,'userdata');

pxdilate=str2num(get(handles.pxdilate,'String'));
pxerode=str2num(get(handles.pxerode,'String'));

limitdens=str2num(get(handles.mindens,'String'));
minnrodetect=str2num(get(handles.mindetect,'String'));
mindiamcluster=str2num(get(handles.minsize,'String'));
maxdiamcluster=str2num(get(handles.maxsize,'String'));
inthresh=str2num(get(handles.intenthresh,'String'));
loc=get(handles.trypushbutton,'userdata');

if handles.voronoi==1 % voronoi
    vorosize=str2num(get(handles.polsize,'String'));
else
    vorosize=[];
end

namefile=handles.namefile;
j=handles.nroroi;
dx=handles.PALMpx;
sigmaloc=dx*2;
listeclucercle=[];
listecluoutcercle=[];
indexcluloc=[];

%convertion area from pixels to µm2
% results=[n size(dataxclu,1) size(datax,1) areaclu size(dataxclu,1)/areaclu typeloc minaxisclu maxaxisclu ];
     % 1 nro cluster - 2 # detect - 3 #detect roi - 4 area cluster - 5 densité - 6 loc- 7 min axis - 8 max axis 

if isempty(resultsclu)==0
    resultsclu(:,4)=resultsclu(:,4).*handles.PALMpx^2; %µm2
    resultsclu(:,5)=resultsclu(:,2)./resultsclu(:,4);
    resultsclu(:,7)=resultsclu(:,7).*handles.PALMpx*1000; %nm
    resultsclu(:,8)=resultsclu(:,8).*handles.PALMpx*1000; %nm
end

I=handles.data;
BW=zeros(size(I,1),size(I,2));
valintthresh=inthresh/100*max(max(I));
for i=1:size(I,1)
    indexj=find(I(i,:)>valintthresh);
    if isempty(indexj)==0
       BW(i,indexj)=1;
    end
end        

if handles.dilate>0
    se = strel('disk',pxdilate,0);
    I2 = imdilate(BW,se);
    I2 = imfill(I2,'holes');
else
    I2=BW;
end
BW=I2;
if handles.erode>0
    se = strel('disk',pxerode,0);
    I2 = imerode(BW,se);
end

% visualization
imageclu=zeros(size(I2,1),size(I2,2));
fond=ones(size(I2,1),size(I2,2));

if isempty(rendxmask)==0 %there are clusters

  for z=1:max(rendxmask(:,4)) %all selected clusters
    index=find(rendxmask(:,4)==z);
    if isempty(index)==0
        for t=1:size(index)
            posx=floor(rendxmask(index(t),1));
            if posx==0; posx=1; end
            posy=floor(rendymask(index(t),1));
            if posy==0; posy=1; end
            imageclu(posy,posx)=1; %new image with pixels that have enough dens
        end
    end
  end
    
else
    auxoutx=handles.datax;
    auxouty=handles.datay;
end % no clusters

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

finalfig=figure;
hold on
    imshow(fond,'InitialMagnification','fit')
    title([handles.namefile,', ROI #',num2str(handles.nroroi)])
    axis off
    
    plot(handles.datax,handles.datay,'.','MarkerSize',5,'Color',[0.5 0.55 0.55]);
    
    if isempty(rendxmask)==0
      plot(rendxmask(:,1),rendymask(:,1),'.','MarkerSize',5,'Color','g'); % all clustered points or clusterd points out loc
    end
if isempty(rendxmask)==0 %there are clusters
    imcontour(imageclu,1)
end

hold off

presentimage2=getframe(gca);  % gets the mask image 
Savename=[namefile,'-roi',num2str(j),'-selection.tif'];
[imagemask2,Map] = frame2im(presentimage2);  
imwrite(imagemask2,Savename,'tif','Compression','none');   
                
% save data detections
xdimmask=size(imagemask,1);
ydimmask=size(imagemask,2);

save([namefile,'-roidata_',num2str(j),'.mat'],'auxoutx','auxouty','imagemask','limitdens','xdimmask', 'ydimmask','rendxmask','rendymask','alphamask','frmask',...
    'minnrodetect','mindiamcluster','maxdiamcluster','dx','sigmaloc','listeclucercle','vorosize','-mat');

figselecmask=get(handles.quitpushbutton,'userdata');

%close(figselecmask) 
close(finalfig)


% analysis not aborted
stopanalysis=0;

% save data
save('auxiliar.mat','imagemask','resultsclu','listeclucercle','inthresh','limitdens','minnrodetect','mindiamcluster','maxdiamcluster','stopanalysis','vorosize','-mat') 

close %MaskRendered

clear  'imagemask' 'mindens'  'xdimmask'  'ydimmask'  'rendxmask'  'rendymask'  'alphamask'  'frmask' 'minnrodetect'  'minsize' 
clear 'maxsize'  'dx'  'sigmaloc' 'presentimage2' 'presentimage' 'listeclucercle'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function quitpushbutton_Callback(hObject, eventdata, handles)
% abort analysis

figselecmask=get(handles.quitpushbutton,'userdata');

% analysis aborted
stopanalysis=1;
save('auxiliar.mat','stopanalysis','-mat') 

close(figselecmask) 
clear threshold handles
close
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
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

function maxsize_Callback(hObject, eventdata, handles)
function maxsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pxdilate_Callback(hObject, eventdata, handles)
function pxdilate_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mindetect_Callback(hObject, eventdata, handles)
function mindetect_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function threshole_Callback(hObject, eventdata, handles)
function threshole_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pxerode_Callback(hObject, eventdata, handles)
function pxerode_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function intenthresh_Callback(hObject, eventdata, handles)
function intenthresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function polsize_Callback(hObject, eventdata, handles)
function polsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function epsilon_Callback(hObject, eventdata, handles)
function epsilon_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function vororadiobutton_Callback(hObject, eventdata, handles)
function nanoradiobutton_Callback(hObject, eventdata, handles)
function skipradiobutton_Callback(hObject, eventdata, handles)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
