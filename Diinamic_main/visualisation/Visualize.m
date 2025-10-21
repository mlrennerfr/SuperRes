function varargout = Visualize(varargin)
%VISUALIZE MATLAB code file for Visualize.fig
% GUI to plot pointillistic and rendered images
%
% Marianne Renner oct2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last Modified by GUIDE v2.5 24-Jun-2022 23:16:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Visualize_OpeningFcn, ...
                   'gui_OutputFcn',  @Visualize_OutputFcn, ...
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
function Visualize_OpeningFcn(hObject, eventdata, handles, varargin)

set(handles.selectfilepushbutton,'userdata',[]);
handles.zoom.x=[];
handles.merge=0;
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%-------------------------------------------------------------------------
function varargout = Visualize_OutputFcn(hObject, eventdata, handles)

varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selectfilepushbutton_Callback(hObject, eventdata, handles)

% load data % simplified version
aux0 = get(handles.selectfilepushbutton,'String');
px_mu = str2num(get(handles.szpx,'String'));
mu=get(handles.muradiobutton,'Value');
handles.zoom.x=[];
set(handles.zoompushbutton,'String','ROI to zoom','FontWeight','Normal');

% data : x,y, fr, alpha, radius
[x,y, fr, alpha, radius,sigma, blink, ratio, z,test1,test2, ffname,filename,control]=openPALMdata(handles,px_mu,0);

set(handles.selectfilepushbutton,'userdata',filename);
handles.ffname = ffname;

disp(' ');
disp('Visualization and analysis');
disp(['File ',ffname,' loaded']);

%disp(max(x))

if exist('Xmatrix','var')
    Nfr = size(Xmatrix,2); 
    Npart = size(Xmatrix,1);
    aux = ~isnan(Xmatrix);
    Npositions = sum(aux(:));
else
    Nfr = max(fr);
    Npart = NaN;
    Npositions = length(x);
end

% Storing
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

clear x y alpha fr radius

% update GUI
set(handles.selectfilepushbutton,'String',aux0);
set(handles.filename1,'String',filename);

% Update handles structure
guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filemergepushbutton_Callback(hObject, eventdata, handles)

% load data merge
px_mu = str2num(get(handles.szpx,'String'));
sigmaloc = str2num(get(handles.widthgaussiansmooth,'String'));
dx = str2num(get(handles.PALMszpx,'String'));
mu=get(handles.muradiobutton,'Value');
handles.ffname2 = [];
handles.x2 = [];
handles.y2 = [];
handles.alphaf2 = [];
handles.fr2 = [];
handles.radius2=[];
handles.sigma2=[];
handles.blink2=[];
handles.ratio2=[];
handles.z2=[];
handles.test12=[];
handles.test22=[];
handles.merge=0;

control=0;
[filename, pathname] = uigetfile('*.*','Select .mat or .tif file for merging');
if filename==0
    set(handles.filename2,'String','')
    set(handles.fireradiobutton,'Enable','on') % no LUT
    set(handles.grayradiobutton,'Enable','on') % no LUT
    return
end
mensaje=msgbox('Please wait');
ffname2 = fullfile(pathname,filename);
[pathstr, name, ext] = fileparts(ffname2);

if strcmp(ext,'.mat')
    % data : x,y, fr, alpha, radius
    [x,y, fr, alpha, radius,sigma, blink, ratio, z,test1,test2, ffname2,filename,control]=openPALMdata(handles,px_mu,1,ffname2);
    control=1;
elseif strcmp(ext,'.tif')
    [~,datamatrix] = tifdataread(ffname2);
    if isstruct(datamatrix)
        datamatrix=double(datamatrix.data);   
    else
        datamatrix=double(datamatrix);  
    end
    control=2;
end
    
handles.ffname2 = ffname2;
disp(' ');
disp('Merging with:');
disp(['File ',ffname2]);

set(handles.filename2,'String',filename)
set(handles.fireradiobutton,'Enable','off') % no LUT
set(handles.grayradiobutton,'Enable','off') % no LUT
handles.merge=1;

if control==1
    % Storing
    handles.x2 = x;
    handles.y2 = y;
    handles.alphaf2 = alpha;
    handles.fr2 = fr;
    handles.radius2=radius;
    handles.sigma2=sigma;
    handles.blink2=blink;
    handles.ratio2=ratio;
    handles.z2=z;
    handles.test12=test1;
    handles.test22=test2;
    handles.datamatrix =[];
    clear x y alpha fr radius
    
elseif control==2 %.tif
    
    % rendered first image to convert size fluo image
    [I,rendx, rendy]=makerendered2(handles,px_mu,sigmaloc,dx,mu);
   % minx=min(rendx)
    %miny=min(rendy)
    
    % verify loc file size, compare to that of rendered image used for ROI
    if size(I,1)==size(datamatrix,1) && size(I,2)==size(datamatrix,2)
        aux=datamatrix;
    else
        % need to change image size
        aux = imresize(datamatrix,[size(I,1)+min(rendy),size(I,2)+min(rendx)]); %considers the "black frame" without detections
    end
    
    %rescale intensity
    Imin=0; Imax=1; %!!!!!!!!!!!!!!!!!!!!!
    aux = double(aux);
    minaux = min(aux(:));
    maxaux = max(aux(:));
    aux = Imin + (aux-minaux)/(maxaux-minaux)*(Imax-Imin);
    
    if size(aux,1)>size(I,1) 
        difx=size(aux,1)-size(I,1);
        dify=size(aux,2)-size(I,2);
        handles.datamatrix=aux(1:size(I,1)+1,dify:(size(aux,2)));
    else
        handles.datamatrix =aux;
    end
    clear datamatrix

else
    handles.merge=0;
    set(handles.filename2,'String','')
    set(handles.fireradiobutton,'Enable','on') % no LUT
    set(handles.grayradiobutton,'Enable','on') % no LUT
    
    clear handles.x2  handles.y2 handles.alphaf2 handles.fr2 handles.radius2
end
    
close(mensaje)

% Update handles structure
guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zoompushbutton_Callback(hObject, eventdata, handles)

mu=get(handles.muradiobutton,'Value');
szpx = str2num(get(handles.szpx,'String'));
dx = str2num(get(handles.PALMszpx,'String'));
nom=get(handles.filename2,'String');
mergeok=1;
if isempty(nom)==1
        mergeok=0;
end

if isempty(handles.zoom.x)==0 % Full image
    set(handles.zoompushbutton,'String','ROI to zoom','FontWeight','Normal');
    handles.zoom.x=[];
else % select ROI
    % image positions
    if ~isfield(handles,'x')
        errordlg('No position data ! Select a file first.');
        return;
    end
    X_mu=handles.x;
    Y_mu=handles.y;
    alpha=handles.alpha;
    fr=handles.fr;
    Z_mu=handles.z;
    
    % remove for visualization
    vect2remove=selectremovealpha(handles,X_mu,Y_mu,alpha);
    if ~isempty(vect2remove)
        X_mu(vect2remove) = [];
        Y_mu(vect2remove) = [];
        alpha(vect2remove) = [];
        fr(vect2remove) = [];
        Z_mu(vect2remove) = [];
    end
    
    if mergeok==1
    
        if isempty(handles.datamatrix)==0 %.tif
            stackmin=min(min(handles.datamatrix));
            stackmax=max(max(handles.datamatrix));
            
            X_mu=X_mu/dx;
            Y_mu=Y_mu/dx;
            xlim([0 max(X_mu)]); ylim([0 max(Y_mu)]);
            
            imshow(handles.datamatrix,[stackmin stackmax],'InitialMagnification','fit');
            color1 = 'b';
            
        else %detections
            
            X_mu2=handles.x2;
            Y_mu2=handles.y2;
            alphaf2=handles.alphaf2;
            fr2=handles.fr2;
            
            % remove
            vect2remove2=selectremovealpha(handles,X_mu2,Y_mu2,alphaf2);
            if ~isempty(vect2remove2)
                [X_mu2,Y_mu2,alphaf2,fr2]=removealpha(vect2remove2,X_mu2,Y_mu2,alphaf2,fr2);
            end

        end
    else
        %color1 = 'k';
    end
    pointsize = str2num(get(handles.dotsize,'String'));
    
    % figure
    fig_points = figure('Name','Select the region to zoom in','Toolbar','figure');
    axis equal;
    hold on
    set(gca,'YDir','reverse')
    xlim([0 max(X_mu)]); ylim([0 max(Y_mu)]);
    xlabel('X (SR pixels)');
    ylabel('Y (SR pixels)');   
    
    I=ones(ceil(max(Y_mu)),ceil(max(X_mu))); 
    imshow(I,'InitialMagnification','fit');
    hold on
    
    if mergeok==1 
        plot(X_mu,Y_mu,'.','MarkerSize',pointsize,'Color','b');
        plot(X_mu2,Y_mu2,'.','MarkerSize',pointsize,'Color','r');
    else
        plot(X_mu,Y_mu,'.','MarkerSize',pointsize,'Color','k');
    end
    
    [BW,x,y]=roipoly;
    x=round(x); %pixelisable
    y=round(y);
    width =max(x)-min(x);
    height=max(y)-min(y);
    
    % pick up ORIGINAL detections and save a new file
    X_mu=handles.x;
    Y_mu=handles.y;
    alpha=handles.alpha;
    fr=handles.fr;
    radius=handles.radius;
    sigma=handles.sigma;
    blink=handles.blink;
    ratio=handles.ratio;
    Z_mu=handles.z;
    test1=handles.test1;
    test2=handles.test2;
    
    in=inpolygon(X_mu,Y_mu,x,y); %index of X_mu,Y_mu points in polygon defined by x,y
    handles.zoom.x=X_mu(in)-min(x); %corrected for the position of the ROI
    handles.zoom.y=Y_mu(in)-min(y); 
    
    handles.zoom.alpha=alpha(in);
    handles.zoom.fr=fr(in);
    handles.zoom.radius=radius(in);
    handles.zoom.sigma=sigma(in);
    handles.zoom.blink=blink(in);
    handles.zoom.ratio=ratio(in);
    handles.zoom.z=Z_mu(in);
    handles.zoom.test1=test1(in);
    handles.zoom.test2=test2(in);
    
    % second image
    if mergeok==1 
        % pick up ORIGINAL detections and save a new file
        X_mu2=handles.x2;
        Y_mu2=handles.y2;
        alpha2=handles.alphaf2;
        fr2=handles.fr2;
        in2=inpolygon(X_mu2,Y_mu2,x,y); %index of X_mu,Y_mu points in polygon defined by x,y
        handles.zoom.x2=X_mu2(in2)-min(x); %corrected for the position of the ROI
        handles.zoom.y2=Y_mu2(in2)-min(y); 
        handles.zoom.alphaf2=alpha2(in2);
        handles.zoom.fr2=fr2(in2);
    end
    
    % plot
    I2 = imcrop(BW, [min(y) min(x) height width] ) ;
    I2=ones(size(I2,2),size(I2,1)); 
    
    % figure
    fig_zoom = figure('Name','Zoomed region','Toolbar','figure');
    set(gca,'YDir','reverse')
    xlim([0 max(handles.zoom.x)]); ylim([0 max(handles.zoom.y)]);
    xlabel('X (SR pixels)');
    ylabel('Y (SR pixels)');   
    axis equal;%axis ij
    hold on
    
    imshow(I2,'InitialMagnification','fit');
    if mergeok==1 
        plot(handles.zoom.x,handles.zoom.y,'.','MarkerSize',pointsize,'Color','b');
        plot(handles.zoom.x2,handles.zoom.y2,'.','MarkerSize',pointsize,'Color','r');
    else
        plot(handles.zoom.x,handles.zoom.y,'.','MarkerSize',pointsize,'Color','k');
    end
    set(handles.zoompushbutton,'String','Full image');
    
end

% Update handles structure
guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savematpushbutton_Callback(hObject, eventdata, handles)

% save .mat file
currentdir=cd;
[pathstr, name, ext] = fileparts(handles.ffname);

if isempty(handles.zoom.x)==1 % Full image
    x=handles.x;
    y=handles.y;
    alpha=handles.alpha;
    fr=handles.fr;
    radius=handles.radius;
    sigma=handles.sigma;
    blink=handles.blink;
    ratio=handles.ratio;
    z=handles.z;
    test1=handles.test1;
    test2=handles.test2;
    def_name=[name,'-clean.mat'];

else %zoom
    
    x=handles.zoom.x;
    y=handles.zoom.y;
    alpha=handles.zoom.alpha;
    fr=handles.zoom.fr;
    radius=handles.zoom.radius;
    sigma=handles.zoom.sigma;
    blink=handles.zoom.blink;
    ratio=handles.zoom.ratio;
    z=handles.zoom.z;
    test1=handles.zoom.test1;
    test2=handles.zoom.test2;
    def_name=[name,'-zoomclean.mat'];
    
end

%% Removing unneeded points from column vectors
vect2remove=selectremovealpha(handles,x,y,alpha);
if ~isempty(vect2remove)
    disp(['Saving new .mat data. Removing ',num2str(length(vect2remove)),' elements from column vectors ...']);
    x(vect2remove) = [];
    y(vect2remove) = [];
    alpha(vect2remove) = [];
    fr(vect2remove) = [];
    radius(vect2remove) = [];
    sigma(vect2remove) = [];
    blink(vect2remove) = [];
    ratio(vect2remove) = [];
    z(vect2remove) = [];
    test1(vect2remove) = [];
    test2(vect2remove) = [];
end

cd(pathstr);
[filename,path] = uiputfile(def_name,'Save image data as :');
cd(path)

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

save(filename, 'matrice_results');

disp('File saved');
cd(currentdir)
clear x y alpha fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pointillistic_Callback(hObject, ~, handles)

nom=get(handles.filename2,'String');
mergeok=1;
if isempty(nom)==1
        mergeok=0;
end
szpx = str2num(get(handles.szpx,'String'));
filename1=get(handles.filename1,'String');
filename2=get(handles.filename2,'String');
mu=get(handles.muradiobutton,'Value');
dx = str2num(get(handles.PALMszpx,'String'));

% image positions
if ~isfield(handles,'x')
    errordlg('No position data ! Select a file first.');
    return;
end
pointsize = str2num(get(handles.dotsize,'String'));

if isempty(handles.zoom.x)==1
    X_mu=handles.x;
    Y_mu=handles.y;
    alpha=handles.alpha;
    fr=handles.fr;
    Z_mu=[];

else
    X_mu=handles.zoom.x;
    Y_mu=handles.zoom.y;
    alpha=handles.zoom.alpha;
    fr=handles.zoom.fr;
    Z_mu=[];
end

% remove
vect2remove=selectremovealpha(handles,X_mu,Y_mu,alpha);
if ~isempty(vect2remove)
    [X_mu,Y_mu,alpha,fr]=removealpha(vect2remove,X_mu,Y_mu,alpha,fr);
end

%pixel size
if mu==1 %already in µm
else
    X_mu=X_mu*szpx;
    Y_mu=Y_mu*szpx;
end

fig_points = figure('Name','Pointillist image','Toolbar','figure');
axis equal;
hold on
set(gca,'YDir','reverse')
xlim([0 max(X_mu)]); ylim([0 max(Y_mu)]);
xlabel('X (mu)');
ylabel('Y (mu)');   

if mergeok==1
    if isempty(handles.datamatrix)==0 %.tif
        stackmin=min(min(handles.datamatrix));
        stackmax=max(max(handles.datamatrix));
        X_mu=X_mu/dx;
        Y_mu=Y_mu/dx;
        xlim([0 max(X_mu)]); ylim([0 max(Y_mu)]);
        xlabel('X (PALM pixels)');
        ylabel('Y (PALM pixels)');   
        imshow(handles.datamatrix,[stackmin stackmax],'InitialMagnification','fit');
        color1 = 'b';
    else %detections
        if isempty(handles.zoom.x)==1
            X_mu2=handles.x2;
            Y_mu2=handles.y2;
            alphaf2=handles.alphaf2;
            fr2=handles.fr2;
        else
            X_mu2=handles.zoom.x2;
            Y_mu2=handles.zoom.y2;
            alphaf2=handles.zoom.alphaf2;
            fr2=handles.zoom.fr2;
        end
        vect2remove2=selectremovealpha(handles,X_mu2,Y_mu2,alphaf2);
        if ~isempty(vect2remove2)
            [X_mu2,Y_mu2,alphaf2,fr2]=removealpha(vect2remove2,X_mu2,Y_mu2,alphaf2,fr2);
        end
        if mu==1 %already in µm
        else
            X_mu2=X_mu2*szpx;
            Y_mu2=Y_mu2*szpx;
        end
        color1 = 'b';
        color2 = 'r';
    end
else
    color1 = 'k';    
end
plot(X_mu,Y_mu,'.','MarkerSize',pointsize,'Color',color1);

if mergeok==1 && isempty(handles.datamatrix)==1 % detections
        plot(X_mu2,Y_mu2,'.','MarkerSize',pointsize,'Color',color2);
        legend([filename1,' in blue'],[filename2,' in red'],'Location','southoutside','Orientation','horizontal');legend('boxoff')
   title([' Selected positions =',num2str(length(X_mu)),'(blue) and ',num2str(length(X_mu2)),'(red)']);
else
   title([' Selected positions =',num2str(length(X_mu))]);
end
hold off

clear x y alpha fr
figure(fig_points);
handles.fig_points = fig_points;

% Update handles structure
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function renderedpushbutton_Callback(hObject, eventdata, handles)

set(handles.saveimagepushbutton,'userdata', []);
szpx = str2num(get(handles.szpx,'String'));
colorcue=get(handles.fireradiobutton,'Value');
nom=get(handles.filename2,'String');
mergeok=1;
if isempty(nom)==1
        mergeok=0;
end

mu=get(handles.muradiobutton,'Value');
sigmaloc = str2num(get(handles.widthgaussiansmooth,'String'));
dx = str2num(get(handles.PALMszpx,'String'));
fr=handles.fr;

if isempty(handles.zoom.x)==1
    X_mu=handles.x;
    Y_mu=handles.y;
    alpha=handles.alpha;
    Z_mu=handles.z;
else
    X_mu=handles.zoom.x;
    Y_mu=handles.zoom.y;
    alpha=handles.zoom.alpha;
    Z_mu=handles.zoom.z;
end

if mu==1 %already in µm
    xdim=ceil(max(X_mu))/dx;
    ydim=ceil(max(Y_mu))/dx;
else
    xdim=ceil(max( X_mu))*szpx/dx;
    ydim=ceil(max( Y_mu))*szpx/dx;
    X_mu =  X_mu*szpx;
    Y_mu =  Y_mu*szpx;
end
alpha(1:10);
gausssmooth = 1;
auxstring = get(handles.renderedpushbutton,'String'); %%%%%%%%%%%%%%%%%

% remove
vect2remove=selectremovealpha(handles,X_mu,Y_mu,alpha);
if ~isempty(vect2remove)
    [X_mu,Y_mu,alpha,fr]=removealpha(vect2remove,X_mu,Y_mu,alpha,fr);
end
    
% Create rendered image
maxx=ceil(max(X_mu/dx));
maxy=ceil(max(Y_mu/dx));
[I,xx1,yy1] = PALM_rendering3(X_mu,Y_mu,alpha,sigmaloc,dx,0,xdim, ydim, gausssmooth);

if size(I,2)>maxx
    I=I(:,1:maxx); %elimination black pixels
end
if size(I,1)>maxy
    I=I(1:maxy,:); %elimination black pixels
end

% Normalize image to maximum intensity of first frame or to user-defined saturation intensity
if isfield(handles,'Imax')
    Imax1 = handles.Imax;
    I = I/Imax1;
end
%---------------------------------------------------------

if mergeok==1
   if isempty(handles.datamatrix)==0 %.tif
      I2=handles.datamatrix;
   else
    
    % Image 2 -----------------------------------------------------------
   
    X_mu2=handles.x2;
    Y_mu2=handles.y2;
    alphaf2=handles.alphaf2;
    fr2=handles.fr2;
    % remove
    vect2remove2=selectremovealpha(handles,X_mu2,Y_mu2,alphaf2);
    if ~isempty(vect2remove2)
        [X_mu2,Y_mu2,alphaf2,fr2]=removealpha(vect2remove2,X_mu2,Y_mu2,alphaf2,fr2);
    end
    
    %pixel size
    if mu==1 %already in µm
    else
        X_mu2=X_mu2*szpx;
        Y_mu2=Y_mu2*szpx;
    end

    % Create rendered image
    maxx2=ceil(max(X_mu2/dx));
    maxy2=ceil(max(Y_mu2/dx));
    
    [I2,xx12,yy12] = PALM_rendering3( X_mu2,Y_mu2,alphaf2,sigmaloc,dx,0,xdim, ydim, gausssmooth);

    if size(I2,2)>maxx
        I2=I2(:,1:maxx); %elimination black pixels
    end
    if size(I2,1)>maxy
        I2=I2(1:maxy,:); %elimination black pixels
    end
   
   end %type of file
end

handles.fig_render{1} = figure('Name','PALM rendered image');

% save data rendered in .rnd file
[pathstr, name, ext] = fileparts(handles.ffname);
currentdir=cd;

%Plot
xlabel('x (mu)');    ylabel('y (mu)');    axis equal xy;
xlim([0 max(X_mu)]);    ylim([0 max(Y_mu)]);    hold on;

if exist('forlegend','var')
    text(xmin,ymin ,forlegend{it},'Color','y','VerticalAlignment','bottom');
end
title(['PALM rendered image with N=',num2str(length(X_mu)),' points']);

if mergeok==0
    image_ch(xx1,yy1,I,[0 1]); % image_ch(xx1,yy1,I);
    colorbar;
    
    if colorcue==1
        colormap(hot);            % colorbar;
    elseif colorcue==0
        colormap(gray)
    end
    
    if exist('I','var') % && ~surfaceplot; %&& ~playsavemovie
        imcontrast(gca);
    end
    filename=[name,'.rnd'];
    
else % merge==1
    
    if size(I2,1)>size(I,1)
        auxiliar=[];
        auxiliar=I2(1:size(I,1),:);
        I2=auxiliar;
        clear auxiliar
    end
    if size(I2,2)>size(I,2)
        auxiliar=[];
        auxiliar=I2(:,1:size(I,2));
        I2=auxiliar;
        clear auxiliar
    end
    
    Imerge=cat(3,I,zeros(size(I,1),size(I,2)),I2) ;
    image_ch(xx1,yy1,Imerge,[0 1]); % image_ch(xx1,yy1,I);
    filename=[name,'-merged.rnd'];
end
     
alphamin = str2num(get(handles.alpha1,'String'));
alphamax = str2num(get(handles.alpha2,'String'));
xdim=size(I,1);
ydim=size(I,2);

% Save rendered data
if isdir('rendered'); else mkdir('rendered'); end 
cd('rendered')

rendx=X_mu/dx;
rendy=Y_mu/dx;
rendy=abs(rendy-(max(rendy)))+1;
handles.rendx=rendx;
handles.rendy=rendy;
save(filename,'xdim', 'ydim','rendx','rendy','alpha','fr','alphamin','alphamax','dx','gausssmooth','sigmaloc','-mat');

cd(currentdir)
set(handles.renderedpushbutton,'String',auxstring);
set(handles.renderedpushbutton,'userdata',I);

clear X_mu Y_mu alpha fr rendx rendy

% Update handles structure
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveimagepushbutton_Callback(hObject, eventdata, handles)

filename=get(handles.filename1,'string');
namefile=strtok(filename,'.');
szpx=str2num(get(handles.szpx,'string'));
PALMszpx=str2num(get(handles.PALMszpx,'string'));
merge= handles.merge;
mu=get(handles.muradiobutton,'Value');

if isfield(handles,'fig_render')
else
   renderedpushbutton_Callback(hObject,eventdata , handles)
end
fig = handles.fig_render{1};
aux= findobj(fig,'Type','Axes');
if length(aux)>1 % makes sure we don't use the colorbar
    if  strcmp(get(aux(1),'Tag'),'Colorbar')
        aux = aux(2);
    else
        aux = aux(1);
    end
end
aux1 = get(aux,'CLim');

if isempty(aux1)==0
    Imax = aux1(2);
    for iframe = 1:length(handles.fig_render)
        aux4 = findobj(handles.fig_render{iframe},'Type','image'); % look for image in figure axes corresponding to frame # iframe
        if length(aux4)>1, % makes sure we take the image and not the colorbar by looking for the largest image data
            for iaux = 1:length(aux4)
                aux5(iaux) = prod(size(get(aux4(iaux),'CData')));
            end
            [aux6 iauxmax] = max(aux5);
            aux4 = aux4(iauxmax);
        end
        I = get(aux4,'CData');
        I = I/Imax; 
        if isempty(handles.zoom.x)==0 %zoom
            [namefile,~] = uiputfile([namefile,'-zoom.tif'],'Save zoomed image as'); 
        end
        if merge==0
            Savename1=[namefile,'-rend.tif'];
        else
            Savename1=[namefile,'-merged-rend.tif'];
        end
        imwrite(I,Savename1,'tif','Compression','none');  % save the image to the disk
        disp(['Rendered image saved as ',Savename1,' normalized to the maximum intensity ',num2str(Imax)]);
    end
else
    msgbox('Create a rendered image first','Error','error')
end

clear iI I

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dotsize_Callback(hObject, eventdata, handles)
function dotsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function widthgaussiansmooth_Callback(hObject, eventdata, handles)
function widthgaussiansmooth_CreateFcn(hObject, eventdata, handles)
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

function muradiobutton_Callback(hObject, eventdata, handles)

function fireradiobutton_Callback(hObject, eventdata, handles)
code=get(hObject, 'value');
if code==1
    set(handles.grayradiobutton,'Value',0);
end

function grayradiobutton_Callback(hObject, eventdata, handles)
code=get(hObject, 'value');
if code==1
    set(handles.fireradiobutton,'Value',0);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%