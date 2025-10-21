function varargout = CreateROI(varargin)
% function varargout = CreateROI(varargin)
%
% independent window to create and visualize ROI
%
% called by Segment.m
%
% Marianne Renner may 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last Modified by GUIDE v2.5 14-May-2021 11:29:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CreateROI_OpeningFcn, ...
                   'gui_OutputFcn',  @CreateROI_OutputFcn, ...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes just before CreateROI is made visible.
function CreateROI_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

set(handles.acceptpushbutton,'Enable','off')
set(handles.clearpushbutton,'Enable','off')
set(handles.transferpushbutton,'Value',0)  % no transfer to other image
set(handles.text5,'Userdata',[]); %to save data about ROIs
set(handles.backgroundpushbutton,'Userdata',[]);

% to fit screen definition/changes in size
h3 = findobj('Type','figure');
txtHand = findall(h3, '-property', 'FontUnits');
set(txtHand, 'FontUnits', 'normalized');
%set(txtHand, 'FontUnits', 'centimeters');

backgroundim=get(handles.backgroundpushbutton,'Userdata');

handles2=cell2mat(varargin{1});

lengthvalue=get(handles2.lengthradiobutton,'value');
set(handles.text7,'value',lengthvalue)

szpx=str2num(get(handles2.szpx,'string'));
dx = str2num(get(handles2.PALMszpx,'string'));
mu=get(handles2.muradiobutton,'Value');
%handles.smallimage=get(handles2.smallradiobutton,'Value');

X_mu=handles2.x;
% change Y coordinates for plotting
% (not for small images: less than 5 µm)
if max(X_mu)<5 
    handles.smallimage=1; 
    Y_mu=handles2.y;
else
    handles.smallimage=0; 
    Y_mu=handles2.y-szpx;
end

alpha = handles2.alpha(:);
fr = handles2.fr(:);
alpha(1:10);
handles.xdim=max(X_mu);
handles.ydim=max(Y_mu);

set(handles.clearpushbutton,'Userdata',dx)
set(handles.text8,'Userdata',mu)
set(handles.text8,'Userdata',szpx)

filename=get(handles2.filename, 'Userdata') ;
[namefile,~]=strtok(filename,'.'); %sin extension
set(handles.filename,'String',filename)
    
%vect2remove=selectremovealpha2(handles2,X_mu,Y_mu,alpha);
%if ~isempty(vect2remove)
%    [X_mu,Y_mu,alpha,fr]=removealpha(vect2remove,X_mu,Y_mu,alpha,fr);
%end
%meanintens=mean(alpha(:));
set(handles.newroipushbutton,'Userdata',[X_mu Y_mu alpha fr])

axes(handles.axes1);  
xlim([0 handles.xdim]); ylim([0 handles.ydim]);
axis equal
set(gca,'YDir','reverse')


xticks([])
yticks([])
xticklabels({})
yticklabels({})

if handles.smallimage==1
    xlim([0 handles.xdim*100]); ylim([0 handles.ydim*100]);
    I=ones(ceil(handles.ydim*100),ceil(handles.xdim*100));
    imshow(I,'InitialMagnification','fit');
    hold on
    plot(X_mu*100,Y_mu*100,'.','MarkerSize',2,'Color','k');
else
    I=ones(ceil(handles.ydim),ceil(handles.xdim)); %!!!!!!!!!!!!!!!!!!
    imshow(I,'InitialMagnification','fit');
    hold on
    plot(X_mu,Y_mu,'.','MarkerSize',2,'Color','k');
end

count=1;
logicalroi=1;
nroorder=0;

%if ROIs, plot ROIs
deleteok=0;
if exist([namefile,'.rgn'],'file')
    savenamedataroi=[namefile,'-ROIinfo.txt'];
    droi=dir(savenamedataroi);
    if isempty(droi)
        alldataroi=[];
    else
        alldataroi=load(savenamedataroi);
    end
    set(handles.text5,'Userdata',alldataroi); 
    data=importdata([namefile,'.rgn']);
    
    if isstruct(data)
        roifile=data.coord;
        dx=data.dx;
        set(handles.roinumber,'String','all')
        set(handles.roinumber,'Userdata',size(roifile,2))
        set(handles.clearpushbutton,'Enable','on')
        if isfield(data,'dist')
         else
            if lengthvalue==1
                answer = questdlg('There is previous ROI data without lenght measurement. Delete and proceed?','Confirm choice','Yes','No','so');
                % Handle response
                switch answer
                    case 'Yes'
                        deleteok=1;
                    case 'No'
                        deleteok=2;
                        string= ['Please close the window'];
                        hbox=msgbox(string);
                end
            end                
        end
    else
        roifile=data;
    end %if sstruct
    
    if deleteok==0
        disp(['Regions from ',namefile,'.rgn loaded'])
        disp(' ')
        for roij=1:size(roifile,2) %all ROIs in the image
            if isempty(roifile{roij}(:,1))==0
                xi=data.x{roij}/dx;
                yi=data.y{roij}/dx;
                if handles.smallimage==1
                    xi=xi*100;
                    yi=yi*100;
                end
                plot(xi,yi,'-r')
                text(mean(xi),mean(yi),num2str(roij),'Color','b','FontSize',16,'FontWeight','Bold')
                hold on
           end
        end
        ROI=data;
        set(handles.acceptpushbutton,'Userdata',ROI);
        set(handles.roinumber,'String',num2str(size(roifile,2)))
        set(handles.roinumber,'Userdata',size(roifile,2))
        if lengthvalue==1
            % measures length
            string= ['You chose to measure the length of ROIs. Use the mouse to draw the line on the rendered image : click to choose each point. Press Backspace or Delete to remove the previously selected point. A double-click ends the selection'];
            hbox=msgbox(string);
        end
    elseif deleteok==1
       delete([namefile,'.rgn'])
       set(handles.acceptpushbutton,'Userdata',[]);
       set(handles.roinumber,'Userdata',0); % number of  ROI
       set(handles.roinumber,'String','0'); % number of  ROI
       disp(['File ',[namefile,'.rgn'],' deleted']);
       if lengthvalue==1
            % measures length
            string= ['You chose to measure the length of ROIs. Use the mouse to draw the line on the rendered image : click to choose each point. Press Backspace or Delete to remove the previously selected point. A double-click ends the selection'];
            hbox=msgbox(string);
        end
    end
else
    set(handles.roinumber,'String','0')
    set(handles.roinumber,'Userdata',0)
end

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------
function varargout = CreateROI_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function backgroundpushbutton_Callback(hObject, eventdata, handles)
% loads and reads background image

disp('Reading image, please wait...')

%mu=get(handles.text8,'Value');
%szpx=get(handles.text8,'Userdata');
dx=get(handles.clearpushbutton,'Userdata');
dataset=get(handles.newroipushbutton,'Userdata'); %[X_mu Y_mu alpha fr]
ROI=get(handles.acceptpushbutton,'Userdata');
count= get(handles.roinumber,'Userdata'); % number of next ROI

[file,path] = uigetfile('*.tif','Load background file');
filename = [path,file];
if path>0
else
    set(handles.backgroundpushbutton,'Userdata',[]);
    set(handles.textbackground,'String','');
    return
end

%[stack_info,datamatrix] = tifdataread(filename);
datamatrix=double(imread(filename));
datamatrix=checkfluoimagesize2(dataset(:,1), dataset(:,2),dataset(:,3),datamatrix,dx);

set(handles.backgroundpushbutton,'Userdata',datamatrix);
set(handles.textbackground,'String',file);

%plot
axes(handles.axes1);  
xlim([0 size(datamatrix,2)]); ylim([0 size(datamatrix,1)]);

%figure
stackmin=min(min(datamatrix));
stackmax=max(max(datamatrix));
imshow(datamatrix,[stackmin stackmax],'InitialMagnification','fit');
hold on
plot(dataset(:,1)/dx,dataset(:,2)/dx,'.','MarkerSize',2,'Color','r');
%axis equal
%set(gca,'YDir','reverse')
if count>0   
    for n=1:count
        plot(ROI.coord{n}(:,1)/dx,ROI.coord{n}(:,2)/dx,'-g');
       % plot((ROI.x{n}/dx)/dx,(ROI.y{n}/dx)/dx,'-b');
        text(mean(ROI.coord{n}(:,1))/dx,mean(ROI.coord{n}(:,2))/dx,num2str(n),'Color','g','fontsize',14,'FontWeight','bold');
    end
end

disp('Done')
disp(' ')

guidata(gcbo, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function transferpushbutton_Callback(hObject, eventdata, handles)
% reads info second image to save ROI on this image too

mu=get(handles.text8,'Value');
szpx=get(handles.text8,'Userdata');

%currentcd=cd;
[file,path] = uigetfile('*.mat*','Load associated file');
if file==0
    set(handles.transferpushbutton,'Value',0);
    return
end

% read data
if mu==0
    handles=loadselectPALMdata(file,szpx,handles);
else
    handles=loadselectPALMdata(file,1,handles); %data already in microns
end
    
set(handles.filetransfer,'Value',1)
set(handles.filetransfer,'String',file); 

X_mu2=handles.x;
Y_mu2=handles.y;
fr2=handles.fr;
alpha2 = handles.alpha(:);
alpha2(1:10);

set(handles.transferpushbutton,'Userdata',[X_mu2 Y_mu2 alpha2 fr2])
%cd(currentcd)

guidata(gcbo, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newroipushbutton_Callback(hObject, eventdata, handles)

transfer=get(handles.filetransfer,'Value');
dx=get(handles.clearpushbutton,'Userdata');
szpx=get(handles.text8,'Userdata');

fullimage=get(handles.fullimageradiobutton,'value');
lengthvalue=get(handles.text7,'value');
ROI=get(handles.acceptpushbutton,'Userdata');
count= get(handles.roinumber,'Userdata'); % number of next ROI
backgroundim=get(handles.backgroundpushbutton,'Userdata');
dataset=get(handles.newroipushbutton,'Userdata'); %[X_mu Y_mu alpha fr]

disp('Preparing tool, please wait...')

aux=get(handles.newroipushbutton,'Userdata');
X_mu=aux(:,1);
Y_mu=aux(:,2);
alpha=aux(:,3);
fr=aux(:,4);
clear aux

xdim=max(X_mu);
ydim=max(Y_mu);

if transfer==1
    ROI2=get(handles.filetransfer,'Userdata');
    aux2=get(handles.transferpushbutton,'Userdata');
    X_mu2=aux2(:,1);
    Y_mu2=aux2(:,2);
    alpha2=aux2(:,3);
    fr2=aux2(:,4);
    clear aux2
end

% refresh image
axes(handles.axes1);  
axis equal
set(gca,'YDir','reverse')

xticks([])
yticks([])
xticklabels({})
yticklabels({})

if isempty(backgroundim)==0
    I=backgroundim;
    %plot
    xlim([0 size(I,2)]); ylim([0 size(I,1)]);
    stackmin=min(min(backgroundim));
    stackmax=max(max(backgroundim));
    imshow(backgroundim,[stackmin stackmax],'InitialMagnification','fit');
    hold on
    plot(dataset(:,1)/dx,dataset(:,2)/dx,'.','MarkerSize',2,'Color','r');
    if count>0   
        for n=1:count
            plot(ROI.coord{n}(:,1)/dx,ROI.coord{n}(:,2)/dx,'-g');
            text(mean(ROI.coord{n}(:,1)/dx),mean(ROI.coord{n}(:,2)/dx),num2str(n),'Color','g','fontsize',14,'FontWeight','bold');
        end
    end

else
    if handles.smallimage==1
        xlim([0 xdim*100]); ylim([0 ydim*100]);
        I=ones(ceil(ydim*100),ceil(xdim*100));
        imshow(I,'InitialMagnification','fit');
        hold on
        plot(X_mu*100,Y_mu*100,'.','MarkerSize',2,'Color','k');
    else
        I=ones(ceil(ydim),ceil(xdim));
        imshow(I,'InitialMagnification','fit');
        hold on
        plot(X_mu,Y_mu,'.','MarkerSize',2,'Color','k');
    end
    if count>0   
        for n=1:count
            if handles.smallimage==1
                plot(ROI.coord{n}(:,1)*100,ROI.coord{n}(:,2)*100,'-r');
                text(mean(ROI.coord{n}(:,1)*100),mean(ROI.coord{n}(:,2)*100),num2str(n),'Color','b','fontsize',14,'FontWeight','bold');
            else
                plot(ROI.coord{n}(:,1),ROI.coord{n}(:,2),'-r');
                text(mean(ROI.coord{n}(:,1)),mean(ROI.coord{n}(:,2)),num2str(n),'Color','b','fontsize',14,'FontWeight','bold');
            end
        end
    end

end

count=count+1;

disp('Done. Please draw the ROI')

if fullimage==0
    [BWROI,x,y]=roipolyold;  
    %[~,x,y]=roipoly;  
    if handles.smallimage==1 % image for ROI 100 times bigger
        x=x/100;
        y=y/100;
    end
else
    x=[0;ceil(xdim);ceil(xdim);0];
    y=[0;0;ceil(ydim);ceil(ydim)];
end

% area ROI
sizexconv=ceil(ydim/dx);    sizeyconv=ceil(xdim/dx);
if isempty(backgroundim)==0
    xconv=x;    yconv=y;
else
    xconv=x/dx;    yconv=y/dx;
end
imagenoire=zeros(sizexconv,sizeyconv);


BW = roipoly(imagenoire,xconv,yconv); %in nm
stats = regionprops(BW,'Area');
surface=stats(1).Area;
surface=surface*(dx^2);

% plot polygon ROI
coord=[x y];
if isempty(backgroundim)==0
    plot(x,y,'-g')
else
    if handles.smallimage==1 
            % image for ROI 100 times bigger
            plot(x*100,y*100,'-r');
    else
        plot(x,y,'-r');
    end

end
hold on  

disp('Done')
%pick points roi
roiselectedx=[];
roiselectedy=[];   
Xcorrected=X_mu;
Ycorrected=Y_mu;
if isempty(backgroundim)==0
    in = inpolygon(Xcorrected,Ycorrected,x*dx,y*dx);
    roiselectedx=[roiselectedx; Xcorrected(in)-ceil(min(x*dx))+1];
    roiselectedy=[roiselectedy; Ycorrected(in)-ceil(min(y*dx))+1];
else
    in = inpolygon(Xcorrected,Ycorrected,x,y);
    roiselectedx=[roiselectedx; Xcorrected(in)-ceil(min(x))+1];
    roiselectedy=[roiselectedy; Ycorrected(in)-ceil(min(y))+1];
end

alpharoi=alpha(in);
frroi=fr(in);

if transfer==1
        %pick points roi image 2
        roiselectedx2=[];
        roiselectedy2=[];  
        Xcorrected2=X_mu2;
        Ycorrected2=Y_mu2;
        if isempty(backgroundim)==0
            in2 = inpolygon(Xcorrected2,Ycorrected2,x*dx,y*dx);
            roiselectedx2=[roiselectedx2; Xcorrected2(in2)-ceil(min(x*dx))+1];
            roiselectedy2=[roiselectedy2; Ycorrected2(in2)-ceil(min(y*dx))+1];
        else
            in2 = inpolygon(Xcorrected2,Ycorrected2,x,y);
            roiselectedx2=[roiselectedx2; Xcorrected2(in2)-ceil(min(x))+1];
            roiselectedy2=[roiselectedy2; Ycorrected2(in2)-ceil(min(y))+1];
        end
        frroi2=fr2(in2);
        alpharoi2=alpha2(in2);
end
xdim2=ceil(max(roiselectedx));
ydim2=ceil(max(roiselectedy));

%figure
figapart =figure('Name',['ROI #', num2str(count)]);

% background???
if handles.smallimage==1
    plot(roiselectedx,roiselectedy,'.k','MarkerSize',5)
else
    Iroi=ones(ceil(ydim2),ceil(xdim2)); %!!!!!!!!!!!!!!!!!!
    imshow(Iroi,'InitialMagnification','fit');
    hold on
    plot(roiselectedx, roiselectedy,'.k','MarkerSize',1);
end

xlabel(['Area: ',num2str(surface),' µm2 - Detections : ',num2str(size(roiselectedx,1)),' - Density: ',num2str(size(roiselectedx,1)/surface),' det/µm2']);
hold off
roidata=[count surface size(roiselectedx,1) size(roiselectedx,1)/surface];
set(handles.text2,'Userdata',roidata);

% mask of the ROI (the same size than the original image)
xdimrend=xdim2/dx;
ydimrend=ydim2/dx;
alpharoi(1:10);
gausssmooth = 1;
[imageroirend,rendx, rendy] = PALM_rendering3( roiselectedx,roiselectedy,alpharoi,dx*2,dx,0,xdimrend, ydimrend, gausssmooth);

figrendered= figure('Name','Rendered image overlaid with detections','Toolbar','figure');
imshow(imageroirend,'InitialMagnification','fit')
hold on
plot(roiselectedx/dx, roiselectedy/dx,'.r','MarkerSize',1);
if lengthvalue==1  %measure of distance
    title('Use the mouse to draw a line to measure the distance')
end

newxdim=ceil(X_mu);
newydim=ceil(Y_mu);
ROI.fr{count}=frroi;
if isempty(backgroundim)==0
    ROI.coord{count}=coord*dx;
    ROI.x{count}=x*dx*dx;
    ROI.y{count}=y*dx*dx;
else
    ROI.coord{count}=coord;
    ROI.x{count}=x*dx;
    ROI.y{count}=y*dx;
end

ROI.in{count}=in; %index of selected points
ROI.alpha{count}=alpharoi;
ROI.xselec{count}=roiselectedx; %index of selected points
ROI.yselec{count}=roiselectedy; %index of selected points
%ROI.imagerend{count}=imageroirend; %crop of rendered image
ROI.mock=imageroirend; %crop of rendered image to be saved as an independant file
ROI.dimx=[1 max(X_mu)];
ROI.dimy=[1 max(Y_mu)];
%ROI.dim=[newxdim newydim];
ROI.dx=dx;

if lengthvalue==1
    disp('Draw a line to measure the distance')
    set(0,'CurrentFigure',figrendered)
    set(figrendered,'Pointer','cross')
    [cx,cy,mock]=improfile;
    plot(cx,cy,'-g');
    clear mock
    dist=0;
    for i=2:size(cx,1)
        dist=dist+ sqrt((cx(i)-cx(i-1))^2+(cy(i)-cy(i-1))^2); %length in nm...???
    end
   % dist=dist/1000; %length in µm
    ROI.dist{count}=dist;
    disp(['Distance in nm: ',num2str(dist)])
    close(figrendered)

end

if transfer==1
        [imageroirend2,rendx2, rendy2] = PALM_rendering3( roiselectedx2,roiselectedy2,alpharoi2,dx*2,dx,0,xdimrend, ydimrend, gausssmooth);
        
        ROI2.in{count}=in2; %index of selected points
        if isempty(backgroundim)==0
            ROI2.coord{count}=coord*dx;
            ROI2.x{count}=x*dx*dx;
            ROI2.y{count}=y*dx*dx;
        else
            ROI2.coord{count}=coord*dx;
            ROI2.x{count}=x*dx;
            ROI2.y{count}=y*dx;
        end
        ROI2.fr{count}=frroi2;
        ROI2.alpha{count}=alpharoi2;
        ROI2.dimx=[1 max(X_mu)];
        ROI2.dimy=[1 max(Y_mu)];
        ROI2.dim=[newxdim newydim];
        ROI2.dx=dx;
        ROI2.xselec{count}=roiselectedx2; 
        ROI2.yselec{count}=roiselectedy2; 
       % ROI2.imagerend{count}=imageroirend2; %crop of rendered image
        ROI2.mock=imageroirend2; %crop of rendered image
        
        if lengthvalue==1
            ROI2.dist{count}=dist;
        end
        set(handles.filetransfer,'Userdata',ROI2)
end

set(handles.acceptpushbutton,'Userdata',ROI)
set(handles.acceptpushbutton,'Enable','on')

clear I x y meanintens xdim ydim ROI ROI2 aux aux2
disp(' ')

% Update handles structure
guidata(gcbo, handles);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clearpushbutton_Callback(hObject, eventdata, handles)

transfer=get(handles.filetransfer,'Value');
ROI=get(handles.acceptpushbutton,'Userdata');
dx=get(handles.clearpushbutton,'Userdata');

count= get(handles.roinumber,'Userdata'); % number of ROI
filename=get(handles.filename,'String');
[namefile,~]=strtok(filename,'.'); %sin extension
savename=[namefile,'.rgn'];
detec=get(handles.newroipushbutton,'Userdata');
X_mu=detec(:,1);
Y_mu=detec(:,2);
alldataroi=get(handles.text5,'Userdata');
savenamedataroi=[namefile,'-ROIinfo.txt'];
backgroundim=get(handles.backgroundpushbutton,'Userdata');
dataset=get(handles.newroipushbutton,'Userdata'); %[X_mu Y_mu alpha fr]

if transfer==1
    filename2=get(handles.filetransfer,'String');
    [namefile2,~]=strtok(filename2,'.'); %sin extension
    savename2=[namefile2,'.rgn'];
    ROI2=get(handles.filetransfer,'Userdata');
end

xdim=max(X_mu);
ydim=max(Y_mu);
nroroi=str2num(get(handles.nroROIclear,'String'));

lengthvalue=get(handles.text7,'value');
clearall=get(handles.allradiobutton,'Value');
order=1;

if clearall==0 
   
    if nroroi<0 || nroroi>count
        msgbox('Error','Please indicate a valid ROI number')
    else
        if count==1 % only one ROI : delete file
            delete(savename)
            delete(savenamedataroi)
            disp(['File ',savename,' deleted']);

            if transfer==1
               delete(savename2)
               disp(['File ',savename2,' deleted']);
               set(handles.filetransfer,'Userdata',[]);
            end
            set(handles.acceptpushbutton,'Userdata',[]);
            set(handles.roinumber,'Userdata',0); % number of  ROI
            set(handles.roinumber,'String','0'); % number of  ROI
            count=0;
            set(handles.acceptpushbutton,'Enable','off')
            set(handles.clearpushbutton,'Enable','off')
            set(handles.text5,'Userdata',[]); 
            set(handles.text2,'Userdata',[]); 
        else
            auxdataroi=[];
            for i=1:count
                if i==nroroi % to be deleted
                else
                    if isempty(auxdataroi)==0
                        auxdataroi(order,:)=alldataroi(i,:);
                        auxdataroi(order,1)=order;
                    else
                        auxdataroi=[];
                    end
                    aux.coord{order}= ROI.coord{i};
                    aux.in{order}= ROI.in{i};
                    aux.xselec{order}= ROI.xselec{i};
                    aux.yselec{order}= ROI.yselec{i};
                  %  aux.imagerend{order}= ROI.imagerend{i};
                    aux.x{order}= ROI.x{i};
                    aux.y{order}= ROI.y{i};
                    aux.dimx= ROI.dimx;
                    aux.dimy= ROI.dimy;
                   % aux.dim= ROI.dim;
                    aux.dx= ROI.dx;
                    if lengthvalue==1
                        aux.dist{order}= ROI.dist{i};
                    end
                    if transfer==1
                        aux2.coord{order}= ROI.coord{i};
                        aux2.x{order}= ROI.x{i};
                        aux2.y{order}= ROI.y{i};
                        aux2.dimx= ROI.dimx;
                        aux2.dimy= ROI.dimy;
                      %  aux2.dim= ROI.dim;
                        aux2.dx= ROI.dx;
                        if lengthvalue==1
                            aux2.dist{order}= ROI.dist{i};
                        end
                        aux2.in{order}= ROI2.in{i};
                        aux2.xselec{order}= ROI2.xselec{i};
                        aux2.yselec{order}= ROI2.yselec{i};
                      %  aux2.imagerend{order}= ROI2.imagerend{i};
                    end
                    order=order+1;
                end % if count
            end % for count
            alldataroi=auxdataroi;
            ROI=aux;
            if transfer==1
                ROI2=aux2;
                set(handles.filetransfer,'Userdata',ROI2);
            end
            count=order-1;
            set(handles.acceptpushbutton,'Userdata',ROI);
            set(handles.text5,'Userdata',alldataroi);
            set(handles.roinumber,'Userdata',count); % number of  ROI
            set(handles.roinumber,'String',num2str(count)); % number of  ROI
            save([namefile,'.rgn'],'ROI','-mat') 
            if isempty(auxdataroi)==0
                save([namefile,'-ROIinfo.txt'],'alldataroi','-ascii') 
            end

        end % count=1
    end % verif nroroi

else % clear all
    delete(savename)
    delete(savenamedataroi)

    disp(['File ',savename,' deleted']);
    if transfer==1
        delete(savename2)
        disp(['File ',savename2,' deleted']);
        set(handles.filetransfer,'Userdata',[]);
    end
    set(handles.acceptpushbutton,'Userdata',[]);
    set(handles.roinumber,'Userdata',0); % number of  ROI
    set(handles.roinumber,'String','0'); % number of  ROI
    count=0;
    set(handles.acceptpushbutton,'Enable','off')
    set(handles.clearpushbutton,'Enable','off')
end

% refresh image
axes(handles.axes1);  
xlim([0 xdim]); ylim([0 ydim]);
axis equal
set(gca,'YDir','reverse')

xticks([])
yticks([])
xticklabels({})
yticklabels({})

if isempty(backgroundim)==0
    I=backgroundim;
    %plot
    xlim([0 size(I,2)]); ylim([0 size(I,1)]);
    stackmin=min(min(backgroundim));
    stackmax=max(max(backgroundim));
    imshow(backgroundim,[stackmin stackmax],'InitialMagnification','fit');

    %imshow(backgroundim,'InitialMagnification','fit');
    hold on
    plot(dataset(:,1)/dx,dataset(:,2)/dx,'.','MarkerSize',2,'Color','r');
    if count>0   
        for n=1:count
            plot(ROI.coord{n}(:,1)/dx,ROI.coord{n}(:,2)/dx,'-g');
            text(mean(ROI.coord{n}(:,1)/dx),mean(ROI.coord{n}(:,2)/dx),num2str(n),'Color','g','fontsize',14,'FontWeight','bold');
        end
    end

else
    if handles.smallimage==1
        xlim([0 xdim*100]); ylim([0 ydim*100]);
        I=ones(ceil(ydim*100),ceil(xdim*100));
        imshow(I,'InitialMagnification','fit');
        hold on
        plot(X_mu*100,Y_mu*100,'.','MarkerSize',2,'Color','k');
    else
        I=ones(ceil(ydim),ceil(xdim));
        imshow(I,'InitialMagnification','fit');
        hold on
        plot(X_mu,Y_mu,'.','MarkerSize',2,'Color','k');
    end
    if count>0   
        for n=1:count
            if handles.smallimage==1
                plot(ROI.coord{n}(:,1)*100,ROI.coord{n}(:,2)*100,'-r');
                text(mean(ROI.coord{n}(:,1)*100),mean(ROI.coord{n}(:,2)*100),num2str(n),'Color','b','fontsize',14,'FontWeight','bold');
            else
                plot(ROI.coord{n}(:,1),ROI.coord{n}(:,2),'-r');
                text(mean(ROI.coord{n}(:,1)),mean(ROI.coord{n}(:,2)),num2str(n),'Color','b','fontsize',14,'FontWeight','bold');
            end
        end
    end

end

disp('Data cleared')
disp(' ')

set(handles.roinumber,'Userdata',count) % number of next ROI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function acceptpushbutton_Callback(hObject, eventdata, handles)

ROI=get(handles.acceptpushbutton,'Userdata');
dx=get(handles.clearpushbutton,'Userdata');
count= get(handles.roinumber,'Userdata') ; % number of next ROI
transfer=get(handles.filetransfer,'Value');
roidata=get(handles.text2,'Userdata');
allroidata=get(handles.text5,'Userdata');
allroidata=[allroidata; roidata];
set(handles.text5,'Userdata',allroidata);
backgroundim=get(handles.backgroundpushbutton,'Userdata');
count=count+1;

filename=get(handles.filename,'String');
[namefile,~]=strtok(filename,'.'); %sin extension
savenameIrend=[namefile,'-',num2str(count),'-Irend.mat'];
Irend=ROI.mock;
save(savenameIrend,'Irend','-mat') 
ROI.mock=[];
clear Irend
save([namefile,'.rgn'],'ROI','-mat') 

save([namefile,'-ROIinfo.txt'],'allroidata','-ascii') 

if transfer==1
    ROI2=get(handles.filetransfer,'Userdata');
    filename2=get(handles.filetransfer,'String');
    [namefile2,~]=strtok(filename2,'.'); %sin extension
    savenameIrend=[namefile2,'-',num2str(count),'-Irend.mat'];
    Irend=ROI2.mock;
    save(savenameIrend,'Irend','-mat') 
    ROI2.mock=[];
    clear Irend
    ROI2.mock=[];
    save([namefile2,'.rgn'],'ROI2','-mat') 
    
end

set(handles.roinumber,'Userdata',count); % number of ROIs
set(handles.roinumber,'String',num2str(count)); % number of ROIs

disp(['ROI #',num2str(count),' saved']); 
axes(handles.axes1); 
hold on

if isempty(backgroundim)==0
    text(mean(ROI.coord{count}(:,1)/dx),mean(ROI.coord{count}(:,2)/dx),num2str(count),'Color','g','fontsize',14,'FontWeight','bold');
else
        if handles.smallimage==1
            text(mean(ROI.coord{count}(:,1)*100),mean(ROI.coord{count}(:,2)*100),num2str(count),'Color','b','fontsize',14,'FontWeight','bold');
        else
            text(mean(ROI.coord{count}(:,1)),mean(ROI.coord{count}(:,2)),num2str(count),'Color','b','fontsize',14,'FontWeight','bold');
        end
end

set(handles.clearpushbutton,'Enable','on')

disp(' ')

clear ROI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function quitpushbutton_Callback(hObject, eventdata, handles)

close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function allradiobutton_Callback(hObject, eventdata, handles)

function nroROIclear_Callback(hObject, eventdata, handles)
function nroROIclear_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fullimageradiobutton_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
