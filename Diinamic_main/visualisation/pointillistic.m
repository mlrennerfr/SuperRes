function pointillistic(handles)
 
X_mu=handles.x-min(handles.x);
Y_mu=handles.y-min(handles.y);
alpha=handles.alpha;
alpha(1:10);

fr=handles.fr;
gausssmooth = 1;

%sigmaloc = str2num(get(handles.widthgaussiansmooth,'String'));
dx = str2num(get(handles.PALMszpx,'String'));
sigmaloc=dx*2;
    
vect2remove=selectremovealpha2(handles,X_mu,Y_mu,alpha); %selection by absolute intensities
if ~isempty(vect2remove)
    [X_mu,Y_mu,alpha,fr]=removealpha(vect2remove,X_mu,Y_mu,alpha,fr);
end
 
% clean data
xmax = max(X_mu);
ymax = max(Y_mu);

% figure
% in nm!!!!
fig_points = figure('Name','Pointillist image','Toolbar','figure');
axis equal;
pointsize = 2;
xlim([0 xmax*1000]); ylim([0 ymax*1000]);
I=ones(ceil(ymax*1000),ceil(xmax*1000)); %!!!!!!!!!!!!!!!!!!
imshow(I,'InitialMagnification','fit');
color = 'k';    
hold on;
areaimage=xmax*ymax; %in µm2
densityimage=size(X_mu,1)/areaimage;

title([' Selected positions =',num2str(length(X_mu)),'. Area image =',num2str(areaimage),...
    ' µm2 (overall density = ',num2str(densityimage),'det/µm2)']);
%xlabel('X (nm)');
%xticks(0:xmax*1000/10:xmax*1000)
%ylabel('Y (nm)');  
%yticks(0:ymax*1000/10:ymax*1000)
handles.pointsplotted =  plot(X_mu*1000,Y_mu*1000,'.','MarkerSize',pointsize,'Color',color);
hold on

figure(fig_points);
handles.fig_points = fig_points;

%rendered
xdim=ceil(max(X_mu)/dx);
ydim=ceil(max(Y_mu)/dx);

%% Create rendered image
[I,xxi,yyi] = PALM_rendering3( X_mu,Y_mu,alpha,dx*2,dx,0,xdim, ydim, gausssmooth);
%aux=I;

%figselected=figure('Name','Select ROIs','Toolbar','figure');
figselected = figure('Name','Rendered image','Toolbar','figure');

hold on;

% rendered in false colors, pixel 0,0 top left
imshow(I,'InitialMagnification','fit')
colormap(gca,'hot')

hold on
 
imwrite(I,'rend.tif')

rendx=X_mu/dx;
rendy=Y_mu/dx;
handles.rendx=rendx;
handles.rendy=rendy;

%set(handles.refreshpushbutton,'userdata',I); %also to recover size of rendered
clear x y X_mu Y_mu alpha fr I

% Update handles structure
guidata(gcbo, handles);
