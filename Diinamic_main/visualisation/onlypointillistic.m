function onlypointillistic(handles)
 
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
fig_points = figure('Name','Pointillist image','Toolbar','figure');
axis equal;
pointsize = 2;
xlim([0 xmax]); ylim([0 ymax]);
I=ones(ceil(ymax),ceil(xmax)); %!!!!!!!!!!!!!!!!!!
imshow(I,'InitialMagnification','fit');
color = 'k';    
hold on;
title([' Selected positions =',num2str(length(X_mu))]);
handles.pointsplotted =  plot(X_mu,Y_mu,'.','MarkerSize',pointsize,'Color',color);
hold on

figure(fig_points);
handles.fig_points = fig_points;

% Update handles structure
guidata(gcbo, handles);
