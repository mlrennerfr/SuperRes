function figurapoints=plotpoints(x,y,fr,handles);
 
figurapoints = figure('Name','Pointillist image','Toolbar','figure');
 
% set boundaries
if isfield(handles,'xlim')
    [xmin,xmax] = deal(handles.xlim(1),handles.xlim(2));
    [ymin,ymax] = deal(handles.ylim(1),handles.ylim(2));
else
    xmin = min(x); xmax = max(x);
    ymin = min(y); ymax = max(y);
end
axis equal;
pointsize = str2num(get(handles.dotsize,'String'));
xlim([xmin xmax]); ylim([ymin ymax]);
 
color = 'k';    
hold on;
 
xlabel('X (mu)');
ylabel('Y (mu)');   
handles.pointsplotted =  plot(x,y,'.','MarkerSize',pointsize,'Color',color);
title([' Total nb of positions =',num2str(length(x))]);
 