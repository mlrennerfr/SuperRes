function fullimageroi(handles)
%function fullimageroi(handles)
%
% Creates a ROI of the full image for the full list of files
%
% Marianne Renner feb 23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

szpx=str2num(get(handles.szpx,'string'));
dx = str2num(get(handles.PALMszpx,'string'));
mu=get(handles.muradiobutton,'Value');

listafiles=get(handles.textlistafiles,'userdata');

for nro=1:size(listafiles,2)
    count=1;
    filename=listafiles{nro}
    [namefile,~]=strtok(filename,'.'); %sin extension
    
    if mu==0
        handles=loadselectPALMdata(filename,szpx,handles);
    else
        handles=loadselectPALMdata(filename,1,handles); %data already in microns
    end
    
    X_mu=handles.x;
    Y_mu=handles.y;
    alpha = handles.alpha(:);
    fr = handles.fr(:);
    alpha(1:10);
    
    xdim=max(X_mu);
    ydim=max(Y_mu);

    x=[0;ceil(xdim);ceil(xdim);0];
    y=[0;0;ceil(ydim);ceil(ydim)];
    coord=[x y];

    %pick points roi
    roiselectedx=[];
    roiselectedy=[];  
    Xcorrected=X_mu;
    Ycorrected=Y_mu;
    
    in = inpolygon(Xcorrected,Ycorrected,x,y);
    roiselectedx=[roiselectedx; X_mu];
    roiselectedy=[roiselectedy; Y_mu];

    
    alpharoi=alpha;
    frroi=fr;
 
    xdim2=ceil(max(roiselectedx));
    ydim2=ceil(max(roiselectedy));

    % mask of the ROI (the same size than the original image)
    xdimrend=xdim2/dx;
    ydimrend=ydim2/dx;
    alpharoi(1:10);
    gausssmooth = 1;
    [imageroirend,rendx, rendy] = PALM_rendering3( roiselectedx,roiselectedy,alpharoi,dx*2,dx,0,xdimrend, ydimrend, gausssmooth);


    newxdim=ceil(X_mu);
    newydim=ceil(Y_mu);
    ROI.fr{count}=frroi;
    ROI.coord{count}=coord;
    ROI.x{count}=x*dx;
    ROI.y{count}=y*dx;

    ROI.in{count}=in; %index of selected points
    ROI.alpha{count}=alpharoi;
    ROI.xselec{count}=roiselectedx; %index of selected points
    ROI.yselec{count}=roiselectedy; %index of selected points
    %ROI.imagerend{count}=imageroirend; %crop of rendered image
    ROI.mock=imageroirend; %crop of rendered image to be saved as an independant file
    ROI.dimx=[1 max(X_mu)];
    ROI.dimy=[1 max(Y_mu)];
    ROI.dx=dx;

    savenameIrend=[namefile,'-',num2str(count),'-Irend.mat'];
    Irend=ROI.mock;
    save(savenameIrend,'Irend','-mat') 
    ROI.mock=[];
    clear Irend
    save([namefile,'.rgn'],'ROI','-mat') 
    disp(' ')
end % loop files

%disp('Done')
