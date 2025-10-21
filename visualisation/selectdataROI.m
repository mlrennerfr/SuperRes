function [handles,resp]=selectdataROI(trcdata, trcdata2, nromol, handles)
%function selectdataROI(trcdata, trcdata2, nromol, handles)
% selects data on ROI (image, traces...) to be represented 
% in a zoomed image
% Marianne Renner 09/09 for SPTrack v4
% Marianne Renner 01/2025 - adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datamatrixs= get(handles.file1,'userdata');

if isstruct(datamatrixs)
    datamatrix=datamatrixs.data;
else
    datamatrix=datamatrixs;
end
    
cum=get(handles.cumradiobutton,'value');

resp=[];

%data en ROI
if isempty(trcdata)==0
    axes(handles.axes1);
    if handles.typetraj==5 % Dinst
       minposx=ceil(min(trcdata(:,3)))-1;
       minposy=ceil(min(trcdata(:,4)))-1;
       maxposx=floor(max(trcdata(:,3)))+1;
       maxposy=floor(max(trcdata(:,4)))+1;
       newtrcdata=trcdata;
       xi(1)=minposx;   xi(2)=maxposx;
       yi(1)=minposy;   yi(2)=maxposy;
    else
       [areaselect,xi,yi]=roipolyold;    %seleccion ROI
       minposx=max(ceil(min(xi)),1);
       minposy=max(ceil(min(yi)),1);
       maxposx=min(floor(max(xi)), handles.param.Xdim);
       maxposy=min(floor(max(yi)), handles.param.Ydim);
       [newtrcdata]=pickpointsfast(areaselect,[minposx maxposx],[minposy maxposy],trcdata,handles.param.Xdim,handles.param.Ydim,handles);
    end
    dimx=maxposx-minposx+1;
    dimy=maxposy-minposy+1;

   % correccion coordenadas
    if isempty(newtrcdata)==0
       newtrcdata(:,3)= newtrcdata(:,3)-minposx+1;
       newtrcdata(:,4)= newtrcdata(:,4)-minposy+1;
    end
    if isempty(trcdata2)==0
       set(handles.filetrc,'userdata',trcdata2); %cambio de variables para pickpoints
       [newtrcdata2]=pickpointsfast(areaselect,xi,yi,trcdata,handles.param.Xdim,handles.param.Ydim,handles);
       set(handles.filetrc,'userdata',trcdata); %trc data original
       if isempty(newtrcdata2)==0
          newtrcdata2(:,3)= newtrcdata2(:,3)-minposx+1;
          newtrcdata2(:,4)= newtrcdata2(:,4)-minposy+1;
       end
    else
        newtrcdata2=[];
    end
else % no trc
   [areaselect,xi,yi]=roipolyold;    %seleccion ROI
    newtrcdata=[];
    newtrcdata2=[];
    minposx=max(ceil(min(xi)),1);
    minposy=max(ceil(min(yi)),1);
    maxposx=min(floor(max(xi)), handles.param.Xdim);
    maxposy=min(floor(max(yi)), handles.param.Ydim);
    dimx=maxposx-minposx+1;
    dimy=maxposy-minposy+1;
end
      
newmovie=[];
msgbox('Creating zoomed image...');

%disp(handles.param.lastimage)

if handles.typefile==4 %merge
      newmovie.gray.data=zeros(dimx,dimy);
      newmovie.red.data=zeros(dimx,dimy);
      newmovie.green.data=zeros(dimx,dimy);
      newmovie.blue.data=zeros(dimx,dimy);
      
      for frame=1:handles.param.lastimage %
        if handles.param.gray.nfram>0 
            if handles.param.gray.nfram>1 
                newmovie.gray(frame).data=handles.cgray.image(frame).data(minposy:maxposy,minposx:maxposx);
            else
                newmovie.gray(1).data=handles.cgray.image(1).data(minposy:maxposy,minposx:maxposx);
            end
            dimx=size(newmovie.gray(1).data,2);   dimy=size(newmovie.gray(1).data,1);
        end
        if handles.param.red.nfram>0
            if handles.param.red.nfram>1 
                newmovie.red(frame).data=handles.cred.image(frame).data(minposy:maxposy,minposx:maxposx);
            else
                newmovie.red(1).data=handles.cred.image(1).data(minposy:maxposy,minposx:maxposx);
            end
            dimx=size(newmovie.red(1).data,2);   dimy=size(newmovie.red(1).data,1);
        end
        if handles.param.green.nfram>0
            if handles.param.green.nfram>1 
                newmovie.green(frame).data=handles.cgreen.image(frame).data(minposy:maxposy,minposx:maxposx);
            else
                newmovie.green(1).data=handles.cgreen.image(1).data(minposy:maxposy,minposx:maxposx);
            end
           dimx=size(newmovie.green(1).data,2);   dimy=size(newmovie.green(1).data,1);
        end
        if handles.param.blue.nfram>0
            if handles.param.blue.nfram>1 
                newmovie.blue(frame).data=handles.cblue.image(frame).data(minposy:maxposy,minposx:maxposx);
            else
                newmovie.blue(1).data=handles.cblue.image(1).data(minposy:maxposy,minposx:maxposx);
            end
           dimx=size(newmovie.blue(1).data,2);   dimy=size(newmovie.blue(1).data,1);
        end
      end %loop frame

else
    if isstruct(datamatrixs)
          for frame=1:handles.param.lastimage
              newmovie(frame).data=datamatrixs(frame).data(minposy:maxposy,minposx:maxposx);
              dimx=size(newmovie(1).data,2);
              dimy=size(newmovie(1).data,1);
          end
    else
           newmovie=datamatrix(minposy:maxposy,minposx:maxposx);
           dimx=maxposy-minposy+1;
           dimy=maxposx-minposx+1;
    end
end

close

if isempty(trcdata)==1 
else
    if isempty(trcdata2)==0
        if handles.frame.lasttrc<max(trcdata2(:,2))
            handles.frame.lasttrc=max(trcdata2(:,2));
        end
    end
end

% data to transfer
varargin{1}=newmovie;
varargin{2}=dimx;
varargin{3}=dimy;
varargin{4}=handles.param;
varargin{5}=newtrcdata;
%varargin{6}=handles.param.lastimage; %frames movie o imagen
varargin{6}=cum; %cumulative plot
varargin{7}=handles.typetraj;
varargin{8}=nromol;
varargin{9}=handles.color;
varargin{10}=newtrcdata2;

% ventana ROI
varargout=zoommovtrack(varargin);
uiwait;
aux=['auxiliar.mat'];
%resp=[0;0;0;0;0];
if length(dir(aux))>0
    det=load(aux);
    detopt = struct2cell(det);
    output=detopt{1} ;
    if isempty(output)==0
        [handles,resp]=docorrectshift(output, trcdata, handles);
    else
        resp=[];
    end
end
delete('auxiliar.mat')

clear trcdata datamatrix newtrcdata newmovie

guidata(gcbo,handles) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
