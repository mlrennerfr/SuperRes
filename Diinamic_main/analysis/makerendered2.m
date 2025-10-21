function [I,rendx, rendy]=makerendered2(handles,szpx,sigmaloc,dx,mu)
%function [I,rendx, rendy]=makerendered2(handles,szpx,sigmaloc,dx,mu)
%creates rendered image I (without plotting) 
%
% input: data.x and data.y : coordinates
%       szpx, dx: pixel size (camera: szpx; super res: dx)
%       sigmaloc: sigma gaussian
%       mu: code for data units (µm or pixels?)
% calls PALMrendering3.m
%
%MR nov 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mu=1;
if mu==1 %already in µm
    xdim=ceil(max(handles.x))/dx;
    ydim=ceil(max(handles.y))/dx;
    X_mu = handles.x(:);
    Y_mu = handles.y(:);
else
    X_mu = handles.x(:)*szpx;
    Y_mu = handles.y(:)*szpx;
    xdim=ceil(max(handles.x))*szpx/dx;
    ydim=ceil(max(handles.y))*szpx/dx;
end


alpha = handles.alpha(:);
fr = handles.fr(:);
alpha(1:10);
gausssmooth = 1;

vect2remove=selectremovealpha2(handles,X_mu,Y_mu,alpha); %selection by absolute intensities

if ~isempty(vect2remove)
    [X_mu,Y_mu,alpha,fr]=removealpha(vect2remove,X_mu,Y_mu,alpha,fr);
end

Nfr = min(max(fr));
iframe = Nfr;
it=iframe;
Nframes1=iframe;

%% Create rendered image
[I,xx1,yy1] = PALM_rendering3( X_mu,Y_mu,alpha,sigmaloc,dx,0,xdim, ydim, gausssmooth);

%figure
%image_ch(xx1,yy1,I,[0 1]); % image_ch(xx1,yy1,I);
%colormap(hot);            % colorbar;

rendx=X_mu/dx;
rendy=Y_mu/dx;

%sizex=size(I,1);
%sizey=size(I,2);

clear X_mu Y_mu alpha fr  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vect2remove=selectremovealpha2(handles,x,y,alpha)

%% Remove points with alphas outside specified intensities range
vect2remove = [];
 
%% Remove points at x=y=0
aux1 = find(x==0);
aux2 = find(y==0);
vect2remove = intersect(aux1,aux2);
if exist('Xmatrix','var')
    matrix2remove = flagmatrixelements(matrix2remove,vect2remove);
end
aux = find(alpha==0);
vect2remove = [vect2remove; aux];
if exist('Xmatrix','var')
    matrix2remove = flagmatrixelements(matrix2remove,vect2remove);
end

alphaminperc = str2num(get(handles.alpha1,'String'));
alphamaxperc = str2num(get(handles.alpha2,'String'));

valmaxalpha=max(handles.alpha);
valminalpha=min(handles.alpha);
difalpha=valmaxalpha-valminalpha;

alphamin=valminalpha+(alphaminperc/100*difalpha);
alphamax=valminalpha+(alphamaxperc/100*difalpha);

aux = union(find(alpha<alphamin),find(alpha>alphamax));

if size(aux,2)>1
    vect2remove = [vect2remove; aux(1)];
else
    vect2remove = [vect2remove; aux];
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y,alpha,fr]=removealpha(vect2remove,x,y,alpha,fr)
%% Removing unneeded points from column vectors
 
    x(vect2remove) = [];
    y(vect2remove) = [];
    if nargin>3
        alpha(vect2remove) = [];
        fr(vect2remove) = [];
    else
        alpha=[];
        fr=[];
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I = image_ch(varargin)
xx1 = varargin{1};
yy1 = varargin{2};
I = varargin{3};
pct_sat = 0;
if ndims(I)==2
    if nargin>3
        clims = varargin{4};
    else
        Imax = prctile(I(:),100-pct_sat);
        Imin = min(I(:));
        clims = [Imin Imax];
    end
    if any(isnan(clims)) || clims(2)<=clims(1)
        warning(['clims (=',num2str(clims),')! Using [0,1] instead !']);
        clims = [0 1];
    end
    I = imagesc(xx1,yy1,I,clims);
else
    image(xx1,yy1,I);
end
colorbar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
