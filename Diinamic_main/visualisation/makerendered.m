function [sizex,sizey]=makerendered(data,szpx,sigmaloc,dx,xdim,ydim)


%X_mu = data.x(:)*szpx; 
%Y_mu = data.y(:)*szpx;
%alpha = data.alpha(:);
%fr = data.fr(:);
alpha(1:10);
gausssmooth = 1;
%xdim=ceil(max(data.x))*szpx/dx;
%ydim=ceil(max(data.y))*szpx/dx;

Nfr = min(max(fr));
iframe = Nfr;
it=iframe;
Nframes1=iframe;

%% Create rendered image
maxx=ceil(max(X_mu/dx));
maxy=ceil(max(Y_mu/dx));

[I,xx1,yy1] = PALM_rendering3( X_mu,Y_mu,alpha,sigmaloc,dx,0,xdim, ydim, gausssmooth);

sizex=size(I,1);
sizey=size(I,2);

clear X_mu Y_mu alpha fr  I


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 