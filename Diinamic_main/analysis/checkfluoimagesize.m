function synimage=checkfluoimagesize(handles, synimage,dx)
%function synimage=checkfluoimagesize(handles, synimage,dx)
%
% resizes fluo image to have the same size than a rendred image
%
% Marianne Renner dec 22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = handles.alpha(:);
alpha(1:10);   

xdim=ceil(max(handles.x)/dx);
ydim=ceil(max(handles.y)/dx);


[I2,xxi,yyi] = PALM_rendering3(handles.x,handles.y,alpha,dx*2,dx,0,xdim, ydim, 1);   
        
if size(I2,1)==size(synimage,1) && size(I2,2)==size(synimage,2)
else
    % need to change image size
    newsynimage = imresize(synimage,[size(I2,1),size(I2,2)]);
    synimage=newsynimage;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
