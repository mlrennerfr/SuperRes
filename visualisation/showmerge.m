function mergedimage=showmerge(handles,frame)
% function mergedimage=showmerge(handles,factor,frame)
% prepares the merged image for showbackground
% Marianne Renner 09/09 for SPTrack v4
% Marianne Renner 05/10 for SPTrack v4
% Marianne Renner 01/2025 - adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off MATLAB:divideByZero
hold off;

if nargin>1
    actualframe=frame;
else
   actualframe=handles.param.actual; %actual frame
end

finalimage=[];
K1=[];

if isfield(handles,'movie') %zoom
    handles.cgray.image=handles.movie.gray;
    handles.cred.image=handles.movie.red;
    handles.cgreen.image=handles.movie.green;
    handles.cblue.image=handles.movie.blue;
    imagered=zeros(handles.Ydim,handles.Xdim);
else
    imagered=zeros(handles.param.Ydim,handles.param.Xdim);
%    set(handles.avibutton,'enable','on');  
   set(handles.plottraj,'enable','on'); 
   set(handles.saveimage,'enable','on');  
end

imagegreen=imagered;
imageblue=imagered;
imagegray=imagered;

if handles.param.gray.nfram==0 
else
    if handles.param.gray.nfram==1 
        imagegray=greyscaleimage(handles.cgray.image(1).data,handles.param.gray.factor,handles.param.gray.factorhigh);
    else 
        imagegray=greyscaleimage(double(handles.cgray.image(actualframe).data),handles.param.gray.factor,handles.param.gray.factorhigh); 
    end
end

if handles.param.red.nfram==0 
else
    if handles.param.red.nfram==1 
        imagered=greyscaleimage(handles.cred.image(1).data,handles.param.red.factor,handles.param.red.factorhigh);
    else 
        imagered=greyscaleimage(double(handles.cred.image(actualframe).data),handles.param.red.factor,handles.param.red.factorhigh); 
    end
end
if handles.param.green.nfram==0 
else
    if handles.param.green.nfram==1 
        imagegreen=greyscaleimage(handles.cgreen.image(1).data,handles.param.green.factor,handles.param.green.factorhigh);
    else 
        imagegreen=greyscaleimage(double(handles.cgreen.image(actualframe).data),handles.param.green.factor,handles.param.green.factorhigh); 
    end
end
if handles.param.blue.nfram==0 
else
    if handles.param.blue.nfram==1 
        imageblue=greyscaleimage(handles.cblue.image(1).data,handles.param.blue.factor,handles.param.blue.factorhigh);
    else 
        imageblue=greyscaleimage(double(handles.cblue.image(actualframe).data),handles.param.blue.factor,handles.param.blue.factorhigh); 
    end
end

if size(imagered,1)==size(imagegreen,1) & size(imagered,1)==size(imageblue,1)
    if size(imagered,2)==size(imagegreen,2) & size(imagered,2)==size(imageblue,2)
        finalimage.data=cat(3,imagered,imagegreen,imageblue) ;
        control=1;
    else
        msgbox('Images must have the same size!','error','error')
        control=0;
    end
else
    msgbox('Images must have the same size!','error','error')
    control=0;
end

if control==1
    if handles.param.gray.nfram>0
        mincolor=min(min(min(finalimage.data)));
        maxcolor=max(max(max(finalimage.data)));
        if isnan(mincolor)==0 & isnan(maxcolor)==0 ;
            dicimage=cat(3,imagegray,imagegray,imagegray) ;
            finalimage.data=imadd(finalimage.data,dicimage);
        else
            dicimage=cat(3,imagegray,imagegray,imagegray) ;
            finalimage.data=dicimage;
        end
    end
    mergedimage=finalimage;
    if isfield(handles,'file1')
        set(handles.file1,'userdata',mergedimage);
    end
    handles.typefile=4;
end

clear imagered imagegreen imageblue imagegray

guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%