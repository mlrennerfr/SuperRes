function newimage=shiftimagePALM(handles,ktev,kteh,image,maxval,code)
% function newimage=shiftimage(handles,ktev,kteh,image,maxval,code)
% shift correction for images
%
% Marianne Renner - 07/10 for movtrack.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newimage=zeros(handles.Ydim,handles.Xdim); %att background

if isempty(image)==0
        if ktev>0 % down 
            if kteh>0 %right
                newimage(1+ktev:handles.Ydim,1+kteh:handles.Xdim)=image(1:handles.Ydim-ktev,1:handles.Xdim-kteh);
            elseif kteh<0 % left
                kteh=-kteh;
                newimage(1+ktev:handles.Ydim,1:handles.Xdim-kteh)=image(1:handles.Ydim-ktev,1+kteh:handles.Xdim);
            elseif kteh==0
                newimage(1+ktev:handles.Ydim,:)=image(1:handles.Ydim-ktev,:);
            end
        elseif ktev<0 % up 
            ktev=-ktev;
            if kteh>0 %right
                newimage(1:handles.Ydim-ktev,1+kteh:handles.Xdim)=image(1+ktev:handles.Ydim,1:handles.Xdim-kteh);
            elseif kteh<0 % left
                kteh=-kteh;
                newimage(1:handles.Ydim-ktev,1:handles.Xdim-kteh)=image(1+ktev:handles.Ydim,1+kteh:handles.Xdim);
            elseif kteh==0
                newimage(1:handles.Ydim-ktev,:)=image(1+ktev:handles.Ydim,:);
            end
        elseif ktev==0
            if kteh>0 %right
                newimage(:,1+kteh:handles.Xdim)=image(:,1:handles.Xdim-kteh);
            elseif kteh<0 % left
                kteh=-kteh;
                newimage(:,1:handles.Xdim-kteh)=image(:,1+kteh:handles.Xdim);
            elseif kteh==0
                newimage=image;
            end
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
