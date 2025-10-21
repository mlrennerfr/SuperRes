function showbackgroundimage(handles)
% function showbackgroundimage(handles)
% display background for .avi creation and shift correction
% Marianne Renner 09/09 for SPTrack v4
% Marianne Renner 01/2025 - adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(handles,'typefile')
else
    set(handles.plottraj,'value',2);
    return
end
typefile=handles.typefile;

XDim=handles.param.Xdim;
YDim=handles.param.Ydim;
Nfram=handles.param.nfram;
actualframe=handles.param.actual; %actual frame
tlag=str2num(get(handles.till,'string'));
name=get(handles.file1,'string');
title=get(handles.showname,'value');
imagedim=[num2str(XDim),' x ',num2str(YDim)];
set(handles.imagedim,'string',imagedim);

if typefile<4                                       % salvo merge: lee matriz
   framematrix=get(handles.file1,'userdata'); %image
 if Nfram>1
   if actualframe==1
      firsty=actualframe;
      lasty=YDim;
   else
      firsty=(actualframe-1)*YDim+1;
      lasty=firsty+YDim-1;
   end
   if typefile==0 
        datamatrix=framematrix(firsty:lasty,:); %spe
   else
        datamatrix=framematrix(actualframe).data; %stk
   end
   datamatrix=double(datamatrix);
 else % one image
     if isstruct(framematrix)
        datamatrix=framematrix.data; % some .tif
     else
        datamatrix=framematrix(:,:); % all the rest
     end
 end
end

%figure
plottraj=get(handles.plottraj,'value');
if plottraj==1 %| plotD==1
   figure;  % for plotting trajectories over one image
else
   axes(handles.axes1);
end
axis([0 XDim 0 YDim]);

if handles.param.lastimage>handles.param.lasttrc
    set(handles.text11,'string',[num2str(actualframe),' (of ',num2str(handles.param.lastimage),')']);
else
    set(handles.text11,'string',[num2str(actualframe),' (of ',num2str(handles.param.lasttrc),')']);
end

if handles.typefile==0
          datamatrix=datamatrix-min(min(datamatrix));
          datamatrix=abs(datamatrix/max(max(datamatrix)));
          imshow(datamatrix,'InitialMagnification','fit');
          hold on
      else
          if handles.typefile== 3
              tiffile=get(handles.file1,'userdata');
              imshow(tiffile,'InitialMagnification','fit');
              hold on;
          elseif handles.typefile== 1 | handles.typefile==2
              stackmin=(min(min(min(datamatrix))));
              stackmax=(max(max(max(datamatrix))));
              imshow((datamatrix(:,:,1)),[stackmin stackmax],'InitialMagnification','fit');
              hold on
         elseif handles.typefile== 4
             datamatrix=showmerge(handles);
             imshow(datamatrix.data,'InitialMagnification','fit');
             hold on
          end
end


resx=handles.param.Xdim/4;
resy=handles.param.Ydim/18;
time=handles.param.actual*tlag;
switch handles.time
       case 1
       text((handles.param.Xdim/20),(handles.param.Ydim-resy),sprintf('Frame : %0.0f',handles.param.actual),'Color',[1 1 1]);
       case 2
       text((handles.param.Xdim/20),(handles.param.Ydim-resy),sprintf('Time : %0.0f',time),'Color',[1 1 1]);
end    
if title==1
     text ((handles.param.Xdim/20), resy, sprintf (name),'Color',[1 1 1]);
end

clear datamatrix

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%