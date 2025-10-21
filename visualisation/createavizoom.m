function createavizoom(handles,firstframe,lastframe,molnro)
% function createavizoom(handles,firstframe,lastframe,molnro)
% Creates avi file from zoomed image
% Marianne Renner 09/09 for SPTrack v4
% Marianne Renner 01/2025 - adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% options
prompt = {'File name ','Frames per second','Start at frame:','End at frame:','Molecule number:', 'Line width:','Compression (0:no, 1: yes)'};
num_lines= 1;
dlg_title = 'Movie file';
def = {'movie','5',num2str(firstframe),num2str(lastframe),molnro,'1.5','1'}; % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
if exit(1) == 0
   return; 
end
savename=answer{1};
speed=str2num(answer{2});
firstframe=str2num(answer{3});
lastframechoose=str2num(answer{4});
molnro=answer{5};
linew=str2num(answer{6});
compress=str2num(answer{7});

if firstframe>lastframe
    msgbox('Wrong values','error','error')
    return
end

if lastframe>lastframechoose
else
    lastframechoose=lastframe;
end

%sliding bar
slidervalue=get(handles.slider1,'Value');
minvalue=get(handles.slider1,'Min');
maxvalue=get(handles.slider1,'Max');
prop=slidervalue/(minvalue+maxvalue);
step=(maxvalue-minvalue)/lastframe;
set(handles.slider1,'SliderStep',[step step]);

count=1;
if size(handles.newtraces,1)==0
    x=handles.traces;
else
    x=handles.newtraces; %shift correction
end
x2=handles.traces2;

if isempty(x)==0
    if firstframe>1
        % correct trajectories
        index=find(x(:,2)>firstframe-1);
        trc=x(index,:);
        if isempty(x2)==0
            index2=find(x2(:,2)>firstframe-1);
            trc2=x2(index,:);
        else
            trc2=[];
        end
    else
        trc=x;
        trc2=x2;
    end
else
    trc=x;
    trc2=x2;
end

  %  showframeROIavi(handles);
 %   hold on

% movie
if firstframe==1
    firstframe=2;
end
if compress ==0
   %movi = avifile(savename,'compression','none','fps',speed,'quality',100)
   movi = VideoWriter(savename,'Uncompressed AVI')
   movi.FrameRate = speed;
else
   %movi = avifile(savename,'compression','Cinepak','fps',speed,'quality',100)
   movi = VideoWriter(savename,'MPEG-4')
   movi.FrameRate = speed;
end
open(movi);
for nroframe=firstframe:lastframechoose
    handles.param.actual=nroframe ;      %actual frame
    set(handles.nroframe,'value',nroframe);
    showframeROIavi(handles);
    hold on
    showtrajavizoom(trc,trc2,handles);
    hold off
    %showtrajavizoom2(trc,trc2,handles);
    %hold on
    %peli(count)=getframe(gca);  % gets the figure for the movie
   % axes(handles.axes1);
   % frame=getframe(gca);  % gets the figure for the movie
    frame=getframe(gcf);  % gets the figure for the movie
    
    writeVideo(movi,frame);

    count=count+1;
    set(handles.framenumber,'string',[num2str(nroframe),' of ',num2str(handles.param.maxfram)]);
    value=step*nroframe;
    set(handles.slider1,'value',value);
    
end
close(movi);
hold off

% avi file
%msgbox(['Saving avi file ',savename])
%set(gca,'xlim',[0 500],'ylim',[0 500],'NextPlot','replace','Visible','off');
%if compress ==0
   %movi = avifile(savename,'compression','none','fps',speed,'quality',100)
%   movi = VideoWriter(savename,'Uncompressed AVI','FrameRate',speed,'quality',100)
%else
   %movi = avifile(savename,'compression','Cinepak','fps',speed,'quality',100)
%   movi = VideoWriter(savename,'FrameRate',speed,'quality',100 )
%end

%movi = addframe(movi,peli);
%movi = close(movi);
%close %msgbox

clear trc movi peli


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%