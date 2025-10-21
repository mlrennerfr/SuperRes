function createavi(handles,firstframe,lastframe,molnro)
% function createavi(handles,firstframe,lastframe,molnro)
% creates .avi file
% Marianne Renner 09/09 for SPtrack v4
% Marianne Renner 01/2025 - adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% options
prompt = {'File name ','Frames per second','Start at frame:','End at frame:','Molecule number:', 'Line width:','Compression (0:no, 1: yes)'};
num_lines= 1;
dlg_title = '.avi file';
def = {'movieavi','5',num2str(firstframe),num2str(lastframe),molnro,'1.5','1'}; % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
if exit(1) == 0;
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

if handles.param.lastimage>handles.param.lasttrc
    lastframe=handles.param.lastimage;
else
    lastframe=handles.param.lasttrc;
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
x=get(handles.filetrc,'userdata'); %trc data

%double labeling
doubletrc=get(handles.filetrc2,'string');
if isempty(doubletrc)==0
    x2=get(handles.filetrc2,'userdata'); %trc data traj 2
else
    x2=[];
end

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

% avi file
%msgbox(['Saving avi file ',savename])
disp('Creating movie')
set(gca,'xlim',[0 500],'ylim',[0 500],'NextPlot','replace','Visible','off');
if compress ==0
   movi = VideoWriter(savename,'Uncompressed AVI');
   movi.FrameRate = speed;
else
   movi = VideoWriter(savename,'MPEG-4');
   movi.FrameRate = speed;
end

open(movi);
for nroframe=firstframe:lastframechoose
    handles.param.actual=nroframe ;      %actual frame
    showbackgroundimage(handles);
    hold on
    showtrajectories(trc,trc2,handles);
    hold off
   % peli(count)=getframe(gca);  % gets the figure for the movie
   
   frame=getframe(gca);  % gets the figure for the movie
    writeVideo(movi,frame);
    count=count+1;
    value=step*nroframe;
    set(handles.slider1,'value',value);
end
close(movi);
hold off
%close msgbox
disp('Done')

clear trc peli movi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%