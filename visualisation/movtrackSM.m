function varargout = movtrackSM(varargin)
%MOVTRACKSM M-file for movtrackSM.fig
% Last Modified by GUIDE v2.5 18-Feb-2019 15:21:30
%
% GUI for visualizing image and trajecotires
%
% Marianne Renner 01/2025 - adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @movtrackSM_OpeningFcn, ...
                   'gui_OutputFcn',  @movtrackSM_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before movtrackSM is made visible.
function movtrackSM_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
% initialize handles 
handles.color.all='r';    %color code trajectories
handles.color.extra='b';
handles.color.peri='y';
handles.color.syn='r';

handles.traject=0;
handles.time=1;
handles.alltraj=0;
handles.trcbutton='1';
handles.trajbutton='0';
set(handles.filetrc,'userdata',[]);
set(handles.filetrc,'value',0);  
set(handles.trcradiobutton,'userdata',1); % perival=1: peri=syn

handles.param.first=1;%first frame
handles.param.actual=1;%actual frame
handles.param.lastimage=1;%last frame
handles.param.lasttrc=1;%last frame with trc

set(handles.plottraj,'value',0);
handles.presentfolder=cd;
handles.presenttrcfolder=cd;
handles.perival=1;
handles.path=cd;
handles.background=' ';
handles.name=' ';

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = movtrackSM_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function linewidth_Callback(hObject, eventdata, handles)
function corrhor_Callback(hObject, eventdata, handles)
function corrvert_Callback(hObject, eventdata, handles)
function delay_Callback(hObject, eventdata, handles)
function molindiv_Callback(hObject, eventdata, handles)
function till_Callback(hObject, eventdata, handles)
function showname_Callback(hObject, eventdata, handles)
function static_Callback(hObject, eventdata, handles)
function fgraylow_Callback(hObject, eventdata, handles)
function fredlow_Callback(hObject, eventdata, handles)
function fgreenlow_Callback(hObject, eventdata, handles)
function fbluelow_Callback(hObject, eventdata, handles)
function fgrayhigh_Callback(hObject, eventdata, handles)
function fredhigh_Callback(hObject, eventdata, handles)
function fgreenhigh_Callback(hObject, eventdata, handles)
function fbluehigh_Callback(hObject, eventdata, handles)

function popupmenu1_Callback(hObject, eventdata, handles)
    handles.time= get(hObject,'Value');
    guidata(hObject, handles);

function trcradiobutton_Callback(hObject, eventdata, handles)
    set(handles.pkradiobutton,'value',0);
    set(handles.accel,'value',0);
    guidata(hObject, handles);

function pkradiobutton_Callback(hObject, eventdata, handles)
    handles.pkbutton=get(hObject,'value');
    set(handles.trcradiobutton,'value',0);
    set(handles.accel,'value',1);
    guidata(hObject, handles);
    
function accel_Callback(hObject, eventdata, handles)
    set(handles.pkradiobutton,'value',1);
    guidata(hObject, handles);    
    
function cumradiobutton_Callback(hObject, eventdata, handles)
function ident_Callback(hObject, eventdata, handles)
    
function loc_Callback(hObject, eventdata, handles)
    set(handles.trajtimecolorradiobutton,'value',0);
    set(handles.rainbowradiobutton,'value',0);
    set(handles.timecolorradiobutton,'value',0);
    guidata(hObject, handles);

function rainbowradiobutton_Callback(hObject, eventdata, handles)
    handles.rainbow=get(hObject,'value');
    set(handles.trajtimecolorradiobutton,'value',0);
    set(handles.loc,'value',0);
    set(handles.timecolorradiobutton,'value',0);
    [colorm]=createcolor(50);
    for i=1:size(colorm,1)
          indcol=1;
          indcol=round(rand(1)*size(colorm,1));
          colorm(i,4)=indcol;
    end
    colorm=sortrows(colorm,4);
    handles.colorm=colorm(:,1:3);
    guidata(hObject, handles);

function timecolorradiobutton_Callback(hObject, eventdata, handles)
    handles.timecolor=get(hObject,'value');
    set(handles.trajtimecolorradiobutton,'value',0);
    set(handles.rainbowradiobutton,'value',0);
    set(handles.loc,'value',0);
    guidata(hObject, handles);
    
function trajtimecolorradiobutton_Callback(hObject, eventdata, handles)
    handles.timecolor=get(hObject,'value');
    set(handles.timecolorradiobutton,'value',0);
    set(handles.rainbowradiobutton,'value',0);
    set(handles.loc,'value',0);
    guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)

slidervalue=get(hObject,'Value');
minvalue=get(hObject,'Min');
maxvalue=get(hObject,'Max');
prop=slidervalue/(minvalue+maxvalue);
if handles.param.lastimage>handles.param.lasttrc
    nfram=handles.param.lastimage;
else
    nfram=handles.param.lasttrc;
end
step=(maxvalue-minvalue)/nfram;
set(handles.slider1,'SliderStep',[step step]);
handles.param.actual=round(prop*nfram);
if handles.param.actual==0
    handles.param.actual=1;
end
if handles.param.actual>nfram
    handles.param.actual=nfram;
end
axes(handles.axes1);
set(handles.text11,'string',[num2str(handles.param.actual),' (of ',num2str(nfram),')']);

showbackground(handles)
hold on
showtraj(handles)
hold off;

guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sets trajectories color
function trajcolor_Callback(hObject, eventdata, handles)

% dialog box to enter codes for trajectories color
prompt = {'Without localization ','Localization out of domains (extra: zero values)','Localization peri-domains (negative values)',...
           'Localization inside domains (positive values)'};
num_lines= 1;
dlg_title = 'Enter codes for colors';
def = {handles.color.all ,handles.color.extra,handles.color.peri,...
    handles.color.syn}; % default values
msgbox('r: red, b: blue, g: green, y:yellow, w: white, m:magenta, c: cyan, k: black','Color codes','help')
answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
if exit(1) == 0;
   return; 
end
handles.color.all= answer{1};    %color code trajectories
handles.color.extra=answer{2}; 
handles.color.peri=answer{3}; 
handles.color.syn=answer{4}; 
close %msgbox

guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% one image or movie
% --- Executes on button press in images.
function images_Callback(hObject, eventdata, handles)

% loads image 
currentcd=cd;
if length(dir(handles.presentfolder))>0
   cd(handles.presentfolder)
end

[file,path] = uigetfile('*.*','Load background file (.stk or .tif)');
filename = [path,file];
if path>0
      handles.presentfolder=path;
end
cd(currentcd)
if filename==0
    return
end
set(handles.file1,'string',file);

% movie
stktrue=0;
msgbox('Reading file...');

answer2=findstr(filename,'.stk');
if isempty(answer2)==1
      answer3=findstr(filename,'.tif');
      if isempty(answer3)==1
         answer4=findstr(filename,'.mat');
         if isempty(answer4)==1
            msgbox('Wrong type of file','Error','error');
            return
         end
      else
         % .tif file       
         info=imfinfo(filename);
         if size(info,2)>1 % movie tif
            stktrue=1;
            [stack_info,datamatrix] = stkdataread(filename);
         else
             [stack_info,datamatrix] = tifdataread(filename);
         end
         
         handles.param.nfram=stack_info.frames;
         handles.param.Xdim=stack_info.x;
         handles.param.Ydim=stack_info.y;
         stktrue=2;
         [fil,col]=size(datamatrix);               
         if col/handles.param.Xdim==3  %rgb
            stktrue=3;
         end
         set(handles.file1,'userdata',datamatrix); %images
      end
   else
      % .stk file       
      [stack_info,stackdata] = stkdataread(filename);
      handles.param.Xdim=stack_info.x;
      handles.param.Ydim=stack_info.y;
      handles.param.nfram=stack_info.frames;
      stktrue=1;
      set(handles.file1,'userdata',stackdata); %images
end

handles.typefile=stktrue;

% positions
handles.param.first=1; %first frame
handles.param.actual=1; %actual frame
handles.param.lastimage=handles.param.nfram; %last frame

if handles.param.lastimage>handles.param.lasttrc
   handles.param.maxfram=handles.param.lastimage;
   set(handles.text11,'string',['1 (of ',num2str(handles.param.lastimage),')']);
else
   handles.param.maxfram=handles.param.lasttrc;
end

close %msgbox

%first image
axes(handles.axes1);
showbackground(handles);
hold off;

%buttons
set(handles.avibutton,'enable','on');  
set(handles.zoompushbutton,'enable','on');  
set(handles.plot3Dpushbutton,'enable','on');  
set(handles.plottraj,'enable','on'); 
set(handles.saveimage,'enable','on');  

%slider
minvalue=get(handles.slider1,'Min');
set(handles.slider1,'value',minvalue);
clear datamatrix

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trajectories file
% --- Executes on button press in traces.
function traces_Callback(hObject, eventdata, handles)

pkfile=get(handles.pkradiobutton,'value');
tipotrc=get(handles.trcradiobutton,'value'); 
if pkfile==1
    tipotrc=2;
end

% positions
if length(dir(handles.presentfolder))>0
   cd(handles.presentfolder)
end
currentcd=cd;

% loads traces file1 
if tipotrc==1
   [trcf,tpath] = uigetfile('*.trc','Load trajectories file'); 
elseif tipotrc==2 %pk
   [trcf,tpath] = uigetfile('*.pk*','Load peaks file'); 
end

if tpath>0
   handles.presenttrcfolder=tpath;
   trcfile = [tpath,trcf]
else
   trcfile = [];
end

if isempty(trcfile)==1
    set(handles.filetrc,'userdata',[]);
    set(handles.filetrc,'value',1);  
    set(handles.filetrc,'string',['']);
    nrotraj=1;
    handles.param.lasttrc=nrotraj; %last frame with trc
    guidata(gcbo,handles) ;
    return
else
    set(handles.filetrc,'string',trcfile);
    handles.traject=get(handles.filetrc,'string');
    set(handles.filetrc,'string',trcf);
    name=get(handles.file1,'string');
    if isempty(name)==0
        set(handles.plottraj,'enable','on'); 
    end

    x=[];
    nrotraj=1;

    if tipotrc==1
        x=load(trcfile);   %.trc
        nrotraj=max(x(:,2));
    elseif tipotrc==2
       pk=load(trcfile);   %.pk
       nrotraj=max(pk(:,1));
       if size(pk,2)>15
            x=[pk(:,1) pk(:,1:3) pk(:,17)]; % add z
       else
          % if size(pk,2)==11
                x=[pk(:,1) pk(:,1:3) pk(:,8)]; %  z: col 8
          % else
         %       x=[pk(:,1) pk(:,1:3)]; % no z
          % end
       end
       accelerate=get(handles.accel,'value');
       if accelerate==1
          fin=0; i=1;
          newpk=[];
          while fin==0
                newpk=[newpk; x(i,:)];
                i=i+10;
                if size(x,1)<i
                   break
                end
          end
          set(handles.accel,'userdata',newpk)
          clear newpk
       end
 
    end

    set(handles.filetrc,'userdata',x);
    set(handles.filetrc,'value',nrotraj)  

    handles.param.lasttrc=nrotraj; %last frame with trc
    if handles.param.lasttrc>handles.param.lastimage
      handles.param.maxfram=handles.param.lasttrc;
      set(handles.text11,'string',['1 (of ',num2str(handles.param.lasttrc),')']);
    else
      set(handles.text11,'string',['1 (of ',num2str(handles.param.lastimage),')']);
      handles.param.maxfram=handles.param.lastimage;
    end
    
end

set(handles.plottraj,'enable','on'); 
set(handles.plot3Dpushbutton,'enable','on'); 

%slider
minvalue=get(handles.slider1,'Min');
set(handles.slider1,'value',minvalue);

clear x
guidata(gcbo,handles) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% files to merge
% --- Executes on button press in gray.
function gray_Callback(hObject, eventdata, handles)

% loads image 
[datamatrix,handles.param.gray,handles.cgray.typefile,file,handles.presentfolder]=readimagemerge('DIC', handles);
if isempty(file)==1
    set(handles.grayname,'string',[]); 
    handles.cgray.typefile=0;
else
    set(handles.grayname,'string',file);
    set(handles.grayname,'userdata',datamatrix); %images
end
clear datamatrix file path
guidata(gcbo,handles) ;

%----------------------------------------------------------
% --- Executes on button press in red.
function red_Callback(hObject, eventdata, handles)

% loads image 
clear handles.moviered

[datamatrix,handles.param.red,handles.cred.typefile,file,handles.presentfolder]=readimagemerge('RED', handles);
if isempty(file)==1
    set(handles.redname,'string',[]); 
    handles.cred.typefile=0;
else
    set(handles.redname,'string',file);
    set(handles.redname,'userdata',datamatrix); %images
end
clear datamatrix file path
guidata(gcbo,handles) ;

%------------------------------------------------------------
% --- Executes on button press in green.
function green_Callback(hObject, eventdata, handles)

% loads image 
[datamatrix,handles.param.green,handles.cgreen.typefile,file,handles.presentfolder]=readimagemerge('GREEN', handles);
if isempty(file)==1
    set(handles.greenname,'string',[]); 
    handles.cgreen.typefile=0;
else
    set(handles.greenname,'string',file);
    set(handles.greenname,'userdata',datamatrix); %images
end
clear datamatrix file path
guidata(hObject,handles) ;

%-----------------------------------------------------------------
% --- Executes on button press in blue.
function blue_Callback(hObject, eventdata, handles)

% loads image 
[datamatrix,handles.param.blue,handles.cblue.typefile,file,handles.presentfolder]=readimagemerge('BLUE', handles);
if isempty(file)==1
    set(handles.bluename,'string',[]); 
    handles.cblue.typefile=0;
else
    set(handles.bluename,'string',file);
    set(handles.bluename,'userdata',datamatrix); %images
end
clear datamatrix file path
guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function corrtrajpushbutton_Callback(hObject, eventdata, handles)

corrhor=str2num(get(handles.corrhor,'string'));
corrvert=str2num(get(handles.corrvert,'string'));
trc=get(handles.filetrc,'userdata'); %trc data
trc(:,3)=trc(:,3)+corrhor;
trc(:,4)=trc(:,4)+corrvert;

pkfile=get(handles.pkradiobutton,'value');
accelerate=get(handles.accel,'value');
if pkfile==1 & accelerate==1
   x=get(handles.accel,'userdata');
   x(:,3)=x(:,3)+corrhor;
   x(:,4)=x(:,4)+corrvert;
   set(handles.accel,'userdata',x);
end

set(handles.filetrc,'userdata',trc); %corr trc data

guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% background image
function showbackground(handles)
% shows stk or tif images (one or more)

if isfield(handles,'typefile')
else % no background image
    x=get(handles.filetrc,'userdata'); %trc data
    if isempty(x)==0
       % black background
       set(handles.file1,'userdata',ones(ceil(max(x(:,4))),ceil(max(x(:,3)))))
       handles.typefile=1;
       handles.param.Xdim=ceil(max(x(:,3)));
       handles.param.Ydim=ceil(max(x(:,4)));
       handles.param.nfram=1;
    else
       set(handles.plottraj,'value',2);
       return
    end
    clear x
end

tlag=str2num(get(handles.till,'string'));
name=get(handles.file1,'string');
title=get(handles.showname,'value');

imagedim=[num2str(handles.param.Xdim),' x ',num2str(handles.param.Ydim)];
set(handles.imagedim,'string',imagedim);

if handles.typefile<4                                       % salvo merge: lee matriz
   framematrix=get(handles.file1,'userdata'); %image
 if handles.param.nfram>1
   if handles.param.actual==1
      firsty=handles.param.actual;
      lasty=handles.param.Ydim;
   else
      firsty=(handles.param.actual-1)*handles.param.Ydim+1;
      lasty=firsty+handles.param.Ydim-1;
   end
   if handles.typefile==0 
        datamatrix=framematrix(firsty:lasty,:); %spe
   else
        datamatrix=framematrix(handles.param.actual).data; %stk
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
axis([0 handles.param.Xdim 0 handles.param.Ydim]);
if handles.param.lastimage>handles.param.lasttrc
    set(handles.text11,'string',[num2str(handles.param.actual),' (of ',num2str(handles.param.lastimage),')']);
else
    set(handles.text11,'string',[num2str(handles.param.actual),' (of ',num2str(handles.param.lasttrc),')']);
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

set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'xtick',[])
set(gca,'ytick',[])

clear datamatrix

guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trajectories
function showtraj(handles)
% for each frame, makes an array with traces of the molecules and plots them

pkfile=get(handles.pkradiobutton,'value');
%doubletrc=get(handles.filetrc2,'string'); %second traj
identify=get(handles.ident,'value');
cum=get(handles.cumradiobutton,'value');
localiz=get(handles.loc,'value');
timecolor=get(handles.timecolorradiobutton,'value');
trajtimecolor=get(handles.trajtimecolorradiobutton,'value');
rainbow=get(handles.rainbowradiobutton,'value');

linew=str2num(get(handles.linewidth,'string'));
molnro=get(handles.molindiv,'string');

x=get(handles.filetrc,'userdata'); %trc data

if isempty(x)==0 % old format: zeros when blinks
    indexnozero=find(x(:,3)>0);
    if isempty(indexnozero)==0
        x=x(indexnozero,:);
    end
end

handles.typetraj=0;
if isfield(handles,'param')
else
    % no background image
    handles.param.Xdim=ceil(max(x(:,3)));
    handles.param.Ydim=ceil(max(x(:,4)));
end

%double labeling
%x2=[];
%if isempty(doubletrc)==0
%    x2=get(handles.filetrc2,'userdata'); %trc data traj 2
%end


if isempty(x)==0  % if the trc file is loaded
  maxmol=max(x(:,1));
  axes(handles.axes1);
  axis([0 handles.param.Xdim 0 handles.param.Ydim]);
  tlag=str2num(get(handles.till,'string'));
  name=get(handles.file1,'string');

  if pkfile==0; %.trc
    codenormal=handles.color.all;

    indexmol=[];
    indexactual=[];
    actualtrc=[];
    colorm=[];
    if timecolor==1 
      [colorm]=createcolor(max(x(:,2)));
    elseif rainbow==1
      colorm=handles.colorm;
    else
      codecol=handles.color.all;
    end
    indcol=1;
    
    k=strfind(molnro,'all');
    if isempty(k)==1
       indexmol=find(x(:,1)==str2num(molnro));
       trc=x(indexmol,:);
    else
       trc=x; 
    end
    
    indexactual=find(trc(:,2)<handles.param.actual+1);
    
    if isempty(indexactual)==0
        actualtrc=trc(indexactual,:);
        for i=1:max(actualtrc(:,1))
            %disp(i)
            indexmol=find(actualtrc(:,1)==i);
            indextrc=find(trc(:,1)==i);
            if cum==1 %cummulative
              condition=1;
            else
                if max(trc(indextrc,2))>handles.param.actual %present before and after
                    condition=1;
                else
                    condition=0;
                end
            end
            if isempty(indexmol)==0 && condition==1 % 
                codecol=handles.color.all;
                %%%%%%%%%%%%%
                if rainbow==1
                    handles.typetraj=3;
                    codecol=colorm(indcol,1:3);
                    indcol=indcol+1;
                    if indcol>50
                        indcol=1;
                    end
                end
                
                auxtrc=actualtrc(indexmol,:);
                indextime=find(auxtrc(:,2)==handles.param.actual);
                if isempty(indextime)==0  
                   % if auxtrc(indextime(1),7)==1 % group 1
                   %     codecol=handles.color.extra; % blue
                   % else
                        codecol=handles.color.all; %red
                   % end
                    
                   % plot(auxtrc(indextime(1),3),auxtrc(indextime(1),4),'Color',codecol,'Marker','.','MarkerSize',20);
                   % hold on
                end
                clear auxtrc
            
            else
                %%%%%%
                if localiz==0 && timecolor==0 && trajtimecolor==0 
                   if rainbow==1
                      handles.typetraj=3;
                      codecol=colorm(indcol,1:3);
                      indcol=indcol+1;
                      if indcol>50
                         indcol=1;
                      end
                   end
                   % normal or rainbow
                   plot(actualtrc(indexmol,3),actualtrc(indexmol,4),'Color',codecol,'LineWidth',linew);
                   hold on
                elseif localiz==0 && timecolor==1  
                   handles.typetraj=4;
                   for g=1:size(indexmol,1)-1
                       codecol=colorm(g,1:3);
                       plot(actualtrc(indexmol(g):indexmol(g+1),3),actualtrc(indexmol(g):indexmol(g+1),4),'Color',codecol,'LineWidth',linew);
                       hold on
                   end
                
                elseif localiz==0 && trajtimecolor==1  
                   
                   
                   
                elseif localiz==1 
                   handles.typetraj=1;
                   if size(trc,2)>5  % trc with localization
                       
                      for f=1:size(indexmol,1)-2
                          if size(trc,2)<7
                                if actualtrc(indexmol(f+1),6)<0 %peri
                                    codecol=handles.color.peri;
                                elseif actualtrc(indexmol(f+1),6)>0 %syn
                                    codecol=handles.color.syn;
                                elseif actualtrc(indexmol(f+1),6)==0 %extra
                                    codecol=handles.color.extra;
                                end 
                          end 
                          plot(actualtrc(indexmol(f):indexmol(f+1),3),actualtrc(indexmol(f):indexmol(f+1),4),'Color',codecol,'LineWidth',linew);
                          hold on
                          
                      end %loop f
                   else
                        plot(actualtrc(indexmol,3),actualtrc(indexmol,4),'Color',codecol,'LineWidth',linew); % no syn
                        hold on
                   end
                end  % options colorcode   
                
                %%%%
                
                
                if identify==1 % numero tray
                  pos=indexmol(size(indexmol,1));
                  text(actualtrc(pos,3)+1,actualtrc(pos,4)+1,sprintf('%0.0f',actualtrc(pos,1)),'Color',[1 1 0]);
                  cifras=num2str(actualtrc(pos,1)); space=size(cifras,2);
                  text(actualtrc(pos,3)+(space*5),actualtrc(pos,4)+1,sprintf('(%0.0f)',size(actualtrc(indexmol),1)),'Color',[1 1 1],'FontSize',7);
               end
               hold on
            %else
            %   plot(handles.param.Xdim,handles.param.Ydim,'.k');  %just no avoid crash!
             %  hold on    
            end % empty indexmol
            
      end % loop actual molecules
    else
      plot(handles.param.Xdim,handles.param.Ydim,'.k');  %just to avoid crash!
      hold on  
    end %empty indexactual

 %   if isempty(x2)==0 % second traj
 %      plottracesframe(x2,handles.param.actual,'g',localiz,linew,handles)
 %   end

  else %pk
    
    accelerate=get(handles.accel,'value');
    if accelerate==1
       x=get(handles.accel,'userdata');
    end

    % frame per frame
    index=find(x(:,2)<handles.param.actual+1);
    if isempty(index)==0
       for i=1:size(index,1)
           if size(x,2)>4
                plot3(x(index(i),3),x(index(i),4),x(index(i),5),'Marker','x','MarkerEdgeColor','r');
           else
                plot(x(index(i),3),x(index(i),4),'Marker','x','MarkerEdgeColor','r');
           end
           hold on
       end
    end
  end %pkfile
  
end % empty x

set(handles.zoompushbutton,'userdata',handles.typetraj);

clear x x2 actualtrc
hold off
guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot in separate window
% --- Executes on button press in plottraj.
function plottraj_Callback(hObject, eventdata, handles)

set(handles.plottraj,'value',1);

% options for representation
linew=str2num(get(handles.linewidth,'string'));
pkchoice=get(handles.pkradiobutton,'value');
doubletrc=get(handles.filetrc2,'string');
molnro=get(handles.molindiv,'string');
identify=get(handles.ident,'value');
%firstbutton=get(handles.firstradiobutton,'value');
localiz=get(handles.loc,'value');
timecolor=get(handles.timecolorradiobutton,'value');
%blink=get(handles.blinking,'value');
rainbow=get(handles.rainbowradiobutton,'value');
%spine=get(handles.spineradiobutton,'value');
%diff=get(handles.Dradiobutton,'value');
%if diff==1 %D
%   diffdata=get(handles.Dradiobutton,'userdata');
%end

if timecolor==1
   [colorm]=createcolor(handles.param.lasttrc);
elseif rainbow==1
   colorm=handles.colorm;
end

x=get(handles.filetrc,'userdata'); %trc data
if isempty(x)==0
    indexnozero=find(x(:,3)>0);
    if isempty(indexnozero)==0
        x=x(indexnozero,:);
    end
end

x2=[];
if isempty(doubletrc)==0
    x2=get(handles.filetrc2,'userdata'); %trc data traj 2
end

%background
showbackground(handles)

if isempty(x)==0  % with trajectories
    
hold on;
indcol=1;
[totfilas, columns] = size (x);
name=get(handles.file1,'string');
option=1; 
trcname=get(handles.filetrc,'string');

if pkchoice==0
    codenormal=handles.color.all;
    indexmol=[];
    indexactual=[];
    actualtrc=[];
    indcol=1;
    
    k=strfind(molnro,'all');
    if isempty(k)==1
       indexmol=find(x(:,1)==str2num(molnro));
       actualtrc=x(indexmol,:);
    else
       actualtrc=x; 
    end

    for i=1:max(actualtrc(:,1))
            indexmol=find(actualtrc(:,1)==i);
            if isempty(indexmol)==0 
                codecol=handles.color.all;
                if localiz==0 & timecolor==0  
                   if rainbow==1
                      codecol=colorm(indcol,1:3);
                      indcol=indcol+1;
                      if indcol>50
                         indcol=1;
                      end
                   end
                   % normal, blinking or rainbow
                   plot(actualtrc(indexmol,3),actualtrc(indexmol,4),'Color',codecol,'LineWidth',linew);
                   hold on
                elseif localiz==0 & timecolor==1  
                   for g=1:size(indexmol,1)-1
                       codecol=colorm(g,1:3);
                       plot(actualtrc(indexmol(g):indexmol(g+1),3),actualtrc(indexmol(g):indexmol(g+1),4),'Color',codecol,'LineWidth',linew);
                       hold on
                   end
                elseif localiz==1 
                   if size(actualtrc,2)>5  % trc with localization, deco or not
                       spine=0; %!!!!!!!!!!!!!!!!!!!!!
                        for f=1:size(indexmol,1)-2
                            if spine==0 | size(actualtrc,2)<7
                                if actualtrc(indexmol(f+1),6)<0 %peri
                                    codecol=handles.color.peri;
                                elseif actualtrc(indexmol(f+1),6)>0 %syn
                                    codecol=handles.color.syn;
                                elseif actualtrc(indexmol(f+1),6)==0 %extra
                                    codecol=handles.color.extra;
                                end 
                            elseif spine==1 
                                if actualtrc(indexmol(f+1),7)>0 %neck
                                       codecol=handles.color.peri;
                                elseif actualtrc(indexmol(f+1),7)<0 %head
                                       codecol=handles.color.syn;
                                elseif actualtrc(indexmol(f+1),7)==0 %extra
                                       codecol=handles.color.extra;
                                end
                            end %spine
                            plot(actualtrc(indexmol(f):indexmol(f+1),3),actualtrc(indexmol(f):indexmol(f+1),4),'Color',codecol,'LineWidth',linew);
                            hold on
                        end
                   else
                        plot(actualtrc(indexmol,3),actualtrc(indexmol,4),'Color',codecol,'LineWidth',linew); % no syn
                        hold on
                   end
                end  % options colorcode   
               if identify==1 % numero tray
                   pos=indexmol(size(indexmol,1));
                  text(actualtrc(pos,3)+1,actualtrc(pos,4)+1,sprintf('%0.0f',actualtrc(pos,1)),'Color',[1 1 0]);
               end
               hold on
            else
            end % empty indexmol
      end % loop actual molecules

  else %pk
    
    accelerate=get(handles.accel,'value');
    if accelerate==1
       x=get(handles.accel,'userdata');
    end
     for i=1:size(x(:,1)) % all points
          if size(x,2)>4
                plot3(x(i,3),x(i,4),x(i,5),'Marker','x','MarkerEdgeColor','r');
          else
                plot(x(i,3),x(i,4),'Marker','x','MarkerEdgeColor','r');
          end
         hold on
     end
      hold off
  end %pkfile

       
if isempty(x2)==0
    plottraces(x2,'g',localiz,linew,handles)
end
hold off

end % empty x

set(handles.plottraj,'value',0);
clear x x2 graph

guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MERGING
% --- Executes on button press in merge.
function merge_Callback(hObject, eventdata, handles)

handles.maxfram=[];
count=1;
grayname=get(handles.grayname,'string');
redname=get(handles.redname,'string');
greenname=get(handles.greenname,'string');
bluename=get(handles.bluename,'string');

set(handles.file1,'string','Merged image');

msgbox('Reading files...');

% image DIC
if isempty(grayname)==0
    handles.maxfram(count)=handles.param.gray.nfram; count=count+1;
    count=count+1;
    if handles.param.gray.nfram>1
       handles.cgray.image=get(handles.grayname,'userdata');
    else
       handles.cgray.image(1).data=get(handles.grayname,'userdata');
    end
    Xdim=handles.param.gray.Xdim;
    Ydim=handles.param.gray.Ydim;
else
    handles.param.gray.nfram=0;
end


% red, green, blue
if isempty(redname)==0
    handles.maxfram(count)=handles.param.red.nfram; count=count+1;
    if handles.param.red.nfram>1
       handles.cred.image=get(handles.redname,'userdata');
    else
       handles.cred.image(1).data=get(handles.redname,'userdata');
    end
    Xdim=handles.param.red.Xdim;
    Ydim=handles.param.red.Ydim;
else
    handles.param.red.nfram=0;
end

if isempty(greenname)==0
    handles.maxfram(count)=handles.param.green.nfram; count=count+1;
    if handles.param.green.nfram>1
       handles.cgreen.image=get(handles.greenname,'userdata');
    else
       handles.cgreen.image(1).data=get(handles.greenname,'userdata');
    end
    Xdim=handles.param.green.Xdim;
    Ydim=handles.param.green.Ydim;
else
    handles.param.green.nfram=0;
end

if isempty(bluename)==0
    handles.maxfram(count)=handles.param.blue.nfram;count=count+1;
    if handles.param.blue.nfram>1
       handles.cblue.image=get(handles.bluename,'userdata');
    else
       handles.cblue.image(1).data=get(handles.bluename,'userdata');
    end
    Xdim=handles.param.blue.Xdim;
    Ydim=handles.param.blue.Ydim;
else
    handles.param.blue.nfram=0;
end

nfram=max(handles.maxfram);
if isempty(nfram)==1
    nfram=1;
end
for r=1:size(handles.maxfram,2)  % ojo si las movies tienen nro frame distinto
    if handles.maxfram(r)>1
        if handles.maxfram(r)<nfram
            nfram=handles.maxfram(r);
        end
    end
end
close %msgbox

% poner todo en handles
handles.param.actual=1;  % frame data
handles.param.lastimage=nfram;  % frame data
if handles.param.lastimage>handles.param.lasttrc
   handles.param.maxfram=handles.param.lastimage;
   set(handles.text11,'string',['1 (of ',num2str(handles.param.lastimage),')']);
else
   handles.param.maxfram=handles.param.lasttrc;
   set(handles.text11,'string',['1 (of ',num2str(handles.param.lasttrc),')']);
end

handles.param.gray.factor=str2num(get(handles.fgraylow,'string'));
handles.param.red.factor=str2num(get(handles.fredlow,'string'));
handles.param.green.factor=str2num(get(handles.fgreenlow,'string'));
handles.param.blue.factor=str2num(get(handles.fbluelow,'string'));
handles.param.gray.factorhigh=str2num(get(handles.fgrayhigh,'string'));
handles.param.red.factorhigh=str2num(get(handles.fredhigh,'string'));
handles.param.green.factorhigh=str2num(get(handles.fgreenhigh,'string'));
handles.param.blue.factorhigh=str2num(get(handles.fbluehigh,'string'));

handles.param.Xdim=Xdim;
handles.param.Ydim=Ydim;
handles.param.nfram=nfram;
handles.typefile=4;

%slider
minvalue=get(handles.slider1,'Min');
set(handles.slider1,'Value',minvalue);
handles.param.actual=1;

showbackground(handles);
set(handles.zoompushbutton,'enable','on');  

guidata(gcbo,handles) ;
clear all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in avibutton.
function avibutton_Callback(hObject, eventdata, handles)

molnro=get(handles.molindiv,'string');
createavi(handles,1,handles.param.maxfram,molnro)

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in zoompushbutton.
function plot3Dpushbutton_Callback(hObject, eventdata, handles)

plottraj3DSM(handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in zoompushbutton.
function zoompushbutton_Callback(hObject, eventdata, handles)

% elije el area a agrandar
warning off all
nromol=get(handles.molindiv,'string');

handles.typetraj=0;
syn=get(handles.loc,'value');
if syn==1
   handles.typetraj=1;
end
rainbow=get(handles.rainbowradiobutton,'value');
if rainbow==1
   handles.typetraj=3;
end
timecolor=get(handles.timecolorradiobutton,'value');
if timecolor==1
   handles.typetraj=4;
end
%firstbutton=get(handles.firstradiobutton,'value');
%if firstbutton==1
%   handles.typetraj=7;
%end

trcdata=get(handles.filetrc,'userdata'); %trc data
trcdata2=get(handles.filetrc2,'userdata'); %trc data double labeling

[handles,resp]=selectdataROI(trcdata, trcdata2,nromol, handles);

if isempty(resp)==0
    
    % dialog box to accept
    qstring=['Accept shift correction?'];
    button = questdlg(qstring); 
    if strcmp(button,'Yes')
        % save new images
        if isdir('newimages'); else; mkdir('newimages');end
        if resp(1)==1
            if isfield(handles,'newtrcdata')
                trcfile=get(handles.filetrc,'string')
                newtrc=handles.newtrcdata;
                save(['newimages\',trcfile],'newtrc','-ascii');
                set(handles.filetrc,'userdata',newtrc); %trc data
            end 
        end
        if resp(2)==1
            handles.cgray.image(1).data=get(handles.grayname,'userdata');
            grayname=get(handles.grayname,'string');
            imwrite(handles.cgray.image(1).data,['newimages\',grayname],'tif','compression','none');
        end
        if resp(3)==1
            handles.cred.image(1).data=get(handles.redname,'userdata');
            redname=get(handles.redname,'string');
            imwrite(handles.cred.image(1).data,['newimages\',redname],'tif','compression','none');
        end
        if resp(4)==1
            handles.cgreen.image(1).data=get(handles.greenname,'userdata');
            greenname=get(handles.greenname,'string');
            imwrite(handles.cgreen.image(1).data,['newimages\',greenname],'tif','compression','none');
        end
        if resp(5)==1
            handles.cblue.image(1).data=get(handles.bluename,'userdata');
            bluename=get(handles.bluename,'string');
            imwrite(handles.cblue.image(1).data,['newimages\',bluename],'tif','compression','none');
        end
    else
        %slider
        minvalue=get(handles.slider1,'Min');
        set(handles.slider1,'Value',minvalue);
        handles.param.actual=1;
        if isfield(handles,'movie')
            handles = rmfield(handles, 'movie'); %deletes handles.movie
        end
        showbackground(handles);
    end

end

if isfield(handles,'movie')
    handles = rmfield(handles, 'movie'); %deletes handles.movie
end


guidata(gcbo,handles) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function zoomDinstpushbutton_Callback(hObject, eventdata, handles)

warning off all
handles.typetraj=5;
nromol=str2num(get(handles.nromolDinst,'string'));
currentdir=cd;
trcdata=[];

% Dinst file
%currentdir=cd;
path=['Dinst\'];
if length(dir(path))>0
    nametrc=get(handles.filetrc,'string');
    if isempty(nametrc)==0
       [namefile,rem]=strtok(nametrc,'.'); %sin extension
       cd(path)
       [file,path] = uigetfile('*.txt*','Load Dinst file');
       if isempty(file)==1
           cd(currentdir)
           return
       end
        %filename = [path,file];
       data=load(file);
       index=find(data(:,1)==nromol);
       if isempty(index)==0
           trcdata=data(index,:);
       else
           msgbox('Trajectory not found','error','error')
           return
       end
    else
        % error
        msgbox('Load a traces file','error','error')
        return
    end
end
cd(currentdir)

trcdata2=[];
if isempty(trcdata)==0
   selectdataROI(trcdata, trcdata2,nromol, handles);
else
    disp('No Dinst data for this molecule')
end
guidata(gcbo,handles) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% trajectories file 2
% --- Executes on button press in traces.
function doubletraces_Callback(hObject, eventdata, handles)

tipotrc=get(handles.trcradiobutton,'value');
% loads traces file1 

%if pkfile==0
if tipotrc==1
    currentcd=cd;
    if length(dir(handles.presentfolder))>0
      cd(handles.presentfolder)
    end
   [trcf,tpath] = uigetfile('*.trc','Load trajectories file'); 
   cd(currentcd)
else
    currentdir=cd;
    %cd('C:\Trajectoires')
   [trcf,tpath] = uigetfile('*.traj','Load trajectories file'); 
   % peri=?
   perival=1;
   cd(currentdir)
end

trcfile = [tpath,trcf];
if trcfile==0
    set(handles.filetrc2,'userdata',[]);
    set(handles.filetrc2,'value',1);  
    set(handles.filetrc2,'string',['']);
    return
end
set(handles.filetrc2,'string',trcfile);
handles.traject2=get(handles.filetrc2,'string');
set(handles.filetrc2,'string',trcf);

%trfile
if isempty(trcfile)==0
    if tipotrc==1
       x=load(trcfile);   %.trc
    else
       [x,a,Te,nb_frames]=trajTRC(trcfile,perival);
    end
    nrotraj=max(x(:,2));
else
    x=[];
    nrotraj=1;
end
%
set(handles.filetrc2,'userdata',x);
set(handles.filetrc2,'value',nrotraj)  
handles.param.lasttrc2=nrotraj; %last frame with trc

%slider
minvalue=get(handles.slider1,'Min');
set(handles.slider1,'value',minvalue);

clear x
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function special_Callback(hObject, eventdata, handles)
%--------------------------------------------------------------------------

% --- Executes on button press in saveimage.
function saveimage_Callback(hObject, eventdata, handles)

[filename,path] = uiputfile('.tif','Save image as') ;
if filename==0
    return
end
presentimage=getframe(gca);  % gets the figure for the movie
[image,Map] = frame2im(presentimage);
imwrite(image,[path,filename],'tif');
clear image presentimage

%----------------------------------------------------------------------------------------------
% plot D in color in separate window
function plotD_Callback(hObject, eventdata, handles)
% plot puntos traj con info D

set(handles.plottraj,'value',1);
S = warning('off', 'all');
trc=get(handles.filetrc,'userdata'); %trc data

% dialog boxs to enter acquisition data
prompt = {'Size sliding window for MSD calculation: ','Pixel size: ','Acquisition time :'};
num_lines= 1;
dlg_title = 'Enter values for:';
def = {'10','167','10'}; % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
   if exit(1) == 0;
       %cd(currentdir);
       return; 
   end
sizewindow=str2num(answer{1}); %
szpx = str2num(answer{2}); 
till = str2num(answer{3}); %
[maxpoints,fil]=size(trc);

%background
showbackground(handles)
        
%definicion colores
categories=6;
[colorstep, cmap]=definecolor(0,[],6); % create color map

if isempty(trc)==0
else
    set(handles.plottraj,'value',0);
    return
end
         
for nro=1:max(trc(:,1)) ;    % cada molecula
    step=[];
    index=find(trc(:,1)==nro);  % puntos de trc de cada mol
             
    if isempty (index)==0
        for u=1:size(trc,1)-sizewindow   % para cada punto (salvo los ultimos!!!!)
            step=[];
            step=trcwindow(trc,index,u,sizewindow); %[step]=trcwindow(trc,index,u,sizewindow)
            if isempty(step)==0 % if there is not blinking longer than maxblink
               if size(step,1)>5
                  [D,b,MSD]=calculMSD(step,szpx,till,4,5);
                  msddata=MSD.rho;
                  if D>0
                    ang= atan2(step(2,3)-step(1,3),step(2,4)-step(1,4)); % arcotangente
                    colorstep=definecolor(D,cmap,categories);
                    distx=[step(1,3); step(1,3)+(0.15*sin(ang))];
                    disty=[step(1,4); step(1,4)+(0.15*cos(ang))];
                    plot (distx, disty,'color',colorstep);   % grafica direccion
                    hold on
                    % plot ((step (:,1))+1, (step (:,2))+1,'color',colorstep,'LineWidth',1.5);   % grafica traces
                    plot ((step (1,3)), (step (1,4)),'color',colorstep,'Marker','.','MarkerSize',5);   % grafica puntos
                    hold on
                 end % D=0
               end
            end % empty step
      end  %puntos
    end  %existe la mol
         
end   %loop

hold off
set(handles.plottraj,'value',0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function util_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
function perisynimage_Callback(hObject, eventdata, handles)
imageperi(handles);

%--------------------------------------------------------------------------
function domainarea_Callback(hObject, eventdata, handles)
measureareasyn

%--------------------------------------------------------------------------
function mask_Callback(hObject, eventdata, handles)

coloc2MIAimages3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in helpbutton.
function help_Callback(hObject, eventdata, handles)

[path]=readfolder; 
open([path,'\help\5-MovieMaker.pdf'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in quit.
function quit_Callback(hObject, eventdata, handles)

clear all;
close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function till_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function firstframe_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function lastframe_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function molindiv_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function nromolDinst_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function delay_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function fgraylow_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function fredlow_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function fgreenlow_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function fbluelow_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% --- Executes during object creation, after setting all properties.
function fgrayhigh_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function fredhigh_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function fgreenhigh_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function fbluehigh_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%% --- Executes during object creation, after setting all properties.
function linewidth_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%% --- Executes during object creation, after setting all properties.
function corrhor_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%% --- Executes during object creation, after setting all properties.
function corrvert_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%function trajtimecolorradiobutton_Callback(hObject, eventdata, handles)
