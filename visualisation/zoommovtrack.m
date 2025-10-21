function varargout = zoommovtrack(varargin)
% ZOOMMOVTRACK M-file for zoommovtrack.fig
% Creates zoomed images from movtrack.m
% possibility to save snapshot or .avi
% displays Dinst if called accordingly
% MR - sep 07 - for SPTrack                                   MatLab7
% MR - mar 09 - for SPTrack v4.0                              MatLab7
% MR - jan 10 - for SPTrack v4.0                              MatLab7
% MR - jun 10 - for SPTrack v4.0                              MatLab7
% MR - jan 22 - SPTrack_v6
% Marianne Renner 01/2025 - adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Last Modified by GUIDE v2.5 12-Jul-2007 23:27:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @zoommovtrack_OpeningFcn, ...
                   'gui_OutputFcn',  @zoommovtrack_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes just before zoommovtrack is made visible.
function zoommovtrack_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

output=[];
set(handles.finishedpushbutton,'userdata',output);
save('auxiliar.mat','output','-mat')

varargout{1} = [];

handles.movie=cell2mat(varargin{1}(1));
handles.Xdim=cell2mat(varargin{1}(2));
handles.Ydim=cell2mat(varargin{1}(3));
handles.param=cell2mat(varargin{1}(4));
handles.traces=cell2mat(varargin{1}(5));
handles.cumulative=cell2mat(varargin{1}(6));
handles.typetraj=cell2mat(varargin{1}(7));
handles.nromolecule=cell2mat(varargin{1}(8));
handles.color=cell2mat(varargin{1}(9));
handles.traces2=cell2mat(varargin{1}(10));
handles.newtraces=handles.traces;

set(handles.framenumber,'string',[num2str(1),' of ',num2str(handles.param.maxfram)]);

if handles.typetraj==5 % Dinst
    set(handles.Dcolor,'enable','on')
    set(handles.Dnumbers,'enable','on')
else
    set(handles.Dcolor,'enable','off')
    set(handles.Dnumbers,'enable','off')
end

if handles.typetraj==3 || handles.typetraj==7
  [colorm]=createcolor(50);
  for i=1:size(colorm,1)
          indcol=round(rand(1)*size(colorm,1));
          colorm(i,4)=indcol;
  end
  colorm=sortrows(colorm,4);
  set(handles.saveimagepushbutton,'userdata',colorm);
end

if handles.param.lastimage>1 % movie
       set(handles.savefilespushbutton,'enable','off');
       set(handles.textadvert,'string','Attention: cannot save .stk');
end
   
if isempty(handles.traces)==0
    set(handles.trajverticalshift,'enable','on')
    set(handles.trajhorizontalshift,'enable','on')
end
handles.newdic.data=[];
handles.newred.data=[];
handles.newgreen.data=[];
handles.newblue.data=[];

if isfield(handles.movie,'gray')
    % not stk!!
    if handles.param.gray.nfram==1
        set(handles.dicverticalshift,'enable','on')
        set(handles.dichorizontalshift,'enable','on')
        handles.newdic.data=handles.movie.gray.data;
    elseif handles.param.gray.nfram==0
        handles.newgray.data=handles.movie.gray.data;
    end
end
if isfield(handles.movie,'red')
    % not stk!!
    if handles.param.red.nfram==1
        set(handles.redverticalshift,'enable','on')
        set(handles.redhorizontalshift,'enable','on')
        handles.newred.data=handles.movie.red.data;
    elseif handles.param.red.nfram==0
        handles.newred.data=handles.movie.red.data;
    end
end
if isfield(handles.movie,'green')
    % not stk!!
    if handles.param.green.nfram==1
        set(handles.greenverticalshift,'enable','on')
        set(handles.greenhorizontalshift,'enable','on')
        handles.newgreen.data=handles.movie.green.data;
    elseif handles.param.green.nfram==0
        handles.newgreen.data=handles.movie.green.data;
    end
end
if isfield(handles.movie,'blue')
    % not stk!!
    if handles.param.blue.nfram==1
        set(handles.blueverticalshift,'enable','on')
        set(handles.bluehorizontalshift,'enable','on')
        handles.newblue.data=handles.movie.blue.data;
    elseif handles.param.blue.nfram==0
        handles.newblue.data=handles.movie.blue.data;
    end
end


set(handles.slider1,'value',0);
clear data dimx dimy traces frames

handles.param.gray.factor=str2num(get(handles.zoomfgraylow,'string'));
handles.param.red.factor=str2num(get(handles.zoomfredlow,'string'));
handles.param.green.factor=str2num(get(handles.zoomfgreenlow,'string'));
handles.param.blue.factor=str2num(get(handles.zoomfbluelow,'string'));
handles.param.gray.factorhigh=str2num(get(handles.zoomfgrayhigh,'string'));
handles.param.red.factorhigh=str2num(get(handles.zoomfredhigh,'string'));
handles.param.green.factorhigh=str2num(get(handles.zoomfgreenhigh,'string'));
handles.param.blue.factorhigh=str2num(get(handles.zoomfbluehigh,'string'));
set(handles.finishedpushbutton,'value',0);
set(handles.panel,'userdata',zeros(4,4));
set(handles.nroframe,'value',1);

axes(handles.axes1);
showframeROI(handles);

guidata(hObject,handles) ;


%--------------------------------------------------------------------
function varargout = zoommovtrack_OutputFcn(hObject, eventdata, handles) 

output=get(handles.finishedpushbutton,'userdata');

varargout{1} = output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Dcolor_Callback(hObject, eventdata, handles)
set(handles.Dnumbers,'value',0);
guidata(hObject, handles);

function Dnumbers_Callback(hObject, eventdata, handles)
set(handles.Dcolor,'value',0);
guidata(hObject, handles);

function zoomfgraylow_Callback(hObject, eventdata, handles)
handles.param.gray.factor=1-str2num(get(hObject,'string'));
showframeROI(handles);
guidata(hObject, handles);

function zoomfredlow_Callback(hObject, eventdata, handles)
handles.param.red.factor=1-str2num(get(hObject,'string'));
showframeROI(handles);
guidata(hObject, handles);

function zoomfgreenlow_Callback(hObject, eventdata, handles)
handles.param.green.factor=1-str2num(get(hObject,'string'));
showframeROI(handles);
guidata(hObject, handles);

function zoomfbluelow_Callback(hObject, eventdata, handles)
handles.param.blue.factor=1-str2num(get(hObject,'string'));
showframeROI(handles);
guidata(hObject, handles);

function zoomfgrayhigh_Callback(hObject, eventdata, handles)
handles.param.gray.factorhigh=1-str2num(get(hObject,'string'));
showframeROI(handles);
guidata(hObject, handles);

function zoomfredhigh_Callback(hObject, eventdata, handles)
handles.param.red.factorhigh=1-str2num(get(hObject,'string'));
showframeROI(handles);
guidata(hObject, handles);

function zoomfgreenhigh_Callback(hObject, eventdata, handles)
handles.param.green.factorhigh=1-str2num(get(hObject,'string'));
showframeROI(handles);
guidata(hObject, handles);

function zoomfbluehigh_Callback(hObject, eventdata, handles)
handles.param.blue.factorhigh=1-str2num(get(hObject,'string'));
showframeROI(handles);
guidata(hObject, handles);

function linewidth_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

function alltrajradiobutton_Callback(hObject, eventdata, handles)
valtraj=get(hObject,'value');
if valtraj==1
    set(handles.Dcolor,'enable','off')
    set(handles.Dnumbers,'enable','off')
else
    set(handles.Dcolor,'enable','on')
    set(handles.Dnumbers,'enable','on')
end
showframeROI(handles);
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zoommergepushbutton_Callback(hObject, eventdata, handles)

handles.param.gray.factor=str2num(get(handles.zoomfgraylow,'string'));
handles.param.red.factor=str2num(get(handles.zoomfredlow,'string'));
handles.param.green.factor=str2num(get(handles.zoomfgreenlow,'string'));
handles.param.blue.factor=str2num(get(handles.zoomfbluelow,'string'));
handles.param.gray.factorhigh=str2num(get(handles.zoomfgrayhigh,'string'));
handles.param.red.factorhigh=str2num(get(handles.zoomfredhigh,'string'));
handles.param.green.factorhigh=str2num(get(handles.zoomfgreenhigh,'string'));
handles.param.blue.factorhigh=str2num(get(handles.zoomfbluehigh,'string'));

showframeROI(handles);
guidata(gcbo,handles) ;

%----------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)

frames=handles.param.maxfram;

slidervalue=get(hObject,'Value');
minvalue=get(hObject,'Min');
maxvalue=get(hObject,'Max');
prop=slidervalue/(minvalue+maxvalue);

frame=round(prop*frames);
step=(maxvalue-minvalue)/frames;
set(handles.slider1,'SliderStep',[step step]);
if frame<1
    frame=1;
end
if frame>frames
    frame=frames;
end
set(handles.nroframe,'value',frame);
set(handles.framenumber,'string',[num2str(frame),' of ',num2str(frames)]);

clear handles.frames frame step prop minvalue maxvalue

showframeROI(handles);

%--------------------------------------------------------------------

function showframeROI(handles)

% display background image
handles.roiclean=[];
actualframe=get(handles.nroframe,'value');

handles.param.gray.factor=str2num(get(handles.zoomfgraylow,'string'));
handles.param.red.factor=str2num(get(handles.zoomfredlow,'string'));
handles.param.green.factor=str2num(get(handles.zoomfgreenlow,'string'));
handles.param.blue.factor=str2num(get(handles.zoomfbluelow,'string'));
handles.param.gray.factorhigh=str2num(get(handles.zoomfgrayhigh,'string'));
handles.param.red.factorhigh=str2num(get(handles.zoomfredhigh,'string'));
handles.param.green.factorhigh=str2num(get(handles.zoomfgreenhigh,'string'));
handles.param.blue.factorhigh=str2num(get(handles.zoomfbluehigh,'string'));

% control sliding bar
if actualframe==0
    actualframe=1;
elseif actualframe>handles.param.maxfram
    actualframe=handles.param.maxfram;
end

% show frame on axesroi
axes(handles.axes1);
set(handles.framenumber,'string',[num2str(actualframe),' of ',num2str(handles.param.maxfram)]);

if isfield(handles.movie,'gray') % merge
    % shift
    val=get(handles.finishedpushbutton,'value');
    if val==1 % show correction
        actualimagen=shownewmerge(handles,actualframe);
    else
        actualimagen=showmerge(handles,actualframe);
    end
    if actualframe==1
        valint(1)=(min(min(min(actualimagen.data))));
        valint(2)=(max(max(max(actualimagen.data))));
        set(handles.zoommergepushbutton,'userdata',valint);
    else
       valint=get(handles.zoommergepushbutton,'userdata');
    end
   imshow(actualimagen.data,[valint(1) valint(2)],'InitialMagnification','fit');
   
else
  if isstruct(handles.movie)
      if size(handles.movie,2)>1
          actualimagen = handles.movie(actualframe).data;
      else
          actualimagen = handles.movie.data;
      end
  else 
      actualimagen = handles.movie;
  end
  if actualframe==1
        valint(1)=(min(min(min(actualimagen))));
        valint(2)=(max(max(max(actualimagen))));
        set(handles.zoommergepushbutton,'userdata',valint);
  else
       valint=get(handles.zoommergepushbutton,'userdata');
  end
  imshow(actualimagen,[valint(1) valint(2)],'InitialMagnification','fit');
end
hold on

clear actualimage roimovie movie

% plot trayectorias
showtraj(handles)

guidata(gcbo,handles) ;

%--------------------------------------------------------------------------
function showtraj(handles)
% for each frame, makes an array with traces of the molecules and plots them

valtraj=get(handles.alltrajradiobutton,'value');

if size(handles.newtraces,1)==0
    x=handles.traces;
else
    x=handles.newtraces; %shift correction
end

x2=handles.traces2;
molnro=handles.nromolecule;
if handles.typetraj==5
    molnro=num2str(molnro);
end

linew=str2num(get(handles.linewidth,'string'));

axes(handles.axes1);
selecD=get(handles.Dcolor,'value');
selecnumber=get(handles.Dnumbers,'value');
if handles.typetraj==5 && selecnumber==1
    handles.typetraj=0;
end

if isempty(x)==0
   actualframe=get(handles.nroframe,'value');
   if actualframe==0
      actualframe=1;
   end
   if size(x,2)>4 %trc
      if handles.typetraj==5 && selecD==1   % Dinst in color
            % definicion colores
            categories=6;
            [colorstep, cmap]=definecolor(0,[],categories); % create color map
            index=find(x(:,2)<actualframe+1);
            if isempty(index)==0
               for j=1: size(index,1)
                   D=x(j,7); % Dinst
                   [colorstep, cmap]=definecolor(D,cmap,0);
                   plot ((x (j,3)), (x (j,4)),'color',colorstep,'Marker','.','MarkerSize',5);   % grafica puntos
                   hold on
               end
            end
            
       else % trajectories
           
          if valtraj==1
              %all traj
              plotall(x,handles)
           
          else

            %maxmol=max(x(:,1));  % numero corregido
            %listmoltotal='';
            %listmol='';
            indexmol=[];
            indexactual=[];
            actualtrc=[];
            colorm=[];
            indcol=1;
            localiz=0;
            identify=0;
            firstbutton=0;
            rainbow=0;
            timecolor=0;
            
            if handles.typetraj==0
            elseif handles.typetraj==1 || handles.typetraj==6
                localiz=1;
            elseif handles.typetraj==2
                blinking=1;
            elseif handles.typetraj==3
                rainbow=1;
            elseif handles.typetraj==4
                timecolor=1;
            elseif handles.typetraj==7
                firstbutton=1;
                rainbow=1;
            end 
            
            if rainbow==1
                 colorm=get(handles.saveimagepushbutton,'userdata');
            end
            
            k=strfind(molnro,'all');
            if isempty(k)==1
                indexmol=find(x(:,1)==str2num(molnro));
                trc=x(indexmol,:);
            else
                trc=x; 
            end

            indexactual=find(trc(:,2)<actualframe+1);
            
            if isempty(indexactual)==0
                actualtrc=trc(indexactual,:);
                
                for i=1:max(actualtrc(:,1))
                    indexmol=find(actualtrc(:,1)==i);
                    indextrc=find(trc(:,1)==i);

                    if handles.cumulative==1 %cumulative
                        condition=1;
                    else
                        if max(trc(indextrc,2))>actualframe %present before and after
                            condition=1;
                        else
                            condition=0;
                        end
                    end
                    
                    if isempty(indexmol)==0 & condition==1 %
                        codecol=handles.color.all;   
                        
                        if firstbutton==1 % primer punto
                            if rainbow==1
                                handles.typetraj=3;
                                codecol=colorm(indcol,1:3);
                                indcol=indcol+1;
                                if indcol>50
                                    indcol=1;
                                end
                            end
                            
                            auxtrc=actualtrc(indexmol,:);
                            indextime=find(auxtrc(:,2)==actualframe);   
                            
                            if isempty(indextime)==0  
                                if auxtrc(indextime(1),6)>0 & auxtrc(indextime(1),6)<1000 %slot
                                    codecol=handles.color.extra; % blue
                                else
                                    codecol=handles.color.all; %red
                                end
                                if auxtrc(indextime(1),7)==1 % group 1
                                    codecol=handles.color.extra; % blue
                                else
                                    codecol=handles.color.all; %red
                                end
                                codecol=handles.color.all;
                                plot(auxtrc(indextime(1),3),auxtrc(indextime(1),4),'Color',codecol,'Marker','.','MarkerSize',20);
                                hold on
                            end % empty(indextime)
                            
                            clear auxtrc
                        else %firstbutton
                            if localiz==0 && handles.typetraj<4
                                if actualtrc(indexmol(size(indexmol,1)),2)<actualframe  && handles.typetraj==2 % blink period
                                    codecol=handles.color.blink;
                                elseif handles.typetraj==3
                                    codecol=colorm(indcol,1:3);
                                    indcol=indcol+1;
                                    if indcol>50
                                        indcol=1;
                                    end
                                end
                                % normal, blinking or raibow
                                plot(actualtrc(indexmol,3),actualtrc(indexmol,4),'Color',codecol,'Linewidth',linew);
                                hold on
                            elseif localiz==0 && handles.typetraj==4  
                               % for g=1:size(indexmol,1)-1
                               %     codecol=colorm(g,1:3)
                               %     plot(actualtrc(indexmol(g):indexmol(g+1),3),actualtrc(indexmol(g):indexmol(g+1),4),'Color',codecol,'Linewidth',linew);
                               %     hold on
                               % end
                            elseif localiz==1
                                if size(trc,2)>5  % trc with localization, deco or not
                                    for f=1:size(indexmol,1)-2
                                        if handles.typetraj<6 || size(trc,2)<7
                                            if actualtrc(indexmol(f+1),6)<0 %peri
                                                codecol=handles.color.peri;
                                            elseif actualtrc(indexmol(f+1),6)>0 %syn
                                                codecol=handles.color.syn;
                                            elseif actualtrc(indexmol(f+1),6)==0 %extra
                                                codecol=handles.color.extra;
                                            end
                                        elseif handles.typetraj==6
                                            if actualtrc(indexmol(f+1),7)>0 %neck
                                                codecol=handles.color.peri;
                                            elseif actualtrc(indexmol(f+1),7)<0 %head
                                                codecol=handles.color.syn;
                                            elseif actualtrc(indexmol(f+1),7)==0 %extra
                                                codecol=handles.color.extra;
                                            end
                                        end %
                                        
                                        plot(actualtrc(indexmol(f):indexmol(f+1),3),actualtrc(indexmol(f):indexmol(f+1),4),'Color',codecol,'Linewidth',linew)
                                        hold on
                                    end %for
                                else
                                    plot(actualtrc(indexmol,3),actualtrc(indexmol,4),'Color',codecol,'Linewidth',linew); % no syn
                                    hold on
                                end %size trc
                            end  % options colorcode  
                        end % firsbutton

                        if identify==1 % numero tray
                            pos=indexmol(size(indexmol,1));
                            text(actualtrc(pos,3)+1,actualtrc(pos,4)+1,sprintf('%0.0f',actualtrc(pos,1)),'Color',[1 1 0]);
                            cifras=num2str(actualtrc(pos,1)); space=size(cifras,2);
                            text(actualtrc(pos,3)+(space*5),actualtrc(pos,4)+1,sprintf('(%0.0f)',size(actualtrc(indexmol),1)),'Color',[1 1 1],'FontSize',7);
                        end
                        
                        if selecnumber==1 % D inst
                            pos=indexmol(size(indexmol,1));
                            cifras=num2str(actualtrc(pos,1)); space=size(cifras,2);
                            text(actualtrc(pos,3)+0.5,actualtrc(pos,4)+0.5,sprintf('%1.4f',actualtrc(pos,7)),'Color',[1 1 1],'FontSize',7);
                        end
                        
                        hold on
                    else
                        plot(handles.param.Xdim,handles.param.Ydim,'.k');  %just no avoid crash!
                        hold on    
                    end % empty indexmol
                end % loop actual molecules
            else
                plot(handles.param.Xdim,handles.param.Ydim,'.k');  %just no avoid crash!
                hold on  
            end %empty indexactual

            if isempty(x2)==0 % second traj
                codecol='r';
                synflag=0;
                plottracesframe(x2,actualframe,codecol,synflag,linew,handles)
            end
          end %all traj
       end% % Dinst

  else %pk
    
    accelerate=get(handles.accel,'value');
    if accelerate==1
       x=get(handles.accel,'userdata');
    end
    % frame per frame
    index=find(x(:,2)<actualframe+1);
    if isempty(index)==0
       for i=1:size(index,1)
           plot(x(index(i),3),x(index(i),4),'Marker','x','MarkerEdgeColor','r');
           hold on
       end
    end
  end %pkfile

end

clear actualtraces x
hold off

guidata(gcbo,handles) ;

%--------------------------------------------------------------------------

% --- Executes on button press in saveimagepushbutton.
function saveimagepushbutton_Callback(hObject, eventdata, handles)

[filename,path] = uiputfile('.tif','Save image as') ;
if filename==0
    return
end
presentimage=getframe(gca);  % gets the figure for the movie
[image,Map] = frame2im(presentimage);
imwrite(image,[path,filename],'tif');
clear image presentimage

%--------------------------------------------------------------------------

function savefilespushbutton_Callback(hObject, eventdata, handles)

movie=handles.movie.data;
x=handles.traces;

if handles.finalimage==1
   [filename,path] = uiputfile('.tif','Save background image as') ;
   if filename==0
       return
   end
   image=movie ;  % 
   imwrite(image,[path,filename],'tif');
   clear image 
else
    msgbox('Cannot save a movie','error','error');
end

if size(x,2)>4 %trc
   [filename,path] = uiputfile('.trc','Save tracking points as:') ;
   if filename==0
      return
   end
   save([path,filename],'x','-ascii');
else
    newx=x(:,2:4);
   [filename,path] = uiputfile('.pk','Save tracking points as:') ;
   if filename==0
      return
   end
   save([path,filename],'newx','-ascii');
end

clear x newx movie


%-----------------------------------------------------------------------
% --- Executes on button press in makeavipushbutton.
function makeavipushbutton_Callback(hObject, eventdata, handles)

createavizoom(handles,1,handles.param.maxfram,'all')
disp('Done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% correct shift
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trajUbutton_Callback(hObject, eventdata, handles)

status=get(handles.trajverticalshift,'enable');
if size(status,2)==2 %on

%UP
handles.thold=get(handles.trajverticalshift,'String');
ktev=str2num(handles.thold);
ktev=ktev-0.2;
set(handles.trajverticalshift,'string',num2str(ktev));
datashift=get(handles.panel,'userdata');
datashift(1,4)=ktev;
set(handles.panel,'userdata',datashift);


handles.newtraces=shifttrc(handles,ktev,datashift(2,4),handles.traces);

indexzero=find(handles.newtraces(:,4)>0);
handles.newtraces=handles.newtraces(indexzero,:);
set(handles.finishedpushbutton,'value',1);
showframeROI(handles)

end %status

guidata(hObject, handles);

%-----------------------------------------------------------------
function trajDbutton_Callback(hObject, eventdata, handles)

status=get(handles.trajverticalshift,'enable');
if size(status,2)==2 %on

%DOWN
handles.thold=get(handles.trajverticalshift,'String');
ktev=str2num(handles.thold);

ktev=ktev+0.2;
set(handles.trajverticalshift,'string',num2str(ktev));
datashift=get(handles.panel,'userdata');
datashift(1,4)=ktev;
set(handles.panel,'userdata',datashift);

handles.newtraces=shifttrc(handles,ktev,datashift(2,4),handles.traces);

indexzero=find(handles.newtraces(:,4)<handles.Ydim);
handles.newtraces=handles.newtraces(indexzero,:);
set(handles.finishedpushbutton,'value',1);

showframeROI(handles)

end

guidata(hObject, handles);

%-----------------------------------------------------------------
function trajLbutton_Callback(hObject, eventdata, handles)

status=get(handles.trajhorizontalshift,'enable');
if size(status,2)==2 %on

%LEFT
handles.thold=get(handles.trajhorizontalshift,'String');
kteh=str2num(handles.thold);
kteh=kteh-0.2;
set(handles.trajhorizontalshift,'string',num2str(kteh));
datashift=get(handles.panel,'userdata');
datashift(2,4)=kteh;
set(handles.panel,'userdata',datashift);

handles.newtraces=shifttrc(handles,datashift(1,4),kteh,handles.traces);
indexzero=find(handles.newtraces(:,3)>0);
handles.newtraces=handles.newtraces(indexzero,:);
set(handles.finishedpushbutton,'value',1);

showframeROI(handles)

end

guidata(hObject, handles);

%-----------------------------------------------------------------
function trajRbutton_Callback(hObject, eventdata, handles)

status=get(handles.trajhorizontalshift,'enable');
if size(status,2)==2 %on

%RIGHT
handles.thold=get(handles.trajhorizontalshift,'String');
kteh=str2num(handles.thold);
kteh=kteh+0.2;
set(handles.trajhorizontalshift,'string',num2str(kteh));
datashift=get(handles.panel,'userdata');
datashift(2,4)=kteh;
set(handles.panel,'userdata',datashift);

handles.newtraces=shifttrc(handles,datashift(1,4),kteh,handles.traces);

indexzero=find(handles.newtraces(:,3)<handles.Xdim);
handles.newtraces=handles.newtraces(indexzero,:);
set(handles.finishedpushbutton,'value',1);

showframeROI(handles)

end

guidata(hObject, handles);

%-----------------------------------------------------------------
function dicUbutton_Callback(hObject, eventdata, handles)
% UP

status=get(handles.dicverticalshift,'enable');
if size(status,2)==2 %on

handles.thold=get(handles.dicverticalshift,'String');
ktev=str2num(handles.thold);
ktev=ktev-1;
set(handles.dicverticalshift,'String',num2str(ktev));
datashift=get(handles.panel,'userdata');
datashift(1,1)=ktev;
set(handles.panel,'userdata',datashift);

handles.newdic.data=shiftimage(handles,ktev,datashift(2,1),handles.movie.gray.data,handles.Ydim,1);
set(handles.finishedpushbutton,'value',1);

showframeROI(handles)

end

guidata(hObject, handles);

%-----------------------------------------------------------------

% --- Executes on button press in dicDbutton.
function dicDbutton_Callback(hObject, eventdata, handles)

status=get(handles.dicverticalshift,'enable');
if size(status,2)==2 %on

% DOWN
handles.thold=get(handles.dicverticalshift,'String');
ktev=str2num(handles.thold);
ktev=ktev+1;
set(handles.dicverticalshift,'String',num2str(ktev));
datashift=get(handles.panel,'userdata');
datashift(1,1)=ktev;
set(handles.panel,'userdata',datashift);

handles.newdic.data=shiftimage(handles,ktev,datashift(2,1),handles.movie.gray.data,handles.Ydim,1);
set(handles.finishedpushbutton,'value',1);

showframeROI(handles)

end

guidata(hObject, handles);

%-----------------------------------------------------------------


% --- Executes on button press in redUbutton.
function redUbutton_Callback(hObject, eventdata, handles)

status=get(handles.redverticalshift,'enable');
if size(status,2)==2 %on

% UP
handles.thold=get(handles.redverticalshift,'String');
ktev=str2num(handles.thold);
ktev=ktev-1;
set(handles.redverticalshift,'String',num2str(ktev));
datashift=get(handles.panel,'userdata');
datashift(1,2)=ktev;
set(handles.panel,'userdata',datashift);

if handles.param.red.nfram==1
    handles.newred.data=shiftimage(handles,ktev,datashift(2,2),handles.movie.red.data,handles.Ydim,1);
    set(handles.finishedpushbutton,'value',1);
    showframeROI(handles)
end

end
guidata(hObject, handles);


%-----------------------------------------------------------------

% --- Executes on button press in redDbutton.
function redDbutton_Callback(hObject, eventdata, handles)
status=get(handles.redverticalshift,'enable');
if size(status,2)==2 %on

% DOWN
handles.thold=get(handles.redverticalshift,'String');
ktev=str2num(handles.thold);
ktev=ktev+1;
set(handles.redverticalshift,'String',num2str(ktev));
datashift=get(handles.panel,'userdata');
datashift(1,2)=ktev;
set(handles.panel,'userdata',datashift);

if handles.param.red.nfram==1
    handles.newred.data=shiftimage(handles,ktev,datashift(2,2),handles.movie.red.data,handles.Ydim,1);
    set(handles.finishedpushbutton,'value',1);
    showframeROI(handles)
end

end

guidata(hObject, handles);


%-----------------------------------------------------------------

% --- Executes on button press in greenUbutton.
function greenUbutton_Callback(hObject, eventdata, handles)
status=get(handles.greenverticalshift,'enable');
if size(status,2)==2 %on

% UP
handles.thold=get(handles.greenverticalshift,'String');
ktev=str2num(handles.thold);
ktev=ktev-1;
set(handles.greenverticalshift,'String',num2str(ktev));
datashift=get(handles.panel,'userdata');
datashift(1,3)=ktev;
set(handles.panel,'userdata',datashift);

if handles.param.green.nfram==1
    handles.newgreen.data=shiftimage(handles,ktev,datashift(2,3),handles.movie.green.data,handles.Ydim,1);
    set(handles.finishedpushbutton,'value',1);
    showframeROI(handles)
end


end

guidata(hObject, handles);

%-----------------------------------------------------------------


% --- Executes on button press in greenDbutton.
function greenDbutton_Callback(hObject, eventdata, handles)
status=get(handles.greenverticalshift,'enable');
if size(status,2)==2 %on

    % DOWN
handles.thold=get(handles.greenverticalshift,'String');
ktev=str2num(handles.thold);
ktev=ktev+1;
set(handles.greenverticalshift,'String',num2str(ktev));
datashift=get(handles.panel,'userdata');
datashift(1,3)=ktev;
set(handles.panel,'userdata',datashift);

if handles.param.green.nfram==1
    handles.newgreen.data=shiftimage(handles,ktev,datashift(2,3),handles.movie.green.data,handles.Ydim,1);
    set(handles.finishedpushbutton,'value',1);
    showframeROI(handles)
end

end



guidata(hObject, handles);

%-----------------------------------------------------------------

% --- Executes on button press in blueUbutton.
function blueUbutton_Callback(hObject, eventdata, handles)

status=get(handles.blueverticalshift,'enable');
if size(status,2)==2 %on

% UP
handles.thold=get(handles.blueverticalshift,'String');
ktev=str2num(handles.thold);
ktev=ktev-1;
set(handles.blueverticalshift,'String',num2str(ktev));
datashift=get(handles.panel,'userdata');
datashift(1,4)=ktev;
set(handles.panel,'userdata',datashift);

if handles.param.blue.nfram==1
    handles.newblue.data=shiftimage(handles,ktev,datashift(2,4),handles.movie.blue.data,handles.Ydim,1);
    set(handles.finishedpushbutton,'value',1);
    showframeROI(handles)
end

end

guidata(hObject, handles);

%-----------------------------------------------------------------

% --- Executes on button press in blueDbutton.
function blueDbutton_Callback(hObject, eventdata, handles)
status=get(handles.blueverticalshift,'enable');
if size(status,2)==2 %on
% DOWN
handles.thold=get(handles.blueverticalshift,'String');
ktev=str2num(handles.thold);
ktev=ktev+1;
set(handles.blueverticalshift,'String',num2str(ktev));
datashift=get(handles.panel,'userdata');
datashift(1,4)=ktev;
set(handles.panel,'userdata',datashift);

if handles.param.blue.nfram==1
    handles.newblue.data=shiftimage(handles,ktev,datashift(2,4),handles.movie.blue.data,handles.Ydim,1);
    set(handles.finishedpushbutton,'value',1);
    showframeROI(handles)
end

end

guidata(hObject, handles);


%-----------------------------------------------------------------

% --- Executes on button press in dicLbutton.
function dicLbutton_Callback(hObject, eventdata, handles)
status=get(handles.dichorizontalshift,'enable');
if size(status,2)==2 %on

% LEFT
handles.thold=get(handles.dichorizontalshift,'String');
kteh=str2num(handles.thold);
kteh=kteh-1;
set(handles.dichorizontalshift,'String',num2str(kteh));
datashift=get(handles.panel,'userdata');
datashift(2,1)=kteh;
set(handles.panel,'userdata',datashift);

if handles.param.gray.nfram==1
    handles.newdic.data=shiftimage(handles,datashift(1,1),kteh,handles.movie.gray.data,handles.Xdim,2);
    set(handles.finishedpushbutton,'value',1);
    showframeROI(handles)
end

end

guidata(hObject, handles);

%-----------------------------------------------------------------

% --- Executes on button press in dicRbutton.
function dicRbutton_Callback(hObject, eventdata, handles)
status=get(handles.dichorizontalshift,'enable');
if size(status,2)==2 %on
% RIGHT
handles.thold=get(handles.dichorizontalshift,'String');
kteh=str2num(handles.thold);
kteh=kteh+1;
set(handles.dichorizontalshift,'String',num2str(kteh));
datashift=get(handles.panel,'userdata');
datashift(2,1)=kteh;
set(handles.panel,'userdata',datashift);

if handles.param.gray.nfram==1
    handles.newdic.data=shiftimage(handles,datashift(1,1),kteh,handles.movie.gray.data,handles.Xdim,2);
    set(handles.finishedpushbutton,'value',1);
    showframeROI(handles)
end
end
guidata(hObject, handles);

%-----------------------------------------------------------------

% --- Executes on button press in redLbutton.
function redLbutton_Callback(hObject, eventdata, handles)
status=get(handles.redhorizontalshift,'enable');
if size(status,2)==2 %on

% LEFT
handles.thold=get(handles.redhorizontalshift,'String');
kteh=str2num(handles.thold);
kteh=kteh-1;
set(handles.redhorizontalshift,'String',num2str(kteh));
datashift=get(handles.panel,'userdata');
datashift(2,2)=kteh;
set(handles.panel,'userdata',datashift);

if handles.param.red.nfram==1
    handles.newred.data=shiftimage(handles,datashift(1,2),kteh,handles.movie.red.data,handles.Xdim,2);
    set(handles.finishedpushbutton,'value',1);
    showframeROI(handles)
end
end
guidata(hObject, handles);

%-----------------------------------------------------------------

% --- Executes on button press in redRbutton.
function redRbutton_Callback(hObject, eventdata, handles)
status=get(handles.redhorizontalshift,'enable');
if size(status,2)==2 %on
% RIGHT
handles.thold=get(handles.redhorizontalshift,'String');
kteh=str2num(handles.thold);
kteh=kteh+1;  
set(handles.redhorizontalshift,'String',num2str(kteh));
datashift=get(handles.panel,'userdata');
datashift(2,2)=kteh;
set(handles.panel,'userdata',datashift);

if handles.param.red.nfram==1
    handles.newred.data=shiftimage(handles,datashift(1,2),kteh,handles.movie.red.data,handles.Xdim,2);
    set(handles.finishedpushbutton,'value',1);
    showframeROI(handles)
end

end
guidata(hObject, handles);

%-----------------------------------------------------------------


% --- Executes on button press in greenLbutton.
function greenLbutton_Callback(hObject, eventdata, handles)
status=get(handles.greenhorizontalshift,'enable');
if size(status,2)==2 %on
% LEFT
handles.thold=get(handles.greenhorizontalshift,'String');
kteh=str2num(handles.thold);
kteh=kteh-1;
set(handles.greenhorizontalshift,'String',num2str(kteh));
datashift=get(handles.panel,'userdata');
datashift(2,3)=kteh;
set(handles.panel,'userdata',datashift);

if handles.param.green.nfram==1
    handles.newgreen.data=shiftimage(handles,datashift(1,3),kteh,handles.movie.green.data,handles.Xdim,2);
    set(handles.finishedpushbutton,'value',1);
    showframeROI(handles)
end
end
guidata(hObject, handles);


%-----------------------------------------------------------------

% --- Executes on button press in greenRbutton.
function greenRbutton_Callback(hObject, eventdata, handles)
status=get(handles.greenhorizontalshift,'enable');
if size(status,2)==2 %on
% RIGHT
handles.thold=get(handles.greenhorizontalshift,'String');
kteh=str2num(handles.thold);
kteh=kteh+1;
set(handles.greenhorizontalshift,'String',num2str(kteh));
datashift=get(handles.panel,'userdata');
datashift(2,3)=kteh;
set(handles.panel,'userdata',datashift);

if handles.param.green.nfram==1
    handles.newgreen.data=shiftimage(handles,datashift(1,3),kteh,handles.movie.green.data,handles.Xdim,2);
    set(handles.finishedpushbutton,'value',1);
    showframeROI(handles)
end
end
guidata(hObject, handles);

%-----------------------------------------------------------------

% --- Executes on button press in blueLbutton.
function blueLbutton_Callback(hObject, eventdata, handles)
status=get(handles.bluehorizontalshift,'enable');
if size(status,2)==2 %on
% LEFT
handles.thold=get(handles.bluehorizontalshift,'String');
kteh=str2num(handles.thold);
kteh=kteh-1;
set(handles.bluehorizontalshift,'String',num2str(kteh));
datashift=get(handles.panel,'userdata');
set(handles.panel,'userdata',datashift);

if handles.param.blue.nfram==1
    handles.newblue.data=shiftimage(handles,datashift(1,4),kteh,handles.movie.blue.data,handles.Xdim,2);
    set(handles.finishedpushbutton,'value',1);
    showframeROI(handles)
end
end
guidata(hObject, handles);

%-----------------------------------------------------------------

% --- Executes on button press in blueRbutton.
function blueRbutton_Callback(hObject, eventdata, handles)
status=get(handles.bluehorizontalshift,'enable');
if size(status,2)==2 %on
% RIGHT
handles.thold=get(handles.bluehorizontalshift,'String');
kteh=str2num(handles.thold);
kteh=kteh+1;
set(handles.bluehorizontalshift,'String',num2str(kteh));
datashift=get(handles.panel,'userdata');
datashift(2,4)=kteh;
set(handles.panel,'userdata',datashift);

if handles.param.blue.nfram==1
    handles.newblue.data=shiftimage(handles,datashift(1,4),kteh,handles.movie.blue.data,handles.Xdim,2);
    set(handles.finishedpushbutton,'value',1);
    showframeROI(handles)
end
end
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in finishedpushbutton.
function finishedpushbutton_Callback(hObject, eventdata, handles)

output(1)=str2num(get(handles.trajverticalshift,'string'));
output(2)=str2num(get(handles.trajhorizontalshift,'string'));
output(3)=str2num(get(handles.dicverticalshift,'string'));
output(4)=str2num(get(handles.dichorizontalshift,'string'));
output(5)=str2num(get(handles.redverticalshift,'string'));
output(6)=str2num(get(handles.redhorizontalshift,'string'));
output(7)=str2num(get(handles.greenverticalshift,'string'));
output(8)=str2num(get(handles.greenhorizontalshift,'string'));
output(9)=str2num(get(handles.blueverticalshift,'string'));
output(10)=str2num(get(handles.bluehorizontalshift,'string'));


%if isfield(handles,'newdic')
if isempty(handles.newdic.data)==0
    handles.movie.gray.data=handles.newdic.data;
end
%if isfield(handles,'newred')
if isempty(handles.newred.data)==0
    handles.movie.red.data=handles.newred.data;
end
%if isfield(handles,'newgreen')
if isempty(handles.newgreen.data)==0
    handles.movie.green.data=handles.newgreen.data;
end
%if isfield(handles,'newblue')
if isempty(handles.newblue.data)==0
    handles.movie.blue.data=handles.newblue.data;
end

if isfield(handles,'newtraces')
    handles.traces=handles.newtraces;
else
    handles.newtraces=handles.traces;
end

set(handles.finishedpushbutton,'userdata',output);
save('auxiliar.mat','output','-mat')
%varargout{1} = output;


guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function zoomfgraylow_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function zoomfredlow_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function zoomfgreenlow_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function zoomfbluelow_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function zoomfgrayhigh_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function zoomfredhigh_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function zoomfgreenhigh_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function zoomfbluehigh_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function linewidth_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function dicverticalshift_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function dichorizontalshift_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function redverticalshift_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function redhorizontalshift_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function greenverticalshift_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function greenhorizontalshift_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function blueverticalshift_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function bluehorizontalshift_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function trajverticalshift_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function trajhorizontalshift_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dichorizontalshift_Callback(hObject, eventdata, handles)
valor=get(handles.hObject,'String');
set(handles.dichorizontalshift,'string',valor);
guidata(hObject, handles);

function dicverticalshift_Callback(hObject, eventdata, handles)
valor=get(handles.hObject,'String');
set(handles.dicverticalshift,'string',valor);
guidata(hObject, handles);


function redhorizontalshift_Callback(hObject, eventdata, handles)
valor=get(handles.hObject,'String');
set(handles.redhorizontalshift,'string',valor);
guidata(hObject, handles);

function redverticalshift_Callback(hObject, eventdata, handles)
valor=get(handles.hObject,'String');
set(handles.redverticalshift,'string',valor);
guidata(hObject, handles);

function greenverticalshift_Callback(hObject, eventdata, handles)
valor=get(handles.hObject,'String');
set(handles.greenverticalshift,'string',valor);
guidata(hObject, handles);

function greenhorizontalshift_Callback(hObject, eventdata, handles)
valor=get(handles.hObject,'String');
set(handles.greenhorizontalshift,'string',valor);
guidata(hObject, handles);

function blueverticalshift_Callback(hObject, eventdata, handles)
valor=get(handles.hObject,'String');
set(handles.blueverticalshift,'string',valor);
guidata(hObject, handles);

function bluehorizontalshift_Callback(hObject, eventdata, handles)
valor=get(handles.hObject,'String');
set(handles.bluehorizontalshift,'string',valor);
guidata(hObject, handles);

function trajverticalshift_Callback(hObject, eventdata, handles)
valor=get(handles.hObject,'String');
set(handles.trajverticalshift,'string',valor);
guidata(hObject, handles);

function trajhorizontalshift_Callback(hObject, eventdata, handles)
valor=get(handles.hObject,'String');
set(handles.trajhorizontalshift,'string',valor);
guidata(hObject, handles);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotall(actualtrc, handles)

localiz=0;
blinking=0;
rainbow=0;
timecolor=0;
firstbutton=0;
linew=str2num(get(handles.linewidth,'string'));
if handles.typetraj==0
            elseif handles.typetraj==1 | handles.typetraj==6
                localiz=1;
            elseif handles.typetraj==2
                blinking=1;
            elseif handles.typetraj==3
                rainbow=1;
            elseif handles.typetraj==4
                timecolor=1;
            elseif handles.typetraj==7
                firstbutton=1;
                rainbow=1;
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
               if firstbutton==1 % primer punto
                   if rainbow==1
                      handles.typetraj=3;
                      codecol=colorm(indcol,1:3);
                      indcol=indcol+1;
                      if indcol>50
                         indcol=1;
                      end
                      plot(actualtrc(indexmol,3),actualtrc(indexmol,4),'Color',codecol,'LineWidth',linew);
                      hold on
                   end
                  plot(actualtrc(indexmol(1),3),actualtrc(indexmol(1),4),'Color',codecol,'Marker','.','MarkerSize',20);
               end
         %  if identify==1 % numero tray
         %%          pos=indexmol(size(indexmol,1));
          %        text(actualtrc(pos,3)+1,actualtrc(pos,4)+1,sprintf('%0.0f',actualtrc(pos,1)),'Color',[1 1 0]);
          %     end
               hold on
            else
            end % empty indexmol
      end % loop actual molecules

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


