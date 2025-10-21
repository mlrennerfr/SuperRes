function varargout = DetectSM(varargin)
%DETECTSM M-file for DetectSM.fig
% Last Modified by GUIDE v2.5 20-Jun-2017 09:49:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DetectSM_OpeningFcn, ...
                   'gui_OutputFcn',  @DetectSM_OutputFcn, ...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes just before DetectSM is made visible.
function DetectSM_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
handles.folder=cd;
%set(handles.datafolder,'string',handles.folder);
alpha=str2num(get(handles.alphavalue,'string'));
set(handles.slider3,'value',alpha/2000); %%%%%%%%!!!!! at max alpha!!!!!!!!!!!

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù
function varargout = DetectSM_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function slider3_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function slider3_Callback(hObject, eventdata, handles)

%alpha_value = str2num(get(handles.alphavalue,'string'));

slidervalue=get(hObject,'Value');
minvalue=get(hObject,'Min');
maxvalue=get(hObject,'Max');
prop=slidervalue/(minvalue+maxvalue);

%!!!!! AT MAX!!!!!!!!!!!

alpha=round(prop*2000);
step=(maxvalue-minvalue)/2000;
set(handles.slider3,'SliderStep',[step step]);

if alpha<1
    alpha=1;
end
if alpha>2000
    alpha=2000;
    %%%%% ATENTION!!!!!!!!!!!!!!!!!
end

%set(handles.slider3, 'value', alpha);
set(handles.alphavalue, 'string', num2str(alpha));
guidata(gcbo,handles) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function loadfiles_Callback(hObject, eventdata, handles)

% folder
path=cd;

%set(handles.datafolder,'string',path);
%handles.folder=get(handles.datafolder,'string');
%handles.folder=[handles.folder,'\'];

st=[]; 

%files
d=dir('*tif*'); % .stk files
st={d.name};
if isempty(st)==1
    d=dir('*nd2*'); % nikon files
    st={d.name};
    if isempty(st)==1
        msgbox(['No files!!'],'','error');
        return
    end
end

%choose data
[files,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',st);
if v==0
     return
end

[f,ultimo]=size(files);
for i=1:ultimo
      listafiles{i}=st{files(i)};
end

handles.file=listafiles{1};
set(handles.moviefile,'userdata',listafiles);

handles.listafiles=get(handles.moviefile,'userdata');  %selected files

filename=handles.file; %first file
if ultimo>1 %batch
  set (handles.moviefile, 'string',['Batch: File ',filename,' (1/',num2str(ultimo),')']) ;
else
  set (handles.moviefile, 'string',handles.file) ;
end

%pushbuttons & radiobuttons
set (handles.calibratepushbutton, 'Enable','on');
set (handles.MTTnotracking, 'Enable','on');

guidata(gcbo,handles) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function calibratepushbutton_Callback(hObject, eventdata, handles)


[options] =calibratePALM (handles);

set(handles.alphavalue,'String',options(1));
set(handles.slider3,'value',options(1)/2000); %%%%%%%%!!!!! at max alpha!!!!!!!!!!!

set(handles.threshold,'String',options(2));
set(handles.windsize,'String',options(3));
set(handles.gaussrad,'String',options(4));
set(handles.defloops,'String',options(5));

guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MTTnotracking_Callback(hObject, eventdata, handles)

currentdir=cd;
add3D=get(handles.add3Dradiobutton,'value');
szpx=str2num(get(handles.pixel,'string'));
if add3D==1
    %ask for calibration curve
   % [file,path] = uigetfile('*calibration3D.mat','Load calibration curve');
    [file,path] = uigetfile('*curvecalib*','Load calibration curve');
    if isempty(file)==1
        return
    end
    cd(path)
    fitresult=load(file); %a: fitresult(1,4) b: fitresult(2,4)
    %fitresult=load(file,'-mat')
    %load(file,'-mat')
end
cd(currentdir)

% loop analysis
handles.listafiles=get(handles.moviefile,'userdata');
%currentdir=cd
%isdir([currentdir,'\pk'])

if isdir([currentdir,'\pk']);else mkdir([currentdir,'\pk']);end


for nromovie=1:size(handles.listafiles,2)
    name_stk=handles.listafiles{nromovie};
    [namefile,rem]=strtok(name_stk,'.'); 
    disp(' ')
    disp('Detection of single molecules using MTT software')
    disp(' ')

    disp(['File ',name_stk]);
    str=['Batch: File ',name_stk,' (',num2str(nromovie),'/',num2str(size(handles.listafiles,2)),')'];
    set (handles.moviefile, 'string', str);

    seuil_alpha = str2num(get(handles.alphavalue,'String'));
    seuil_detec_1vue = str2num(get(handles.threshold,'String'));
    wn = str2num(get(handles.windsize, 'String'));
    r0 = str2num(get(handles.gaussrad, 'String'));
    nb_defl = str2num(get(handles.defloops,'String'));
    activation_sig_fit=get(handles.gaussfitradiobutton, 'value');
    
    info=imfinfo(name_stk);

  %  k=strfind(name_stk,'nd2');
  %  if k>0
   %     data = bfopen(name_stk);
   %     stack = data{1, 1};
   %     Nb_image_ds_stack = size(stack, 1);
   % else
        stack=[];
        %Nb_image_ds_stack =length(info);
   % end
    
    Nb_image_ds_stack =length(info);
    output_dir=[name_stk, '_', 'output'];    
    
    %[fr,x,y,alpha,radius, radiusmean,sommevalid] = detectionMTT2(stack,Nb_image_ds_stack, seuil_detec_1vue, wn, r0, nb_defl, seuil_alpha,activation_sig_fit);  
    
    [fr,x,y,alpha,radius, radiusmean,sigma, blink,sommevalid] = detectionMTTSR(name_stk , stack, Nb_image_ds_stack, seuil_detec_1vue, wn,r0, nb_defl, seuil_alpha, activation_sig_fit);
    
     disp(['Mean number of valid particles: ',num2str(sommevalid/Nb_image_ds_stack)]);
     %disp(matrice_results)
     
     % clean zero values
     %% Remove points at x=y=0
     aux1 = find(x==0);
     aux2 = find(y==0);
     vect2remove = intersect(aux1,aux2);
     if ~isempty(vect2remove)
         x(vect2remove) = [];
         y(vect2remove) = [];
         alpha(vect2remove) = [];
         fr(vect2remove) = [];
         radius(vect2remove) = [];
         sigma(vect2remove) = [];
         blink(vect2remove) = [];
     end
     
     clear pk
     pk(:,1)=fr; 
     pk(:,2)=x;  
     pk(:,3)=y;    
     pk(:,4)=radius*2; 
     pk(:,5)=alpha; 
     pk(:,6)=blink;
     pk(:,7)=sigma; 
     pk(:,8)= zeros(size(fr,1)); %ratio widths
     pk(:,9)= zeros(size(fr,1)); %z
     pk(:,10)= zeros(size(fr,1)); %test?
     pk(:,11)= zeros(size(fr,1)); %test?

     %3D?
     if add3D==1
         if isdir(['pk3']); else; mkdir(['pk3']); end %VOIR

         disp('Adding z coordinate...')  
         %waitbarhandle3D=waitbar( 0,'Please wait...','Name',['Detecting peaks in 3D on ',file]);
         %  pk=go3D(file,pk,Image,ImagePar,curve,szpx)
         waitbarhandle=waitbar( 0,'Please wait...','Name','Adding z coordinate') ;

        [ImagePar,Image] = stkdataread(name_stk);
      %   pk=go3D(pk,stack,Nb_image_ds_stack,curve,wn);
      
         [datawidth,pk]=go3D(pk,Image,Nb_image_ds_stack,fitresult,wn,seuil_alpha,szpx,waitbarhandle);
         clear Image ImagePar
         close(waitbarhandle);
         
         % prepare file for visp
         %pkvisp=[pk(:,2) pk(:,3) pk(:,8) pk(:,5) pk(:,1)];         
         %[x y z alpha fr];
        % save([namefile, '.3d'], 'pkvisp','-ascii');
         %save([namefile, '-datawidth.txt'], 'datawidth','-ascii');
         
       %  if isempty(totalpk)==0
       %      index=find(totalpk(:,17)<1000 & totalpk(:,17)~=NaN);
       %      save (['pk3\',namefile,'.pk3'], 'totalpk','-ascii');
       %      save (['pk3\',namefile,'.cal'], 'datawidth','-ascii');
     %    end
     
     
         if isempty(pk)==0
            % totalpk=pk(find(pk(:,8)<1000),:);  %ATT limites z d'apres courbe calib
             indexpk=find(pk(:,8)<1000 & pk(:,10)<0.99);  %ATT limites z d'apres courbe calib
              % totalpk1=pk(indexpk,:);
           
           %  indexpk=find(totalpk1(:,10)<0.99);  %ATT limites z d'apres courbe calib
          %   totalpk=totalpk1(indexpk,:);
             
            % indexpk2=find(pk(indexpk,10)<1);  %ATT limites z d'apres courbe calib
             totalpk=pk(indexpk,:);
             
           %  save (['pk3\',namefile,'.pk'], 'pk','-ascii');
            % save (['pk3\',namefile,'.pk31'], 'totalpk1','-ascii');
             save (['pk3\',namefile,'.pk3'], 'totalpk','-ascii');
             save (['pk3\',namefile,'.cal'], 'datawidth','-ascii');
         end

     else %no 3D
         %z=zeros(size(x,2),1);
         totalpk=pk;
         totalpk(:,8)=zeros(size(totalpk,1),1);
         totalpk(:,9)=zeros(size(totalpk,1),1);
         totalpk(:,10)=zeros(size(totalpk,1),1);
         totalpk(:,11)=zeros(size(totalpk,1),1);
         
         path=[currentdir,'/pk'];
         cd(path);
         save([namefile, '.pk'], 'pk','-ascii');
         cd(currentdir)

     end

     clear matrice_results
     % matrice_results(1,:)=fr; %pk(:,1)
     %  matrice_results(2,:)=x;   %pk(:,2);
     %  matrice_results(3,:)=y;    %pk(:,3);
     %  matrice_results(4,:)=alpha;  %pk(:,5);
     %  matrice_results(5,:)=radius; %pk(:,4)/2;
     %  matrice_results(6,:)=sigma;  %pk(:,7);
     %  matrice_results(7,:)=blink; %pk(:,6);
     %   matrice_results(8,:)=pk(:,8); %ratio
     %   matrice_results(9,:)=pk(:,9); %z
     %   matrice_results(10,:)=pk(:,10); %test
     %   matrice_results(11,:)=pk(:,11); %test
     
     matrice_results(1,:)=totalpk(:,1);
     matrice_results(2,:)=totalpk(:,2);
     matrice_results(3,:)=totalpk(:,3);
     matrice_results(4,:)=totalpk(:,5);
     matrice_results(5,:)=totalpk(:,4)/2;
     matrice_results(6,:)=totalpk(:,7);
     matrice_results(7,:)=totalpk(:,6);
     matrice_results(8,:)=totalpk(:,8); %ratio
     matrice_results(9,:)=totalpk(:,9); %z
     matrice_results(10,:)=totalpk(:,10); %test
     matrice_results(11,:)=totalpk(:,11); %test
     
     %name=[namefile,'_','peaks'];
     
     save([namefile, '.mat'], 'matrice_results');

     if add3D==1
         disp(['Results saved in ', namefile,'.mat and ', namefile,'.pk3']);
     else
         disp(['Results saved in ', namefile,'.mat and ', namefile,'.pk']);
     end
     
     disp(' ')
    clear global matrice_results pk totalpk
end

if activation_sig_fit & radiusmean~=0
    set(handles.gaussrad, 'String', num2str(radiusmean))
elseif activation_sig_fit==0
    set(handles.gaussrad, 'String', num2str(r0))
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tracking
function trackpushbutton_Callback(hObject, eventdata, handles)

buildTracksSM


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Menus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------------
function detection_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------

function peaktest_Callback(hObject, eventdata, handles)
testspeakPALM

% --------------------------------------------------------------------
function fuse_Callback(hObject, eventdata, handles)
    
% open pontillistic data of two files, fuse them (add frame number for the
% second)

[x,y, fr, alpha, radius,sigma, blink, ffname1,filename1]=readPALMdata2(handles);

[x2,y2, fr2, alpha2, radius2,sigma2, blink2,ffname2,filename2]=readPALMdata2(handles);

disp(max(x))
disp(max(x2))

qstring=['Fuse files ',filename1,' and ',filename2,'? (only 2D data)'];
button = questdlg(qstring); 
if strcmp(button,'Yes')
    x=[x; x2];
    y=[y; y2];
    lastfr=max(fr);
    fr=[fr; fr2+lastfr];
    alpha=[alpha; alpha2];
    radius=[radius; radius2];
    sigma=[sigma; sigma2];
    blink=[blink; blink2];
    if isfield(handles,'Xmatrix')
        Xmatrix = x1;
        Ymatrix = y1;
        alphamatrix = alpha1;
        frmatrix = fr1;
        uisave({'x','y','alpha','fr','Xmatrix','Ymatrix','alphamatrix','frmatrix'},filename1);
    else
        uisave({'x','y','alpha','fr','radius','sigma','blink'},filename1);
    end
    disp('File saved');
end


% --------------------------------------------------------------------
function corrstage_Callback(hObject, eventdata, handles)

CorrectDrift

% --------------------------------------------------------------------
function corrcolor_Callback(hObject, eventdata, handles)

CorrectColorShift 

% --------------------------------------------------------------------
function corrmultiple_Callback(hObject, eventdata, handles)

% dialog box to enter recognition criteria for images
prompt = {'Max. distance for multiple detections (nm):','Pixel size (nm):','Max. blinking time (frames):'};
num_lines= 1;
dlg_title = 'Correction of multiple detections';
def = {'50','160','350'};
answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
if exit(1) == 0;
    return;
end
dist=str2num(answer{1}); 
szpx=str2num(answer{2}); 
blinktime=str2num(answer{3});
%dist=dist/str2num(answer{2}); %px

%szpx=szpx/1000; %µm

%files
d=dir('*.mat*'); % .stk files
st={d.name};
if isempty(st)==1
    msgbox(['No files!!'],'','error');
    return
end

%choose data
[listafiles,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',st);
if v==0
     return
end

for nromovie=1:size(listafiles,2)
    matfile=st{listafiles(nromovie)} 
   disp=['Correction of multiple detections: File ',matfile];
   disp(' ');
    [filename,rem]=strtok(matfile,'.');
    correctedfile=[filename,'-corr.mat']  
    
    aux=correctdetections3(matfile,dist,szpx,blinktime);
    
    corrpeak=aux(find(isnan(aux(:,2))==0),:);
    clear aux   
    
    current=cd;
    if isdir('pk');else mkdir('pk');end
    cd('pk')
    save([filename, '-corr.pk'], 'corrpeak','-ascii')
    cd(current)
    
    % save new data
    if isempty(corrpeak)==0
        matrice_results(1,:)=corrpeak(:,1);
        matrice_results(2,:)=corrpeak(:,3);
        matrice_results(3,:)=corrpeak(:,2);
        matrice_results(4,:)=corrpeak(:,4);
        matrice_results(5,:)=corrpeak(:,5);
        matrice_results(6,:)=corrpeak(:,7);
        matrice_results(7,:)=corrpeak(:,6);
        
        savename=[filename, '-corr.mat'];
        save(savename, 'matrice_results'); % 
        clear matrice_results
    end
end %loop

% --------------------------------------------------------------------
function docurve_Callback(hObject, eventdata, handles)

curvecalibration3D

% --------------------------------------------------------------------
function groupcalib_Callback(hObject, eventdata, handles)

groupcurves2

% --------------------------------------------------------------------
% --------------------------------------------------------------------
function trajectories_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------
% --------------------------------------------------------------------
function MTTresults_Callback(hObject, eventdata, handles)
% convert MTT files

convertMTTtrcnotrack

% --------------------------------------------------------------------
function manrec_Callback(hObject, eventdata, handles)

 
%------------------------------------------------------------------------
% --- Executes on button press in autorec.
function autorec_Callback(hObject, eventdata, handles)
% automatic reconnection

%maxblink, distmax, minpoints
% dialog box to enter parameters 2
%prompt = {'Max blinking time (frames):','Max distance during blinking (pixels):','Minimum length of trajectories:'};
%num_lines= 1;
%dlg_title = 'Automatic Reconnection (SPTrack routine)';
%def = {'3','1','3'}; % default values
%answer  = inputdlg(prompt,dlg_title,num_lines,def);
%exit=size(answer);
%if exit(1) == 0;
%       return; 
%end
%maxblink = str2num(answer{1});%#ok<ST2NM>
%distmax = str2num(answer{2});%#ok<ST2NM>   %blinking
%minpoints = str2num(answer{3});

% dialog box to enter parameters 2
prompt = {'Max blinking time (frames):','Max distance during blinking (pixels):','Maximum distance for tlag=1 (pixels):','Minimum length of trajectories:'};
num_lines= 1;
dlg_title = 'Automatic Reconnection (SPTrack routine)';
def = {'3','1','NaN','3'}; % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
if exit(1) == 0;
       return; 
end

blinktime = str2num(answer{1});%#ok<ST2NM>
blinkdist = str2num(answer{2});%#ok<ST2NM>   %blinking
limdist = str2num(answer{3});
mintrace = str2num(answer{4});


% files
currentdir=cd;
path=[cd,'\trc'];
cd(path);
controlf=1;

d=dir('*trc*');
st = {d.name};

if isempty(st)==1
     msgbox(['No trajectory files to reconnect!!'],'','error');
     cd(currentdir)
     return
end
%choose data
[listafiles,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',st);
if v==0
    cd(currentdir)
     return
end
[f,ultimo]=size(listafiles);


for cont=1:ultimo   
    trcfile=st{listafiles(cont)}; disp(['Reconnecting trajectories in file ',trcfile]); disp(' '); % con extension
    [namefile,rem]=strtok(trcfile,'.');
    pause(0.001) % to show waitbar

   % rectrc=elongatetrack(trcfile, detoptions(19), detoptions(20), detoptions(21), handles);  
%file, maxblink, distmax, minpoints
    %[namefile,rem]=strtok(handles.file,'.');
    
    
        %reconnection
    
    trc=load(trcfile);
    reconnectedfile=trc;

    
    if size(blinktime,2)>1;            
            %series of reconnections
            number=size(blinktime,2);
               series=1;   
               while series<number
                     blink=blinktime(series);
                     distmax=blinkdist(series);
                     disp('  ');
                     disp('Reconnecting trajectories...');
                     waitbarhandle=waitbar( 0,'Please wait...','Name',['Reconnecting trajectories']);
                     reconnectedfile=reconnectfasttrc(reconnectedfile,blink,distmax,3,limdist,waitbarhandle);
                     series=series+1;
                     close(waitbarhandle);
               end
               % last one: uses minpoints from the window
               blink=blinktime(series);
               distmax=blinkdist(series);
               disp('  ');
               disp('Reconnecting trajectories...');
               waitbarhandle=waitbar( 0,'Please wait...','Name',['Reconnecting trajectories']);
               reconnectedfile=reconnectfasttrc(reconnectedfile,blink,distmax,mintrace,limdist,waitbarhandle);
               close(waitbarhandle);
     else
            %one reconnection
            disp('  ');
            disp('Reconnecting trajectories...');
            waitbarhandle=waitbar( 0,'Please wait...','Name',['Reconnecting trajectories']);
            reconnectedfile=reconnectfasttrc(trc,blinktime,blinkdist,mintrace,limdist,waitbarhandle);
            close(waitbarhandle);
     end
    
    if isempty(reconnectedfile)==0 
        %save everything
        save([namefile,'.con.trc'],'reconnectedfile','-ascii');
    else
        disp('Empty .con.trc file');
    end

end % loop files

clear rectrc
cd(currentdir); % comes back

guidata(gcbo,handles) ;


% --------------------------------------------------------------------
function clean_Callback(hObject, eventdata, handles)
%Clean trajectories
cleantrcSM


% --------------------------------------------------------------------
function locdomains_Callback(hObject, eventdata, handles)

onlylocalizationSM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function add3Dpushbutton_Callback(hObject, eventdata, handles)

add3Dpk(handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------------
function diffusion_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------

function calculdiffindiv_Callback(hObject, eventdata, handles)

trajindivSM

% --------------------------------------------------------------------
function summarydiff_Callback(hObject, eventdata, handles)

compileindivSM

% --------------------------------------------------------------------

function stats_Callback(hObject, eventdata, handles)

analizeresults
% --------------------------------------------------------------------

%function helpsm_Callback(hObject, eventdata, handles)

% helpPALM



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y, fr, alpha, radius,sigma,blink,ffname,filename]=readPALMdata2(handles)


[filename, pathname, filterindex] = uigetfile('*.mat','Select .mat file containing detection results');
if filename==0
    return
end

ffname = fullfile(pathname,filename);
[pathstr, name, ext] = fileparts(ffname);
%set(handles.selectfilepushbutton,'FontWeight','normal');

    if strcmp(ext,'.mat')
        S = load(ffname);
        if isfield(S,'matrice_results')
            %disp('File without tracking information');
            
            aux = S.matrice_results;
            clear S;
            fr = aux(1,:); fr = fr(:);
            
           % disp(aux(3,:))
           % disp(px_mu)
            
          %  x = aux(3,:) * px_mu;
          %  x = x(:);
           % y = aux(2,:) * px_mu; 
           % y = y(:);
           
            x = aux(3,:);
            x = x(:);
            y = aux(2,:); 
            y = y(:);
            alpha = aux(4,:); alpha = alpha(:);
            if size(aux,1)>4
                radius=(aux(5,:));
            else
                radius=zeros(size(x,2));
            end
            
            sigma= aux(6,:);
            blink= aux(7,:);
           
        elseif isfield(S,'alpha') && isfield(S,'fr') && isfield(S,'x') && isfield(S,'y')
            
           % disp('coucou2')
            
            x = S.x; % * px_mu;
            x = x(:);
            y = S.y; % * px_mu;
            y = y(:);
            alpha = S.alpha; alpha = alpha(:);
            fr = S.fr; fr = fr(:);
            if isfield(S,'radius')
                radius=S.radius;
            else
                radius=zeros(size(x,2));
            end
            %validfile = 1;
        else
            %validfile = 0;
            %             errordlg('This file is not valid !','Warning','modal');
            disp('Invalid file!');
            return
        end
   end
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pixel_Callback(hObject, eventdata, handles)
function pixel_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Te_Callback(hObject, eventdata, handles)
function Te_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function threshold_Callback(hObject, eventdata, handles)
function threshold_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function windsize_Callback(hObject, eventdata, handles)
function windsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function gaussrad_Callback(hObject, eventdata, handles)
function gaussrad_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function defloops_Callback(hObject, eventdata, handles)
function defloops_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function alphavalue_Callback(hObject, eventdata, handles)
alpha=str2num(get(hObject, 'string'));
set(handles.slider3,'value',alpha/2000); %%%%%%%%!!!!! at max alpha!!!!!!!!!!!
guidata(gcbo,handles) ;

function alphavalue_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function add3Dradiobutton_Callback(hObject, eventdata, handles)

function gaussfitradiobutton_Callback(hObject, eventdata, handles)
sig_fit=get(hObject, 'value');
if sig_fit==1;
    set(handles.gaussrad,'enable','off');
else
    set(handles.gaussrad,'enable','on');
end

function moviefile_Callback(hObject, eventdata, handles)
function moviefile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
