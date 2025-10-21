function [options]=calibratePALM (handles)
% function [options] = calibratePALM (handles)
% displays SubImage and detects peaks
% allows changing parameters
% in blue, all the detected peaks,
% in red, the peaks left after cutoffs
% 
% Marianne Renner - avril 09 for SPTrack_v4.m   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

listafiles=get(handles.moviefile,'userdata');


file=listafiles{1};
disp(' ')
disp(['Calibration on file ',file])
disp(' ')

stack=[];
info=imfinfo(file);
if length(info)<2
    Nb_image_ds_stack=floor(info.FileSize/info.StripByteCounts);
else
    Nb_image_ds_stack=length(info);
end

detoptions.file=file;
detoptions.stack=stack;
detoptions.lastframe=Nb_image_ds_stack;

detoptions.seuil_alpha = str2num(get(handles.alphavalue,'String'));
detoptions.seuil_detec_1vue = str2num(get(handles.threshold,'String'));
detoptions.wn = str2num(get(handles.windsize, 'String'));
detoptions.r0 = str2num(get(handles.gaussrad, 'String'));
detoptions.nb_defl = str2num(get(handles.defloops,'String'));
options(1)=detoptions.seuil_alpha  ;
options(2)=detoptions.seuil_detec_1vue;         % threshold
options(3)=detoptions.wn;                       % window size
options(4)=detoptions.r0;                       % gaussian radius
options(5)=detoptions.nb_defl;                  % number deflation loops

%gui
varargout=calibrationSuperRes(detoptions);
uiwait;

%new values
pathdet='detecoptions.mat';
det=load(pathdet);
detopt = struct2cell(det);
detoptions=detopt{1};

options(1)=detoptions.seuil_alpha  ;
options(2)=detoptions.seuil_detec_1vue;         % threshold
options(3)=detoptions.wn;                       % window size
options(4)=detoptions.r0;                       % gaussian radius
options(5)=detoptions.nb_defl;                  % number deflation loops
   
%clear detoptions;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% end of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

