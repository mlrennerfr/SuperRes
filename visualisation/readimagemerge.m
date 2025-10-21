function [datamatrix,param,stktrue,file,path]=readimagemerge(colorname,handles)
% function [datamatrix,param,stktrue,file,path]=readimagemerge(colorname,handles)
% Reads images for merging
% Marianne Renner 09/09 for SPTrack v4
% Marianne Renner 01/2025 - adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datamatrix=[];
param.Xdim=0;
param.Ydim=0;
param.nfram=0;
stktrue=0;

% loads image 
currentcd=cd;
if length(dir(handles.presentfolder))>0
   cd(handles.presentfolder)
end

[file,path] = uigetfile('*.*',['Load ',colorname,' file (.spe, .stk or .tif)']);
filename = [path,file];
handles.presentfolder=path;
cd(currentcd)
if filename==0
    file=[];
    path=currentcd;
else
  [datamatrix,param.Xdim,param.Ydim,param.nfram,stktrue]=readimage(filename);
end
%param=[Xdim, Ydim, nfram];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%