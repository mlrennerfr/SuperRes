function [select]=pickpointsfast(BW,xi,yi,data,maxx,maxy,handles);
%function [select]=pickpointsfast(BW,xi,yi,data,maxx,maxy,handles);
%
% selects tracking points inside a region - fast version
%
% MR fév 08 - SPTrack pack                 MatLab 7
% Marianne Renner 01/2025 - adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3
   data=get(handles.filetrc,'userdata');
   maxx=handles.param(1);  
   maxy=handles.param(2);
end
%newselectrace=[];
%selectrace=[];
%count=1;

% preseleccion trc
index1=find(data(:,3)>min(xi) & data(:,3)<max(xi));
aux=data(index1,:);
index2=find(aux(:,4)>min(yi) & aux(:,4)<max(yi));
select=aux(index2,:);
clear aux index1 index2 

%for it=1:size(data2,1)
%    ty=round(data2(it,3)); % x
%    tx=round(data2(it,4)); %y
%    ojo las di vuelta
%    if tx>maxy;          tx=maxy;      end
%    if ty>maxx;         ty=maxx;     end
%    if tx<1;          tx=1;      end
%    if ty<1;         ty=1;     end
%    if BW(tx,ty)>0          
%        newselectrace(count,:)=data2(it,:);   % puntos dentro de la ROI
%        count=count+1;
%    end
%end

clear data BW 

%eof