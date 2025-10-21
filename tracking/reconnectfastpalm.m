function [res]=reconnectfastpalm(path,file,maxblink,distmax,mintrace,handles,blink, waitbarhandle)
% function [res]=reconnectfastpalm(path,file,maxblink,distmax,mintrace,handles,blink, waitbarhandle)
% connects trajectories between the limits of maxblink and maxdistance
% path: data folder
% file: .trc file
% maxblink : maximum off time (frames)
% distmax : maximum distance between the trajectories (end point of the first one to start
% point of the second one)
% mintrace : minimum number of point to keep the trajectory
% blink: max npeiod of blinking allowed
% waitbarhandle: to actualize wait bar
%
% Marianne Renner mar 2011 -  for SPTpalm v1
% Marianne Renner 01/2025 - adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load file
[namefile,rem]=strtok(file,'.');
if handles.usecontrc==0
   trcfile=[path,'\trc\',namefile,'.trc'];
elseif handles.usecontrc==1
   trcfile=[path,'\trc\',namefile,'.con.trc'];
end
if length(dir(trcfile))>0		
      trc=load(trcfile);
      disp(['File ' ,trcfile, ' loaded.']);
   else
      disp(['Couldn''t find .trc file ',trcfile]);
      res=[];
      return
end
if isempty(trc)==1 
      res=[];
      disp('Empty .trc file');
      return
end
newtrc=trc(:,1:5);    % no toma en cuanta la loc!!!!!
maxspots=max(trc(:,1));
totalspot=maxspots;

if maxblink>0 & distmax>0 % =0 mock

% vectors with first and last points
startpoint(1,:)=trc(1,:); %first point first spot
endpoint(maxspots,:)=trc(size(trc,1),:); % last point last spot
for i=2:maxspots
    index=find(trc(:,1)==i);
    if isempty(index)==0
        endpoint(i-1,:)=trc(index(1)-1,:); % last point previous spot
        startpoint(i,:)=trc(index(1),:); % first point actual spot
    end
end

% loops reconnection

reconnect=1;
count=1;
control=0;
reconnect=1;
iteration=1;
%nroiter=1;

while reconnect==1 & iteration<500 % iteration while it is possible do reconnections or number of iterations < 500
  control=0; 
  %nroiter=nroiter+1;
  for i=1:maxspots
       candidats=[];
       candidats=find(startpoint(:,2)>endpoint(i,2) & startpoint(:,2)< (endpoint(i,2)+ maxblink) ); % all points after up to maxblink
       if isempty(candidats)==0
           dist=[];
           for j=1:size(candidats,1)
               % vector of distances (space and time)
               dist(j,1)=sqrt((startpoint(candidats(j),3)-endpoint(i,3))^2+(startpoint(candidats(j),4)-endpoint(i,4))^2);
               dist(j,2)=startpoint(candidats(j),2)-endpoint(i,2); % dif time
               dist(j,3)=j; % order
           end
           gooddist=[];
           gooddist=find(dist(:,1)<distmax); % at good distance
           if isempty(gooddist)==0
              goodlista=dist(gooddist,:); 
              goodlista=sortrows(goodlista,2); % in ascending order by difference in time
              index=[];
              %trajectory to incorporate: the one that is nearer in space and time
              index=find(newtrc(:,1)==startpoint(candidats(goodlista(1,3)),1)); 
              if isempty(index)==0
                 newtrc(index,1)=i; % renumber the trajectory to incorporate
                 control=1;
              end
           end
       end
  end % loop spots
  
  if control==0 % no new reconnection
      break
  end
  
  % renumbers trajectories without missing number
  newtrc=sortrows(newtrc,1); % sort by trajectory number
  maxspots=max(newtrc(:,1));
  count=1;
  for i=1:maxspots
      index=[];
      index=find(newtrc(:,1)==i);
      if isempty(index)==0
         newtrc(index(:),1)=count;
         count=count+1;
      end
  end
  % vectors with first and last points
  maxspots=max(newtrc(:,1));
  startpoint=[];
  endpoint=[];
  startpoint(1,:)=newtrc(1,:); %first point first spot
  endpoint(maxspots,:)=newtrc(size(newtrc,1),:); % last point last spot
  for i=2:maxspots
    index=find(newtrc(:,1)==i);
    endpoint(i-1,:)=newtrc(index(1)-1,:); % last point previous spot
    startpoint(i,:)=newtrc(index(1),:); % first point actual spot
  end

  % actualizes waitbar
  if exist('waitbarhandle')
     waitbar(iteration/maxspots,waitbarhandle,['Iteration # ',num2str(iteration)]);
  end
  iteration=iteration+1;

end % while reconnect

else % mock
    
    newtrc=trc;
end
          
%renumbers trajectories and keeps the ones with more than mintrace points 
% and with blinking periods of less than blink number of points
newmaxspots=max(newtrc(:,1));
index=[];
first=newtrc(1,1);
nro=1;
finaltrc=[];
cuenta=1;
for i=first:newmaxspots
    index=find(newtrc(:,1)==i);
    if isempty(index)==0
        aux=[];
        cuenta=cuenta+1;
        %countblink=0;
        %disp(mintrace)
        nropoints=size(index,1);
        
       if nropoints>mintrace                % mintrace
          % j=2;
          % while j<size(index,1)+1
          %     countblink=newtrc(index(j),2)-newtrc(index(j-1),2);
          %     if countblink > blink
          %         j=size(index,1)+1;
                   %countblink=0;
          %         break
          %     end
          %     j=j+1;
          % end
          % if countblink < blink + 1                %blink
               aux(:,:)=newtrc(index,:);
               aux(:,1)=nro;
               finaltrc=[ finaltrc ; aux];
               nro=nro+1;
          % end
       end
    end
end

   

%results
res=finaltrc;

if isempty(res)==0
   tracefin=max(finaltrc(:,1));
   disp([Num2str(totalspot) ' trajectories at the beginning.']);
   disp([Num2str(cuenta-1) ' trajectories after reconnection.']);
   %disp([Num2str(tracefin) ' trajectories after filtering the short ones (less than ', Num2str(mintrace), ' points and less than ', Num2str(blink), ' frames of blinking).']);
   disp([Num2str(tracefin) ' trajectories after filtering the short ones (less than ', Num2str(mintrace), ' points.']);
   disp(' ');
   %report
   text=[Num2str(cuenta) ' trajectories after reconnection and ',Num2str(tracefin) ' trajectories after filtering the short ones (less than ' Num2str(mintrace) ' points).'];
   updatereport(handles,text)
else
   disp('No trajectory left with enough number of points.');
   %report
   text=['No trajectory left with enough number of points.'];
   updatereport(handles,text,1)
   res=[];
end

clear newtrc finaltrc trc

% end of file