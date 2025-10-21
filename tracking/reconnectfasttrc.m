function [res]=reconnectfasttrc(trc,maxblink,distmax,mintrace,limdist,waitbarhandle)
% function [res]=reconnectfasttrc(trc,maxblink,distmax,mintrace)
% connects trajectories between the limits of maxblink and maxdistance
% maxblink : maximum off time (frames)
% distmax : maximum distance between the trajectories (end point of the first one to start
% point of the second one)
% mintrace : minimum number of point to keep the trajectory
% limdist : maximum distance allowed for a displacement of one tlag
%
% Marianne Renner oct 2007 -  for SPTrack.m
% MR august 09 for SPTrack v4
% MR dec 2012 for BuildTracks.m
% Marianne Renner 01/2025 - adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
          
%control distance and renumbering trajectories, it keeps the ones with more than mintrace points    
newmaxspots=max(newtrc(:,1));
index=[];
first=newtrc(1,1);
nro=1;
finaltrc=[];
cuenta=1;
for i=first:newmaxspots
    index=find(newtrc(:,1)==i);
    if isempty(index)==0
        %aux=newtrc(index(1),:);
        aux=newtrc(index,:);
        cuenta=cuenta+1;
        %for j=2:size(index,1)
        %    if newtrc(index(j),2)-newtrc(index(j-1),2)==1 %one tlag
        %        dist=sqrt((newtrc(index(j),3)-newtrc(index(j-1),3))^2+(newtrc(index(j),4)-newtrc(index(j-1),4))^2)
        %        if dist<limdist
        %            aux=[aux; newtrc(index(j),:)];
        %        end
        %    else
        %        aux=[aux; newtrc(index(j),:)];
        %    end
        %end
       if size(aux,1)>mintrace
          aux(:,1)=nro;            
          finaltrc=[ finaltrc ; aux];
          nro=nro+1;
       end
    end
end

   

%results
res=finaltrc;
if isempty(res)==0
   tracefin=max(finaltrc(:,1));
   disp([num2str(totalspot) ' trajectories at the beginning.']);
   disp([num2str(cuenta-1) ' trajectories after reconnection.']);
   disp([num2str(tracefin) ' trajectories after filtering the short ones (less than ' num2str(mintrace) ' points).']);
   disp(' ');
   %report
else
   disp('No trajectory left with enough number of points.');
   %report
   res=[];
end

clear newtrc finaltrc trc

% end of file