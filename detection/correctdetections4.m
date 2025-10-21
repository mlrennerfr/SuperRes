function corrpeak=correctdetections4(matfile,dist,maxblink)
% function corrpeak=correctdetections4(matfile,dist,maxblink)
% runs the correction for multiple detection
% all detections within dist distance during maxblink period
% are considered as coming from a single molecule.
% the position of this molecule is kept, as the mean position of all the
% detections
%
% MR 08/2019 SuperRes_v2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

validfile = 0;
corrpeak=[];
%dist=dist*2; % side of a square

%read data
S=load(matfile) ;
if isfield(S,'matrice_results')
            aux = S.matrice_results;
            clear S;
            peak=[aux(1,:)' aux(3,:)' aux(2,:)'  aux(4,:)' aux(5,:)'*2 ];    
            validfile = 1;
else
            %validfile = 0;
            disp('Invalid file!');
            return
end

if validfile==0
    return
end

% group peaks same molecule
nSeq=max(peak(:,1));
peak=[peak zeros(size(peak,1),1)]; %last column: number of individual mol
waitbarhandle=waitbar( 0,'Please wait...','Name',['Selecting detections of the same molecule in ',matfile]) ;
    
count=0;
corrpeak=[];    
colmol=size(peak,2);


for i=1:size(peak,1) %loop all peaks
    if exist('waitbarhandle')
        waitbar(i/size(peak,1),waitbarhandle,['Peak # ',num2str(i)]);
    end

    if peak(i,colmol)==0 % not identified yet
        
        truecandidates=[];
        count=count+1;
        peak(i,colmol)= count; %VOIR
        
        x1=peak(i,2)-dist;         x2=peak(i,2)+dist;
        y1=peak(i,3)-dist;         y2=peak(i,3)+dist;

        candidates=inpolygon(peak(:,2),peak(:,3),[x1 x2],[y1 y2]); % in a square around the present peak !!!!!!!!!!!!!!!!
        indexpeak=find(candidates(:)>0);
        
      %  disp(count)
      %  disp('indexpeak')
      %  disp(size(indexpeak))
        
        if isempty(indexpeak)==0   
            
            truecandidates(1)=indexpeak(1); % first one included by default

            for j=2:size(indexpeak,1)
                
                % consecutives or between blink limit
                if peak(indexpeak(j),1)-peak(indexpeak(j-1),1)<maxblink
                    truecandidates=[truecandidates; indexpeak(j)] ;
                end
                
            end % loop candidates
            
           % disp('truecandidates')
           % disp(size(truecandidates))

            peak(truecandidates(:),colmol)= count;
        else
            disp('empty square...')
            
        end 
    end % not identified yet
        
end   % loop peaks

if exist('waitbarhandle')
    close(waitbarhandle);
end

% mean of multiple detections
waitbarhandle=waitbar( 0,'Please wait...','Name',['Correcting detections in ',matfile]) ;

% VOIR if zeros!!!!!!!!!!!!
%indexzero=find(peak(:,colmol)==0);
%disp('peaks')
%disp(size(peak))
%disp('zeros')
%disp(size(indexzero))
%corrpeak=peak(indexzero,:);
%disp(max(peak(:,colmol)))

for i=1:max(peak(:,colmol)) %loop all peaks
    
    if exist('waitbarhandle')
        waitbar(i/max(peak(:,colmol)),waitbarhandle,['Detection # ',num2str(i)]);
    end
    
   % disp(i)
    index=find(peak(:,colmol)==i);
    
    if isempty(index)==0
        %mean position for all the peaks with the same number
      %  corrpeak=[corrpeak; peak(index(1),1) mean(peak(index(:),2)) mean(peak(index(:),3)) peak(index(1),4:5)]; 
     %   peak=[aux(1,:)' aux(3,:)' aux(2,:)'  aux(4,:)' aux(5,:)'*2 ];  
     
    % disp(size(index))

        corrpeak=[corrpeak; aux(1,index(1))' mean(aux(3,index(:))) mean(aux(2,index(:))) aux(4:size(aux,1),index(1))']; 
    end
    
end

disp([num2str(size(peak,1)),'  detections in the original file, ',num2str(size(corrpeak,1)),'  detections after correction.'])

if exist('waitbarhandle')
    close(waitbarhandle);
end

save('peak.txt','peak','-ascii')
save('corrpeak.txt','corrpeak','-ascii')

%disp('Done')
disp(' ')

%plotcontrol=0;
%if plotcontrol
          figure
       Im=ones(ceil(max(peak(:,3))),ceil(max(peak(:,2)))); %!!!!!!!!!!!!!!!!!!
       imshow(Im,'InitialMagnification','fit');
       hold on
       plot(peak(:,2),peak(:,3),'o','MarkerSize',5,'Color','b');
       hold on
       plot(corrpeak(:,2),corrpeak(:,3),'x','MarkerSize',5,'Color','r');
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    
    
    