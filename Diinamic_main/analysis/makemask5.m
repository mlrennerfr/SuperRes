function [rendxmask, rendymask,frmask,alphamask, results,listeboundary,countclusters,auxoutx,auxouty,rendxmasknano, rendymasknano, results nano, countclustersnano] =makemask5(labeled,datax,datay,fr,alpha, mindens, mindetec,minsize,maxsize, minpointsnano)
%function [rendxmask, rendymask,frmask,alphamask, results,listeboundary,countclusters,auxoutx,auxouty] =makemask5(labeled,datax,datay,fr,alpha, mindens, mindetec,minsize,maxsize,minpointsnano)
%
% Main script for finding candidate clsters and for selecting them
%
% Marianne Renner oct2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<11
    labeledloc=[];
end
waitbarhandle=waitbar( 0,'Please wait...','Name','Selecting clusters');
BW=zeros(size(labeled,1),size(labeled,2));

rendxmask=[];
rendymask=[];
frmask=[];
alphamask=[];
results=[];
resultsnano=[];
listeboundary=[];
countclusters=1;
countnanoclusters=1;
auxoutx=[datax zeros(size(datax,1),1)];
auxouty=[datay zeros(size(datay,1),1)];
rendxmasknano=[];
rendymasknano=[];

stats = regionprops(labeled,'PixelList','Area') ; 
nb_clust=size(stats,1);
listpoints=[];

for n=1:nb_clust % chaque cluster
    
        clust=stats(n).PixelList;
        
       for cc=1:size(clust,1)
           index=find(round(datax(:))==clust(cc,1));
           if isempty(index)==0
               index2=find(round(datay(index))==clust(cc,2));
               if isempty(index2)==0
                   if size(index2,1)> mindens % minimum density per pixel
                       BW(clust(cc,2), clust(cc,1))=1; %new image with pixels that have enough dens
                       listpoints=[listpoints index(index2)'];
                   end
               end
           end
       end
       
       % area as convex hull????
       % boundary....
       
         % for DBSCAN:  AreaC = polyarea(Contour(:,1),Contour(:,2)); but
         % needs vertices (Contour)
         % [K,area] = convhull(datax(listpoints),datay(listpoints))

       
end

se = strel('disk',1,0);
BW = imdilate(BW,se);
BW = imfill(BW,'holes');
se = strel('disk',1,0);
BW = imerode(BW,se);
[labeled2,nb_clust] = bwlabel(BW,4);

stats = regionprops(labeled2,'PixelList','Area','MajorAxisLength','MinorAxisLength') ; %voir!!!!!!!!!!!!!!!!!

countclusters=1;
for n=1:nb_clust % each cluster
    
        clust=stats(n).PixelList;
        areaclu=stats(n).Area;
        minaxisclu=stats(n).MinorAxisLength;
        maxaxisclu=stats(n).MajorAxisLength;
        
        dataxclu=[];
        datayclu=[];
        frclu=[];
        alphaclu=[];
        listeclu=[];
        
       for cc=1:size(clust,1)
           index=find(round(datax(:))==clust(cc,1));
           if isempty(index)==0
               index2=find(round(datay(index))==clust(cc,2));
               if isempty(index2)==0
                       listeclu=[listeclu; index(index2)]; %index of detections in each pixel of the cluster
               end
           end
       end

       if isempty(listeclu)==0 &&  size(listeclu,1)> mindetec  
           
             disp(['Cluster # ',num2str(n),', ',num2str(size(listeclu,1)),' detections']);

            indextotal=[];
            aux=[];
             controlloc=0;
             
            % count detections for each pixel and select those with dens>mindens
           
            for t=1:size(listeclu,1)
                 
               if exist('waitbarhandle')
                   waitbar(t/size(listeclu,1),waitbarhandle,['Cluster # ',num2str(n),', ',num2str(floor(t/size(listeclu,1)*100)),'%']);
               end

                posx=round(datax(listeclu(t))); if posx==0; posx=1; end 
                %if posx>max(datax(listeclu(t))); posx=posx-1; end
                if posx>size(labeled,2); posx=size(labeled,2); end
                
                posy=round(datay(listeclu(t))); if posy==0; posy=1; end 
                %if posy>max(datay(listeclu(t))); posy=posy-1; end
                if posy>size(labeled,1); posy=size(labeled,1); end
               
                if labeled(posy,posx)>0
                  if posx==1
                      disp('!!!!')
                  end
                    aux=[aux; listeclu(t) datax(listeclu(t)) datay(listeclu(t))];
                end

           end
            
           indextotal=[indextotal; aux n*ones(size(aux,1),1)];
           
            % total in the cluster
            if isempty(indextotal)==0
                    dataxclu=[datax(indextotal(:,1)) indextotal(:,2) indextotal(:,4)];
                    datayclu=[datay(indextotal(:,1)) indextotal(:,3) indextotal(:,4)];
                    frclu=[fr(indextotal(:,1)) indextotal(:,4)];
                    alphaclu=[alpha(indextotal(:,1)) indextotal(:,4)];
            end
            
            %normalization : maximum number of detections  !!!!!!!!!  
          %  if size(dataxclu,1)>maxnumberdet
           %     dataxclu=dataxclu(1:maxnumberdet,:);
            %    datayclu=datayclu(1:maxnumberdet,:);
            %    frclu=frclu(1:maxnumberdet,:);
            %    alphaclu=alphaclu(1:maxnumberdet,:)';
            %    disp('Max number of detections reached.')
            %end
            
            % density and size
            if size(dataxclu,1)/areaclu > mindens && areaclu > minsize && areaclu < maxsize
              
                % output
                %rendxmask=[rendxmask; dataxclu nb_clust*ones(size(dataxclu,1),1)];  % att nro cluster if nanoclusters!!!
                %rendymask=[rendymask; datayclu nb_clust*ones(size(datayclu,1),1)];
                
                %boundary -> inpolygon 
                [polyin,areacluster]=boundary(dataxclu(:,1),datayclu(:,2));
                in = inpolygon(dataxclu(:,1),datayclu(:,1),dataxclu(polyin,1), datayclu(polyin,1));
                
                rendxmask=[rendxmask; dataxclu(in,:) countclusters*ones(size(dataxclu(in,:),1),1)];  
                rendymask=[rendymask; datayclu(in,:) countclusters*ones(size(dataxclu(in,:),1),1)];
                
                %%%% att!!!!!!!
                frmask=[frmask; frclu];
                alphamask=[alphamask; alphaclu]; 
                countclusters=countclusters+1;

               
                if minpointsnano>0  %checks for nanodomains inside the cluster
                    % use Delaunay triangulation to check distances between neighbours. 
                    % Thresholding by analysing mean distance within points
                    
                    % pts=[dataxclu(:,1),datayclu(:,1)];
                    [goodpoints, meandist]=delaunaygraph(dataxclu(:,1),datayclu(:,1)); 
                    
                    %DBSCAN using meandist as epsilon                    
                    idx = dbscan(goodpoints,meandist*3,minpointsnano); %idx = dbscan(X,epsilon,minpts)
                    
                    idx=idx+2; %to avoid negative values
                    % figure
                    % gscatter(goodpoints(:,1),goodpoints(:,2),idx);
                    
                    for nronano=min(idx):max(idx)
                        newdataxclu=goodpoints(find(idx(:)==nronano),1);
                        newdatayclu=goodpoints(find(idx(:)==nronano),2);
                        
                        if isempty(newdataxclu)==0
                            
                            %boundary -> inpolygon
                            [polyin,areanano]=boundary(newdataxclu(:),newdatayclu(:));
                            if isempty(polyin)==0
                                in = inpolygon(dataxclu(:,1),datayclu(:,1),newdataxclu(polyin), newdatayclu(polyin)); %includes all detection within the boundary
                                    
                                rendxmasknano=[rendxmasknano; dataxclu(in,:) n*ones(size(dataxclu(in,:),1),1) countnanoclusters*ones(size(dataxclu(in,:),1),1)];  
                                rendymasknano=[rendymasknano; datayclu(in,:) n*ones(size(dataxclu(in,:),1),1) countnanoclusters*ones(size(dataxclu(in,:),1),1)];

                                resultsnano=[resultsnano; n countnanoclusters size(in,1) size(dataxclu,1) areanano size(in,1)/areanano];
                                        % nro cluster - nro nanocluster- # detect - #detect cluster - area nanocluster - densité
                               
                                countnanoclusters=countnanoclusters+1;
                            end %polyin
                        end %data inside
                    end %all nanodomains
                end %check for nanodomains
                
                if isempty(indextotal)==0
                    auxoutx(indextotal(:,1))=1; % data about localization (in and out)
                    auxouty(indextotal(:,1))=1;
                end

                typeloc=0;
                
                if max(dataxclu(:,size(dataxclu,2)))>0 %loc in
                     typeloc=1;
                end
                
                results=[results; n size(dataxclu,1) size(datax,1) areacluster size(dataxclu,1)/areacluster];
                    % nro cluster - # detect - #detect roi - area cluster - densité
                
            end %check size
                            
        end % empty list and mindetec
end            

close(waitbarhandle);


output.rendxmask=rendxmask;
output.rendymask=rendymask;
output.frmask=frmask;
output.alphamask=alphamask;
output.results=results;
output.listeboundary=listeboundary;
output.countclusters=countclusters;
output.auxoutx=auxoutx;
output.auxouty=auxouty;
output.rendxmasknano=rendxmasknano;
output.rendymasknano=rendymasknano;
output.resultsnano=resultsnano;
output.countclustersnano=countclustersnano;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function linesprofile=definelines(statsbound,nb_clust)

angle= -pi:0.01:pi;

for t=1:nb_clust
    maxdiam=statsbound(t).MajorAxisLength;
    dist=maxdiam/2 + maxdiam/10;
    posx=statsbound(t).Centroid(1);
    posy=statsbound(t).Centroid(2);
    for a=1:size(angle)
        x1(a)=round(posx-dist*cos(angle(a)));
        y1(a)=round(posy-dist*sin(angle(a)));
        x2(a)=round(posx+dist*cos(angle(a)));
        y2(a)=round(posy+dist*sin(angle(a)));
    end
    linesprofile{t}=[x1 y1 x2 y2] ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%