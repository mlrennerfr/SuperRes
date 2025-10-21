function [output,epsilon, labeled2] =selectclusters(labeled,datax,datay,fr,alpha, mindens,mindenspx,mindetec,minsize,maxsize, minpointsnano,epsilon)
% function [output,epsilon, labeled2] =selectclusters(labeled,datax,datay,fr,alpha, mindens, mindetec,minsize,maxsize, minpointsnano,epsilon)
% 
% selects clusters from a selection made by segmentation (intensity threshold or
% Voronoi tesselation) on a grid
% keeps pixels with more than mindens density and mindetect number of
% detections, with an area between minsize and maxsize
% Marianne Renner oct 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


waitbarhandle=waitbar( 0,'Please wait...','Name','Selecting clusters');
BW=zeros(size(labeled,1),size(labeled,2));

output.rendxmask=[];
rendxmask=[];
rendymask=[];
frmask=[];
alphamask=[];
results=[];
resultsnano=[];
listeboundary=[];
countnanoclusters=0;
countclusters=0;

auxoutx=[datax zeros(size(datax,1),1)];
auxouty=[datay zeros(size(datay,1),1)];
rendxmasknano=[];
rendymasknano=[];

stats = regionprops(labeled,'PixelList','Area') ; 
nb_clust=size(stats,1);
listpoints=[];

if mindenspx>0 % only intensity threshold
    for n=1:nb_clust % chaque cluster
    
       clust=stats(n).PixelList;
       
       for cc=1:size(clust,1)
           index2 = intersect(find(round(datax(:))==clust(cc,1)),find(round(datay(:))==clust(cc,2)));
           if isempty(index2)==0
               if size(index2,1)> mindenspx % minimum density per pixel
                   BW(clust(cc,2), clust(cc,1))=1; %new image with pixels that have enough dens
                   listpoints=[listpoints index2'];
               end
           end
       end
    end
else
    BW=labeled;
end

se = strel('disk',1,0);
BW = imdilate(BW,se);
BW = imfill(BW,'holes');
se = strel('disk',1,0);
BW = imerode(BW,se);
[labeled2,nb_clust] = bwlabel(BW,4);


stats = regionprops(labeled2,'PixelList','Area','MajorAxisLength','MinorAxisLength') ; %voir!!!!!!!!!!!!!!!!!

for n=1:nb_clust % each cluster
    
        clust=stats(n).PixelList;
        areaclu=stats(n).Area;
       
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
           if exist('waitbarhandle')
                waitbar(n/nb_clust,waitbarhandle,['Cluster # ',num2str(n)]);
           end
           for t=1:size(listeclu,1)
                posx=round(datax(listeclu(t))); if posx==0; posx=1; end 
                if posx>size(labeled,2); posx=size(labeled,2); end
                posy=round(datay(listeclu(t))); if posy==0; posy=1; end 
                if posy>size(labeled,1); posy=size(labeled,1); end
                if labeled(posy,posx)>0
                  if posx==1
                      %disp('!!!!')
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

           [polyin,areacluster]=boundary(dataxclu(:,1),datayclu(:,2));

           % density and size
          % if size(dataxclu,1)/areaclu > mindens && areaclu > minsize && areaclu < maxsize
           if size(dataxclu,1)/areacluster > mindens && areaclu > minsize && areaclu < maxsize
              
                %boundary -> inpolygon 
                in = inpolygon(dataxclu(:,1),datayclu(:,1),dataxclu(polyin,1), datayclu(polyin,1));
                countclusters=countclusters+1;

                rendxmask=[rendxmask; dataxclu(in,1) n*ones(size(dataxclu(in,1),1),1) countclusters*ones(size(dataxclu(in,1),1),1)];  
                rendymask=[rendymask; datayclu(in,1) n*ones(size(dataxclu(in,1),1),1) countclusters*ones(size(dataxclu(in,1),1),1)];
                
                frmask=[frmask; frclu];
                alphamask=[alphamask; alphaclu]; 
               
                if minpointsnano>0  %checks for nanodomains inside the cluster
                    % use Delaunay triangulation to check distances between neighbours. 
                    % Thresholding by analysing mean distance within points
                    
                    indataxclu=dataxclu(in,1);
                    indatayclu=datayclu(in,1);
                        
                    [goodpoints, meandist]=delaunaygraph(indataxclu,indatayclu); 

                    if epsilon==0   %calculation auto from distance between points
                        epsilon=meandist*3;  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        epsilon=round(epsilon*100)/100;
                    end
                    
                    %DBSCAN                     
                    idx = dbscan(goodpoints,epsilon,minpointsnano); %idx = dbscan(X,epsilon,minpts)
                    
                    % figure
                    % gscatter(goodpoints(:,1),goodpoints(:,2),idx);
                    
                    for nronano=1:max(idx)
                        newdataxclu=goodpoints(find(idx(:)==nronano),1);
                        newdatayclu=goodpoints(find(idx(:)==nronano),2);
                        
                        if isempty(newdataxclu)==0
                            
                            %boundary -> inpolygon
                            [polyin,areanano]=boundary(newdataxclu(:),newdatayclu(:));
                            if isempty(polyin)==0
                                in = inpolygon(dataxclu(:,1),datayclu(:,1),newdataxclu(polyin), newdatayclu(polyin)); %includes all detection within the boundary
                                rendxmasknano=[rendxmasknano; dataxclu(in,1) n*ones(size(dataxclu(in,1),1),1) countnanoclusters*ones(size(dataxclu(in,1),1),1)];  
                                rendymasknano=[rendymasknano; datayclu(in,1) n*ones(size(dataxclu(in,1),1),1) countnanoclusters*ones(size(dataxclu(in,1),1),1)];
                                countnanoclusters=countnanoclusters+1;

                                resultsnano=[resultsnano; n countnanoclusters numel(dataxclu(in)) size(dataxclu,1) areanano size(in,1)/areanano];
                                        % nro cluster - nro nanocluster- #detect - #detect cluster - area nanocluster - density
                               
                                
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
output.countclustersnano=countnanoclusters;

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