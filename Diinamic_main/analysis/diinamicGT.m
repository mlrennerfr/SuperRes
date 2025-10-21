function [clusterdata, dataclustersD, dataclustersGT,elapsedTime]=diinamicGT(handles,Irend, namefile,closecode)
% launches analysis

% results
clusterdata=[];
dataclustersD=[];
dataclustersGT=[];

pxdilate=1;
pxerode=1;
mindens=handles.mindens;
mindenspx=handles.mindenspx;
mindetect=handles.mindetect;
minsize=handles.minsize;
inthresh=handles.inthresh;
maxsize=100000;

handles.dim=[size(Irend,1) size(Irend,2)];

% points ROI
roiselectedx=handles.x; 
roiselectedy=handles.y; 
handles.datax=roiselectedx /handles.dx;
handles.datay=roiselectedy /handles.dx;

facteur=handles.dx*1000;

gtdatax=handles.GTdata(:,1)/facteur;
gtdatay=handles.GTdata(:,2)/facteur;

valintthresh=inthresh/100*max(max(Irend));
BW=zeros(size(Irend,1),size(Irend,2));
for i=1:size(Irend,1)
    indexj=find(Irend(i,:)>valintthresh);
    if isempty(indexj)==0
        BW(i,indexj)=1;
    end
end
if pxdilate>0
    se = strel('disk',pxdilate,0);
    I2 = imdilate(BW,se);
    I2 = imfill(I2,'holes');
else
    I2=BW;
end
BW=I2;
if pxerode>0
    se = strel('disk',pxerode,0);
    I2 = imerode(BW,se);
end
       
[labeledall,~] = bwlabel(I2,4);

maxsize=(((maxsize/2)^2)*pi)/(handles.dx*1000)^2 ;
minsize=(((minsize/2)^2)*pi)/(handles.dx*1000)^2;

minpointsnano=0;
epsilon=1;

tic;

% keep clusters with more than mindens density of points
[output,epsilon,labeled2] =selectclusters(labeledall,handles.datax ,handles.datay,handles.fr,handles.alpha,mindens,mindenspx,mindetect,minsize,maxsize,minpointsnano,epsilon);

elapsedTime = toc;

countclusters=output.countclusters;
disp([num2str(countclusters),' clusters found'])
disp(' ')

if closecode==0
    figurerend= figure('Name','Segmented rendered image','Toolbar','figure');
    imshow(labeled2,'InitialMagnification','fit')
end

rendxmask=output.rendxmask;
rendymask=output.rendymask;
countclusters=output.countclusters;


if isempty(rendxmask)==0
    
  % visualization
  imageclu=zeros(size(I2,1),size(I2,2));
  fond=ones(size(I2,1),size(I2,2));

  for z=1:max(rendxmask(:,3)) %all selected clusters
    index=find(rendxmask(:,3)==z);
    if isempty(index)==0
        for t=1:size(index)
            posx=floor(rendxmask(index(t),1));
            if posx==0; posx=1; end
            posy=floor(rendymask(index(t),1));
            if posy==0; posy=1; end
            imageclu(posy,posx)=1; %new image with pixels that have enough dens
        end
    end
  end
  
  se = strel('ball',5,5);
  dilatedI = imdilate(imageclu,se);
  imageclu=imerode(dilatedI,se);
  for z=1:size(imageclu,1)
    index=find(imageclu(z,:)>0);
    if isempty(index)==0
        imageclu(z,index)=1;
    end
  end
  imageclu = imfill(imageclu,'holes');
  
  % figure clusters
  figuremask=figure;
  hold on
  imshow(fond,'InitialMagnification','fit')
  axis off
  
  plot(handles.datax,handles.datay,'o','MarkerSize',5,'Color',[0.5 0.55 0.55]);
  
  if isempty(rendxmask)==0
      
        %possibility to plot with different colors!!
        colorcode=(jet(max(rendxmask(:,3))+1));
        for ll=1:size(colorcode)
            order=ceil(rand*size(colorcode,1));
            newcolorcode(ll,:)=colorcode(order,:);
        end
        count=1;
        for jj=1:max(rendxmask(:,3))
            index=find(rendxmask(:,3)==jj);
            if isempty(index)==0
               plot(rendxmask(index,1),rendymask(index,1),'o','MarkerSize',3,'Color',newcolorcode(jj,:)); % all clustered points or clusterd points out loc
               hold on
               [polyin, areaclu]=boundary(rendxmask(index,1),rendymask(index,1));
               plot(rendxmask(index(polyin),1),rendymask(index(polyin),1),'Color',newcolorcode(jj,:)); % all clustered points or clusterd points out loc
               
               % save data
               clusterdata=[clusterdata; rendxmask(index,1) rendymask(index,1) count*ones(size(index,1),1)]; % x - y - clu number
               dataclustersD=[dataclustersD; count size(index,1) areaclu]; %#clu #det area
               count=count+1;

            end
        end
  end
  count=1;
  for pp=1:max(handles.GTdata(:,3))
      index=find(handles.GTdata(:,3)==pp);
      if isempty(index)==0
          plot(gtdatax(index),gtdatay(index),'.k'); 
          hold on
          [polyin, areaclu]=boundary(gtdatax(index),gtdatay(index));
          plot(gtdatax(index(polyin)),gtdatay(index(polyin)),'k'); 
          dataclustersGT=[dataclustersGT; count size(index,1) areaclu]; %#clu #det area
          count=count+1;

      end
  end
  
  title(namefile,'Interpreter','none')
  hold off
  
  savename=[namefile,'-plot.fig'];
  savefig(savename)
  
  if closecode==1
      close(figuremask)
  end
  
else
    disp('No clusters left')
end

    
disp('Done')



