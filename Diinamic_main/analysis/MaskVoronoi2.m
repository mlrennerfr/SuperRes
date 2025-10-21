function [surfacearea, mask]=MaskVoronoi2(fondo,roiselectedx,roiselectedy,limitarea)
%function [surfacearea, mask]=MaskVoronoi2(fondo,roiselectedx,roiselectedy,limitarea)
%
% Calculates a mask from Voronoi tessels (criteria: minimum surface of
% polygons)
% Calculates area
%
% Marianne Renner version for Clustering.m 06/16
% Marianne Renner version for ClusterDens.m 05/20
% Marianne Renner version for Diinamic.m 01/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mensaje=msgbox('Creating Voronoi tesselation','Please wait');
positions=[roiselectedx,roiselectedy];
[v,c] = voronoin(positions);

figuremask=figure;
fondob=zeros(size(fondo,1),size(fondo,2));
imshow(fondob,'InitialMagnification','fit')
hold on
%plot(roiselectedx,roiselectedy,'.','MarkerSize',5,'Color','r');
%hold on
surfacearea=[];   

for gg=1:length(c)
    patch(v(c{gg},1),v(c{gg},2),gg);
    surfacearea=[surfacearea; gg polyarea(v(c{gg},1),v(c{gg},2))];
    hold on
end
plot(roiselectedx,roiselectedy,'.','MarkerSize',5,'Color','r');
hold off     

indexlimit=find(surfacearea(:,2)<limitarea);
if isempty(indexlimit)==0
    selecregions=surfacearea(indexlimit,:);
end

aux=fondob;
for kk=1:length(selecregions)
    gg=selecregions(kk,1);
    roi=roipoly(fondob,v(c{gg},1),v(c{gg},2));
    aux=aux+roi;
end
mask=aux;

%figure
%imshow(aux,'InitialMagnification','fit')
%hold on
%plot(roiselectedx,roiselectedy,'.','MarkerSize',5,'Color','r');

if exist('mensaje') %#ok<EXIST>
    close(mensaje)
end

%save(['areasROI', num2str(j),'.txt'],'surfacearea','-ascii');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

