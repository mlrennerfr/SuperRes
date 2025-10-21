function [goodpoints, meandist]=delaunaygraph (datax,datay)
%function [goodpoints, meandist]=delaunaygraph (datax,datay)
%
% used to calculate epsilon for DBSCAN in sub-domain detection
%
% Marianne Renner dec22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pts=[datax datay];
DT = delaunayTriangulation(pts); 
e=DT.edges;
allpoints=DT.Points;

% Distances of the edges for weights
plist1 = pts(e(:,1),:); %start points
plist2 = pts(e(:,2),:);  %end points
diffs = plist2-plist1; %distance between all the points forming edges
if diffs<0
    diffs=diffs*-1;
end
dists = vecnorm(diffs');
meandist=mean(dists);

limitdist=meandist/2;

%eliminate conexions if distance to high to be in the same nanodomain
indexok=find(dists<limitdist);
newe=e(indexok,:);
allpossible=[newe(:,1); newe(:,2)];

newpoints=allpossible(1); 
for i=1:size(allpossible,1) %all the possible points
    index=find(newpoints(:)==allpossible(i));
    if isempty(index)==1
        newpoints=[newpoints;allpossible(i)];
    end
end

newpoints=sortrows(newpoints,1);

goodpoints=allpoints(newpoints,:);


% create Graph
%G = graph(e(:,1), e(:,2), dists);

%figure
%plot(allpoints(:,1),allpoints(:,2),'.b');
%figure
%plot(goodpoints(:,1),goodpoints(:,2),'.g');

%p = plot(G,'EdgeLabel',G.Edges.Weight);

% only short distances
%newG=graph(e(indexok,1), e(indexok,2), dists(indexok));

%figure
%p = plot(newG,'EdgeLabel',newG.Edges.Weight);

%analyze centrality
%figure
%p = plot(G);
%wcc = centrality(G,'closeness','Cost',G.Edges.Weight);
%p.NodeCData = wcc;
%title('Closeness Centrality Scores - Weighted')

%figure
%p = plot(newG);

%wcc = centrality(newG,'closeness','Cost',newG.Edges.Weight);
%p.NodeCData = wcc;
%title('Closeness Centrality Scores - Weighted')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
