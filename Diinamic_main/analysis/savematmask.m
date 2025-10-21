function savematmask(handles,indexin,mu, rendxin,rendyin,Savename)

indexminx=find(handles.x==min(handles.x));
indexminy=find(handles.y==min(handles.y));

indexmaxx=find(handles.x==max(handles.x));
indexmaxy=find(handles.y==max(handles.y));


matrice_results(1,:)=handles.fr';
matrice_results(2,:)=handles.y';
matrice_results(3,:)=handles.x';
matrice_results(4,:)=handles.alpha';
matrice_results(5,:)=handles.radius';
if size(handles.sigma)>0
    matrice_results(6,:)=handles.sigma';
end
if isfield(handles, 'blink')
    matrice_results(7,:)=handles.blink';
else
     matrice_results(7,:)=zeros(1,size(matrice_results,2));
end

if isfield(handles,'ratio') %3D
    matrice_results(8,:)=handles.ratio';
    matrice_results(9,:)=handles.z';
    matrice_results(10,:)=handles.test1';
    matrice_results(11,:)=handles.test2';
else
    matrice_results(8,:)=zeros(1,size(matrice_results,2));
    matrice_results(9,:)=zeros(1,size(matrice_results,2));
    matrice_results(10,:)=zeros(1,size(matrice_results,2));
    matrice_results(11,:)=zeros(1,size(matrice_results,2));
    
end

aux=matrice_results(:,indexin);
if mu==0
    aux(3,:)=rendxin';
    aux(2,:)=rendyin';
end

%corners
alreadyin=find(indexin(:)==indexminx);
if isempty(alreadyin)
    %disp(size(matrice_results(:,indexminx)))
    aux=[aux matrice_results(:,indexminx)]; % minx
end
alreadyin=find(indexin(:)==indexminy);
if isempty(alreadyin)
    aux=[aux matrice_results(:,indexminy)]; % minx
end
alreadyin=find(indexin(:)==indexmaxx);
if isempty(alreadyin)
    aux=[aux matrice_results(:,indexmaxx)]; % minx
end
alreadyin=find(indexin(:)==indexmaxy);
if isempty(alreadyin)
    aux=[aux matrice_results(:,indexmaxy)]; % minx
end

matrice_results=aux;

save(Savename, 'matrice_results');       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            
            
