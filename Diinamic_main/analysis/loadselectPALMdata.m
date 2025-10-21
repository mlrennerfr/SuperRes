function handles=loadselectPALMdata(filename,szpx,handles)
%function handles=loadselectPALMdata(filename,szpx,handles)
%
% reads .mat file containing detection results
% called by different script in Diinamic
% szpx: camera pixel size
%
% Marianne Renner 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

points=[];
S = load(filename);
sigma=[];

if isfield(S,'Xmatrix')
           % disp(' File with tracking information');
            Xmatrix = S.Xmatrix * szpx;
            Ymatrix = S.Ymatrix * szpx;
            %Xmatrix = S.Xmatrix;
           % Ymatrix = S.Ymatrix;
            alphamatrix = S.alphamatrix;
            frmatrix = S.frmatrix;
            if isfield(S,'radiusmatrix')
                radiusmatrix = S.radiusmatrix;
            else
                radiusmatrix=zeros(size(Xmatrix,1),1);
            end
            x = Xmatrix(:);
            y = Ymatrix(:);
            alpha = alphamatrix(:);
            fr = frmatrix(:);
            radius = radiusmatrix(:); %????
            clear S;
            
elseif isfield(S,'matrice_results')
            aux = S.matrice_results;
            clear S;
            fr = aux(1,:); fr = fr(:);
            x = aux(3,:) * szpx;
           % x = aux(3,:);
            x = x(:);
            y = aux(2,:) * szpx; 
           % y = aux(2,:); 
            y = y(:);
            alpha = aux(4,:); alpha = alpha(:);
            radius=(aux(5,:));
            if size(aux,1)>5
                sigma=(aux(7,:));
            end
            
elseif isfield(S,'alpha') && isfield(S,'fr') && isfield(S,'x') && isfield(S,'y')
            x = S.x * szpx;
            y = S.y * szpx;
          %  x = S.x ;
            x = x(:);
           % y = S.y ;
            y = y(:);
            alpha = S.alpha; alpha = alpha(:);
            fr = S.fr; fr = fr(:);
            if isfield(S,'radius')
                radius=S.radius;
            else
                radius=zeros(size(x,1),1);
            end
            radius=radius(:);
else
            %validfile = 0;
            disp('Invalid file!');
            return
end


clear S


%% Remove points at x=y=0
aux1 = find(x==0);
aux2 = find(y==0);
vect2remove = intersect(aux1,aux2);
if exist('Xmatrix','var')
    matrix2remove = flagmatrixelements(matrix2remove,vect2remove);
end
 
if ~isempty(vect2remove)
    x(vect2remove) = [];
    y(vect2remove) = [];
    alpha(vect2remove) = [];
    fr(vect2remove) = [];
    radius(vect2remove) = [];
end
 
% Storing
handles.x = x;
handles.y = y;
handles.alpha = alpha;
handles.fr = fr;
handles.radius=radius;
handles.sigma=sigma;


clear x y alpha fr radius
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%