function vect2remove=selectremovealpha(handles,x,y,alpha)
%% Remove points with alphas outside specified percentile range

vect2remove = [];

%% Remove points at x=y=0
aux1 = find(x==0);
aux2 = find(y==0);
vect2remove = intersect(aux1,aux2);
if exist('Xmatrix','var')
    matrix2remove = flagmatrixelements(matrix2remove,vect2remove);
end

%% Remove points with alpha=0
aux = find(alpha==0);
vect2remove = [vect2remove; aux];
if exist('Xmatrix','var')
    matrix2remove = flagmatrixelements(matrix2remove,vect2remove);
end

pctmin = str2num(get(handles.alpha1,'String'));
pctmax = str2num(get(handles.alpha2,'String'));
if exist('Xmatrix','var')
    aux = alphamatrix(:);
else
    aux = alpha;
end
aux = aux(aux>0);
alphamin = prctile(aux,pctmin);
alphamax = prctile(aux,pctmax);

aux = union(find(alpha<alphamin),find(alpha>alphamax));
vect2remove = [vect2remove; aux];