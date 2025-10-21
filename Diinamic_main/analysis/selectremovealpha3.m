function vect2remove=selectremovealpha2(handles,x,y,alpha)

%% Remove points with alphas outside specified intensities range
vect2remove = [];
 
%% Remove points at x=y=0
aux1 = find(x==0);
aux2 = find(y==0);
vect2remove = intersect(aux1,aux2);
if exist('Xmatrix','var')
    matrix2remove = flagmatrixelements(matrix2remove,vect2remove);
end
aux = find(alpha==0);
vect2remove = [vect2remove; aux];
if exist('Xmatrix','var')
    matrix2remove = flagmatrixelements(matrix2remove,vect2remove);
end

percalphamin = handles.alpha1;
percalphamax = handles.alpha2;

valminalpha=min(handles.alpha(:));
valmaxalpha=max(handles.alpha(:));
difalpha=valmaxalpha-valminalpha;

alphamin=valminalpha+(percalphamin/100*difalpha);
alphamax=valminalpha+(percalphamax/100*difalpha);

aux = union(find(alpha<alphamin),find(alpha>alphamax));

if size(aux,2)>1
    vect2remove = [vect2remove; aux(1)];
else
    vect2remove = [vect2remove; aux];
end
   