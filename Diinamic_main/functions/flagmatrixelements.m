function M = flagmatrixelements(M,vectorindices)
siz = size(M);
aux = M(:);
aux(vectorindices) = 1;
M = reshape(aux,siz);
 