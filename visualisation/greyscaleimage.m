function C=greyscaleimage(Im,low,high)
% function C=greyscaleimage(Im,low,high)
% transforms a double type image 'Im' in an other image with grey values 
% between 0 and 1. 'low' and 'high' arguments allow to modify the grey 
% scale of the image (values in percent)
%_____________________________________________________________________


l1=(max(Im(:))-min(Im(:)))*low+min(Im(:));
l2=(max(Im(:))-min(Im(:)))*high+min(Im(:));
if l2-l1==0
    C=Im;
else
    a=1/(l2-l1);
    b=-l1*a;
    f=@(x) x*a+b;
    C=f(Im);
end
C(C>1)=1;
C(C<0)=0;
clear l1 l2 a b f

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
