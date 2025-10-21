function I = image_ch(varargin)
xx1 = varargin{1};
yy1 = varargin{2};
I = varargin{3};
pct_sat = 0;
if ndims(I)==2
    if nargin>3
        clims = varargin{4};
    else
        Imax = prctile(I(:),100-pct_sat);
        Imin = min(I(:));
        clims = [Imin Imax];
    end
    if any(isnan(clims)) || clims(2)<=clims(1)
        warning(['clims (=',num2str(clims),')! Using [0,1] instead !']);
        clims = [0 1];
    end
    I = imagesc(xx1,yy1,I,clims);
else
    image(xx1,yy1,I);
end
colorbar;

