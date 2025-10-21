
%% Render PALM image
function [I,xx1,yy1] = PALM_rendering3(varargin)


x = varargin{1};
y = varargin{2};
alpha = varargin{3};
sigmas = varargin{4};
dx = varargin{5};
removeduplicates = varargin{6};
xdim = varargin{7};
ydim = varargin{8};

x = x(:);
y = y(:);
dy = dx;

% axes
xmin=1;
ymin=1;
xmax=max(x);
ymax=max(y);

xx1 = 1*dx:dx:xdim*dx;
yy1 = 1*dx:dy:ydim*dx;
Nx = length(xx1)-1;
Ny = length(yy1)-1;

Ipd1 = zeros(Ny,Nx); % initializes the image

iii = find(~isnan(x));
iii = iii(:)';
x = x(iii);
y = y(iii);
alpha = alpha(iii);
N = length(iii);


%% Create a histogram of positions (optionally weighted by the intensity)
pos = [x,y];
edges{1} = xx1;
edges{2} = yy1; 

aux = hist3(pos,'Edges',edges)';
Ih = aux(1:end-1,1:end-1);

%% convolve the histogram with Gaussian kernel using its separability
%if smooth
    q = 3; % half-size of the kernel in terms of number of standard deviations
    Nk = ceil(q*sigmas/dx)*2+1;
    aux1 = linspace(-q*sigmas, q*sigmas,Nk);     % x^2
    Ik = exp(-aux1.^2/(2*sigmas^2));  % exp( - (x^2+y^2)/(2*sigma^2) )
    if 1==1 % normalize to unit energy by area
        Ik = Ik/sum(Ik(:));
    end
    hcol = Ik(:);
    hrow = hcol';
    Ipd1 = conv2(hcol,hrow,Ih,'same');
%else
%    Ipd1 = Ih;
%end

% normalize maximum intensity to 1
if 1==0
    aux  = max(Ipd1(:));
    Ipd1 = Ipd1/aux;
end


I = Ipd1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

