function matrice_results =correctshiftstage3(x,y, fr, alpha, radius,handles)
%function matrice_results =correctshiftstage3(x,y, fr, alpha, radius,handles)
%
% from PALMvis, correction of  stage drift
%
% Marianne Renner, for SuperRes programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xmax = max(handles.x);
ymax = max(handles.y);
pointsize = 3;

Im=ones(ceil(xmax),ceil(ymax));
figure
xlabel('X (mu)');
ylabel('Y (mu)');   

imshow(Im,'initialmagnification','fit')
hold on
plot(handles.x,handles.y,'.','MarkerSize',pointsize,'Color','k');
hold on

title([' Total nb of positions =',num2str(length(handles.x))]);


frames = 1:max(handles.fr);
% Let user draw polygons around beads
ip = 1;
stop = 0;

while ~stop
    
    title(['Please draw a polygon around bead #',num2str(ip)]);
    % roipoly!!        %h{ip} = imfreehand(gca);
    
    [~,x,y]=roipolyold;  % ver version vieja!!!!
    %[areaselect,xi,yi]=roipolyold;    %seleccion ROI

    coord{ip}=[x y];
    plot(x,y,'-r');
    hold on      
    button = questdlg('Select another bead ?','','Yes','No','Yes');
    stop = strcmp(button,'No');
    ip = ip+1;
end

handles.totalbeads=[];
%Nbeads = length(h);
ip=ip-1;
Nbeads = ip;
disp([num2str(Nbeads),' beads selected']);
% create matrices from bead position data
xbeads = NaN(ip,length(frames));
ybeads = NaN(ip,length(frames));
frbeads = NaN(ip,length(frames));
alphabeads = NaN(ip,length(frames));   

for ip= 1:Nbeads
    hpos=coord{ip} ;   % roipoly!!   %hpos = getPosition(h{ip})   
    % find positions inside polygon
    aux = find(inpolygon(handles.x,handles.y,hpos(:,1),hpos(:,2))) ;
    handles.totalbeads=[handles.totalbeads;aux];
    if isempty(aux)
        warndlg('No positions inside this polygon !');
    end
    xbead = handles.x(aux);
    ybead = handles.y(aux);
    frbead = handles.fr(aux);
    alphabead = handles.alpha(aux);
    % statistics on positions inside polygon
    xbeads(ip,frbead) = xbead;
    ybeads(ip,frbead) = ybead;
    alphabeads(ip,frbead) = alphabead;
   % xbeadsmed(ip) = nanmedian(xbead);
   % ybeadsmed(ip) = nanmedian(ybead);
    %text(max(xbead),max(ybead),['Landmark #',num2str(ip),': sigma(x) =',num2str(std(xbead)*1000),' nm, sigma(y) =',num2str(std(ybead)*1000),'nm']);
    plot(xbeads',ybeads','.r');
    hold on
end

hold off
handles.xbeads = xbeads;
handles.ybeads = ybeads;
handles.alphabeads = alphabeads;

    %xbeads = handles.xbeads;
    %ybeads = handles.ybeads;
    
    Nbeads = size(xbeads,1);
    
    message=msgbox('Caclulating drift, please wait');

    [xshift,yshift,xbeadsdev0,ybeadsdev0,xbeadsdev1,ybeadsdev1] = estimate_drift(frames,xbeads,ybeads,handles);
    
    close(message)
    
   % if isfield(handles,'fig_bead')
   %     figure(handles.fig_bead);
   % else
        handles.fig_bead = figure('Name','Bead trajectory information');
   % end
    
    % Plot the histograms of bead displacements along x and y
    subplot(3,2,1);
    [n,xout] = hist(xbeadsdev1',100);
    stairs(xout,n); hold on;
    xlabel('x-median(x) (mu)');
    ylabel('#');
    %     title(['std. dev. (x) = ',num2str(std(xbead)),' mu']);
    subplot(3,2,2);
    [n,xout] = hist(ybeadsdev1',100);
    stairs(xout,n); hold on;
    xlabel('y-median(y) (mu)');
    ylabel('#');

    % Plot the estimated drift as function of time
    subplot(3,1,2);
    plot(frames,xbeadsdev1); hold on;
    plot(frames,xshift,'g-','LineWidth',2);

    xlabel('frame #');
    ylabel('x-median(x) (mu)');
    subplot(3,1,3);
    plot(frames,ybeadsdev1); hold on;
    %     plot(frames,ysmoothdev,'r-'); hold on;
    plot(frames,yshift,'g-','LineWidth',2);
    xlabel('frame #');
    ylabel('y-median(y) (mu)');
    handles.xshift = xshift;
    handles.yshift = yshift;

%set(handles.estimatedriftpushbutton,'String',aux1);
%close(handles.figurapoints)


x = handles.x(:); 
y = handles.y(:);
alpha = handles.alpha(:);
fr = handles.fr(:);
%radius=handles.radius(:);

px_mu = str2num(get(handles.szpx,'String'));

%xdim=ceil(max(handles.x));
%ydim=ceil(max(handles.y));

xcorr = NaN(size(x));
ycorr = NaN(size(y));
frmax = max(fr);

% Corrected position coordinates
Npos = length(x);
if ~isfield(handles,'xshift')
    warndlg('You first need to estimate the drift !');
end
xshift = handles.xshift;
yshift = handles.yshift;
%if any(isnan(xshift)) || any(isnan(yshift))
 %   warndlg('There are NaNs in the computed shifts !');
%end

fr2 = fr(1:Npos);
xshift2 = xshift(fr2);
yshift2 = yshift(fr2);
xcorr = x(:) - xshift2(:);
ycorr = y(:) - yshift2(:);

fdrift=figure('Name','Pointillist image before and after drift correction');
 set(fdrift,'Position',get(0,'ScreenSize')*0.9 + 10); % maximize figure size
 
sp(1) = subplot(2,1,1);
plot(x,y,'k.','MarkerSize',2);
axis equal tight;
title('before drift correction');
sp(2) = subplot(2,1,2);
%plot(newxcorr,newycorr,'k.','MarkerSize',2);
plot(xcorr,ycorr,'k.','MarkerSize',2);
axis equal tight;
title('after drift correction');

linkaxes(sp,'xy');

button = questdlg('Do you want to accept the corrected positions ?','Please answer','Yes','No','Yes');
clear newxcorr newycorr

if strcmp(button,'Yes')
    disp('Drift corrected');
    
    handles.x = xcorr;
    handles.y = ycorr;
  % disp('done');
  
  %VOIR!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  
    % clean beads
    qstring='Delete signals from beads?';
    button = questdlg(qstring); 
    if strcmp(button,'Yes')
        totalbeads=handles.totalbeads;
        handles.x(totalbeads)=[];
        handles.y(totalbeads)=[];
        handles.alpha(totalbeads)=[];
        handles.fr(totalbeads)=[];
        handles.radius(totalbeads)=[];
    end
    
    % OJO!!!!!!!!!!!!
   % vect2remove=selectremovealpha(handles,xcorr,xcorr,handles.alpha);
   % if ~isempty(vect2remove)
    %    [x,y,alpha,fr]=removealpha(vect2remove,handles.x,handles.y,handles.alpha,handles.fr);
    %end
    x=handles.x;
    y=handles.y;
    fr=handles.fr;
    radius=handles.radius;
%    figurapoints=plotpoints(x,y,fr,handles);
   
    matrice_results(1,:)=fr';
    matrice_results(2,:)=y';
    matrice_results(3,:)=x';
    matrice_results(4,:)=alpha';
    matrice_results(5,:)=radius';
    
else
    disp('Positions left unchanged.');
end

close(fdrift)
close(handles.fig_bead)
%close(figurapoints)
clear x y alpha fr totalbeads xcorr ycorr xshift2 yshift2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliar functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Estimate drift from bead trajectories
function [xshift,yshift,xbeadsdev0,ybeadsdev0,xbeadsdev1,ybeadsdev1] = estimate_drift(frames,xbeads,ybeads,handles)

span = handles.smoothframes;

Nbeads = size(xbeads,1);
for ip=1:Nbeads
    xbead = xbeads(ip,:);
    ybead = ybeads(ip,:);
    xbeadsmed(ip) = nanmedian(xbead);
    ybeadsmed(ip) = nanmedian(ybead);
    xsmooth(ip,:) = smooth_ch(xbead,span,10);
    ysmooth(ip,:) = smooth_ch(ybead,span,10);
end
% subtract starting value
xbeadsdev0 = xbeads-repmat(xbeads(:,1),1,size(xbeads,2));
ybeadsdev0 = ybeads-repmat(ybeads(:,1),1,size(ybeads,2));

% subtract median value
xbeadsdev1 = xbeads-repmat(xbeadsmed(:),1,size(xbeads,2));
ybeadsdev1 = ybeads-repmat(ybeadsmed(:),1,size(ybeads,2));

% subtract starting value from x- and y-shifts
xsmoothdev = xsmooth-repmat(xsmooth(:,1),1,size(xsmooth,2));
ysmoothdev = ysmooth-repmat(ysmooth(:,1),1,size(ysmooth,2));

% Estimate the drift
xshift = nanmedian(xsmoothdev,1);
yshift = nanmedian(ysmoothdev,1);

% Interpolate drift across gap regions
frxshift0 = frames(~isnan(xshift));
xshift0 = xshift(~isnan(xshift));
xshift = interp1(frxshift0,xshift0,frames,'linear');
fr0 = frames(~isnan(yshift));
yshift0 = yshift(~isnan(yshift));
yshift = interp1(fr0,yshift0,frames,'linear');

% Extrapolate drift where needed
aux = find(isnan(xshift));
if ~isempty(aux)
    if aux(end)-aux(1)+1 ~= length(aux)
        warndlg('PROBLEM: contact Christophe !');
    else
        xshift(aux) = xshift(aux(1)-1);
    end
end
aux = find(isnan(yshift));
if ~isempty(aux)
    if aux(end)-aux(1)+1 ~= length(aux)
        warndlg('PROBLEM: contact Christophe !');
    else
        yshift(aux) = yshift(aux(1)-1);
    end
end
if any(isnan(xshift)) || any(isnan(yshift))
    warndlg('There are NaNs in the computed shift ! contact Christophe !');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute distance between pairs of beads and output into a format suitable for plotting
function [dxmat,dymat,dmat,pair] = dist_beads(xbeads,ybeads)
Nbeads = size(xbeads,1);
Nfr = size(xbeads,2);
%disp(['Computing distances between all pairs of the ',num2str(Nbeads),' beads..']);
dxmat = NaN(Nbeads,Nbeads,Nfr);
% TODO: do this without the for loop to speed things up
for ifr = 1:Nfr
    aux = repmat(xbeads(:,ifr),1,Nbeads);
    dxmat1(:,:,ifr) = aux-aux';
    aux = repmat(ybeads(:,ifr),1,Nbeads);
    dymat1(:,:,ifr) = aux-aux';
end
d2mat1 = dxmat1.^2+dymat1.^2;
dmat1 = sqrt(d2mat1);

Npairs = Nbeads*(Nbeads-1);
dxmat = NaN(Npairs,Nfr);
dymat = NaN(Npairs,Nfr);
dmat = NaN(Npairs,Nfr);
ip = 1;
for ib = 1:Nbeads-1
    for jb = ib+1:Nbeads
        %     [iaux,jaux] = ind2sub([Nbeads Nbeads],ip);
        dxmat(ip,:) = dxmat1(ib,jb,:);
        dymat(ip,:) = dymat1(ib,jb,:);
        dmat(ip,:) = dmat1(ib,jb,:);
        pair(ip,:) = [ib jb];
        ip = ip+1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Smoothing
function ys = smooth_ch(y,span,Nmin)
ys = movingwindowfilter(y,span,Nmin,@median);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filter over moving windows
function ys = movingwindowfilter(y,span,Nmin,fhandle)

ys = NaN(size(y));
N = length(y);
left = floor((span-1)/2);
right = left;
for i=1:N
    ileft = max(1,i-left);
    iright = min(i+right,N);
    window = ileft:iright;
    aux = y(window);
    aux1 = aux(~isnan(aux));
    if length(aux1)>Nmin
        ys(i) = fhandle(aux1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = flagmatrixelements(M,vectorindices)
siz = size(M);
aux = M(:);
aux(vectorindices) = 1;
M = reshape(aux,siz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
