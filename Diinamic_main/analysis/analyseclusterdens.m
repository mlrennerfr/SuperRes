function [auxdens, Densparam, stopanalysis]=analyseclusterdens(handles)
% function [auxdens, Densparam, stopanalysis]=analyseclusterdens(handles)
%
% Analysis of clustering for Diinamic. One file each time.
%
% Marianne Renner 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% analysis code : parameters for density analysis
pxdilate=str2num(get(handles.pxdilate,'string'));
pxerode=str2num(get(handles.pxerode,'string'));
limitdens=str2num(get(handles.mindens,'string'));
minnrodetect=str2num(get(handles.minnrodetec,'string'));
mindiamcluster=str2num(get(handles.mindiam,'string')); %/(dx*1000) %!!!!!!!
maxdiamcluster=str2num(get(handles.limitdiam,'string')); %/(dx*1000)
maxnumberdet=str2num(get(handles.maxnumberdet,'string'));
percinthresh=str2num(get(handles.intensthresh,'string'));% intensity thershold (%)
polsize=str2num(get(handles.polsize,'string'));
voro=get(handles.vororadiobutton,'Value');
stopanalysis=0;
filename=handles.file;
namefile=get(handles.namefile,'String');

if size(namefile,2)==0
    [namefile,~]=strtok(filename,'.') %sin extension
    set(handles.namefile,'String',namefile);
end
dx = str2num(get(handles.PALMszpx,'String'));
sigmaloc = dx*2;
szpx = str2num(get(handles.szpx,'String'));

%if ROIs, plot ROIs
if exist([namefile,'.rgn'])==2
    data=importdata([namefile,'.rgn']);
    if isstruct(data)
        roifile=data.coord;
        indata=data.in; %indexes points
        Irend=data.imagerend; %rendered of the ROI
        dx=data.dx;
        if isfield(data,'dist')
            % distance
            distance=data.dist;
        end
    else
        roifile=data;
    end %if sstruct
    disp(['Loading regions from ',namefile,'.rgn'])
else
    % roifile= whole image
    % not implemented    
end % of exist    

% localization with respect to a second image
localizcode=0;
%--------------------------------------------------------------------
disp(' ');
disp('Analysing clustering by density');

resultsout=[];
resdistance=[];
distance=[];

% ROIs
disp([num2str(size(roifile,2)),' rois'])
    
for roij=1:size(roifile,2) %all ROIs in the image
    
    if stopanalysis==0
            disp(' ')
            disp(['Mask in ROI # ', num2str(roij)])  

            %select data ROI
            listeclucercle=[];
            resultsclu=[];
            roiselectedx=[];
            roiselectedy=[];
            frselected=[];
            alphaselected =[];
            
            % crop image by ROI
            BW=zeros(size(handles.I,1),size(handles.I,2));
            
           % disp(roifile)
            
            xi=roifile{roij}(:,1)/dx; % coordinates of the ROI in PALM pixels size
            yi=roifile{roij}(:,2)/dx; 
            
            in= indata{roij};
            
            minxi=floor(min(xi));
            maxxi=ceil(max(xi));
            minyi=floor(min(yi));
            maxyi=ceil(max(yi));  
            minyi=max(1,floor(min(yi)));
            maxyi=min(size(handles.I,1),ceil(max(yi)));

            if minxi<1; minxi=1; end % to avoid points out of the image size
            if maxxi>size(handles.I,2); maxxi=size(handles.I,2);  end
            if minyi<1; minyi=1; end
            if maxyi>size(handles.I,1); maxyi=size(handles.I,1) ; end  
            
            roi=roipolyold(BW,xi,yi);   % mask of the ROI (the same size than the rendered image)
            imageroi=immultiply(roi,handles.I); % picks up pixels of the rendered image that corresponds to the ROI
            imageroi=imageroi(minyi:maxyi,minxi:maxxi); % crops to the size of the ROI
            
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VOIR!!!!!!!!!
           %  ROI of the localization image
            synroi=[];
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %rendered ROI
            imageroirend=Irend{roij}; 

            %pick points ROI
            roiselectedx=data.xselec{roij}; %index of selected points
            roiselectedy=data.yselec{roij}; %index of selected points
            
            % frames and intensities
            frselected=[frselected; handles.fr(in)];             
            alphaselected=[alphaselected;handles.alpha(in)];  
            
            if isempty(roiselectedx)==1
                disp('No data')
                return
            else
                varargin{1} = handles.I;
                varargin{2} = []; %roifile;
                varargin{3} = handles.fr;
                varargin{4} = handles.alpha;
                varargin{5} = filename;  
                varargin{6} = dx;
                varargin{7} = szpx;
                
                if length(dir('auxiliar.mat'))>0  
                    delete auxiliar.mat        %clears previous results
                end
                
                varargout=DetectClusters(varargin); 
                uiwait;   
              
                %read results (in auxiliar.mat) 
                control=0;
                aux='auxiliar.mat'; 

                if length(dir(aux))>0  
                    load(aux,'-mat')
                    
                        resultsout=[resultsout; roij*ones(size(resultsclu,1),1) resultsclu];
                    
                        % # clusters per distance in ROI
                         if isempty(distance)==0     % ATT distance en pixels!!!!!!!!!!!!
                            resdistance=[resdistance; roij size(resultsclu,1) distance{roij} size(resultsclu,1)/distance{roij} ];
                                 %  nro roi - number of clusters - length - number of clusters/length
                         end
                    
                        imagename=[namefile,'-roi',num2str(roij),'-maskrend.mat'];
                        save(imagename,'imagemask','-mat') 
                    
                
                     delete 'auxiliar.mat'
                end % empty results from maskrendered (aux)
            end % empty roiselected
    end %stopanalysis
end % loop roij
 
% parameters

auxdens{1}=resultsout;
auxdens{2}=resdistance;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newBWlocmatrix=loadloc(locfilename,sizerend)
% load localization file 
 
newBWlocmatrix=[];  

if isempty(dir(locfilename))==1
    disp('No localization file found');
    return
end

locmatrix=double(imread(locfilename));

%binarize
level = graythresh(locmatrix);
BWlocmatrix = im2bw(locmatrix,level);

%convert to rendered size
newBWlocmatrix = imresize(BWlocmatrix,[sizerend(1),sizerend(2)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
