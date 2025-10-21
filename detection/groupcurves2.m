function groupcurves2
% function groupcurves2
% groups together calibration curves
% fits the resulting group of points
% saves the total list of points and the fit
%
% MR juin 17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% choose folder
dialog_title=['Select data folder'];
path = uigetdir(cd,dialog_title);
if path==0
    return
end
cd(path)

%choose data files
d = dir('*curvecalib.txt*');
st = {d.name};
if isempty(st)==1
   msgbox(['No files!!'],'Select files','error')
   return
end
[listafiles,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',st);
if v==0
   return
end

[f,last]=size(listafiles);
totaldata=[];

%loop files
for cont=1:last
    file=st{listafiles(cont)};
    [namefile,rem]=strtok(st{listafiles(cont)},'.');
    
    data=load(file);
    totaldata=[totaldata; data];
end
totaldata=sortrows(totaldata,1); %sort by pos z

plot(totaldata(:,2),totaldata(:,1),'*r')

% options
prompt = {'Keep (min max): ','Pixel size (nm):'};
num_lines= 1;
dlg_title = 'Calibration curves';
def = {'-400 400','160'}; % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
if exit(1) == 0;
   return; 
end
limits=sort(str2num(answer{1}));
szpx=str2num(answer{2});

newtotaldata=[];
indexmin=find(totaldata(:,1)>=limits(1));
if isempty(indexmin)==0
    indexmax=find(totaldata(indexmin,1)<=limits(2));
    if isempty(indexmax)==0   
        newtotaldata=totaldata(indexmin(indexmax),:);
    end
end

%fit
if isempty(newtotaldata)==0
    
    %cftool (newtotaldata(:,2),newtotaldata(:,1))

    [fitresult, gof] = createFitcalib(newtotaldata(:,2), newtotaldata(:,1));
    
    %y = feval(fitresult,2.1554)

    result(:,1)=newtotaldata(:,2);
   % result(:,2)=newtotaldata(:,1)/1000; % in µm
    result(:,2)=newtotaldata(:,1)/szpx; % in pk
    
    save('datacalibration3D.txt','result','-ascii')
    save('calibration3D.mat','fitresult','-mat')
    disp('Files saved')
    
else
    disp('No data left')
end

%EOF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%EOF


