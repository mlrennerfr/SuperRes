function compareGT

prompt = {'Pixel size (nm):','Localization accuracy:','mindens cluster','mindens px',...
    'min nro detec','minsize cluster','Int threshold'};
num_lines= 1;
dlg_title = 'Simulation parameters 1';
def = {'160','10','0.1','0','4','2','25'};


% multiple detections
%def = {'160','10','0.1','0.01','20','25','20'};  test scenario 2
%def = {'160','10','0.1','0.1','10','5','0'};  test scenario 2 geom3.2 très bizarre
%def = {'160','10','1','0.1','30','60','25'};  test scenario 3
% def = {'160','10','0.1','0','10','70','20'}; test scenario 4 : mauvais!
%def = {'160','10','0.01','0.5','30','50','20'};  test scenario 5 situation bizarre
%def = {'160','10','0.01','0','30','70','20'}; test scenario 6
%def = {'160','10','0.01','0','50','120','8'};  test scenario 7 moche et bizarre
%def = {'160','10','0.01','0','10','50','20'}; test scenario 8 pas terrible
% def = {'160','10','0.8','0.7','10','75','15'}; test scenario 9 bof
% def = {'160','10','0.3','0.3','10','65','15'}; test scenario 10


% ground truth
%def = {'160','10','0.1','0','5','25','20'};  test scenario 2
%def = {'160','10','0.1','0','10','25','30'};  test scenario 3
%def = {'160','10','0.1','0','4','2','25'};  test scenario 4
%def = {'160','10','0.1','0','7','7','23'}; test scenario 5
%def = {'160','10','0.01','0','30','70','20'};  test scenario 6
%def = {'160','10','0.01','0','5','5','8'}; test scenario 7
%def = {'160','10','0.12','0','5','5','10'}; test scenario 8
%def = {'160','10','0.05','0','10','20','15'}; test scenario 9
%def = {'160','10','0.01','0','5','5','25'};  test scenario 10


answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
if exit(1) == 0
       return; 
end
handles.szpx=str2num(answer{1})/1000  ; 
handles.dx=str2num(answer{2})/1000  ; 

handles.mindens=str2num(answer{3});
handles.mindenspx=str2num(answer{4});
handles.mindetect=str2num(answer{5});
handles.minsize=str2num(answer{6});
handles.inthresh=str2num(answer{7});

% comparison to ground truth
GTdata=[];
handles.GTdata=[];

%load GT data
[fileGT,path] = uigetfile('*.csv','Load ground truth(.csv)');
if path>0
    filenameGT=[path,fileGT]
    [namefile,rem]=strtok(filenameGT,'.');
    handles.GTdata=csvread(filenameGT,1,0);
end

handles.GTdata(:,1)=handles.GTdata(:,1)+abs(min(handles.GTdata(:,1)));
handles.GTdata(:,2)=handles.GTdata(:,2)+abs(min(handles.GTdata(:,2)));

 %load test data
[file,path] = uigetfile('*.mat','Load detection data (.mat)');
if path>0
    filename=[path,file];
    handles=loadselectPALMdata(filename,1,handles); %data already in microns
end

X_mu=handles.x ;
Y_mu=handles.y ;
alpha = handles.alpha;
alpha(1:10);   
xdim=ceil(max(X_mu)/handles.dx);
ydim=ceil(max(Y_mu)/handles.dx);
        
% Create rendered image
[Irend,xxi,yyi] = PALM_rendering3(X_mu,Y_mu,alpha,handles.dx*2,handles.dx,0,xdim, ydim, 1);

%figure
%imshow(Irend,'InitialMagnification', 'fit')
%hold on
%plot(X_mu,Y_mu,'.r')

[clusterdata, dataclustersD, dataclustersGT,elapsedtime]=diinamicGT(handles,Irend, namefile,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





