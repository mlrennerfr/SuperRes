function convertcsvTStomat
%function convertcsvTStomat
%
%converts ThunderSTORM-type data (.csv) into .mat format usable by Diinamic
% 
%input: 
%.csv file with the following data in column number:
%    col1=id
%    col2=fr
%    col3=x
%    col4=y
%    col5=z    
%    col6=sigma1
%    col7=sigma2
%    col8:intensity
%    col9:offset
%    col10
%    col11:chi
%    col12: uncertainty
%
% output:
%.mat file with the structure (all sizes in µm)
%
%    matrice_results(1,:)=fr;
%    matrice_results(3,:)= x; 
%    matrice_results(2,:)= y;
%    matrice_results(4,:)=alpha;
%    matrice_results(5,:)=radius;  
%    matrice_results(6,:)=sigma;  
%    matrice_results(7,:)=blink;  
%    matrice_results(8,:)=ratio;  
%    matrice_results(9,:)=z;  
%    matrice_results(10,:)=test1;  
%    matrice_results(11,:)=test2;
%
% Marianne Renner 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


d = dir('*.csv'); %
st={d.name};
if isempty(st)==1
   msgbox('No files!!','Select files','error');
   return
end
[files,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',st);
if v==0
   return
end

for i=1:size(files,2)   
    file=st{files(i)}
    [namefile,rem]=strtok(file,'.');
    savename=[namefile,'.mat'];
    
    data = csvread(file,1,0);

    x=data(:,3)/1000; % in µm
    y=data(:,4)/1000; % in µm
    alpha=data(:,8);
    fr=data(:,2);
    radius=data(:,6)/2;
    sigma=data(:,12);
    blink=ones(size(data,1),1); 
    ratio=ones(size(data,1),1);  %data(:,6)/data(:,7);
    z=data(:,5)/1000; % in µm
    test1=data(:,11);
    test2=zeros(size(data,1),1);

    matrice_results=zeros(11,size(data,1));
    matrice_results(1,:)=fr;
    matrice_results(3,:)= x; 
    matrice_results(2,:)= y;
    matrice_results(4,:)=alpha;
    matrice_results(5,:)=radius;  
    matrice_results(6,:)=sigma;  
    matrice_results(7,:)=blink;  
    matrice_results(8,:)=ratio;  
    matrice_results(9,:)=z;  
    matrice_results(10,:)=test1;  
    matrice_results(11,:)=test2;

    save(savename, 'matrice_results');

    clear data
    clear matrice_results

end

clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

