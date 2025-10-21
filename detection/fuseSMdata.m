function fuseSMdata
% function fuseSMdata(handles)
% open pontillistic data files, fuse them (change frame number for the
% next one)
%
% MR sept 19 for SuperRes_v2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matrice_results=[];
currentdir=cd;

d = dir('*mat*'); % movie
st={d.name};
if isempty(st)==1
   msgbox(['No files!!'],'Select files','error');
   return
end
[files,v] = listdlg('PromptString','Select files to concatenate:','SelectionMode','multiple','ListString',st);
if v==0
   return
end

%first file
disp('Concatenating files:')
file=st{files(1)}
[namefile,rem]=strtok(file,'.'); %sin extension
dataall=loadSMdata(file);
lastframe=max(dataall.fr);
names=[];

for i=2:size(files,2)
        
    file2=st{files(i)}  
    names=[names;file2];

    dataSM=loadSMdata(file2);
    dataall.x=[dataall.x; dataSM.x];     
    dataall.y=[dataall.y; dataSM.y];
    dataall.alpha=[dataall.alpha; dataSM.alpha];
    dataall.fr=[dataall.fr; dataSM.fr + lastframe ];
    
    dataall.radius=[dataall.radius; dataSM.radius];
    dataall.sigma=[dataall.sigma; dataSM.sigma];
    dataall.blink=[dataall.blink; dataSM.blink];
    dataall.ratio=[dataall.ratio; dataSM.ratio];
    dataall.z=[dataall.z; dataSM.z];
    dataall.test1=[dataall.test1; dataSM.test1];
    dataall.test2=[dataall.test2; dataSM.test2];

    lastframe=max(dataall.fr);

    clear dataSM

end %loop files

qstring=['Confirm?'];
button = questdlg(qstring); 

if strcmp(button,'Yes')
    % save data
    matrice_results(1,:)=dataall.fr';
    matrice_results(3,:)=dataall.x';
    matrice_results(2,:)=dataall.y';
    matrice_results(4,:)=dataall.alpha';
    matrice_results(5,:)=dataall.radius';
    matrice_results(6,:)=dataall.sigma';
    matrice_results(7,:)=dataall.blink';
    matrice_results(8,:)=dataall.ratio';
    matrice_results(9,:)=dataall.z';
    matrice_results(10,:)=dataall.test1';
    matrice_results(11,:)=dataall.test2';
    
    %pk
    pk(:,1)=matrice_results(1,:); 
  %  pk(:,2)=matrice_results(2,:);  
  %  pk(:,3)=matrice_results(3,:);    
    pk(:,2)=matrice_results(3,:);  
    pk(:,3)=matrice_results(2,:);    
    
    pk(:,5)= matrice_results(4,:); 
    pk(:,4)=matrice_results(5,:)*2; 
    pk(:,7)=matrice_results(6,:);  
    pk(:,6)=matrice_results(7,:);
    
    if size(matrice_results,1)>7
        pk(:,8)=matrice_results(8,:); %ratio
        pk(:,9)=matrice_results(9,:); %z
        pk(:,10)=matrice_results(10,:); %test
        pk(:,11)=matrice_results(11,:); %test
        
        if max(pk(:,9))==0 && min(pk(:,9))==0 % no values for z
            if isdir('pk'); else mkdir('pk');end
            pkpath=['pk'];
            savename=[namefile, '.pk'];
        else
            if isdir('pk3'); else mkdir('pk3');end
            pkpath=['pk3'];
            savename=[namefile, '.pk3'];
        end
    else
        if isdir('pk'); else mkdir('pk');end
        pkpath=['pk'];
        savename=[namefile, '.pk'];
    end
    
    save([namefile, '-concat.mat'], 'matrice_results'); % one image, corrected for stage drift
    
    cd(pkpath);
    save([namefile, '-concat.pk'], 'pk','-ascii');
    cd(currentdir)
    
    disp(['Saved as : ',namefile,'-concat.mat and ',namefile, '-concat.pk'])
    disp(' ')
end % confirm



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




