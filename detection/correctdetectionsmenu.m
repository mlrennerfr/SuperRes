function correctdetectionsmenu
% function correctdetectionsmenu
% launchs correction for multiple detections
% with space and time criteria
%
% MR 2019 - SuperRes_v2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currentdir=cd;

% dialog box to enter recognition criteria for images
prompt = {'Max. distance between detections (SMLM pixels):','Maximum blinking period (frames):'};
num_lines= 1;
dlg_title = 'Correction of multiple detections';
def = {'0.5','300'};
answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
if exit(1) == 0
    return;
end
dist=str2num(answer{1});
maxblink=str2num(answer{2});

%files
d=dir('*.mat*'); % .stk files
st={d.name};
if isempty(st)==1
    msgbox(['No files!!'],'','error');
    return
end

%choose data
[listafiles,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',st);
if v==0
     return
end


for nromovie=1:size(listafiles,2)
    matfile=st{listafiles(nromovie)} 
   disp=['Correction of multiple detections: File ',matfile];
   disp(' ');
    [filename,rem]=strtok(matfile,'.');
    correctedfile=[filename,'-corr.mat']  ;
    
    aux=correctdetections4(matfile,dist,maxblink);
    
    corrpeak=aux(find(isnan(aux(:,2))==0),:);
    clear aux   
    
    % save new data
    if isempty(corrpeak)==0
        matrice_results(1,:)=corrpeak(:,1);
        matrice_results(2,:)=corrpeak(:,3);
        matrice_results(3,:)=corrpeak(:,2);
        matrice_results(4,:)=corrpeak(:,4);
        matrice_results(5,:)=corrpeak(:,5);
        matrice_results(6,:)=corrpeak(:,6);
        matrice_results(7,:)=corrpeak(:,7);
        matrice_results(8,:)=corrpeak(:,8);
        matrice_results(9,:)=corrpeak(:,9);
        matrice_results(10,:)=corrpeak(:,10);
        matrice_results(11,:)=corrpeak(:,11);
        savename=[filename, '-corr.mat'];
        save(savename, 'matrice_results'); % 
        
        if isdir('pk'); else; mkdir('pk'); end
        cd('pk')
        save([filename,'.corr.pk'], 'corrpeak','-ascii'); % 
        cd(currentdir)
        
        disp=['Corrected files ',filename,'-corr.mat and -corr.pk saved'];
        disp(' ');

        clear matrice_results corrpeak
    end
    

end %loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

