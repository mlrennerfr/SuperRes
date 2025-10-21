function convertMTTpeaks
%function convertMTTpeaks
%
% Marianne Renner 01/2025 - adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
w=1.3;
opt(15)=5;
%Dini=detoptions(16); % Dpred
%a=detoptions(18); % size pixel
%persistance=detoptions(15);
%Te=detoptions(17);

currentdir=cd;
dialog_title=['Select data folder'];
path = uigetdir(cd,dialog_title);
if path==0
    return
end

handles.folder=path;
cd(path)

% dialog box to select files
d=dir('*.mat*');
st = {d.name};
[listafiles,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',st);
if v==0
    return
end
[f,ultimo]=size(listafiles); 

if isdir ('pk'); else; mkdir ('pk'); end

% loop over files
for cont=1:ultimo   
    % file names
    file=[st{listafiles(cont)}]
    [namefile,rem]=strtok(file,'.');  %raiz 
    structure=load(file)
    
    if isfield(structure,'matrice_results')
        peaks=structure.matrice_results';
    else
        peaks=[structure.fr  structure.x  structure.y  structure.alpha];
    end
    
    %disp(peaks)
    
    disp([num2str(size(peaks,1)),' peaks were detected in total']);
    save (['pk\',namefile,'.pk'], 'peaks','-ascii'); 
    
end

clear peaks 

%%%%%%%%%%%%%%%%%%%%%%
                
            
        
