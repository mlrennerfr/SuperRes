function convertMTTtrc
% function convertMTTtrc
%
% Convertion from MTT formatto SuperRes format
%
% Marianne Renner - adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


currentdir=cd;



dialog_title=['Select data folder'];
path = uigetdir(cd,dialog_title);
if path==0
    return
end

cd(path)

% dialog box to select files
d=dir('*.mat*');
st = {d.name};
[listafiles,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',st);
if v==0
    return
end
[f,ultimo]=size(listafiles); 

if isdir ('trc'); else; mkdir ('trc'); end
if isdir ('traj'); else; mkdir ('traj'); end
if isdir ('pk'); else; mkdir ('pk'); end

% loop over files
for cont=1:ultimo   
    % file names
    file=[st{listafiles(cont)}]
    [namefile,rem]=strtok(file,'.')  %raiz 
    structure=load(file)
        %    save([stack,'.mat'],'x','y','alpha','fr','Xmatrix','Ymatrix','alphamatrix','frmatrix');
    
        frames=structure.fr
    trc=[frames  structure.x structure.y structure.alpha]
    for i=1:size(trc,1)
        trc=[i trc(i,:)];
    end
    
    disp(trc)
    
    save(['trc\',namefile,'.con.trc'],'reconnectedfile','-ascii');

    [fit]=creastruct(trcfilename);
    
    writetraj(fn, fit, frames, typefile, detoptions, handles)
    
end

clear datamatrix stackdata stack_info results resultscorr

%%%%%%%%%%%%%%%%%%%%%%
                
            
        
