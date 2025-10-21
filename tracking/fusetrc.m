function fusetrc
% function fusetrc
% fuses files of trajectories
%
% MR - jul 16 - 
% Marianne Renner 01/2025 - adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currentdir=cd;
start_path=[cd,'\trc'];
dialog_title=['Select data folder'];
directory_name = uigetdir(start_path,dialog_title);
if directory_name==0
    return
end
trcpath=directory_name;
k=strfind(trcpath,'trc');
path=(trcpath(1:k-1));

cd(trcpath)
%choose data
d = dir('*.con.trc*');
st = {d.name};
saveextension='.con.trc';

if isempty(st)==1
    d = dir('*.trc*');
    saveextension='.trc';
    st = {d.name};
    if isempty(st)==1
       msgbox(['No files!!'],'Select files','error')
       return
    end
end
[listafiles,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',st);
if v==0
   return
end
[f,ultimo]=size(listafiles);
cd(currentdir)  

for cont=1:ultimo   % toda la lista de archivos
    
  %trc  
  cd(trcpath)
  file=st{listafiles(cont)};
  
  [namefile,rem]=strtok(st{listafiles(cont)},'.');
  trc =load(file);                                        % load trc (x)
  
  if cont==1
      trctotal=trc;
  else
      maxtrc=max(trctotal(:,1))
      maxframe=max(trctotal(:,2))
      trctotal=[trctotal; trc(:,1)+maxtrc trc(:,2)+maxframe trc(:,3:size(trc,2))];
  end

end

[filename,path] = uiputfile([namefile,'.con.trc'],'Save file as :');
if isequal(filename,0) | isequal(path,0)
else
    save(filename,'trctotal','-ascii');
end
disp(['File ',filename,' saved'])

cd(currentdir)  
clear trc trctotal

% end of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

