function convertMTTtrcnotrack
%function convertMTTtrcnotrack
%
% Marianne Renner - adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
w=1.3;
opt(15)=5;
%Dini=detoptions(16); % Dpred
%a=detoptions(18); % size pixel
%persistance=detoptions(15);
%Te=detoptions(17);

currentdir=cd;

% dialog box to enter parameters 
prompt = {'Dini:','Pixel size (nm):','Time between images (ms):','NA:','Wavelength:'};
num_lines= 1;
dlg_title = 'Dwell time analysis';
def = {'0.01','160','50','1.45','647'}; % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
if exit(1) == 0;
       return; 
end
opt(16)=sort(str2num(answer{1}));
opt(18)=str2num(answer{2});
opt(17)=str2num(answer{3});
opt(23)=str2num(answer{4});
opt(24)=str2num(answer{5});

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

if isdir ('trc'); else; mkdir ('trc'); end
if isdir ('traj'); else; mkdir ('traj'); end
if isdir ('pk'); else; mkdir ('pk'); end

% loop over files
for cont=1:ultimo   
    % file names
    file=[st{listafiles(cont)}]
    [namefile,rem]=strtok(file,'.');  %raiz 
    structure=load(file);
    
    if isfield(structure,'matrice_results')
        peaks=structure.matrice_results';
    else
        peaks=[structure.fr  structure.x  structure.y  structure.alpha];
    end
    
    %disp(peaks)
    
    disp([num2str(size(peaks,1)),' peaks were detected in total']);
    save (['pk\',namefile,'.pk'], 'peaks','-ascii'); 
    nroframes=max(peaks(:,1));
    
    for m=1:nroframes
        newpk= [];
       index=find(peaks(:,1)==m);
       if isempty(index)==0
           newpk=peaks(index,:);
           j=1;
           for i=1:size(newpk,1)
               objet(j).centre=newpk(i,2:3);
               j=j+1;
           end 
           plan(m).objet=objet;
           plan(m).Nb_objets=size(newpk,1);
           clear objet
       else % no peaks
           plan(m).objet=[];
           plan(m).Nb_objets=0;
       end
   end

   % tracking
   traj=initialtrackerMTT(file,plan,nroframes, opt, handles);
   
   save(['trc\',namefile,'.trc'],'traj','-ascii'); %
   
   if (length(traj)>0) 
      nrotrc=max(traj(:,1));
   else
      nrotrc=0;
   end
   disp(['and the initial tracking constructed ' Num2str(nrotrc) ' trajectories.']); disp('  ');

    
end

clear peaks traj

%%%%%%%%%%%%%%%%%%%%%%
                
            
        
