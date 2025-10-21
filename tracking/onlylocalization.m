function onlylocalization(handles)
% function onlylocalization(handles)
% assign localization to trajectories (Scripts MVE)
%
% Marianne Renner mar 09 for SPTrack.m v4.0                     MatLab 7.00
% Marianne Renner jan 2022 for SPTrack_v6
%
% Marianne Renner 01/2025- adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currentdir=cd;
%filedefnames=get(handles.filedefin,'userdata');

% dialog boxs to enter acquisition data
prompt = {'Identifier for localization files: '};
num_lines= 1;
dlg_title = 'Enter';
def = {'-loc'}; % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
   if exit(1) == 0
       return; 
   end
domaindef=answer{1}; %

%difpar=get(handles.difparameters,'userdata');  % other values than default
%perisyn=str2num(difpar{4});

% carga files y loop general 
dialog_title='Select data folder';
path = uigetdir(cd,dialog_title);
if path==0
    return
end
datapath=[path,'\traj'];
cd(datapath)

%choose data
d = dir('*traj*');
st = {d.name};
if isempty(st)==1
   msgbox('No .traj files!!','Select files','error');
   return
end
[listafiles,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',st);
if v==0
   return
end

[f,ultimo]=size(listafiles);

for cont=1:ultimo
    cd(datapath)
    %pn=st{listafiles(cont)};
    [namefile,rem]=strtok(st{listafiles(cont)},'.');
    cd(path);
    domainfile=[namefile,domaindef,'.tif'];
    if length(dir(domainfile))>0
       % reads files
       disp(['Domain file: ',domainfile]);
       Image=imread(domainfile);
       Image=double(Image);
       handles.synimage=Image;
       localization(domainfile,path,namefile,handles);                      %localization MVE
    else
       disp(['Domain file: ',domainfile]);
       disp('Domain file not found. Without localization');
    end
end

cd(currentdir);
clear Image

% end of file

