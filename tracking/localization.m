function localization(domainfile,pn,fn,handles); 
% function localization(domainfile,pn,fn,handles); 
% assign localization to trajectories (Scripts MVE)
%
% Marianne Renner mar 09 for SPTrack.m v4.0                     MatLab 7.00
% Marianne Renner 01/2025 - adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

waitbarhandle=waitbar( 0,'Please wait...','Name',['Localizing trajectories over ',fn]);

%Image=handles.synimage;
currentdir = cd;
trajfolder=[pn,'\traj'];
cd(trajfolder);
fntraj=[fn,'.traj'];
%load(fntraj,'-mat');

fit = loadfit(fntraj);

%current=fit(1).new_spot;
%if isfield(fit,'new_spot')
   current=fit.new_spot;
%else
 %   current=fit.spot;
%end
nb_frames=recadrage.Nz ;
if isfield(handles,'difparameters')
    difpar=get(handles.difparameters,'userdata');  % other values than default
    perisyn=str2num(difpar{4});
else
    perisyn=handles.perisyn;
end

[img_zone,nb_syn]=numerotesyn(handles.synimage,perisyn);
current=rajoutzone(current,nb_frames,img_zone,perisyn,waitbarhandle);

%fit(1).new_spot=current;
fit.new_spot=current;
%save(fntraj,'source','parametres','fit','-mat');
save(fntraj,'source','information','recadrage','parametres','fond','nb_fits','nofit','fit','-mat');

%disp(' '); 
disp('Localization done'); disp(' ');

cd(currentdir);
close(waitbarhandle);
clear current fit img_zone

% end of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%