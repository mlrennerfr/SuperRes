function localizeSM(handles)
% function localizeSM(handles)
% assign localization to trajectories 
%
% Marianne Renner mar 09 for SPTrack.m v4.0                     MatLab 7.00
% Marianne Renner jan 2022 for SPTrack_v6
% Marianne Renner for SuperRes_v4_app feb 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currentdir=cd;

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

% carga files y loop general 
%dialog_title='Select data folder';
%path = uigetdir(cd,dialog_title);
%if path==0
%    return
%end
%datapath=[path,'\trc'];
datapath=[currentdir,'\trc'];
cd(datapath)

%choose data
d = dir('*trc*');
st = {d.name};
if isempty(st)==1
   msgbox('No .trc files!!','Select files','error');
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
    file=st{listafiles(cont)};
    [namefile,rem]=strtok(file,'.');
    trc=load(file);

    cd(currentdir);
    domainfile=[namefile,domaindef,'.tif']
    waitbarhandle=waitbar( 0,'Please wait...','Name',['Localizing trajectories over ',file]);

    if length(dir(domainfile))>0
        % reads files
        disp(['Domain file: ',domainfile]);
        Image=imread(domainfile);
        Image=double(Image);
        handles.synimage=Image;
       
        [img_zone,nb_syn]=numerotesyn(handles.synimage,1);
       % current=rajoutzone(current,nb_frames,img_zone,perisyn,waitbarhandle);

        addcol=6-size(trc,2);
        if addcol>0
            for j=1:addcol
                trc=[trc zeros(size(trc,1),1)];
            end
        end

        for point=1:size(trc,1)

            posx=round(trc(point,3));
            posy=round(trc(point,4));
            if posx<1; posx=1; end
            if posy<1; posy=1; end
            if posx>size(img_zone,1); posx=size(img_zone,1); end
            if posy>size(img_zone,2); posy=size(img_zone,2); end

          %  trc(point,6)=img_zone(posx,posy);
            trc(point,6)=handles.synimage(posy,posx);
        end %for



        % SMOOTH !!!!!



        

        cd(datapath)
        save([namefile,'.loc.trc'],'trc','-ascii')
        disp('Localization done'); disp(' ');

        cd(currentdir)

    else
       disp('Localization not done, file not found'); disp(' ');

    end % with domain


    cd(currentdir);
    close(waitbarhandle);
    clear img_zone

end %loop files

cd(currentdir);
clear Image

% end of file

