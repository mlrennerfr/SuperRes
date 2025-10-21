function onlylocalizationSM
% function onlylocalizationSM
% assign localization to trajectories (Scripts MVE)
%
% Marianne Renner jul 2016
% Marianne Renner 01/2025 - adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currentdir=cd;
% dialog box 
prompt = {'Identifier for localisation file','Size perisynaptic ring (pixels):'};
num_lines= 1;
dlg_title = 'Enter';
def = {'-binary','2'}; % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
if exit(1) == 0;
   return; 
end
domaindef=answer{1};
perisyn=str2num(answer{2});

% carga files y loop general 
dialog_title=['Select data folder'];
path = uigetdir(cd,dialog_title);
if path==0
    return
end
datapath=[path,'\trc'];
cd(datapath)

%choose data
d = dir('*.trc*');
st = {d.name};
if isempty(st)==1
   msgbox(['No files!!'],'Select files','error');
   return
end
[listafiles,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',st);
if v==0
   return
end

[f,ultimo]=size(listafiles);

for cont=1:ultimo
    cd(datapath)
    file=st{listafiles(cont)};
    [namefile,rem]=strtok(st{listafiles(cont)},'.');
    cd(path);
    domainfile=[namefile,domaindef,'.tif']
    cd(datapath)
    traces=load(file);
    if size(traces,2)<6
        traces=[traces zeros(size(traces,1),1)];
    end
    cd(path);
    if length(dir(domainfile))>0
       % reads files
       disp(['Domain file: ',domainfile]);
       Image=imread(domainfile);
       Image=double(Image);
       
              figure
              stackmin=min(min(Image));
              stackmax=max(max(Image));
       imshow(Image,[stackmin stackmax],'InitialMagnification','fit')

       %numbering
       level1 = graythresh(Image);
       ImageBW = im2bw(Image,level1);
       
      % figure
      % imshow(ImageBW,'InitialMagnification','fit')
       
       [img_dist,pos_syn_proche]=bwdist(ImageBW);
       [labeled,nb_dom] = bwlabel(ImageBW,4);
       statimage=regionprops(labeled,'Centroid');

       img_zone=labeled;
       y=find(img_dist<=perisyn & img_dist~=0);
       img_zone(y)=-labeled(pos_syn_proche(y)); %perisyn: negative value
       
       labeled=img_zone; %with perisyn
       
       %save labeled image!!!!!!
       save([namefile,domaindef,'-num.mat'],'labeled','-mat');

       
       for i=1:size(traces,1)
           posx=ceil(traces(i,3));
           if posx==0; posx=1; end
           if posx>size(labeled,2); posx=size(labeled,2); end
           posy=ceil(traces(i,4));
           if posy==0; posy=1; end
           if posy>size(labeled,1); posy=size(labeled,1); end
           
           if labeled(posy,posx)>0
               nrodom=labeled(posy,posx);
               traces(i,6)=nrodom;
           %    areadom(nrodom,2)=statimage(nrodom).Area;
            %   areadom(nrodom,3)=areadom(nrodom,3)+1; % one more traj in this domain
           end
       end
       
       figure
       imshow(labeled,'InitialMagnification','fit')
       hold on
       for i=1:nb_dom
           posx=statimage(i).Centroid(1);
           posy=statimage(i).Centroid(2);
           text(posx,posy,num2str(i),'Color',[0 0 1])
       end
       cd(datapath)
       save([namefile,'.loc.trc'],'traces','-ascii');
       
    else
       disp(['Domain file not found. Without localization']);
    end
end

cd(currentdir);
clear Image

disp('Done')

% end of file

