function cleantrcSM
% function cleantrcSM
% deletes selected trajectories
%
% MR - mar 06 - v 1.0       for trackdiffusion.m            MatLab6p5p1
% MR - mar 09 - for SPTrack v4.0                             MatLab7.00
% MR SuperRes jul 2016
% MR add the posibility to use a previously created mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

controlf=1;
currentdir=cd;
cut=1;

% dialog box 
prompt = {'Use saved mask? (0: no, 1: yes)','Localize over domains? (0: no, 1: yes)','Identifier for background file (domain image to do the localization)','Size perisynaptic ring (pixels):'};
num_lines= 1;
dlg_title = 'Cleaning options';
def = {'0','0','-binary','2'}; % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
if exit(1) == 0;
   return; 
end
previousmask=str2num(answer{1});
localize=answer{2};
identif=answer{3};
perisyn=str2num(answer{4});

%folder and data
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


%--------------------------------------------------------------------------
for cont=1:ultimo   % toda la lista de archivos
    
  %trc  
  cd(trcpath)
  file=st{listafiles(cont)};
  [namefile,rem]=strtok(st{listafiles(cont)},'.');
  namechar=(size(namefile,2));
  x =load(file);                                        % load trc (x)
  disp(['File ' ,file, ' loaded.']);
  cd(currentdir)  
  control=1;

  dicfile=[currentdir,'\',namefile,identif,'.tif']   %looks for background image
  if length(dir(dicfile))==0
      % black image
      datamatrix=zeros(ceil(max(x(:,4)))-floor(min(x(:,4))), ceil(max(x(:,3)))-floor(min(x(:,3))));
      stack_info.x=ceil(max(x(:,4)))-floor(min(x(:,4)));
      stack_info.y=ceil(max(x(:,3)))-floor(min(x(:,3)));
  else
      disp(['Background image:', dicfile]);
   %  [stack_info,datamatrix] = tifdatareadclean(dicfile);
     [stack_info,datamatrix] = tifdataread(dicfile);
     if isfield(datamatrix,'data')
         datamatrix=datamatrix.data;
     end
  end
  
  mask=zeros(size(datamatrix,1),size(datamatrix,2));
  figclear=figure;
  axis ([0 stack_info.y 0 stack_info.x]);
  otra=1;
  firstime=1;
  control=1;
  count = 1;
  fila = 1;
  newtrc = x;  % archivo traces: trabajo sobre el auxiliar hasta ultimo momento
  if isdir('clean\trc\');else; mkdir ('clean\trc\');end

  if previousmask==0
      
      while otra==1;   %loop general limpieza
    
        stackmin=(min(min(min(datamatrix))));
        stackmax=(max(max(max(datamatrix))));
        imshow((datamatrix(:,:,1)),[stackmin stackmax],'InitialMagnification','fit');
        hold on
        aux=[];
        del=[];
  
        if isempty(newtrc)==0
           for m=1:max(newtrc(:,1))
               indice=find(newtrc(:,1)==m);
               if isempty(indice)==0
                  for i = 1:size(indice,1)
                      graph(i,:)=newtrc(indice(i),:);   % archivo auxiliar con los puntos de cada trayectoria
                  end
                  plot ((graph (:,3)), (graph (:,4)), 'b-');   % grafica traces 
                  hold on
               end
               graph=[];
           end
         end %newtrcempty

        % dialog box to enter new data
        if firstime==0;
           qstring=['more areas?'];
           button = questdlg(qstring); 
           if strcmp(button,'Yes')
              otra=1;
           else 
              otra=0;
              break
           end
        end
        firstime=0;

        if otra==1
            
            %rutina limpieza
            [BW,xi,yi]=roipolyold;    %seleccion ROI
            % crea un nuevo archivo trc sin las moleculas que estan dentro del area
            % seleccionada
            count=1;
            ind=1;
            del=[];
            selectrace=[];
            newselectrace=[];
            vectormol=[];
        
            [newselectrace,aux]=pickpointsclear2(xi,yi,newtrc);
        
            if isempty(newselectrace)==0
                for t=1:max(newselectrace(:,1))
                    indexsel=[];
                    indexsel=find(newselectrace(:,1)==t);
                    if isempty(indexsel)==0
                        plot (newselectrace(indexsel(:),3), newselectrace(indexsel(:),4), 'r-');   % grafica traces seleccionadas
                        hold on;
                    end
                end
                 % dialog box to confirm
                 qstring=['Confirm deleting?'];
                button = questdlg(qstring); 
                if strcmp(button,'Yes')
                     newtrc=aux;
                     mask=mask+BW; %
                end
                clear newselectrace aux
            else
                % no select
                newtrc=aux; %all
             end % hay mol selec
        end %otra
      
      end %while
      
      %save mask     MODIFIER!!!! save coordinates to use inpolygon
      %afterwards!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      cd(path)
      cd('clean\');
      save([namefile,'-cleanmask.mat'],'mask','-mat');
      cd(path)
      
  else %previous mask
      cd(path)
      cd(['clean\']);
      load([namefile,'-cleanmask.mat']);
      
      cd(path)
      xi=[1,size(mask,2)]; %VOIR
      yi=[1,size(mask,1)];
      
      [newselectrace,aux]=pickpointsclear(mask,xi,yi,newtrc,stack_info.x,stack_info.y,[]); %mask instead of polygon coordinates
      
      figure
      stackmin=(min(min(min(datamatrix))));
      stackmax=(max(max(max(datamatrix))));
      imshow((datamatrix(:,:,1)),[stackmin stackmax],'InitialMagnification','fit');
      hold on
      if isempty(newtrc)==0
           for m=1:max(newtrc(:,1))
               indice=find(newtrc(:,1)==m);
               if isempty(indice)==0
                  for i = 1:size(indice,1)
                      graph(i,:)=newtrc(indice(i),:);   % archivo auxiliar con los puntos de cada trayectoria
                  end
                  plot ((graph (:,3)), (graph (:,4)), 'b-');   % grafica traces 
                  hold on
               end
               graph=[];
           end
         end %newtrcempty

      
      if isempty(newselectrace)==0
          for t=1:max(newselectrace(:,1))
              indexsel=[];
              indexsel=find(newselectrace(:,1)==t);
              if isempty(indexsel)==0
                  plot (newselectrace(indexsel(:),3), newselectrace(indexsel(:),4), 'r-');   % grafica traces seleccionadas
                  hold on;
              end
          end
          % dialog box to confirm
          qstring=['Confirm deleting?'];
          button = questdlg(qstring); 
          if strcmp(button,'Yes')
              newtrc=aux;
              %mask=mask+BW; %
          end
          clear newselectrace aux
      else
          % no select
          newtrc=aux; %all
      end % hay mol selec
  end

  if isempty(newtrc)==0
      
    hm=msgbox('Please wait','Re-numbering trajectories','help');
   % renumerotacion para tener nromol consecutivas
    auxtrc=[];
    contaux=1;
    nromol=[];
    conseq=1;
    totalmol=max(newtrc(:,1));
    for q=1:totalmol
       indexmol=find(newtrc(:,1)==q);
       if isempty(indexmol)==0
          nromol(contaux,1)=q;
          nromol(contaux,2)=conseq;
          for k=1:(size(indexmol,1))
             auxtrc(contaux,:)=newtrc(indexmol(k),:);
             auxtrc(contaux,1)=conseq;
             contaux=contaux+1;
          end
          conseq=conseq+1;
       end
    end
   
   
    domainfile=[namefile,identif,'.tif'];
   
    if length(dir(domainfile))>0 && localize==1 %only if wanted
       % reads files
       disp(['Domain file: ',domainfile]);
       Image=imread(domainfile);
       Image=double(Image);
       
       %numbering
       level1 = graythresh(Image);
       ImageBW = im2bw(Image,level1);
       [img_dist,pos_syn_proche]=bwdist(ImageBW);
       [labeled,nb_dom] = bwlabel(ImageBW,4);
       statimage=regionprops(labeled,'Centroid');

       img_zone=labeled;
       y=find(img_dist<=perisyn & img_dist~=0);
       img_zone(y)=-labeled(pos_syn_proche(y)); %perisyn: negative value
       labeled=img_zone; %with perisyn
       
       for i=1:size(auxtrc,1)
           posx=ceil(auxtrc(i,3));
           if posx==0; posx=1; end
           if posx>size(labeled,2); posx=size(labeled,2); end
           posy=ceil(auxtrc(i,4));
           if posy==0; posy=1; end
           if posy>size(labeled,1); posy=size(labeled,1); end
           
           if labeled(posy,posx)>0
               nrodom=labeled(posy,posx);
               auxtrc(i,6)=nrodom;
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
       
    else
       disp(['Domain file not found. Without localization']);
    end

    
   % guarda todo en carpeta clean
   cd(path)
   cd('clean\');
   cd('trc\');
   
   if length(dir(domainfile))>0 && localize==1 %only if wanted
       save([namefile,'.loc.trc'],'auxtrc','-ascii');
   else
       save([namefile,'.con.trc'],'auxtrc','-ascii');
   end

   disp('Files saved')
   cd(currentdir)

 end   
   disp('  ');
   hold off
   close(figclear)
   close(hm)
end

%close(hm)
cd(currentdir)  

% end of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

