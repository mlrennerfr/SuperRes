function convertsizefluo


prompt = {'Microscope pixel size (µm):','PALM pixel width (µm):','Sigma gaussian:','Identifier for the rendered image:','Identifier for the fluo image:'};
num_lines= 1;
dlg_title = 'Parameters for rendered image';
def = {'0.16','0.01','0.01','_kcc2','_geph'}; % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
if exit(1) == 0;
       return; 
end
szpx=str2num(answer{1})  ;   
dx=str2num(answer{2})  ;   
sigma=str2num(answer{3})  ;     
ident1=answer{4}  ;     
ident2=answer{5}  ;     

% carga files y loop general 
dialog_title=['Select data folder'];
directory_name = uigetdir(cd,dialog_title);
if directory_name==0
    return
end
cd(directory_name);

d = dir(['*',ident1,'.mat']); % movie
st={d.name};
if isempty(st)==1
   msgbox(['No files!!'],'Select files','error');
   return
end
[listafiles,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',st);
if v==0
   return
end

for nromovie=1:size(listafiles,2)
    file=[st{listafiles(nromovie)}]
     [namefile,rem]=strtok(file,'.');

    data=loadPALMdatanew(file,0); 
    
    %name second file
    [namefile,rem]=strtok(file,'.');
    k = findstr(namefile, ident1);
    rootname=namefile(1:k-1);
    file2=[rootname,ident2,'.tif']
    
    control=0;
    if isempty(dir(file2))==0
        control=1;  
    end
    
    if control==1
        [~,datamatrix] = tifdataread(file2);
        datamatrix=double(datamatrix.data);
        
        % rendered
        I=makerendered2(data,szpx,sigma,dx) ;   
        
        disp(size(I))
        
        disp(size(datamatrix))
        
        data2 = imresize(datamatrix,[size(I,1),size(I,2)]);
        
        
        
        
    end
    
    

end % option color


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dataread=loadPALMdatanew(file,dx)

points=[];
S = load(file);
dataread.x=[];
dataread.y=[];
           
if isfield(S,'matrice_results')
    aux = S.matrice_results;
    clear S;
    fr = aux(1,:); fr = fr(:);
    % x = aux(3,:) * px_mu;
    x = aux(3,:);
    x = x(:);
    % y = aux(2,:) * px_mu; 
    y = aux(2,:); 
    y = y(:);
    alpha = aux(4,:); alpha = alpha(:);
    radius=(aux(5,:));
    sigma=(aux(6,:));
    blink=(aux(7,:));
    if size(aux,1)>7
        ratio=aux(8,:);
        z=aux(9,:);
        test1=aux(10,:);
        test2=aux(11,:);
    else
        ratio=[];
        z=[];
        test1=[];
        test2=[];
    end
    
else
    disp('Invalid file!');
    return
end


if dx>0
    x=x/dx;
    y=y/dx;
end

dataread.x = x;
dataread.y = y;
dataread.alpha = alpha;
dataread.fr = fr;
dataread.radius=radius;
dataread.sigma=sigma;
dataread.blink=blink;
dataread.ratio=ratio;
dataread.z=z;
dataread.test1=test1;
dataread.test2=test2;

clear aux


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 