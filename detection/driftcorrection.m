function driftcorrection
% function driftcorrection
% correct drift images
% Marianne Renner SuperRes programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% dialog box to enter recognition criteria for images
prompt = {'Pixel size (nm):','PALM pixel size (nm):'};
num_lines= 1;
dlg_title = 'Correction of multiple detections';
def = {'107','10'};
answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
if exit(1) == 0;
    return;
end
px_mu=str2num(answer{1})/1000;
dx=str2num(answer{2})/1000;

%files
d=dir('*.mat*'); % .stk files
st={d.name};
if isempty(st)==1
    msgbox(['No files!!'],'','error');
    return
end

%choose data
[listafiles,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',st);
if v==0
     return
end



for nromovie=1:size(listafiles,2)
    matfile=st{listafiles(nromovie)} ;
   disp=['Correction of stage drift'];
   disp(' ');
    [filename,rem]=strtok(matfile,'.');
    
    %handles=loadPALMdata(file,px_mu,0,handles);
    [x,y, fr, alpha, radius,ffname,filename]=PALMdata(matfile,px_mu);
    
    if dx>0
        x=x/dx;
        y=y/dx;
    end

    
     matrice_results =correctshiftstage3(x,y, fr, alpha, radius);  

     xdim=ceil(max(matrice_results(2,:))); %x?????
     ydim=ceil(max(matrice_results(3,:)));
     
     figure
     imshow(zeros(xdim,ydim),'InitialMagnification','fit');
     hold on
     plot(matrice_results(2,:),matrice_results(3,:),'.','MarkerSize',2,'Color','y');

    
    % att dx!!!!!!
    

    
    currentdir=cd;
    k = strfind(currentdir, 'driftcorrected');
    if isempty(k)==1
        if isdir('driftcorrected'); else mkdir('driftcorrected'); end
        cd('driftcorrected')  
    end

    % save new data
    if isempty(matrice_results)==0
        auxy= matrice_results(2,:);
        auxx= matrice_results(3,:); 
        
        matrice_results(3,:)=auxy/px_mu; %reconversion
        matrice_results(2,:)=auxx/px_mu;
        
        xdim=ceil(max(matrice_results(2,:)));
        ydim=ceil(max(matrice_results(3,:)));   
        
        %correction coordinates!!
        % att: µm!
        minx1=min(handles.x)/px_mu;
        miny1=min(handles.y)/px_mu;
        maxx1=max(handles.x)/px_mu;
        maxy1=max(handles.y)/px_mu;   
        
        indexneg=find(matrice_results(2,:)<minx1);
        if isempty(indexneg)==0
            matrice_results(:,indexneg)=NaN;
        end
        indexneg=find(matrice_results(3,:)<miny1);
        if isempty(indexneg)==0
            matrice_results(:,indexneg)=NaN;
        end
        indexextra=find(matrice_results(3,:)>maxx1);
        if isempty(indexextra)==0
            matrice_results(:,indexextra)=NaN;
        end
        indexextra=find(matrice_results(2,:)>maxy1);
        if isempty(indexextra)==0
            matrice_results(:,indexextra)=NaN;
        end
        
        % false points for size
        last=size(matrice_results,2);
        aux= [matrice_results(1,last) miny1 minx1  matrice_results(4,last)  matrice_results(5,last)];
        matrice_results=[matrice_results, aux'];
        aux=[matrice_results(1,last) maxy1 maxx1  matrice_results(4,last)  matrice_results(5,last)];
        matrice_results=[matrice_results, aux'];
        
        savename=[namefile, '.mat'];
        save(savename, 'matrice_results'); % one image, corrected for stage drift
    end % empty matrice
    
end %loop files
        

cd(currentdir)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
