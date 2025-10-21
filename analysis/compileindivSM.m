function compileindivSM
%function compileindivSM
% compiles diffusion analysis data (trajindivSM)
% Marianne Renner 07/16
% MR 07/17 includes Pc and trajectory area
% Marianne Renner - verified on SuperRes_v4, 09/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


currentdir=cd;

count=1;
data=[];
Dtotalextra=[];
Dtotalsyn=[];
msdtotalextra=[];
msdtotalsyn=[];
results=[];
logical=1;
meanmsdextra=[];
meanmsdsyn=[];
%distlength=[];
resPc=[];
areatrajsyn=[];
areatrajextra=[];
resPcsyn=[];
resPcextra=[];

 % dialog box to enter recognition criteria for images
prompt = {'Maximum length for trajectories (frames): '};
num_lines= 1;
dlg_title = 'Length thresholding';
def = {'1000'};
answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
if exit(1) == 0;
    return;
end

maxlength=str2num(answer{1});

dialog_title=['Select data folder'];
path = uigetdir(cd,dialog_title);
if path==0
    return
end

first=1;

while logical % allows entering data from different folders

   if length(dir([path,'\diff\']))==0
        msgbox('No folder \diff\ with analysis of diffusion','error','error')
        return
    end
    cd([path,'\diff\'])

    d=dir('*.tns*');
    st = {d.name};
    if isempty(st)==1
        msgbox('Wrong folder','error','error')
        return
    end
    
    [files,v] = listdlg('PromptString','Select files:','SelectionMode','multiple','ListString',st);
    if v==0
        return
    end
    
    [f,ultimo]=size(files);
    
    disp(' ')
    disp('Summary of diffusion data')
    
    for i=1:ultimo
        
        file=st{files(i)}
       
        load(file,'-mat'); %tns!!
        nrotraj=data(1).nrotraj;
        count=0;
        till=data(1).till;
%        perival=data(1).peri; %1:p=s 2:p=e
      %%  if perival==1 && i==1
        %    disp('Perisynaptic taken as synaptic')
        %else
        %    disp('Perisynaptic taken as extrasynaptic')
        %end
        
        maxtlagmsd=10; %!!!!!!!!!!!!!!!!
        
        
        cont=0;

        if first==1 % only first file!!  
            msdtotalextra=zeros(maxtlagmsd+1,1); %maxmsd
            msdtotalextra(1)=0;
            msdtotalextra(2)=0;
            for t=3:maxtlagmsd+2
                msdtotalextra(t)=till*(t-3);
            end
            msdtotalsyn=msdtotalextra;
           % msdtotalsyn(:,1)=[0; msdtotalsyn(1:size(msdtotalsyn,1)-1,1)];
           % msdtotalextra(:,1)=[0; msdtotalextra(1:size(msdtotalextra,1)-1,1)];
            first=0;
        end
        
        for j=1:nrotraj
            
            % disp(' ')
            %  disp(['Trajectory # ',num2str(j)]);
            if exist('waitbarhandle')
                waitbar(j/nrotraj,waitbarhandle,['Trajectory # ',num2str(j)]);
            end
            traj=data(j).traj;
            %peritype = data(1).peri;
            
            Dsyn=[];
            Dextra=[];
 
            msdextra=zeros(maxtlagmsd+2,1); %maxmsd
            for t=3:maxtlagmsd+2
                msdextra(t)=till*(t-2);
            end
            msdsyn=msdextra;

          if isempty(traj)==0
              
              if traj.length<maxlength
              
              %distlength=[distlength; traj.length];
              resPc=[resPc; j*ones(size(traj.fill,1),1) traj.fill];
              areatrajsyn=[areatrajsyn; traj.areatrajsyn];
              areatrajextra=[areatrajextra; traj.areatrajextra];
                
              if isfield(traj,'segm')
                
                for k=1:traj.nrosegm
                    
                    lengthsegm=size(traj.segm(k).data,1);
                    loc=traj.segm(k).data(1,6);
                    
                    if loc==0
                        Dextra=[Dextra; i j k traj.segm(k).D loc lengthsegm];
                        %disp(traj.segm(k))
                        %disp(traj.segm(k).fill)
                        if isempty(traj.segm(k).fill)==0
                            resPcextra=[resPcextra; i j k traj.segm(k).fill(5) loc lengthsegm];
                        else
                            resPcextra=[resPcextra; i j k NaN loc lengthsegm];
                        end
                        
                        new=zeros(maxtlagmsd+2,1);
                        new(1)= loc; %loc
                        new(2)= lengthsegm; %length
                        
                        lim=min(size(traj.segm(k).msd,1),maxtlagmsd);
                        for t=2:lim+1
                            new(t+1)=traj.segm(k).msd(t-1,2);
                        end
                        if size(new,1)==1
                            new=new';
                        end
                        
                        if size(new,1)>maxtlagmsd
                            %lim=maxtlagmsd+1;
                           % if new(lim)==0 | new(lim)==msdextra(lim,size(msdextra,2))
                           % else
                                msdextra=[msdextra new(1:maxtlagmsd+2)];
                           % end
                        else
                             msdextra=[msdextra new];
                        end
                        
                    else %syn?
                        
                       % if perival==1 %p=s
                            if loc==0
                                control=0; %extra
                            else
                                control=1; %syn
                            end
                        %else %p=e
                        %    if loc>0
                        %        control=1; %syn
                        %    else
                        %        control=0; %extra
                        %    end
                        %end
                        
                        if control==0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% segment extrasyn
                            Dextra=[Dextra; i j k traj.segm(k).D loc lengthsegm];
                             if isempty(traj.segm(k).fill)==0
                                 resPcextra=[resPcextra; i j k traj.segm(k).fill(5) loc lengthsegm];
                             else
                                 resPcextra=[resPcextra; i j k NaN loc lengthsegm];
                             end

                            new=zeros(maxtlagmsd+2,1);
                            new(1)= loc; %loc
                            new(2)= lengthsegm; %length

                            lim=min(size(traj.segm(k).msd,1),maxtlagmsd); % minimum # of time lags
                           
                            for t=2:lim+1
                                new(t+1)=traj.segm(k).msd(t-1,2);
                            end
                            
                            if size(new,1)>maxtlagmsd % enough time lags
                               % lim=maxtlagmsd+1;
                              %  if new(lim)==0 | new(lim)==msdextra(lim,size(msdextra,2))
                              %  else
                                    msdextra=[msdextra new(1:maxtlagmsd+2)];
                               % end
                            else
                                msdextra=[msdextra new];
                            end
                        else %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% segment synaptic

                            new=zeros(size(msdsyn,1)+2,1);                                 %!!!!!!!!!!!
                            Dsyn=[Dsyn; i j k traj.segm(k).D loc lengthsegm];
                             if isempty(traj.segm(k).fill)==0
                                 resPcsyn=[resPcsyn; i j k traj.segm(k).fill(5) loc lengthsegm];
                             else
                                 resPcsyn=[resPcsyn; i j k NaN loc lengthsegm];
                             end

                            new(1)= loc; %loc       
                             new(2)= lengthsegm; %length
                           
                            lim=min(size(traj.segm(k).msd,1),maxtlagmsd);
                           % disp(traj.segm(k).msd)
                            
                            for t=2:lim+1
                               new(t+1)=traj.segm(k).msd(t-1,2);
                            end
                            
                           % disp(new)
                            
                            if size(new,1)>maxtlagmsd
                               % lim=maxtlagmsd+1
                               % if new(lim)==0 | new(lim)==msdsyn(lim,size(msdsyn,2))
                               % else
                                    msdsyn=[msdsyn new(1:maxtlagmsd+2)];
                               % end
                            else
                                msdsyn=[msdsyn new];
                            end
                           % disp(msdsyn)
                        end % control

                    end %loc
                end %nrosegm
            end %segm           
    
            warning off all
            if size(msdsyn,2)>1
              %  disp(msdsyn)
                indexzero=find(msdsyn(4:size(msdsyn,1),2)==0);
                if isempty(indexzero)==0
                    msdsyn(indexzero+3,2)=NaN;
                end
                %disp(msdsyn)
                msdsyn(:,1)=[0; msdsyn(1:size(msdsyn,1)-1,1)];
                msdtotalsyn=[msdtotalsyn msdsyn(:,2:size(msdsyn,2))];
            end
            if size(msdextra,2)>1
                indexzero=find(msdextra(4:size(msdextra,1),2)==0);
                if isempty(indexzero)==0
                    msdextra(indexzero+3,2)=NaN;
                end
                
               % disp(size(msdextra))
               % disp(size(msdtotalextra))
 
                msdextra(:,1)=[0; msdextra(1:size(msdextra,1)-1,1)];
                msdtotalextra=[msdtotalextra msdextra(:,2:size(msdextra,2))];
            end


            Dtotalextra=[Dtotalextra; Dextra];
            Dtotalsyn=[Dtotalsyn; Dsyn];
          
            Dextra=[];
            Dsyn=[];
            
            end % maxlength
            
          end %emptytraj
          
        end %nrotraj
    end % loop files

    cd(currentdir)

    % dialog box to enter new data from another folder
    qstring=['more data folders?'];
    button = questdlg(qstring); 
    if strcmp(button,'Yes')
        logical=1;
        dialog_title=['Select data folder (with \diff\stab folder)'];
        path = uigetdir(cd,dialog_title);
        if path==0
            return
        end
    else
        break    
        logical=0
    end

end % while

% disp(msdtotalsyn)

for i=3:size(msdtotalextra,1)
    meanmsdextra(i-2,1)=msdtotalextra(i,1);
    indexnan=[];
    data=msdtotalextra(i,2:size(msdtotalextra,2));
    indexnan=isnan(data);
    newdata=data(find(indexnan==0));
    meanmsdextra(i-2,2)=mean(data(find(indexnan==0))); %MSD
    meanmsdextra(i-2,3)=std(data(find(indexnan==0))); %SD
    meanmsdextra(i-2,4)=meanmsdextra(i-2,3)/sqrt(size(msdtotalextra,2)-1); %sem
    meanmsdextra(i-2,5)=median(data(find(indexnan==0))); %median
    if i>3
        meanmsdextra(i-2,6)=size(newdata,2); %number of steps
    else
        meanmsdextra(i-2,6)=0;%first value
    end
end

if size(msdtotalsyn,2>1)
    for i=3:size(msdtotalsyn,1)
     meanmsdsyn(i-2,1)=msdtotalsyn(i,1);
     indexnan=[];
     data=msdtotalsyn(i,2:size(msdtotalsyn,2));
     indexnan=isnan(data);
     %disp(i)
     newdata=data(find(indexnan==0));
     meanmsdsyn(i-2,2)=mean(data(find(indexnan==0))); %MSD
     meanmsdsyn(i-2,3)=std(data(find(indexnan==0))); %SD
     meanmsdsyn(i-2,4)=meanmsdsyn(i-2,3)/sqrt(size(msdtotalsyn,2)-1); %sem
     meanmsdsyn(i-2,5)=median(data(find(indexnan==0))); %median
     if i>3
         meanmsdsyn(i-2,6)=size(newdata,2);
     else
         meanmsdsyn(i-2,6)=0; %first value
     end
    end
end

Pcextra=[];
Pcsyn=[];
if size(resPc,2)>5 %localization
    indexextra=find(resPc(:,7)==0);
    if isempty(indexextra)==0
        Pcextra=resPc(indexextra,:);
    end
    indexsyn=find(resPc(:,7)>0);
    if isempty(indexsyn)==0
        Pcsyn=resPc(indexsyn,:);
    end
    resPc=Pcextra;
else
    Pcextra=resPc;
end


start_path=currentdir;
dialog_title=['Save data in'];
sn = uigetdir(start_path,dialog_title);
if sn==0
        return
end
cd(sn)
def_name=['savename.txt'];
[savename,path] = uiputfile(def_name,'Savename:');
if isequal(savename,0) | isequal(path,0)
else
    [namesavefile,rem]=strtok(savename,'.'); 
    save([namesavefile,'-Dtotalout.txt'],'Dtotalextra','-ascii');
    if isempty(Dtotalsyn)==0
        save([namesavefile,'-Dtotalin.txt'],'Dtotalsyn','-ascii');
    end
    save([namesavefile,'-Pcsegmout.txt'],'resPcextra','-ascii');
    if isempty(Dtotalsyn)==0
        save([namesavefile,'-Pcsegmin.txt'],'resPcsyn','-ascii');
    end
    save([namesavefile,'-MSDtotalout.txt'],'msdtotalextra','-ascii');
    if isempty(msdtotalsyn)==0
        save([namesavefile,'-MSDtotalin.txt'],'msdtotalsyn','-ascii');
    end
    save([namesavefile,'-averagedMSDout.txt'],'meanmsdextra','-ascii');
    if isempty(meanmsdsyn)==0
        save([namesavefile,'-averagedMSDin.txt'],'meanmsdsyn','-ascii');
    end
    if isempty(Pcsyn)==0
        save([namesavefile,'-Pcin.txt'],'Pcsyn','-ascii');
    end
    save([namesavefile,'-Pcout.txt'],'resPc','-ascii');
    save([namesavefile,'-areatrajout.txt'],'areatrajextra','-ascii');
    if isempty(areatrajsyn)==0
        save([namesavefile,'-areatrajin.txt'],'areatrajsyn','-ascii');
    end
   % save([namesavefile,'-lengths.txt'],'distlength','-ascii');

end
        
cd(currentdir)
    
clear all

%eof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    