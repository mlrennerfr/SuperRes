function savePALMdata(namefile, matrice_results, px_mu, auxx,auxy, handles)
% function savePALMdata(namefile, matrice_results, px_mu, auxx,auxy, handles)
%
% Marianne Renner, SuperRes programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
currentdir=cd;
%auxy= matrice_results(2,:);
%auxx= matrice_results(3,:);  
matrice_results(3,:)=auxy/px_mu; %reconversion
matrice_results(2,:)=auxx/px_mu;

%xdim=ceil(max(matrice_results(2,:)));
%ydim=ceil(max(matrice_results(3,:)));

%correction coordinates!!
% att: µm!
minx1=min(handles.x)/px_mu;
miny1=min(handles.y)/px_mu;
maxx1=max(handles.x)/px_mu;
maxy1=max(handles.y)/px_mu;
%minx1=min(auxx)/px_mu;
%miny1=min(auxy)/px_mu;
%maxx1=max(auxx)/px_mu;
%maxy1=max(auxy)/px_mu;

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
aux= [matrice_results(1,last) miny1 minx1  matrice_results(4,last)  matrice_results(5,last)...
        matrice_results(6,last) matrice_results(7,last) matrice_results(8,last) matrice_results(9,last) matrice_results(10,last) matrice_results(11,last)];
matrice_results=[matrice_results, aux'];
aux=[matrice_results(1,last) maxy1 maxx1  matrice_results(4,last)  matrice_results(5,last)...
                matrice_results(6,last) matrice_results(7,last) matrice_results(8,last) matrice_results(9,last) matrice_results(10,last) matrice_results(11,last)];
matrice_results=[matrice_results, aux'];
    
clear pk
pk(:,1)=matrice_results(1,:); 
pk(:,2)=matrice_results(2,:);  
pk(:,3)=matrice_results(3,:);    
pk(:,5)= matrice_results(4,:); 
pk(:,4)=matrice_results(5,:)*2; 
pk(:,7)=matrice_results(6,:);  
pk(:,6)=matrice_results(7,:);

if size(matrice_results,1)>7
    pk(:,8)=matrice_results(8,:); %ratio
    pk(:,9)=matrice_results(9,:); %z
    pk(:,10)=matrice_results(10,:); %test
    pk(:,11)=matrice_results(11,:); %test
    
    if max(pk(:,9))==0 && min(pk(:,9))==0 % no values for z
        if isdir('pk'); else mkdir('pk');end
        pkpath='pk';
        pksavename=[namefile, '.pk'];
    else
        if isdir('pk3'); else mkdir('pk3');end
        pkpath='pk3';
        pksavename=[namefile, '.pk3'];
        % prepare file for visp
               % indexoui=find(pk(:,8)<1000);
              %  disp('Pixel size for Visp : 160 nm')
              %  szpxum=160/1000;
               % if isempty(indexoui)==0
                %    pkvisp=[pk(indexoui,2)*szpxum pk(indexoui,3)*szpxum pk(indexoui,9)*szpxum pk(indexoui,5) pk(indexoui,1)];    %[x y z alpha fr];
                %    save([namefile, '.3d'], 'pkvisp','-ascii');
                %    clear pkvisp indexoui
              %  end
    end
else
    if isdir('pk'); else mkdir('pk');end
    pkpath='pk';
    pksavename=[namefile, '.pk'];
end

savename=[namefile,'.mat'];
save(savename, 'matrice_results'); % second image, corrected for color drift

cd(pkpath);
save(pksavename, 'pk','-ascii');
clear x y alpha fr radius matrice_results pk

cd(currentdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
