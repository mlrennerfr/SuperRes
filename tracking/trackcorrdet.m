function newpk = trackcorrdet(cleanpeaks, Dini, Te,a,blinkdist,blinktime)
% function newpk = trackcorrdet(cleanpeaks, Dini, Te,a,blinkdist,blinktime)
%
% builds trajectories using data of .pk files from DetectSM to correct 
% for multiple detections
% all detections belonging to a "trajectory" (belonging to the same
% molecule) are averaged to keep only one value
% output: pk file corrected (newpk)
%
% Modified from initialtracking 
% Marianne Renner june 2016
% Marianne Renner 01/2025 - adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off all

%cleanpeaks=load(file);

newpk=[];
w=1.3;
persistance=5;

MSDparam.nb_fit = 5;
MSDparam.nb_points_fraction_traj = 25;
MSDparam.nbMSD = MSDparam.nb_points_fraction_traj;
MSDparam.Dini = Dini;
MSDparam.alpha = 25;
MSDparam.a = a;
MSDparam.Te = Te;
Nb_spots=0;

Nz=max(cleanpeaks(:,1));

waitbarhandle = waitbar( 0,'Please wait...','Name',['Finding positions of the same molecule']) ;
    
   
for m=1:Nz
    
    index=find(cleanpeaks(:,1)==m);
    if isempty(index)==0
        newpk=cleanpeaks(index,:);
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
    %cd(currentdir);
    
    %traj=initialtracker(handles.file,plan, Nz,detoptions, handles);
    
for z=1:Nz
    if exist('waitbarhandle')
       waitbar(z/Nz,waitbarhandle,['Frame # ',num2str(z)]);
    end
    if (Nb_spots == 0)% Recherche du premier plan dans lequel on détecte un objet.
                if plan(z).Nb_objets
                    for q=1:plan(z).Nb_objets
                        Nb_spots = Nb_spots + 1;
                        plan(z).objet(q).antecedent = [0,0];
                        plan(z).objet(q).spot = Nb_spots;
                        plan(z).objet(q).D = Dini;
                        s_s(Nb_spots,:)=[z q plan(z).objet(q).D];
                        buffer(Nb_spots).coordinates = [plan(z).objet(q).centre(1), plan(z).objet(q).centre(2), z];
                        buffer(Nb_spots).D = Dini;
                    end;
                end;
     else
                % Construction des trajectoires.
                if (plan(z).Nb_objets> 0)
                    for q=1:plan(z).Nb_objets
                        plan(z).objet(q).antecedent = [0,0];
                    end;
                    s_s_p = s_s(z-s_s(:,1) <= persistance,:);
                    [Nb_s_s_p,v]=size(s_s_p);
                    clear v q;
                    dist_centres = zeros(plan(z).Nb_objets,Nb_s_s_p);
                    r_zone = zeros(1,Nb_s_s_p);
                    for j=1:Nb_s_s_p
                        
                        %s_s_p(j,3)= dmax
                       % if s_s_p(j,3)>Dini*10
                        %    s_s_p(j,3)=Dini*10;
                       % end
                        %disp(s_s_p(j,3))
                        r_zone(j)=f_dist_max(w,s_s_p(j,3),a,Te,abs(z-s_s_p(j,1)));
                        if r_zone(j)>Dini*150
                            r_zone(j)=Dini*150;
                        end
                        for i=1:plan(z).Nb_objets
                            dist_centres(i,j)=norm(plan(z).objet(i).centre-plan(s_s_p(j,1)).objet(s_s_p(j,2)).centre,2);
                            % Pénalisation si incompatibilité avec le coefficient de diffusion.
                            if dist_centres(i,j) > r_zone(j)
                                dist_centres(i,j)= 512*2; %%%%%%%%%%%%
                            end;
                        end;
                    end;
                    antecedent=zeros(1,plan(z).Nb_objets);
                    successeur=zeros(1,Nb_s_s_p);
                    cond_d = 1;
                    while ((sum(successeur)<Nb_s_s_p) & (sum(antecedent)<plan(z).Nb_objets) & cond_d)
                        sans_antecedent = find(antecedent~=1);
                        sans_successeur = find(successeur~=1);
                        [tri_min,tri_ind] = sort(min(dist_centres(sans_antecedent,sans_successeur)));
                        cond_d= tri_min(1) <= r_zone(tri_ind(1));
                        if cond_d
                            [i_p_c,i_p_s_s]=find(dist_centres(sans_antecedent,sans_successeur)==tri_min(1));
                            % Mise à jour de 'antecedent' et 'successeur'.
                            antecedent(sans_antecedent(i_p_c))=1;
                            successeur(sans_successeur(i_p_s_s))=1;
                            % Mise à jour de 'buffer.coordinates'.
                            traj = plan(s_s_p(sans_successeur(i_p_s_s),1)).objet(s_s_p(sans_successeur(i_p_s_s),2)).spot;
                            local_coordinates  = [plan(s_s_p(sans_successeur(i_p_s_s),1)).objet(s_s_p(sans_successeur(i_p_s_s),2)).centre(1),...
                                     plan(s_s_p(sans_successeur(i_p_s_s),1)).objet(s_s_p(sans_successeur(i_p_s_s),2)).centre(2),...
                                     z];
                            buffer(traj).coordinates = [ buffer(traj).coordinates ; local_coordinates];
                            clear local_coordinates;
                            % Mise à jour de 'plan'.
                            plan(z).objet(sans_antecedent(i_p_c)).antecedent = [s_s_p(sans_successeur(i_p_s_s),1:2)];
                            plan(z).objet(sans_antecedent(i_p_c)).spot = traj;
                            plan(z).objet(sans_antecedent(i_p_c)).D = comp_MSD(buffer(traj).coordinates,MSDparam);
                            % Mise à jour de 'buffer.D'.
                            buffer(traj).D = [ buffer(traj).D ; plan(z).objet(sans_antecedent(i_p_c)).D];
                            % Mise à jour de 's_s'.
                            s_s(traj,:) = [z sans_antecedent(i_p_c) plan(z).objet(sans_antecedent(i_p_c)).D];
                            clear indice traj;
                        end
                    end
                    sans_antecedent = find(antecedent~=1);
                    sans_successeur = find(successeur~=1);
                    
                    
                    for i=1:length(sans_antecedent)
                        Nb_spots=Nb_spots+1;
                        % Mise à jour de 'plan'.
                        plan(z).objet(sans_antecedent(i)).antecedent = [0,0];
                        plan(z).objet(sans_antecedent(i)).spot = Nb_spots;
                        plan(z).objet(sans_antecedent(i)).D = Dini;
                        % Mise à jour de 's_s'.
                        s_s(Nb_spots,:)=[z sans_antecedent(i) plan(z).objet(sans_antecedent(i)).D];
                        % Mise à jour de 'buffer'.
                        buffer(Nb_spots).coordinates = [plan(z).objet(sans_antecedent(i)).centre(1), plan(z).objet(sans_antecedent(i)).centre(2), z];
                        buffer(Nb_spots).D = Dini;
                    end     
                    clear s_s_p Nb_s_s_p sans_antecedent sans_successeur antecedent cond_d;
                    clear successeur tri_ind tri_min i_p_c i_p_s_s dist_centres r_zone;
                end
      end;
end;

% Définition de la structure 'spot'.
if exist('waitbarhandle')
       waitbar(Nz/Nz,waitbarhandle,['Saving results...']);
end
ind_fit=1;
fit.nb_spots = 0;
for ind_spot = 1:Nb_spots
    tampon(ind_spot).nb_points = 0;
    tampon(ind_spot).nb_segments = 0;
    for z=1:Nz
        for ind_obj = 1:plan(z).Nb_objets
            if (plan(z).objet(ind_obj).spot == ind_spot)
                tampon(ind_spot).nb_points = tampon(ind_spot).nb_points + 1;
                if tampon(ind_spot).nb_points == 1
                    tampon(ind_spot).nb_segments = 1;
                    tampon(ind_spot).segment(tampon(ind_spot).nb_segments).length = 1;
                    tampon(ind_spot).segment(tampon(ind_spot).nb_segments).coordinates = [plan(z).objet(ind_obj).centre(1), plan(z).objet(ind_obj).centre(2), z];
                else
                    if (z-1) == tampon(ind_spot).segment(tampon(ind_spot).nb_segments).coordinates(tampon(ind_spot).segment(tampon(ind_spot).nb_segments).length,3)
                        tampon(ind_spot).segment(tampon(ind_spot).nb_segments).length = 1 + tampon(ind_spot).segment(tampon(ind_spot).nb_segments).length;
                        tampon(ind_spot).segment(tampon(ind_spot).nb_segments ).coordinates(tampon(ind_spot).segment(tampon(ind_spot).nb_segments).length,:) = [plan(z).objet(ind_obj).centre(1), plan(z).objet(ind_obj).centre(2), z];
                    else
                        tampon(ind_spot).nb_segments = tampon(ind_spot).nb_segments + 1;
                        tampon(ind_spot).segment(tampon(ind_spot).nb_segments).length = 1;
                        tampon(ind_spot).segment(tampon(ind_spot).nb_segments ).coordinates(tampon(ind_spot).segment(tampon(ind_spot).nb_segments).length,:) = [plan(z).objet(ind_obj).centre(1), plan(z).objet(ind_obj).centre(2), z];
                    end;
                end;
            end;
        end;
    end;
    % On éjecte les trajectoires ne contenant qu'un seul point.
  % if tampon(ind_spot).nb_points>1
        fit.nb_spots = fit.nb_spots + 1;
        fit.spot(fit.nb_spots)=tampon(ind_spot);
   % end;

end;
clear Nb_spots tampon ind_obj plan;



% create trc file
trc=[];
count=1;
for ind_spot=1:fit.nb_spots
    nb_segments=fit.spot(ind_spot).nb_segments;
    for ind_seg=1:nb_segments
        coord=fit.spot(ind_spot).segment(ind_seg).coordinates;
        for j=1:size(coord,1)
           if coord(j,1)>0
              trc(count,1)=ind_spot;
              trc(count,2)=coord(j,3); %frame
              trc(count,3)=coord(j,1); %x
              trc(count,4)=coord(j,2); %y
              trc(count,5)=0;
              count=count+1;
           end
        end
    end
end
    
%close(waitbarhandle);
newpk=[];
    
%reconnection to correct blinking
reconnectedfile=reconnectfastcorrdet(trc,blinktime,blinkdist,1,blinkdist,waitbarhandle); %(trc,blinktime,blinkdist,mintrace,limdist,waitbarhandle);
    
if isempty(reconnectedfile)==0 
    
       % disp(reconnectedfile)
        
        pkaux = [cleanpeaks zeros(size(cleanpeaks,1),1)];
        
        for i=1:size(reconnectedfile,1)
            
            index=find(pkaux(:,2)==reconnectedfile(i,3));
            if isempty(index)==0
                index2=find(pkaux(index,3)==reconnectedfile(i,4));
                 if isempty(index2)==0
                     pkaux(index(index2),size(pkaux,2))=reconnectedfile(i,1);
                 end
            end
        end
        
      %  disp(pkaux)
      %  disp(size(pkaux))


        for i=1:max(pkaux(:,size(pkaux,2)))+1
            indexi=find(pkaux(:,size(pkaux,2))==i-1);
            if isempty(indexi)==0
                if i==0 %not tracked
                    newpk=[newpk; pkaux(indexi,1:size(pkaux,2)-1)];
                else
                    %index=find(pkaux(:,size(pkaux,2))== pkaux(i,size(pkaux,2)));
                    if size(pkaux,2)>7
                        newpk=[newpk; pkaux(indexi(1),1) mean(pkaux(indexi,2)) mean(pkaux(indexi,3))  mean(pkaux(indexi,4)) mean(pkaux(indexi,5)) mean(pkaux(indexi,6)) mean(pkaux(indexi,7))...
                            mean(pkaux(indexi,8)) mean(pkaux(indexi,9)) mean(pkaux(indexi,10)) mean(pkaux(indexi,11))];
                    else
                        newpk=[newpk; pkaux(indexi(1),1) mean(pkaux(indexi,2)) mean(pkaux(indexi,3))  mean(pkaux(indexi,4)) mean(pkaux(indexi,5)) mean(pkaux(indexi,6)) mean(pkaux(indexi,7))];
                    end
                end
            end
        end
       
end %empty reconnected

close(waitbarhandle);

%disp(size(newpk))

plotcontrol=0;
if plotcontrol
          figure
       Im=ones(ceil(max(pkaux(:,3))),ceil(max(pkaux(:,2)))); %!!!!!!!!!!!!!!!!!!
       imshow(Im,'InitialMagnification','fit');
       hold on
       plot(pkaux(:,2),pkaux(:,3),'o','MarkerSize',5,'Color','b');
       hold on
       plot(newpk(:,2),newpk(:,3),'x','MarkerSize',5,'Color','r');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
