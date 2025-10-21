function current=rajoutzone(current,nb_frames,img_zone,perisyn,waitbarhandle)
%function current=rajoutzone(current,nb_frames,img_zone,perisyn,waitbarhandle)
%cette fonction rajoute une matrice (nb de frames,4) current.spot(ind_spot).localisation.coord
%qui contient les x,y, z et le type de zone dans laquelle se trouve chaque spot 
%d'un film donné (x et y sont mis à zeros lors d'un blink):
%0 pour une zone extrasyn
%un nombre pair pour chaque zone synaptique
%(eventuellement un nombre impair pour chaque zone peri correspondant à une
%synapse : cf numerote_syn_peri)
%-1 pour un moment de blink
%elle rajoute aussi le paramètre perisyn(2 pixels à priori) dans current.spot(ind_spot).localisation(ind_rajout).perisyn
% MVE 07
% Mod by MR 09 for SPTrack v4
% Marianne Renner 01/2025 - adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


s=size(img_zone);
if current.nb_spots<length(current.spot)
    current.nb_spots=length(current.spot);
end
%nb_spots= current.nb_spots;

for ind_spot=1:current.nb_spots
    if exist('waitbarhandle')
       waitbar(ind_spot/current.nb_spots,waitbarhandle,['Spot # ',num2str(ind_spot)]);
    end

    nb_points=0;
    temp=zeros(nb_frames,4);
    temp(:,4)=-1.*ones(nb_frames,1);
    temp(:,3)=(1:1+nb_frames-1)';
    nb_segments=current.spot(ind_spot).nb_segments;
    for ind_seg=1:nb_segments
        %L=current.spot(ind_spot).segment(ind_seg).length;
        L=length(current.spot(ind_spot).segment(ind_seg).coordinates(:,1));
        nb_points=nb_points+L;
        for p=1:L
            i=round(current.spot(ind_spot).segment(ind_seg).coordinates(p,2));
            j=round(current.spot(ind_spot).segment(ind_seg).coordinates(p,1));
            k=current.spot(ind_spot).segment(ind_seg).coordinates(p,3)-1+1;
            if i>s(1)
                i=s(1);
            elseif i<1
                i=1;
            end
            if j>s(2)
                j=s(2);
            elseif j<1
                j=1;
            end
            temp(k,:)=[current.spot(ind_spot).segment(ind_seg).coordinates(p,1:3) img_zone(i,j)];
            clear i j k
        end
        clear L p
    end
    clear nb_segments
    current.spot(ind_spot).localisation(1).coord=temp;
    current.spot(ind_spot).localisation(1).perisyn=perisyn;
    current.spot(ind_spot).nb_points=nb_points;
    clear temp nb_points
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
