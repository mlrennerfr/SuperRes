function [img_zone,nb_syn]=numerotesyn(imgFM_seuil,perisyn)
%function [img_zone,nb_syn]=numerotesyn(imgFM_seuil,perisyn)
%� partir de l'image FM binaire 'imgFM_seuil', on
%attribut un num�ro pair � chaque zone synapse et un numero impair � chaque
%zone p�risynaptique correspondante (ex: syn n�150, zone peri n�149)
%on consid�re qu'une zone perisynaptique est situ�e � moins d'une distance
% 'perisyn' d'une synapse (exprim� en pixel)
% cette fonction retourne alors l'image 'img_zone' correspondant
% MVE 07
% Marianne Renner 01/2025 - adapted to SuperRes_v4
%==========================================================================

%obtention d'une carte des distances et d'une carte donnant la position de
%la synapse la plus proche
[img_dist,pos_syn_proche]=bwdist(imgFM_seuil);

%attribution d'un numero pair � chaque synapse
[img_syn,nb_syn] = bwlabel(imgFM_seuil,4);
img_syn=2.*img_syn;

%meme numero pair que la synapse la plus proche � chaque zone p�risynaptique d�fini par une distance <ou = au 
%param�tre "perisyn" par rapport � la synapse la plus proche
img_zone=img_syn;
y=find(img_dist<=perisyn & img_dist~=0);
img_zone(y)=img_syn(pos_syn_proche(y))-1;

clear  y img_dist pos_syn_proche img_syn 