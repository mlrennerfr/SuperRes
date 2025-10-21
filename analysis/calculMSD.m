function [D,b,MSD]=calculMSD(trc,szpx,till,nb_fit,maxtime)
% function [D,b,MSD]=calculMSD(trc,szpx,till,nb_fit,maxtime)
% calculates MSD over the trajectory trc
% size pixel szpx in nm
% time between images till in ms
% number of points to fit for D: nb_fit
% MSD.rho=msddata;
% MSD.N_denominateur=nro points averaged 
% MSD.time_lag=time interval
% MSD.sigma=error
%
% MR 07/08 - rev 12/08. From calcul_global_MSD_coeff_gui.m
% mod for SPTrack.m
% Marianne Renner 09/2025 : verified for SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nb_points = length(trc(:,1));

if nargin<5
    nbMSD=nb_points;
else
    if nb_points >maxtime
       nbMSD=maxtime;
    else
        nbMSD=nb_points;
    end
end

warning off MATLAB:divideByZero;
t_lag = 1;
MSD.time_lag(t_lag,1) = 0;
MSD.rho(t_lag,1) = 0;
MSD.N_denominateur(t_lag,1) = nb_points;

for m=1:nbMSD-1
    nb_points_concern = 0;
    S = 0;
    for i=1:nb_points
        clear indice;
        indice = find(trc(:,2)-trc(i,2)==m);
        if not(isempty(indice))
            %indice
            nb_points_concern = nb_points_concern + 1;
            %disp(norm(trc(i,3:4)-trc(indice,3:4),2))
            S = S + ((szpx/1000)^2)*norm(trc(i,3:4)-trc(indice(1),3:4),2).^2;
        end;
    end;

    t_lag = t_lag+1;
    MSD.N_denominateur(t_lag,1) = nb_points_concern;
    MSD.time_lag(t_lag,1)=m.*(till/1000);
    if nb_points_concern > 0
        MSD.rho(t_lag,1)=S/nb_points_concern;
    else
        MSD.rho(t_lag,1)=NaN;
    end;
end;
clear S i indice S nb_points_concern m t_lag;


%===========================================
% calcul de la précision sur la MSD, sigma = ecart type issu du calcul de
% Quian
%=======================================================
sigma=[];
Ft=[];
Ntotal=length(MSD.rho);

%disp(MSD.rho)

for n=1:Ntotal-1 %taille des segments en nombre d'intervalle
    ind_MSD=n+1;% puisque les segments de 2 points correspondent au 2ème pt de la MSD
    Na_theorique=Ntotal-n;
    Na=MSD.N_denominateur(ind_MSD);
    if Na==0
        F=0;
    elseif Na>=n
        F=(4*n*n*Na+2*Na+n-n*n*n)/(6*n*Na*Na);
    elseif Na<n
        F=1+(Na*Na*Na-4*n*Na*Na+4*n-Na)/(6*n*n*Na);
    end
    Ft= [Ft; F];
    sigma=[sigma ; sqrt(F*MSD.rho(ind_MSD)*MSD.rho(ind_MSD))];
end
MSD.sigma=[sigma(1) ; sigma];

%===============================================
% calcul du coeff de diffusion
%=================================================

y=MSD.rho(2:nb_fit);
x=MSD.time_lag(2:nb_fit);
s=MSD.sigma(2:nb_fit);
[pente,b,sqrtChi2sNm2,corr]=RegLinPondQ(x(~isnan(y)),y(~isnan(y)),s(~isnan(s)));
D=pente/4; %%2D!
if b<0 | b>(0.03*0.03) | D<0
    b=0.005*0.005;
    D=(x(~isnan(y))\(y(~isnan(y))-b)/4);
end
if D<0
    D=0;b=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
