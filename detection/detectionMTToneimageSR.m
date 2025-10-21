function [matrice_results, radiusmean,sommevalid] = detectionMTToneimageSR(name_stk,stack, ind , Nb_image_ds_stack, seuil_detec_1vue, wn,...
    r0, nb_defl, seuil_alpha, activation_sig_fit)
%function [matrice_results, radiusmean,sommevalid] = detectionMTToneimageSR(name_stk,stack, ind , Nb_image_ds_stack, seuil_detec_1vue, wn,...
%    r0, nb_defl, seuil_alpha, activation_sig_fit)
%
% MTT algorythm
% adapted to SuperRes
% fev 2020: adapted to read large .tiff files
%
% Marianne Renner fev 20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



pfa=seuil_detec_1vue;
n_deflt=nb_defl;
demi_wn = ceil(wn/2);
x=[];
y=[];
alpha=[];
fr=[];
radius=[];
sommevalid=0;
matrice_results=[];

    if isempty(stack)==0
        input = stack{ind, 1};
        input=double(input);
    else
       % input = imread(name_stk,ind);
        input = imread_big(name_stk,ind);
        input = double(input);
    end

    [idim, jdim] = size(input) ;
    
    [lest,ldec,dfin] = detect_et_estime_part_1vue (input, wn, r0, pfa, activation_sig_fit) ;
    
    input_deflt = deflat_part_est(input, lest, wn);
    lestime1 = lest;
    if n_deflt==0
        lestime= lestime1;
    else
        for n=1:n_deflt
            
            [l,ld,d,N] = detect_et_estime_part_1vue (input_deflt, wn, r0, pfa, activation_sig_fit) ;
            
            lestime = [lestime1 ; l];
            dfin = dfin | d ;
            
            input_deflt = deflat_part_est(input_deflt, l, wn);
            
        end%for
    end
    test_all = lestime(:,7) & (lestime(:,4)>seuil_alpha) & ...
        (lestime(:,2)>demi_wn) & (lestime(:,2)<idim-demi_wn) & ...
        (lestime(:,3)>demi_wn) & (lestime(:,3)<jdim-demi_wn) ;

    [ind_valid, tmp] = find(test_all);
    
    if isempty(ind_valid)
        %sprintf('No particule')
        %nb_valid=0;
        tab_param = NaN(8, 8) ;
        tab_var = NaN(8, 8) ;
   else
        nb_valid = max(size(ind_valid));
        sommevalid=sommevalid + nb_valid;
      %  fprintf(stderr, 'nb validated particles: %d\n', nb_valid) ;
        %==========================================================================
        %% on alloue les tableaux de sorties
        tab_param = zeros(nb_valid, 8) ;
        tab_var = zeros(nb_valid, 8) ;
        %tab_moy = zeros(nb_valid, 1+7*(T-T_off+1)) ;
        %==========================================================================
        %% on initialise les tableaux de sorties
        tab_param(:,1) = (1:nb_valid)' ;
        tab_param(:,2:8) = [ind.*ones(nb_valid,1), ...
            lestime(ind_valid,[2,3,4,6]), ... %% i j alpha rayon
            lestime(ind_valid,5), ...
            Nb_image_ds_stack*ones(nb_valid,1)];
        %tab_param(:,2+7*(5+15)) = 2*ones(nb_valid,1) ;

        %% valeurs initiales (par default)
        %% dans calcul_reference(traj, t, param, T)
        tab_var(:,1)   = (1:nb_valid)' ;
        tab_var(:,2:8) = [ones(nb_valid,1), ...
            zeros(nb_valid,4), ...
            lestime(ind_valid,5), ...  %% sig2_b
            Nb_image_ds_stack*ones(nb_valid,1) ];
        %tab_var(:,2+7*(5+15)) = 2*ones(nb_valid,1) ;
        tab_moy = tab_var ;

        %==========================================================================
        param((1:nb_valid),1)=tab_param((1:nb_valid),1);
        param((1:nb_valid),(2:8))=tab_param((1:nb_valid),2:8);
        var((1:nb_valid),1)=tab_var((1:nb_valid),1);
        var((1:nb_valid),(2:8))=tab_var((1:nb_valid),2:8);
        moy((1:nb_valid),1)=tab_moy((1:nb_valid),1);
        moy((1:nb_valid),(2:8))=tab_moy((1:nb_valid),2:8);

        Structure(ind)=struct('param', {param}, 'var', {var}, 'moy', {moy});
        %==========================================================================

        x=[x, Structure(ind).param(:,3)'];
        y=[y, Structure(ind).param(:,4)'];
        alpha=[alpha, Structure(ind).param(:,5)'];
        fr=[fr, Structure(ind).param(:,2)'];
        radius=[radius, Structure(ind).param(:,6)'];
        
        clear param var moy tab_var tab_moy tab_param
    end
%end

radiusmean = median(radius);
if isempty(fr)==0
    matrice_results(1,:)=fr;
    matrice_results(2,:)=x;
    matrice_results(3,:)=y;
    matrice_results(4,:)=alpha;
    matrice_results(5,:)=radius;
    %control=1;
   % if Nb_image_ds_stack==1
   %     figure
   %     imagesc(input), colormap(gray)
   %     hold on
   %     plot(y+1, x+1, 'r*', 'MarkerSize', 12)
   %     axis equal
   %     axis image
   %     axis tight
   %     radiusmean = r0;
   % else
    %    if isempty(x) && isempty(y) && isempty(alpha) && isempty(fr)
         %   warndlg('No particule detected')
    %    end
   % end
end

%t1=cputime;
%h=floor((t1-t0)/3600);
%m=floor(rem((t1-t0),3600)/60);
%s=round(rem(rem((t1-t0),3600),60));

%close(waitbarhandle);

%calculation_time=sprintf('%.2dh%.2dm%.2ds\n', h, m, s)
%cd(repertoire)

%clear x y fr param var moy tab_var tab_moy tab_param Structure sommevalid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

