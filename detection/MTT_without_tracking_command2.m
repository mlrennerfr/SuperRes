function [seuil_detec_1vue, wn, r0, nb_defl, seuil_alpha, filtering_threshold, radiusmean] = MTT_without_tracking_command2(name_stk , Nb_image_ds_stack, seuil_detec_1vue, wn,...
    r0, nb_defl, seuil_alpha, output_dir, filtering_threshold, activation_sig_fit);
%
% from MTT
%
% Modified by Marianne Renner for SuperRes programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%if strcmp(VERSION, 'MATLAB')
%    stderr = 1 ;
%else
 %   clear stderr ;
%end %if

%chemin_init=cd;

%nbre_image=Nb_STK;

wn=wn;
r0=r0;
pfa=seuil_detec_1vue;
n_deflt=nb_defl;
demi_wn = ceil(wn/2);
seuil_alpha = seuil_alpha ;
%name_stk = stack;
%pos=[];
x=[];
y=[];
alpha=[];
fr=[];
radius=[];
sommevalid=0;

waitbarhandle=waitbar( 0,'Please wait...','Name',['Peak detection in ',name_stk]) ;


t0=cputime;

for ind=1:Nb_image_ds_stack
    
    disp(ind)
    
    if exist('waitbarhandle')
       waitbar(ind/Nb_image_ds_stack,waitbarhandle,['Frame # ',num2str(ind)]);
    end

    %disp(name_stk);
    %fprintf(stderr, 'image n° %s\n', num2str(ind)) ;
    input = imread(name_stk,ind);
    input = double(input);
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
end



if Nb_image_ds_stack==1
    figure
    imagesc(input), colormap(gray)
    hold on
    plot(y+1, x+1, 'r*', 'MarkerSize', 12)
    axis equal
    axis image
    axis tight
    radiusmean = r0;
else
    if isempty(x) && isempty(y) && isempty(alpha) && isempty(fr)
        warndlg('No particule detected')
    else

        %mkdir(output_dir)
        %cd(output_dir)
        radiusmean = median(radius)
        matrice_results(1,:)=fr;
        matrice_results(2,:)=x;
        matrice_results(3,:)=y;
        matrice_results(4,:)=alpha;
        matrice_results(5,:)=radius;

        name=[name_stk,'_','peaks'];
        save([name, '.mat'], 'matrice_results')
    end
end

%t1=cputime;
%h=floor((t1-t0)/3600);
%m=floor(rem((t1-t0),3600)/60);
%s=round(rem(rem((t1-t0),3600),60));

disp(['Mean number of valid particles: ',num2str(sommevalid/Nb_image_ds_stack)]);

close(waitbar)

%calculation_time=sprintf('%.2dh%.2dm%.2ds\n', h, m, s)
%cd(repertoire)

%clear x y fr param var moy tab_var tab_moy tab_param Structure sommevalid



