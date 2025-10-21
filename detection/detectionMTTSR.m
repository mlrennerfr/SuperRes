function [fr,x,y,alpha,radius, radiusmean,sigma, blink, sommevalid] = detectionMTTSR(name_stk , stack, Nb_image_ds_stack, seuil_detec_1vue, wn,...
    r0, nb_defl, seuil_alpha, activation_sig_fit)
%function [fr,x,y,alpha,radius, radiusmean,sommevalid] = detectionMTTSR(name_stk , stack, Nb_image_ds_stack, seuil_detec_1vue, wn,...
%    r0, nb_defl, seuil_alpha, activation_sig_fit)
%
% Detection from MTT programs
% Adapted for SuperRes (2016)
% fev 2020: adapted to read large .tiff files
%
% Marianne Renner fev 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%k=strfind(name_stk,'nd2');

waitbarhandle=waitbar( 0,'Please wait...','Name',['Peak detection in ',name_stk]) ;


for ind=1:Nb_image_ds_stack
   
    if exist('waitbarhandle')
       waitbar(ind/Nb_image_ds_stack,waitbarhandle,['Frame # ',num2str(ind)]);
    end
    
    %if k==1
    if isempty(stack)==0
        input = stack{ind, 1};
        input=double(input);
    else
       input = imread(name_stk,ind);
      %  input = imread_big(name_stk,[ind ind]);
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
        tab_param = NaN(8, 8) ;
        tab_var = NaN(8, 8) ;
   else
        nb_valid = max(size(ind_valid));
        sommevalid=sommevalid + nb_valid;
        %==========================================================================
        %% on alloue les tableaux de sorties
        tab_param = zeros(nb_valid, 8) ;
        tab_var = zeros(nb_valid, 8) ;
        %==========================================================================
        %% on initialise les tableaux de sorties
        tab_param(:,1) = (1:nb_valid)' ;
        tab_param(:,2:8) = [ind.*ones(nb_valid,1), ...
            lestime(ind_valid,[2,3,4,6]), ... %% i j alpha rayon
            lestime(ind_valid,5), ... %sig2_b
            Nb_image_ds_stack*ones(nb_valid,1)];

        %% valeurs initiales (par default)
        %% dans calcul_reference(traj, t, param, T)
        tab_var(:,1)   = (1:nb_valid)' ;
        tab_var(:,2:8) = [ones(nb_valid,1), ...
            zeros(nb_valid,4), ...
            lestime(ind_valid,5), ...  %% sig2_b
            Nb_image_ds_stack*ones(nb_valid,1) ];
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
        sigma2=[radius, Structure(ind).param(:,7)'];
        blink2=[radius, Structure(ind).param(:,8)'];
        
        controllength=size(x,2);
        sigma=sigma2(1:controllength);
        blink=blink2(1:controllength);
        
      %  disp(size(blink))
        
        clear param var moy tab_var tab_moy tab_param
    end
end

radiusmean = median(radius);


close(waitbarhandle);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



