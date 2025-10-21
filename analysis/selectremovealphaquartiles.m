function [x,y]=selectremovealphaquartiles(pctmin,pctmax,x,y,alpha)
% function [x,y]=selectremovealphaquartiles(pctmin,pctmax,x,y,alpha)
%
% from PALMvis
%
% Marianne Renner 01/2025 - adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Remove points with alphas outside specified percentile range
vect2remove = [];
%pctmin = str2num(get(handles.alpha1,'String'));
%pctmax = str2num(get(handles.alpha2,'String'));
%if exist('Xmatrix','var')
%    aux = alphamatrix(:);
%else
    aux = alpha;
%end
aux = aux(aux>0);
alphamin = prctile(aux,pctmin);
alphamax = prctile(aux,pctmax);  
aux = union(find(alpha<alphamin),find(alpha>alphamax));
vect2remove = [vect2remove; aux];

if ~isempty(vect2remove)
    %% Removing unneeded points from column vectors
    x(vect2remove) = [];
    y(vect2remove) = [];
    %if nargin>3
     %   alpha(vect2remove) = [];
     %   fr(vect2remove) = [];
    %else
    %    alpha=[];
    %    fr=[];
    %end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

