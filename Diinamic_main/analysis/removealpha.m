function [x,y,alpha,fr]=removealpha(vect2remove,x,y,alpha,fr)
%% Removing unneeded points from column vectors
 
x(vect2remove) = [];
y(vect2remove) = [];
if nargin>3
        alpha(vect2remove) = [];
        fr(vect2remove) = [];
else
        alpha=[];
        fr=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    