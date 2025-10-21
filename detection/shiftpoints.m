function newtrc=shiftpoints(handles,ktev,kteh,trc)
% function newtrc=shifttrc(handles,ktev,kteh,trc)
% shift correction for trc
%
% Marianne Renner - 07/10 for movtrack.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newtrc=trc;

if isempty(trc)==0
        if ktev>0 % down 
            if kteh>0 %right
                newtrc(:,1)=newtrc(:,1)+kteh;
            elseif kteh<0 % left
                kteh=-kteh;
                newtrc(:,1)=newtrc(:,1)-kteh;
            elseif kteh==0
            end
            newtrc(:,2)=newtrc(:,2)+ktev;
        elseif ktev<0 % up 
            ktev=-ktev;
            if kteh>0 %right
                newtrc(:,1)=newtrc(:,1)+kteh;
            elseif kteh<0 % left
                kteh=-kteh;
                newtrc(:,1)=newtrc(:,1)-kteh;
            elseif kteh==0
            end
             newtrc(:,2)=newtrc(:,2)-ktev;
        elseif ktev==0
            if kteh>0 %right
                newtrc(:,1)=newtrc(:,1)+kteh;
            elseif kteh<0 % left
                kteh=-kteh;
                newtrc(:,1)=newtrc(:,1)-kteh;
            elseif kteh==0
            end
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
