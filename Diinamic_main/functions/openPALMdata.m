function [x,y, fr, alpha, radius,sigma, blink, ratio, z,test1,test2, ffname,filename,control]=openPALMdata(handles,px_mu,file2,ffname)

control=1;
x=[];
y=[];
fr=[];
alpha=[];
radius=[];
sigma=[];
blink=[];
ratio=[];
z=[];
test1=[];
test2=[];

message=msgbox('Please wait','Data','Warn');

if nargin<4
    [filename, pathname] = uigetfile('*.mat','Select .mat file containing detection results');
    if filename==0
        if file2==1
            control=0;
        end
        ffname=[];
        filename=[];
        return
    end
    ffname = fullfile(pathname,filename);
else
end

[pathstr, filename, ext] = fileparts(ffname);

    if strcmp(ext,'.mat')
        S = load(ffname);
        if isfield(S,'matrice_results')
            aux = S.matrice_results;
            clear S;
            fr = aux(1,:); fr = fr(:);
            x = aux(3,:);
            x = x(:);
            y = aux(2,:); 
            y = y(:);
            alpha = aux(4,:); alpha = alpha(:);
            if size(aux,1)>4
                radius=(aux(5,:));
            else
                radius=zeros(size(x,2),1);
            end
            if size(aux,1)>5
                sigma=aux(6,:)';
                blink=aux(7,:)';
                ratio=aux(8,:)';
                z=aux(9,:)';
                test1=aux(10,:)';
                test2=aux(11,:)';
            else
                sigma=zeros(size(x,2),1);
                 blink=zeros(size(x,2),1);
                ratio=zeros(size(x,2),1);
                z=zeros(size(x,2),1);
                test1=zeros(size(x,2),1);
                test2=zeros(size(x,2),1);
            end

        elseif isfield(S,'alpha') && isfield(S,'fr') && isfield(S,'x') && isfield(S,'y')
            
            x = S.x; % * px_mu;
            x = x(:);
            y = S.y; % * px_mu;
            y = y(:);
            alpha = S.alpha; alpha = alpha(:);
            fr = S.fr; fr = fr(:);
            if isfield(S,'radius')
                radius=S.radius;
            else
                radius=zeros(size(x,2));
            end
            sigma=zeros(size(x,2));
            blink=zeros(size(x,2));
            ratio=zeros(size(x,2));
            z=zeros(size(x,2));
            test1=zeros(size(x,2));
            test2=zeros(size(x,2));
           
        else
            disp('No data loaded');
            if file2==1
                control=0;
            end
            return
        end
    end
%end
close(message)

guidata(gcbo,handles) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
