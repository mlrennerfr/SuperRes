function dataSM=loadSMdata(file)

S = load(file);

if isfield(S,'matrice_results')
            aux = S.matrice_results;
            clear S;
            fr = aux(1,:); fr = fr(:);
            x = aux(3,:);
            x = x(:);
            y = aux(2,:); 
            y = y(:);
            alpha = aux(4,:); alpha = alpha(:);
            radius=(aux(5,:)); radius = radius(:);
            sigma=(aux(6,:)); sigma=sigma(:);
            blink=(aux(7,:)); blink=blink(:);
            if size(aux,1)>7
                ratio=aux(8,:); ratio=ratio(:);
                z=aux(9,:);z=z(:);
                test1=aux(10,:);test1=test1(:);
                test2=aux(11,:);test2=test2(:);
            else
                ratio=[];
                z=[];
                test1=[];
                test2=[];
            end

else
            disp('Invalid file!');
            return
end

clear S

dataSM.x = x;
dataSM.y = y;
dataSM.alpha = alpha;
dataSM.fr = fr;
dataSM.radius=radius;
dataSM.sigma=sigma;
dataSM.blink=blink;
dataSM.ratio=ratio;
dataSM.z=z;
dataSM.test1=test1;
dataSM.test2=test2;