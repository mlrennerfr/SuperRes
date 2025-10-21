function [datamatrix,Xdim,Ydim,nfram,stktrue]=readimage(filename)
% function [datamatrix,Xdim,Ydim,nfram,stktrue]=readimage(filename)
% reads .stk, .tif, .spe...
% Marianne Renner 09/09 for SPTrack v4
% Marianne Renner 01/2025 - adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stktrue=0;
%msgbox('Reading file...');
answer=findstr(filename,'.spe'); answerb=findstr(filename,'.SPE');
if isempty(answer)==1 && isempty(answerb)==1
            answer2=findstr(filename,'.stk');
            if isempty(answer2)==1
                answer3=findstr(filename,'.tif');
                if isempty(answer3)==1
                    answer4=findstr(filename,'.mat');
                    if isempty(answer4)==1
                       msgbox('Wrong type of file','Error','error');
                       return
                    end
                else
                    
                    % .tif file       
                    info=imfinfo(filename);
                    if size(info,1)>1 % movie tif
                        stktrue=1;
                        [stack_info,datamatrix] = stkdataread(filename);
                        nfram=stack_info.frames;
                    else % une image
                        [stack_info,datamatrix] = tifdataread(filename);
                        stktrue=2;
                        nfram=1;
                    end
                    
                    Xdim=stack_info.x;
                    Ydim=stack_info.y;
                    if isfield(datamatrix,'data') && nfram==1
                        image=datamatrix.data;
                       % disp(size(image))
                        clear datamatrix
                        datamatrix=image;
                        clear image
                    end
                    [fil,col]=size(datamatrix);              
                    if col/Xdim==3  %rgb
                        stktrue=3;
                    else
                        if isstruct(datamatrix)
                        else
                            datamatrix=double(datamatrix);
                        end
                    end
                    
                end
            else
                 % .stk file       
                 [stack_info,datamatrix] = stkdataread(filename);
                 Xdim=stack_info.x;
                 Ydim=stack_info.y;
                 nfram=stack_info.frames;
                 stktrue=1;
             end
        else
            % .spe file
            [datamatrix p]= spedataread (filename);
            Xdim=p(1);
            Ydim=p(2)/p(4);
            nfram=p(4);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

