function convertpktolama
% function convertpktolama
% prepares data file to be read and analysed by LaMA
%
% Marianne Renner 01/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input data
[file,path] = uigetfile('*.pk*','Load detections data (.pk)');
filename = [path,file];
if file==0
    return
end

[namefile,rem]=strtok(file,'.');
pk=load(file);

newdata=[pk(:,2) pk(:,3) pk(:,1) pk(:,4)];

save([namefile,'-lama.txt'],'newdata','-ascii')

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
