function reportclustering(Densparam,filenames)
%function reportclustering(Densparam,filenames)
%
% creates .txt file with parameters values of DetectClusters analysis
%
% Marianne Renner oct 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str = date;
name=['Parameters_clustering_',str,'.txt'];

fi = fopen(name,'wt');
if fi<3
    error('File not found or readerror.');
end

report{1}= 'Parameters for clustering analysis';
report{2}= ' ';

if Densparam.voro==1
    report{3}= ['Analysis of density after Voronoi tesselation. Minimum polygon size: ',num2str(Densparam.vorosize)];
    report{4}=' ';
else
    report{3}=['Analysis of density after segmentation with rendered image. Intensity threshold: ',num2str(Densparam.inthresh),'%. Density threshold: ',num2str(Densparam.mindenspx)];
    report{4}= ['Segmented image improved with dilation: ',num2str(Densparam.pxdilate),' px and erosion: ',num2str(Densparam.pxerode),' px.'];
end
report{5}=' ';
report{6}= 'Cluster selection criteria:';
report{7}= ['Minimum density: ',num2str(Densparam.limitdens), '. Minimum number of detections: ',num2str(Densparam.minnrodetect)];
report{8}= ['Minimum cluster size (pixels): ',num2str(Densparam.mindiamcluster),'. Maximum cluster size (pixels): ',num2str(Densparam.maxdiamcluster)];
report{9}= ' ';

if Densparam.nano==1
    report{10}= ['Including nanodomain detection by DBSCAN. Minimum number of detections in the nanodomain : ',num2str(Densparam.minpointsnano)];
    if Densparam.autoepsilon==1
        report{11}= ['Search radius calculated automatically : ',num2str(Densparam.epsilon)];
    else
        report{11}= ['Search radius : ',num2str(Densparam.epsilon)];
    end
else
    report{10}='No nanodomain detection';
    report{11}=' ';
end

report{12}= 'Files:';

for celda=1:12
    fseek(fi,20,0);
    fprintf(fi,'%20s\n',report{celda});
end

for celda=1:size(filenames,2)
    fseek(fi,0,0);
    fprintf(fi,'%20s\n',filenames{celda});
end
fclose(fi);

