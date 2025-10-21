function readnd2image(file)

data = bfopen(file);

seriesCount = size(data, 1);
series1 = data{1, 1};

series1_planeCount = size(series1, 1);
series1_plane1 = series1{1, 1};
series1_label1 = series1{1, 2};

stackmin=min(min(series1_plane1));
stackmax=max(max(series1_plane1));
figure
imshow(series1_plane1,[stackmin stackmax], 'InitialMagnification','fit')