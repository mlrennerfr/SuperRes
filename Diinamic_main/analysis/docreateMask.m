function docreateMask(handles)
%function docreateMask(handles)
%
% create mask from rendered image
% called from CreateROI.m
%
% Marianne Renner 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%mask form rendered image
savecode=get(handles.savematradiobutton,'Value');
szpx=str2num(get(handles.szpx,'string'));
dx = str2num(get(handles.PALMszpx,'string'));
percinthresh = str2num(get(handles.intensthresh,'string'));

sigmaloc = 2*dx;
minsize = str2num(get(handles.minsizeclumask,'String')); %in nm
areanm=(minsize/2)^2*pi;
areapx=(dx*1000)^2 ; %nm^2
minarea=areanm/areapx;

filename=handles.ffname;
[namefile,rem]=strtok(filename,'.') ;
mu=get(handles.muradiobutton,'Value');

% rendered
[I,rendx, rendy]=makerendered2(handles,szpx,sigmaloc,dx,mu);

%disp('Size (SMLM px)')
%disp(size(I,1))
%disp(size(I,2))

valintthresh=percinthresh/100*max(max(I));
BW=zeros(size(I,1),size(I,2));
for i=1:size(I,1)
    indexj=find(I(i,:)>valintthresh);
    if isempty(indexj)==0
        BW(i,indexj)=1;
    end
end
       
aux=zeros(size(BW,1),size(BW,2));
I2=aux;

level = graythresh(BW);
BW1 = im2bw(BW,level); %second image: first mask
for g=1:size(BW,2)
    for t=1:size(BW,1)
        aux(t,g)=BW(abs(t-size(BW,1))+1,g);
    end
end
stats=regionprops(BW1,'PixelList','Area');
nroclu=size(stats,1);

for n=1:nroclu
    D=size(stats(n,1).PixelList,1);
    areaclu=stats(n,1).Area ;
    if areaclu>=minarea %size area in pixels...
        v1=0;
        m=1;
        while m<=D
            px=stats(n,1).PixelList(m,:);
            I2(px(2),px(1))=1;
            m=m+1;
        end
    end
end % nroclu
    
maskfig=figure;
imshow(I2,'InitialMagnification','fit')
hold on

plot(rendx+1,rendy+1,'.r','MarkerSize',2);
%plot(handles.x+1),handles.y+1,'.r','MarkerSize',2);
axis off
hold off

%save
qstring='Accept?';
button = questdlg(qstring);
if strcmp(button,'Yes')
    Savename=[namefile,'-mask.tif'];
    imwrite(I2,Savename,'tif','Compression','none');  % save the image to the disk
    disp(['Mask image saved as ',Savename]);
    
    if savecode ==1
        % extract detections colocalizing with the mask and save a new .mat
        % file
        Savename=[namefile,'-mask.mat'];
        
        indexin=[];
        
        for i=1:size(rendx,1)
                if floor(rendx(i))<1; rendx(i)=1; end
                if floor(rendy(i))<1; rendy(i)=1; end
               % disp(rendx(i))
               % disp(rendy(i))
                if ~isnan(rendx(i)) && ~isnan(rendy(i))
                    if I2(floor(rendy(i)),floor(rendx(i)))==1 % coloc with mask
                        indexin=[indexin;i];
                    end
                end
        end
        
       % figure
        %imshow(I2,'InitialMagnification','fit')
        %hold on
        %plot(rendx(indexin),rendy(indexin),'.r')
        
        savematmask(handles,indexin, mu,rendx(indexin),rendy(indexin), Savename)
        disp(['File with detections co-localizing with the mask saved as ',Savename]);

        
    end

end

clear X Y alpha fr I I2 BW

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
