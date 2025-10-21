function cut=cutbyloctrc(data, minpoints, smooth,option,codeperi);
% function cut=cutbyloctrc(data, minpoints, smooth);
% separates data (trc) into segments with the same localization
% if segments are more than minpoints long, it calculates mean and median
% values
% smooth: in reserve
%
% Marianne Renner sep 09 SPTrack_v4
% Marianne Renner 01/2025 - adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resumen=[];
order=1;
count=1;
ini=1;

if nargin<4
    option=0;
    codeperi=0;
end

if option==0
   loc=data(1,6); %syn
   col=6;
elseif option>0
   loc=data(1,7); %spine
   col=7;
end

% smooth!!!!!!
for pos=2:size(data,1)-5
    if data(pos,col)~=loc
        index=find(data(pos+1:pos+5,col)==loc);
        if size(index,1)>3 %!!!!!!!!!!!!!!!!!!!!!!!!!
            data(pos,col)=loc;
        else
            loc=data(pos,col);
        end
    end
end
% last
for pos=size(data,1)-5:size(data,1)
    data(pos,col)=data(size(data,1)-5,col);
end

loc=data(1,col);


for pos=2:size(data,1)-1
    pres=data(pos,col);
    if pres==loc
    else % change loc
        if pos>ini+5 %!!!!!!!!!!!!!!!!! changes of at least 5 frames
            aux.segment(order).data=data(ini:pos-1,:);
            ini=pos;
            loc=pres; 
            locmin=min(aux.segment(order).data(:,col));
            locmax=max(aux.segment(order).data(:,col));
            if locmin~=locmax % otro smoothing
                indexm=find(aux.segment(order).data(:,col)==locmin);
                if size(indexm,1)>size(aux.segment(order).data(:,col))/2
                    aux.segment(order).data(:,col)=locmin;
                else
                    aux.segment(order).data(:,col)=locmax;
                end
            end
            
            if order>1
                finanterior=size(aux.segment(order-1).data,1);
                if aux.segment(order).data(1,col)==aux.segment(order-1).data(finanterior,col)
                    diftiempo=aux.segment(order).data(1,2)-aux.segment(order-1).data(finanterior,2);
                    if diftiempo<15 %!!!!!!!!!!! fusiona segmentos con la misma loc
                        aux.segment(order-1).data=[aux.segment(order-1).data; aux.segment(order).data];
                        aux.segment(order).data=[];
                    else
                        order=order+1;
                    end %diftiempo
                else
                    order=order+1;
                end %same loc
            else
                order=order+1;
            end %order
        end
    end
end

% last one
if pos>ini
    if order>1
    finanterior=size(aux.segment(order-1).data,1);
                if data(ini,col)==aux.segment(order-1).data(finanterior,col)
                    diftiempo=data(ini,2)-aux.segment(order-1).data(finanterior,2);
                    if diftiempo<15 %!!!!!!!!!!! fusiona segmentos con la misma loc
                        aux.segment(order-1).data=[aux.segment(order-1).data; data(ini:size(data,1),:)];
                        order=order-1;
                    else
                        aux.segment(order).data=data(ini:size(data,1),:);
                    end %diftiempo 
                else
                    aux.segment(order).data=data(ini:size(data,1),:);
                end
    else
        aux.segment(order).data=data(ini:size(data,1),:);
    end
else
    order=order-1;
end
aux.nrosegm=order;

if col==7
    disp('Nro segments (loc by morphology):')
    disp(order)
end

if option==2
   nrospine=order;
   order=1;
   % loc spine + syn
   %col=6;
   for nro=1:nrospine
       clear data
       data=aux.segment(nro).data;
       ini=1;
     if size(data,1)>0
       loc=data(1,6); %syn
       for pos=2:size(data,1)-1
           pres=data(pos,6);
           if codeperi==0
               if pres<0
                   pres==0; %peri=extra
               end
           else
               if pres<0
                   pres=-pres; %peri=syn
               end
           end
           if pres==loc
           else % change loc
               if pos>ini+5 %!!!!!!!!!!!!!!!!! changes of at least 5 frames
                   aux2.segment(order).data=data(ini:pos-1,:);
                   ini=pos+1;
                   loc=pres; 
                   locmin=min(aux2.segment(order).data(:,6));
                   locmax=max(aux2.segment(order).data(:,6));
                   if locmin~=locmax % otro smoothing
                       indexm=find(aux2.segment(order).data(:,6)==locmin);
                       if size(indexm,1)>size(aux2.segment(order).data(:,6))/2
                           aux2.segment(order).data(:,6)=locmin;
                       else
                           aux2.segment(order).data(:,6)=locmax;
                       end
                   end
                   cut.segment(order).data=aux2.segment(order).data;
                   ini=pos+1;
                   loc=pres; 
                   order=order+1;
               end %pos
           end %change loc
       end %loop data
       
       % last
       cut.segment(order).data=data(ini:pos-1,:);
       cut.nrosegm=order;
       order=order+1;
       
     end %emptytrc

   end %loop segments spine
   
   disp('Nro segments (loc by morphology and domains)')
   disp(order-1)

else
  cut=aux;
end

cut.data=data;


clear aux

%%%%%%%%%%%%%%%%%%%%%%%%%%%


            
