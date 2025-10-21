function showtrajectories(x,x2,handles)
%function showtrajectories(x,x2,handles)
% for each frame, makes an array with traces of the molecules and plots them
% called by shift correction
%
% Marianne Renner 07/10 for movtrack
% Marianne Renner 01/2025 - adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pkfile=get(handles.pkradiobutton,'value');
double=get(handles.filetrc2,'string'); %second traj
doubletrc=[];
if isempty(double)==0 % two trajectory data
    doubletrc=get(handles.filetrc2,'userdata');
end
identify=get(handles.ident,'value');
%firstbutton=get(handles.firstradiobutton,'value');
firstbutton=0;

localiz=get(handles.loc,'value');
%timecolor=get(handles.timecolorradiobutton,'value');
timecolor=0;
%blink=get(handles.blinking,'value');
blink=0;

rainbow=get(handles.rainbowradiobutton,'value');
%spine=get(handles.spineradiobutton,'value');
linew=str2num(get(handles.linewidth,'string'));
molnro=get(handles.molindiv,'string');
handles.typetraj=0;

if isempty(x)==0  % if the trc file is loaded
  maxmol=max(x(:,1));
  axes(handles.axes1);
  axis([0 handles.param.Xdim 0 handles.param.Ydim]);
  tlag=str2num(get(handles.till,'string'));
  name=get(handles.file1,'string');
  actualframe=handles.param.actual; %actual frame
  lastframe=handles.param.lasttrc;

  if pkfile==0; %.trc
    codenormal=handles.color.all;

    indexmol=[];
    indexactual=[];
    actualtrc=[];
    colorm=[];
   % if timecolor==1 
   %   [colorm]=createcolor(lastframe);
   % else
    if rainbow==1
     % colorm=handles.colorm;
    else
      codecol=handles.color.all;
    end
    indcol=1;
    
    k=strfind(molnro,'all');
    if isempty(k)==1
       indexmol=find(x(:,1)==str2num(molnro));
       trc=x(indexmol,:);
    else
       trc=x; 
    end

    indexactual=find(trc(:,2)<actualframe+1);
    if isempty(indexactual)==0
        actualtrc=trc(indexactual,:);
        for i=1:max(actualtrc(:,1))
            indexmol=find(actualtrc(:,1)==i);
            indextrc=find(trc(:,1)==i);
            if isempty(indexmol)==0 & max(trc(indextrc,2))>actualframe % present before and after
                codecol=handles.color.all;
                
                %%%%%%%%%%%%%
                if firstbutton==1 % primer punto
                   if rainbow==1
                      handles.typetraj=3;
                      codecol=colorm(indcol,1:3);
                      indcol=indcol+1;
                      if indcol>50
                         indcol=1;
                      end

                      plot(actualtrc(indexmol,3),actualtrc(indexmol,4),'Color',codecol,'Linewidth',linew);
                      hold on
                   end
                   auxtrc=actualtrc(indexmol,:);
                   indextime=find(auxtrc(:,2)==handles.param.actual);
                   %disp(actualtrc(indexmol,2))
                   
                   if isempty(indextime)==0  
                        %disp(auxtrc(indextime(1),:))
                      %  if auxtrc(indextime(1),6)>0 & auxtrc(indextime(1),6)<1000 %slot
                      %    codecol=handles.color.extra; % blue
                      %  else
                      %    codecol=handles.color.all; %red
                      %  end
                        if auxtrc(indextime(1),7)==1 % group 1
                          codecol=handles.color.extra; % blue
                        else
                          codecol=handles.color.all; %red
                        end
                        plot(auxtrc(indextime(1),3),auxtrc(indextime(1),4),'Color',codecol,'Marker','.','MarkerSize',20);
                        hold on
                    end
                    clear auxtrc
                    
                else %not first point
                    %%%%%%
                    
                    if localiz==0 & timecolor==0  
                        if actualtrc(indexmol(size(indexmol,1)),2)<actualframe  & blink==1 % blink period
                            codecol=handles.color.blink;
                            handles.typetraj=2;
                        elseif rainbow==1
                            handles.typetraj=3;
                            codecol=colorm(indcol,1:3);
                            indcol=indcol+1;
                            if indcol>50
                                indcol=1;
                            end
                        end
                        % normal, blinking or raibow
                        plot(actualtrc(indexmol,3),actualtrc(indexmol,4),'Color',codecol,'Linewidth',linew);
                        hold on
                    elseif localiz==0 & timecolor==1  
                        handles.typetraj=4;
                        for g=1:size(indexmol,1)-1
                            codecol=colorm(g,1:3);
                            plot(actualtrc(indexmol(g):indexmol(g+1),3),actualtrc(indexmol(g):indexmol(g+1),4),'Color',codecol,'Linewidth',linew);
                            hold on
                        end
                    elseif localiz==1
                        handles.typetraj=1;
                        if size(trc,2)>5  % trc with localization
                            for f=1:size(indexmol,1)-2
                               % if spine==0 | 
                               if size(trc,2)<7
                                    if actualtrc(indexmol(f+1),6)<0 %peri
                                        codecol=handles.color.peri;
                                    elseif actualtrc(indexmol(f+1),6)>0 %syn
                                        codecol=handles.color.syn;
                                    elseif actualtrc(indexmol(f+1),6)==0 %extra
                                        codecol=handles.color.extra;
                                    end
                               else   %if spine==1
                                    if actualtrc(indexmol(f+1),7)>0 %neck
                                        codecol=handles.color.peri;
                                    elseif actualtrc(indexmol(f+1),7)<0 %head
                                        codecol=handles.color.syn;
                                    elseif actualtrc(indexmol(f+1),7)==0 %extra
                                        codecol=handles.color.extra;
                                    end
                                end %spine
                                plot(actualtrc(indexmol(f):indexmol(f+1),3),actualtrc(indexmol(f):indexmol(f+1),4),'Color',codecol,'Linewidth',linew);
                                hold on
                            end %loop f
                        else
                            plot(actualtrc(indexmol,3),actualtrc(indexmol,4),'Color',codecol,'Linewidth',linew); % no syn
                            hold on
                        end
                    end  % options colorcode
                    %%%%
                end
                
                %% if firstbutton==1 % primer punto
                % plot(actualtrc(indexmol(1),3),actualtrc(indexmol(1),4),'Color',codecol,'Marker','o','MarkerSize',10);
                % end
                if identify==1 % numero tray
                    pos=indexmol(size(indexmol,1));
                    text(actualtrc(pos,3)+1,actualtrc(pos,4)+1,sprintf('%0.0f',actualtrc(pos,1)),'Color',[1 1 0]);
                    cifras=num2str(actualtrc(pos,1)); space=size(cifras,2);
                    text(actualtrc(pos,3)+(space*5),actualtrc(pos,4)+1,sprintf('(%0.0f)',size(actualtrc(indexmol),1)),'Color',[1 1 1],'FontSize',7);
                end
                hold on
            else
                plot(handles.param.Xdim,handles.param.Ydim,'.k');  %just no avoid crash!
                hold on    
            end % empty indexmol
        end % loop actual molecules
    else
        plot(handles.param.Xdim,handles.param.Ydim,'.k');  %just no avoid crash!
        hold on  
    end %empty indexactual

    if isempty(doubletrc)==0 % second traj
       plottracesframe(doubletrc,actualframe,'g',localiz,linew,handles)
    end

  else %pk
    
    accelerate=get(handles.accel,'value');
    if accelerate==1
       x=get(handles.accel,'userdata');
    end

    % frame per frame
    index=find(x(:,2)<actualframe+1);
    if isempty(index)==0
       for i=1:size(index,1)
           if size(x,2)>4
                plot3(x(index(i),3),x(index(i),4),x(index(i),5),'Marker','x','MarkerEdgeColor','r');
           else
                plot(x(index(i),3),x(index(i),4),'Marker','x','MarkerEdgeColor','r');
           end
           hold on
       end
    end
  end %pkfile
  
end % empty x

set(handles.zoompushbutton,'userdata',handles.typetraj);

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%