function showtrajavizoom(x,x2,handles)
%function showtrajavizoom(x,x2,handles)
% for each frame, makes an array with traces of the molecules and plots them
% Marianne Renner 01/2025 - adapted to SuperRes_v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

molnro=handles.nromolecule;
if handles.typetraj==5
    molnro=num2str(molnro);
end
linew=str2num(get(handles.linewidth,'string'));

axes(handles.axes1);
selecD=get(handles.Dcolor,'value');
selecnumber=get(handles.Dnumbers,'value');
if handles.typetraj==5 && selecnumber==1
    handles.typetraj=0;
end



if isempty(x)==0
    
   lastframe=handles.param.lasttrc;
   actualframe=get(handles.nroframe,'value');
   if actualframe==0;
      actualframe=1;
   end

   if size(x,2)>4 %trc
       
       if handles.typetraj==5 & selecD==1   % Dinst in color
           
            % c%definicion colores
            categories=6;
            [colorstep, cmap]=definecolor(0,[],categories); % create color map
            index=find(x(:,2)<actualframe+1);
            if isempty(index)==0
               for j=1: size(index,1)
                   D=x(j,7); % Dinst
                   [colorstep, cmap]=definecolor(D,cmap,0);
                   plot ((x (j,3)), (x (j,4)),'color',colorstep,'Marker','.','MarkerSize',5);   % grafica puntos
                   hold on
               end
            end
            
       else % trajectories
           
            maxmol=max(x(:,1));  % numero corregido
            indexmol=[];
            indexactual=[];
            actualtrc=[];
            colorm=[];
            indcol=1;
            localiz=0;
            identify=0;
            firstbutton=0;
            rainbow=0;
      
            if handles.typetraj==1 | handles.typetraj==6
                localiz=1;
            elseif handles.typetraj==4 %timecolor
                [colorm]=createcolor(lastframes);
            elseif handles.typetraj==3 % rainbow
                colorm=get(handles.saveimagepushbutton,'userdata');
            elseif handles.typetraj==7 
                firstbutton=1;
                colorm=get(handles.saveimagepushbutton,'userdata');%!!!!!!!!!!!!!!!!!!
                rainbow=1;
            end
            
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

                    %  plot(actualtrc(indexmol,3),actualtrc(indexmol,4),'Color',codecol,'Linewidth',linew);
                   %   hold on
                   end
                    auxtrc=actualtrc(indexmol,:);
                    indextime=find(auxtrc(:,2)==handles.param.actual);
                    %disp(actualtrc(indexmol,2))
                    if isempty(indextime)==0  
                        %disp(auxtrc(indextime(1),:))
                       % if auxtrc(indextime(1),6)>0 & auxtrc(indextime(1),6)<1000 %slot
                       %   codecol=handles.color.extra; % blue
                      %  else
                      %    codecol=handles.color.all; %red
                      %  end 
                        if auxtrc(indextime(1),7)==1 % group1
                          codecol=handles.color.extra; % blue
                        else % group 2
                          codecol=handles.color.all; %red
                        end 

                        plot(auxtrc(indextime(1),3),auxtrc(indextime(1),4),'Color',codecol,'Marker','.','MarkerSize',20);
                        hold on
                    end
                    clear auxtrc
                else
                    %%%%%%

                        if localiz==0 & handles.typetraj<4
                            if actualtrc(indexmol(size(indexmol,1)),2)<actualframe  & handles.typetraj==2 % blink period
                                codecol=handles.color.blink;
                            elseif handles.typetraj==3
                                codecol=colorm(indcol,1:3);
                                indcol=indcol+1;
                                if indcol>50
                                    indcol=1;
                                end
                            end
                            % normal, blinking or raibow
                            plot(actualtrc(indexmol,3),actualtrc(indexmol,4),'Color',codecol,'Linewidth',linew);
                            hold on
                        elseif localiz==0 & handles.typetraj==4  
                            for g=1:size(indexmol,1)-1
                                %codecol=colorm(g,1:3);
                                codecol=colorm(i,1:3);
                                plot(actualtrc(indexmol(g):indexmol(g+1),3),actualtrc(indexmol(g):indexmol(g+1),4),'Color',codecol,'Linewidth',linew);
                                hold on
                            end
                        elseif localiz==1 
                            if size(trc,2)>5  % trc with localization, deco or not
                                for f=1:size(indexmol,1)-2
                                    if handles.typetraj<6 | size(trc,2)<7
                                        if actualtrc(indexmol(f+1),6)<0 %peri
                                            codecol=handles.color.peri;
                                        elseif actualtrc(indexmol(f+1),6)>0 %syn
                                            codecol=handles.color.syn;
                                        elseif actualtrc(indexmol(f+1),6)==0 %extra
                                            codecol=handles.color.extra;
                                        end 
                                    elseif handles.typetraj==6
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
                                end
                            else
                                plot(actualtrc(indexmol,3),actualtrc(indexmol,4),'Color',codecol,'Linewidth',linew); % no syn
                                hold on
                            end
                        end  % options colorcode   
                                        %%%%

                end


                        if identify==1 % numero tray
                            pos=indexmol(size(indexmol,1));
                            text(actualtrc(pos,3)+1,actualtrc(pos,4)+1,sprintf('%0.0f',actualtrc(pos,1)),'Color',[1 1 0]);
                            cifras=num2str(actualtrc(pos,1)); space=size(cifras,2);
                            text(actualtrc(pos,3)+(space*5),actualtrc(pos,4)+1,sprintf('(%0.0f)',size(actualtrc(indexmol),1)),'Color',[1 1 1],'FontSize',7);
                        end
                        
                        if selecnumber==1 % D inst
                            pos=indexmol(size(indexmol,1));
                            cifras=num2str(actualtrc(pos,1)); space=size(cifras,2);
                            text(actualtrc(pos,3)+0.5,actualtrc(pos,4)+0.5,sprintf('%1.4f',actualtrc(pos,7)),'Color',[1 1 1],'FontSize',7);
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
            


            if isempty(x2)==0 % second traj
                codecol='r';
                synflag=0;
                plottracesframe(x2,actualframe,codecol,synflag,linew,handles)
            end
       end% % Dinst

  else %pk
    
    accelerate=get(handles.accel,'value');
    if accelerate==1
       x=get(handles.accel,'userdata');
    end
    % frame per frame
    index=find(x(:,2)<actualframe+1);
    if isempty(index)==0
       for i=1:size(index,1)
           plot(x(index(i),3),x(index(i),4),'Marker','x','MarkerEdgeColor','r');
           hold on
       end
    end
  end %pkfile

  plot(1,1,'.k');  %just no avoid crash!
  hold on    

end


clear actualtraces x
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%