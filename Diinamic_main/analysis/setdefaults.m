function handles=setdefaults(Densparam, handles)

if Densparam.maxdiamcluster==0
    % no previous parameters
    inthresh=1;
    limitdens=0.1;
    mindenspx=0.1;
    pxdilate=1;
    pxerode=1;
    minnrodetect=30;
    mindiamcluster=10;
    maxdiamcluster=10000;
    polsize=10;
    voro=0;
    epsilon=1;
    autoep=1;
    minpointsnano=50;
    donano=0;
else
    inthresh=Densparam.inthresh;
    limitdens=Densparam.limitdens;
    mindenspx=Densparam.mindenspx;
    pxdilate=Densparam.pxdilate;
    pxerode=Densparam.pxerode;
    minnrodetect=Densparam.minnrodetect;
    mindiamcluster=Densparam.mindiamcluster;
    maxdiamcluster=Densparam.maxdiamcluster;
    polsize=Densparam.vorosize;
    voro=Densparam.voro;
    epsilon=Densparam.epsilon;
    autoep=Densparam.autoepsilon;
    minpointsnano=Densparam.minpointsnano;
    donano=Densparam.nano;
end
    
set(handles.pxdilate,'String',pxdilate);
set(handles.pxerode,'String',pxerode);
set(handles.mindens,'String',limitdens);
if voro==0
    set(handles.mindenspx,'String',mindenspx);
    set(handles.mindenspx2,'String','0');
else
    set(handles.mindenspx,'String','0');
    set(handles.mindenspx2,'String',mindenspx);
end
set(handles.mindetect,'String',minnrodetect);
set(handles.minsize,'String',mindiamcluster);
set(handles.maxsize,'String',maxdiamcluster);
set(handles.intenthresh,'String',inthresh);
if isfield(handles,'autoepsilonradiobutton')
    set(handles.autoepsilonradiobutton,'Value',autoep);
    if autoep==0
        epsilon=0;
    end
    set(handles.epsilon,'String',epsilon);
    set(handles.nanoradiobutton,'Value',donano);
    set(handles.mindetecnano,'String',minpointsnano);
end
set(handles.vororadiobutton,'Value',voro);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
