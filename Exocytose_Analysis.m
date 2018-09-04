    function Exocytose_Analysis
    %%
    if size(findobj('name','Exocytose Analysis'),1)>0
        close(findobj('name','Exocytose Analysis'))
    end
    clear all
    close all

    hexo=figure('name','Exocytose Analysis');
    %     set(hexo,'MenuBar','OFF')
    screensize = get( 0, 'Screensize' );
%     set(hexo,'Position',[33 88 560 225]);
    set(hexo,'Position',[0.0179*screensize(3) 0.0857*screensize(4) 0.3333*screensize(3)  0.2143*screensize(4)]); 
    uicontrol('Style','PushButton','Position',[10 10 100 40],...
        'String','Get Data','Callback',@GetData)
    uicontrol('Style','PushButton','Position',[110 10 100 40],...
        'String','Play Movie','Callback',@PlayMovie)

%     uicontrol('Style','PushButton','Position',[10 90 100 40],...
%         'String','Get Synapses','Callback',@ExtractSynapses2)
%     uicontrol('Style','PushButton','Position',[110 90 100 40],...
%         'String','Get DeltaF','Callback',@GetDeltaF)

    uicontrol('Style','PushButton','Position',[110 50 100 40],...
        'String','Show Spots','Callback',@DisplaySpots)
    uicontrol('Style','PushButton','Position',[10 50 100 40],...
        'String','Get Events','Callback',@DoAnalysis)

%     uicontrol('Style','PushButton','Position',[10 130 100 40],...
%         'String','Get Bleach','Callback',@GetBleach)
%     uicontrol('Style','PushButton','Position',[110 130 100 40],...
%         'String','Get Decay/Ampli','Callback',@GetDecayAmplitude)   
%     uicontrol('Style','PushButton','Position',[10 170 100 40],...
%         'String','Do Final Analysis','Callback',@Final_Analysis)

    % Keep buttons
    uicontrol('Style','ToggleButton','Position',[350 10 50 40],...
        'String','Synaptic','Tag','KeepbuttonS','enable','off')
    uicontrol('Style','ToggleButton','Position',[350 50 50 40],...
        'String','Dendtric','Tag','KeepbuttonD','enable','off')
    uicontrol('Style','ToggleButton','Position',[350 90 50 40],...
        'String','NaN','Tag','Keepbutton','enable','off')
    %     uicontrol('Style','ToggleButton','Position',[350 10 50 40],...
    %         'String','Keep','Tag','Keepbutton','enable','off')

    uicontrol('Style','ToggleButton','Position',[400 10 50 40],...
        'String','Delete','Tag','Deletebutton','enable','off')
    %     uicontrol('Style','ToggleButton','Position',[300 10 50 40],...
    %         'String','Delete','Tag','Deletebutton','enable','off')
    
    
    % Decay Categorie
    h = uibuttongroup('visible','on','Position',[0.47 .05 .13 .55],...
        'Backgroundcolor',[0.8 0.8 0.8],'title','Decays');
    uicontrol('Style','radiobutton','Position',[5 5 60 20],...
        'String','Type 1','Tag','Type1button','parent',h,'enable','off')
    uicontrol('Style','radiobutton','Position',[5 30 60 20],...
        'String','Type 2','Tag','Type2button','parent',h,'enable','off')
     uicontrol('Style','radiobutton','Position',[5 55 60 20],...
        'String','Type 3','Tag','Type3button','parent',h,'enable','off')
    

    end


    function GetData(hexo,eventdata)

    % Get Name,Spots,MaxProjection
    [StreamName,dirname] = uigetfile('Stream_GluA1-SEP_001.tif','tif');
    cd(dirname)
    [RedName] = uigetfile('Homer-DsRed_001.tif','tif');
    
    load([StreamName(1:end-4) '_Spots_wDeltaF.mat'])

    % Get I (Stream.tif)
    info = imfinfo(StreamName);
    h = waitbar(0,'Please wait reading Stream...');
    for i=1:numel(info)
        I(i).data = imread(StreamName,i);
        waitbar(i/numel(info))
    end
    close(h)
    I_Max = f_StackProjections(I,'Maximum');

%     if exist(RedName)==1
        % Get  Red Projection
        RedChannelName = RedName;
        %         RedChannelName = 'Homer-DsRed_001.tif';
        %         RedChannelName = uigetfile('SEP_1_Norm_8bits_Homer.tif','tif');
        info = imfinfo(RedChannelName);
        h = waitbar(0,'Please wait reading Homer-DsRed_001_Projection.tif...');
        for i=1:numel(info)
            Red(i).data = imread(RedChannelName,i);
            waitbar(i/numel(info))
        end
        close(h)
        Synapses = f_StackProjections(Red,'Maximum');

        gdata.Red_Channel = Synapses;
        gdata.Synapses = Synapses;
        figure(1211)
        imshow(gdata.Synapses,[])
%     else
%         gdata.Red_Channel = I_Max;
%         gdata.Synapses = I_Max;
%     end
    % Get Keepers if exist
    p = dir(cd);
    K = dir('*_Keepers.txt');
    if size(K,1) > 0
        Keepers = load(K(1).name);
        gdata.Keepers = Keepers ;
    end
    % Get DeltaF if exist
    K = dir('*_DeltaF.txt');
    if size(K,1) > 0
        DeltaF = load(K(1).name);
        gdata.DeltaF = DeltaF ;
    end
    gdata.Green_Channel = I_Max;
    guidata(hexo,gdata)

    %%%%%%%%%%%%%%%%%%%%%%

%     figure('name','Make Mask and Skeleton')
%     %         [cmap] = makecolormaps(gdata.Synapses, 'Red');
%     %         Homer_RGB=ind2rgb(gdata.Synapses, cmap);
%     [cmap] = makecolormaps(Synapses, 'Green');
%     SEP_RGB=ind2rgb(I_Max,cmap);
%     IRGB = SEP_RGB;% + Homer_RGB;
%     imshow(IRGB)
%     text(10,10,'Overlay Align','color','w')
%     hold on
%     title('Draw Mask & Skeleton')
%     uicontrol('Style','ToggleButton','Position',[5 5 50 50],'String','Done','Value',0,'tag','DoneRegions')
% 
%     MaskImage = dir(['*_Projection_Mask.tif']);
%     % if exist([fname(1:end-4) '_Skel.tif'])==2
%     if size(MaskImage,1)
%         BW_Mask = double(imread(MaskImage.name));
%     else
%         BW_Mask = zeros(size(IRGB,1),size(IRGB,2));
%     end
% 
%     % ISHH Si existe ca marche pas
%     SkelImage = dir(['*_Skel.tif']);
%     % if exist([fname(1:end-4) '_Skel.tif'])==2
%     if size(SkelImage,1)
%         BW_Skel = imread(SkelImage(end).name);
%         IRGB = IRGB+ ind2rgb(BW_Skel,[0 0 0 ; 1 1 0]);
%         BW_Skel = im2bw(BW_Skel,0);
%     else
%         BW_Skel = zeros(size(IRGB,1),size(IRGB,2));
%     end
%     while get(findobj('tag','DoneRegions'),'Value')== 0
%         % Local Translocation
%         imshow(IRGB);
%         text(10,15,'Select ROI Region','color', 'w')
%         h = imfreehand(gca);
%         api = iptgetapi(h);
%         position= api.getPosition();
%         BW = poly2mask (position(:,1), position(:,2), size(IRGB,1), size(IRGB,2));
%         BW_Mask = BW_Mask + BW;
%         BWS = bwmorph(BW,'thin',Inf);
%         BW_Skel = BW_Skel + BWS;
%         BWRGB = ind2rgb(BWS,[0 0 0 ; 1 1 0]);
%         IRGB = BWRGB + IRGB;
%     end
%     imwrite(BW_Mask,'Homer-DsRed_001_Projection_Mask.tif','tif','Compression','none')
%     imwrite(BW_Skel,'Homer-DsRed_001_Projection_Skel.tif','tif','Compression','none')
%     close(findobj('name','Make Mask and Skeleton'))




    gdata.Spots = Spots;
    gdata.I_Max = I_Max;
    gdata.I = I;
    gdata.StreamName = StreamName;
%     gdata.BW_Mask = BW_Mask;
%     gdata.BW_Skel = BW_Skel;
    guidata(hexo,gdata)

    PlayMovie(hexo)
    DisplaySpots(hexo)
    end

    function DoAnalysis(hexo,eventdata)
    gdata = guidata(hexo);
    Spots2 = gdata.Spots2;
    set(findobj('Tag','Keepbutton'),'enable','on')
    set(findobj('Tag','KeepbuttonS'),'enable','on')
    set(findobj('Tag','KeepbuttonD'),'enable','on')
    set(findobj('Tag','Deletebutton'),'enable','on')
    
    set(findobj('Tag','Type1button'),'enable','on')
    set(findobj('Tag','Type2button'),'enable','on')
    set(findobj('Tag','Type3button'),'enable','on')
    
    Keepers = [];
    IdedDeltaF = [];
    k=1;
    i=1;
    refitValue = 0;
    refitPeakValue =0;
    resetend=0;
    while i <= size(Spots2,1)
%     for i = 1:size(Spots2,1)
        Spots2(i,:)
         if refitValue == 0 || resetend==1;
            % Get DeltaF first approx.
            framebefore = 5;
            frameafter = 10;
            if Spots2(i,3)+frameafter > numel(gdata.I)
                DeltaFrame = (Spots2(i,3)-framebefore:numel(gdata.I))';
            elseif Spots2(i,3)-framebefore <= 0
                DeltaFrame = (1:Spots2(i,3)+frameafter)';
            else
                DeltaFrame = (Spots2(i,3)-framebefore:Spots2(i,3)+frameafter)';
            end
            %   Ided(:,StartFrame:EndFrame) = [x y frame]
            Ided =  [ones(size(DeltaFrame,1),1).* Spots2(i,1)  ones(size(DeltaFrame,1),1).*Spots2(i,2) DeltaFrame];
            % DeltaF
            [DeltaF]  = f_DeltaF_SquareD(Ided, gdata.I,2, 3);
            Ided = [Ided  DeltaF];

            % trouve l'évènement d'insertion le plus grand
            variation = DeltaF-circshift(DeltaF,1);
            peak = find(variation==max(variation));

            %Get Full DeltaF
            framebefore = 15;
            if Ided(peak,3)-framebefore <= 0
                DeltaFrame = (1:numel(gdata.I))';
            else
                DeltaFrame = (Ided(peak,3)-framebefore:numel(gdata.I))';
            end
            %   Ided(:,StartFrame:EndFrame) = [x y frame]
            Ided2 =  [ones(size(DeltaFrame,1),1).* Ided(1,1)  ones(size(DeltaFrame,1),1).*Ided(1,2) DeltaFrame ];
            [DeltaF2]  = f_DeltaF_SquareD(Ided2, gdata.I,2,3);
            Ided2 = [Ided2  DeltaF2];
            
            if resetend==0
                Ided2 = [Ided2 bwlabel(Ided2(1:end,4)>0)];
                eventnum = Ided2(Ided2(:,4) == max(Ided2(:,4)),5);
                zeroframe =  Ided2(Ided2(:,5) == eventnum(1)+1,3); 
                if size(zeroframe,1) > 0
                    Ided2 = Ided2(Ided2(:,3)<=zeroframe(end),:);
                end
            end
            if refitValue > 0 
                Ided2 = Ided2(Ided2(:,3) <= refitValue,:);
            end
         else
            Ided2 = Ided2(Ided2(:,3) <= refitValue,:);
         end
        
        %Fit for decay DeltaF
        variation = Ided2(:,4)-circshift(Ided2(:,4),1);
        if  refitPeakValue==0
            peak = find(variation==max(variation));
            DeltaF3 = Ided2(peak:end,4);
        else
            peak = find(Ided2(:,3)==refitPeakValue);
            DeltaF3 = Ided2(peak:end,4);
%             Ided2 = Ided2(peak-framebefore:end,:);
        end
        rr =0;
        exp2 = 1;
        try
            fun = @(c,x) c(1)*exp(-x/c(2))+c(3)*exp(-x/c(4));
            b0 = [max(DeltaF3) 1 max(DeltaF3) 1];
            X = 1:size(DeltaF3,1);
            [b,r] = nlinfit(X',DeltaF3,fun,b0);
            fitted = fun(b,X); 
            ssresid = sum(r.^2);
            sstotal = sum(DeltaF3);
            ssreg =sstotal- ssresid;
            rr = ssreg/sstotal;
        catch ME
            exp2 = 0;
            disp('Double Exponential failed')
        end
        try
        if exp2==0
            fun = @(c,x) c(1)*exp(-x/c(2));
            b0 = [max(DeltaF3) 1];
            X = 1:size(DeltaF3,1);
            [b,r] = nlinfit(X',DeltaF3,fun,b0);
            b(3)=0;b(4)=0;
            fitted = fun(b,X);
            ssresid = sum(r.^2); %http://office.microsoft.com/en-ca/excel-help/linest-function-HP010069838.aspx
            sstotal = sum(DeltaF3);
            ssreg =sstotal- ssresid;
            rr = ssreg/sstotal;
        end
        catch ME
            disp('Single Exponential failed too')
            fitted = zeros(1,size(Ided2(:,3),1)-size([zeros(1,peak-1)],2));
            b(1)=0; b(2)=0; b(3)=0; b(4)=0;
        end
        
        if size(DeltaF3,1) < 10
            Amplitude = max(DeltaF3(1:size(DeltaF3,1)));
        else
            Amplitude = max(DeltaF3(1:10));
        end
%             DeltaF = IdedDeltaF(IdedDeltaF(:,5)==i,:);
            % Keeper = [x y frame particule localisation Amplitude Amp1 tau1 Amp2 tau2]
%             Keepers = [Keepers ;DeltaF(peak,1:3) DeltaF(peak,5:6)  Amplitude b(1) b(2) b(3) b(4) ];
% % 
        if size(findobj('name','DeltaF Event'),1)>0
            close(findobj('name','DeltaF Event'))
        end
        hDeltaF = figure('Name','DeltaF Event');
        hold on      
%         set(hDeltaF,'Position',[1168 92 512 420])
        screensize = get( 0, 'Screensize' );
        set(hDeltaF,'Position',[0.7*screensize(3) 0.0857*screensize(4) 0.3333*screensize(3)  0.4*screensize(4)]); 
        
        uicontrol('Style','ToggleButton','Position',[5 5 50 50],'String','Refitend','Value',0,'tag','RefitEndTag')
        uicontrol('Style','ToggleButton','Position',[5 105 50 50],'String','Refitpeak','Value',0,'tag','RefitPeakTag')
        uicontrol('Style','ToggleButton','Position',[5 55 50 50],'String','ResetEnd','Value',0,'tag','ResetEndTag')
        uicontrol('Style','ToggleButton','Position',[5 155 50 50],'String','Back','Value',0,'tag','BackTag')
% %        
        plot(Ided2(:,3),Ided2(:,4),'-k')
        plot(Ided2(:,3),[zeros(1,peak-1) fitted],'r')
        l= ['Amplitude=' num2str(Amplitude) '\nAmp1=' num2str(b(1)) ' tau1=' num2str(b(2)) '\nAmp2=' num2str(b(3)) ' tau2=' num2str(b(4))...
            '\nRatio=' num2str(b(2)/b(4)) '\nR2=' num2str(rr)];
        hlegend = legend(sprintf(l));
        if b(2) > 1000 || b(4) > 1000  
            set(hlegend, 'Color', [1 0 0])
        elseif b(2) < 0 || b(4) < 0
            set(hlegend, 'Color', [1 1 0])
        else
            set(hlegend, 'Color', 'None')
        end
        axis([Ided2(1,3) Ided2(end,3) -0.5 2])
        hold off
        
        % Suggest Decay type
        Ratio=round(b(2)/b(4));
        if Ratio ==1 || Ratio == Inf
            set(findobj('Tag','Type1button'),'Value',1)
        else
            set(findobj('Tag','Type2button'),'Value',1)
        end
        
        
        %% display evnets in movie
        squareD=10;
        Roi= [];
        p=1;
        for time=min(Ided2(:,3)):max(Ided2(:,3))
            Roi(p).data = f_cut_square_on_Image(Ided2(Ided2(:,3)==time,1:2), gdata.I(time).data, squareD);
            Roi(p).data = imresize(Roi(p).data,2);
            p=p+1;
        end
        f_StackDisplay('Event',Roi,'data') 

        % Montage SEP
        Montage = f_MakeMontage(Ided, gdata.I,20);
        Montage = imresize(Montage,2);
        if size(findobj('name','Montage Event'),1)>0
            close(findobj('name','Montage Event'))
        end
        hmontage = figure('Name','Montage Event');

%         set(hmontage,'Position',[200 768 1476 168])
        set(hmontage,'Position',[0.1190*screensize(3) 0.73*screensize(4) 0.8786*screensize(3)  0.16*screensize(4)]); 
        imshow(Montage,[],'initialMagnification',200)
        text(10,10,[num2str(i) '/' num2str(size(Spots2,1))],'color','w')

        % Overlay Synapses
        if size(findobj('name','Synapses Montage'),1)>0
            close(findobj('name','Synapses Montage'))
        end
        hmontage2 = figure('Name','Synapses Montage');
%         set(hmontage2,'Position',[34   318   488   410])
        set(hmontage2,'Position',[0.0202*screensize(3) 0.3029*screensize(4) 0.2905*screensize(3)  0.16*screensize(4)]); 
        RoiHomer = cut_square_on_Image(Ided, gdata.Synapses, 40);
        % Make subtraction from two first images
        id = Ided(Ided(:,4) == max(Ided(:,4)),1:3);
        RoiSEP = cut_square_on_Image(id, gdata.I(id(3)).data, 40);
        RoiSEP0 = cut_square_on_Image(Ided, gdata.I(Ided(1,3)).data, 40);
        RoiSEP1 = cut_square_on_Image(Ided, gdata.I(Ided(2,3)).data, 40);
        RoiSEP0 = imadd(RoiSEP0,RoiSEP1)/2;
        RoiSEP=RoiSEP-RoiSEP0;
        Homer_RGB=ind2rgb(RoiHomer, makecolormaps(RoiHomer,'Red'));
        SEP_RGB=ind2rgb(RoiSEP,makecolormaps(RoiSEP, 'Green'));
        imshow(Homer_RGB+SEP_RGB,'InitialMagnification',400)

        % Set events on player
        set(findobj('tag','Movieeditslider'),'String',num2str(Spots2(i,3)))
        set(findobj('tag','Movieslider'),'Value',Spots2(i,3))
        figure(findobj('name','Movie Player'))
        hold on
        plot(Spots2(i,1),Spots2(i,2),'sr','Markersize',6,'tag','eventPlot')
        hold off

        while sum([get(findobj('Tag','Keepbutton'),'Value') get(findobj('Tag','Deletebutton'),'Value') get(findobj('Tag','KeepbuttonS'),'Value') get(findobj('Tag','KeepbuttonD'),'Value')...
            get(findobj('tag','RefitEndTag'),'Value') get(findobj('tag','RefitPeakTag'),'Value')  get(findobj('tag','ResetEndTag'),'Value') get(findobj('tag','BackTag'),'Value')]) == 0
            waitforbuttonpress;
            pause(0.2)
        end
        
        if get(findobj('tag','RefitEndTag'),'Value')==1
            [x]=ginput(1);
           refitValue = floor(x(1));
           set(findobj('tag','RefitEndTag','Value',0))
        end
        
        if get(findobj('tag','RefitPeakTag'),'Value')==1
            [x]=ginput(1);
            refitPeakValue = floor(x(1));
            set(findobj('tag','RefitPeakTag','Value',0))
        end
        
        if get(findobj('tag','ResetEndTag'),'Value')==1
            refitValue =0;
            refitPeakValue = 0;
            resetend=1;
            set(findobj('tag','ResetEndTag','Value',0))
        end
        
        if get(findobj('tag','BackTag'),'Value')==1
            resetend=0;
            refitValue =0;
            refitPeakValue =0;
            if i>1
                i=i-1;
            end
            set(findobj('tag','BackTag','Value',0))
        end
        
        if get(findobj('Tag','Keepbutton'),'Value')==1 || get(findobj('Tag','KeepbuttonS'),'Value')==1 || get(findobj('Tag','KeepbuttonD'),'Value')==1
            % Get Localisation
            if get(findobj('Tag','Keepbutton'),'Value')
                localization = 0; %'NaN'
            elseif get(findobj('Tag','KeepbuttonS'),'Value')
                localization = 1; % Synaptic
            elseif get(findobj('Tag','KeepbuttonD'),'Value')
                localization = 2; %Dendritique
            end
            %Get Deacy type
            if get(findobj('Tag','Type1button'),'Value')
                Decaytype = 1; %'NaN'
            elseif get(findobj('Tag','Type2button'),'Value')
                Decaytype = 2; % Synaptic
            elseif get(findobj('Tag','Type3button'),'Value')
                Decaytype = 3; %Dendritique
            end
            
            % Keeper = [x y frame particule localisation Decaytype Amplitude Amp1 tau1 Amp2 tau2]
            Keepers = [Keepers ;Ided2(peak,1:3) k localization Decaytype Amplitude b(1) b(2) b(3) b(4) ];
            IdedDeltaF = [IdedDeltaF ; Ided2(:,1:3) k*ones(size(Ided2,1),1) Ided2(:,4) localization*ones(size(Ided2,1),1)];
            k=k+1;
            delete(findobj('tag','eventPlot'))
            figure(findobj('name','Movie Player'))
            hold on
            plot(Spots2(i,1),Spots2(i,2),'sy','Markersize',6)
            text(Spots2(i,1)+5,Spots2(i,2)+5,['E ' num2str(i)  'F ' num2str(Spots2(i,3))],'color','r','fontsize',6)
            hold off
            % save images
            p=cd;
            a=dir('Montages');
            if size(a,1)==0
                mkdir('Montages')
            end
            cd([p '\Montages'])
            imwrite(uint16((Montage*65335)/max(Montage(:))),['Montage ' num2str(i) '.tif'],'Compression','none')
            imwrite(Homer_RGB+SEP_RGB,['Overlay ' num2str(i) '.tif'],'Compression','none')
            cd(p)
            i=i+1;
            refitValue = 0;
            refitPeakValue=0;
            resetend=0;
        end

        if get(findobj('Tag','Deletebutton'),'Value')==1
            delete(findobj('tag','eventPlot'))
            i=i+1;
            refitValue = 0;
            refitPeakValue=0;
        end
        set(findobj('Tag','Keepbutton'),'Value',0)
        set(findobj('Tag','KeepbuttonS'),'Value',0)
        set(findobj('Tag','KeepbuttonD'),'Value',0)
        set(findobj('Tag','Deletebutton'),'Value',0)
    end

%           Keepers = [x y frame amplitude]
    if isempty(Keepers)==1
        Keepers = [];
        save([gdata.StreamName(1:end-4) '_Keepers.txt'],'Keepers','-ascii')
        Freq = 0;
        save([gdata.StreamName(1:end-4) '_Freq.txt'], 'Freq', '-ascii', '-double')
        IdedDeltaF = [];
        save([gdata.StreamName(1:end-4) '_DeltaF.txt'], 'IdedDeltaF', '-ascii', '-double')
    else
%         Keepers = removeRecEvents2(Keepers,10,10);
        save([gdata.StreamName(1:end-4) '_Keepers.txt'],'Keepers','-ascii')
        %Calcule fréquence
        numevents = size(Keepers,1);
%         lenghtD =  sum(gdata.BW_Skel(:)) * 0.16; % in um
        time = (numel(gdata.I)) * 0.1; %in sec
        Freq = (numevents/time)
%         FreqNorm =  (numevents/time)/lenghtD;
        save([gdata.StreamName(1:end-4) '_Freq.txt'], 'Freq', '-ascii', '-double')
%         save([gdata.StreamName(1:end-4) '_FreqNorm.txt'], 'FreqNorm', '-ascii', '-double')
        save([gdata.StreamName(1:end-4) '_DeltaF.txt'], 'IdedDeltaF', '-ascii', '-double')
    end

    gdata.DeltaF = IdedDeltaF;
    gdata.Keepers = Keepers;
    gdata.Freq = Freq;
    guidata(hexo,gdata)

    try close(4); close(5); close(6);
    catch ME
    end


    
    GetDeltaF(hexo)
    end

    function DisplaySpots(hexo,eventdata)
    
    gdata = guidata(hexo);
    Spots = gdata.Spots;
    if size(findobj('name','Display Spots'),1)>0
        close(findobj('name','Display Spots'))
    end

    hspots = figure('name','Display Spots');
%     set(hspots,'Position',[35 280 560 420])
    screensize = get( 0, 'Screensize' );
    set(hspots,'Position',[0.0202*screensize(3) 0.2667*screensize(4) 0.3333*screensize(3)  0.4*screensize(4)]); 
    
    uicontrol('Style','text','Position',[10 10 120 20],'String','Max area under curve')
    hslider = uicontrol('Style','slider','Position',[200 10 200 20],...
        'Min',0,'Max',max(Spots(:,4)),'Value',0.65,...
        'tag','sliderspots','Callback',@DisplayEvents);
    uicontrol('Style','edit','Position',[135 10 50 20],...
        'tag','editslider','String',2)

    uicontrol('Style','text','Position',[10 40 80 30],'String','Distance Threshold')
    uicontrol('Style','Edit','Position',[100 40 30 30],...
        'String','10','tag','Tag_Dist_th','Callback',@DisplayEvents)
    uicontrol('Style','text','Position',[10 80 80 30],'String','Time Threshold')
    uicontrol('Style','Edit','Position',[100 80 30 30],...
        'String','10','tag','Tag_Time_th','Callback',@DisplayEvents)
    %         uicontrol('Style','ToggleButton','Position',[10 35 40 40],...
    %             'String','Done','tag','DoneSpotsDisplay','Callback',@DisplayEvents);

        function DisplayEvents(hspots,eventdata)
            gdata = guidata(hexo);
            Spots = gdata.Spots;
            imshow(gdata.I_Max,[])
            Spots2 = sortrows(Spots,-4);
            set(findobj('tag','editslider'),'String',num2str(get(findobj('tag','sliderspots'),'Value')))
            Spots2 = Spots(Spots(:,4)>get(findobj('tag','sliderspots'),'Value'),:);
            Spots2 = removeRecEvents2(Spots2,str2num(get(findobj('tag','Tag_Dist_th'),'String')),str2num(get(findobj('tag','Tag_Time_th'),'String')));
            hold on
            if size(Spots2,1)<500
                for i=1:size(Spots2,1)
                    plot(Spots2(i,1),Spots2(i,2),'sr','Markersize',6)
                    %text(Spots2(i,1)+5,Spots2(i,2)+5,['Event ' num2str(i)  'Frame' num2str(Spots2(i,3))],'color','r')
                    text(10,10,['Selected ' num2str(size(Spots2,1)) ' Spots'],'color','w')
                end
                hold off
            end
            gdata.Spots2 = Spots2;
            guidata(hexo,gdata)
        end
    end

    function PlayMovie(hexo,eventdata)
    gdata = guidata(hexo);
    if size(findobj('name','Movie Player'),1)>0
        close(findobj('name','Movie Player'))
    end

    hplayer=figure('name','Movie Player');
%     set(hplayer,'Position',[600 90 560 420]);
    screensize = get( 0, 'Screensize' );
    set(hplayer,'Position',[0.3571*screensize(3) 0.0857*screensize(4) 0.3333*screensize(3)  0.4*screensize(4)]); 
    
    set(hplayer,'Toolbar','Figure')
    uicontrol('Style','ToggleButton','Position',[10 10 60 40],...
        'String','Play','Callback',@Player,'tag','PlayTag');
    uicontrol('Style','edit','Position',[135 10 50 20],...
        'tag','Movieeditslider','String',1)
    hslider = uicontrol('Style','slider','Position',[200 10 200 20],...
        'Min',1,'Max',numel(gdata.I),'Value',1,'Sliderstep', [0.001 0.005],...
        'tag','Movieslider','Callback',@SliderPlayer);
    uicontrol('Style','checkbox','Position',[425 10 80 20],...
        'tag','AutoScale','String','AutoScale','Value',0)
    uicontrol('Style','checkbox','Position',[10 60 60 20],...
        'tag','OverlayTag','String','Overlay','Value',0)
    [cmap] = makecolormaps(gdata.Synapses, 'Red');
    Homer_RGB=ind2rgb(gdata.Synapses, cmap);
    [cmap] = makecolormaps(gdata.I(1).data, 'Green');


    axesPlayer = axes('parent',findobj('name','Movie Player'),'tag','MovieAxes');
    image = imshow(gdata.I(1).data,[],'InitialMagnification', 'fit', 'parent', axesPlayer);
    caxis auto




        function Player(hplayer,eventdata)
            i=0;
            set(findobj('tag','PlayTag'),'String','Stop')
            while get(findobj('tag','PlayTag'),'Value')==1
                i=i+1;
                if i>numel(gdata.I)
                    i=1;
                end
                if get(findobj('tag','AutoScale'),'Value') ==0
                    caxis manual
                else
                    caxis auto
                end
                if get(findobj('tag','OverlayTag'),'Value')==0
                    set(image,'CData',gdata.I(i).data);
                else
                    SEP_RGB=ind2rgb(gdata.I(i).data, cmap);
                    set(image,'CData',SEP_RGB+Homer_RGB);
                end

                drawnow
                set(findobj('tag','Movieslider'),'Value',i)
                set(findobj('tag','Movieeditslider'),'String',num2str(i))
            end
            set(findobj('tag','PlayTag'),'String','Play')
        end

        function SliderPlayer(hplayer,eventdata)
            i = round(get(findobj('tag','Movieslider'),'Value'));
            set(findobj('tag','Movieeditslider'),'String',num2str(i))
            set(image,'CData',gdata.I(i).data);
        end

    end

    function ExtractSynapses2(hexo,~)
    gdata = guidata(hexo);
    Keepers = gdata.Keepers;
    I = gdata.Red_Channel.data;
    I_RGB = ind2rgb(I,makecolormaps(I,'Trans'));
    [BW]= f_feature_extract(I,7,2);
    BW_RGB = ind2rgb(BW,[0 0 0 ; 0.5 0 0]);
    figure; imshow(I_RGB + BW_RGB);
    %         hold on
    %%
    BW_Dilate = imdilate(BW ,ones(4));

    I_peri = logical(BW_Dilate-BW);
    I_peri_RGB = ind2rgb(I_peri,[0 0 0 ; 0.5 0.5 0]);
    figure; imshow(I_RGB + BW_RGB+I_peri_RGB);
    %
    hold on
    for i=1:size(Keepers,1)
        plot(Keepers(i,1),Keepers(i,2),'.g')
    end
    hold off

    % Keeps tracks under BW_Dilate
    [r,c]= find(BW_Dilate>0);
    xy = [c r];
    [tf] = ismember(Keepers(:,1:2),xy,'rows');
    Keepers_Synapses= Keepers(tf,:);


    end

    function ExtractSynapses(hexo,~)
    gdata = guidata(hexo);

    f = figure('name','Segmentation','numbertitle','on','Position',[100 100 500 150],'tag', 'Segemntaion parameters');
    set(findobj('tag', 'Segemntaion parameters'),'Position',[82  100   500   150])

    % Window W uicontrols
    Adaph   =uicontrol('Style', 'text', 'String', 'Adaptive threshold parameters',...
        'Position', [10 125 160 20], 'parent', f);
    slideW  =uicontrol('Style','slider','Position', [75,90,100,25],...
        'Value',12,'Min',1,'Max',30,'SliderStep',[0.05 0.5],...
        'Callback',@Segment_synapses , 'parent', f, 'BusyAction', 'cancel');
    Wtext   =uicontrol('Style', 'text', 'String', 'W',...
        'Position', [10 90 20 25], 'parent', f);
    Wh      =uicontrol('Style', 'text', 'String',...
        round(get(slideW,'Value')),'Position', [40 90 20 25], 'parent', f);

    % kSD uicontrols
    slidekSD =uicontrol('Style','slider','Position', [75,50,100,25],...
        'Value',5,'Min',0,'Max',20,'SliderStep',[0.01 0.25],...
        'Callback',@Segment_synapses, 'parent', f, 'BusyAction', 'cancel');
    kSDtext =uicontrol('Style', 'text', 'String', 'kSD',...
        'Position', [10 50 20 25], 'parent', f);
    kSDh    =uicontrol('Style', 'text', 'String', round(get(slidekSD,'Value')),...
        'Position', [40 50 25 25], 'parent', f);

    %clusters size uicontrols
    sizeh   =uicontrol('Style', 'text', 'String', 'Cluster properties parameters',...
        'Position', [200 125 160 20], 'parent', f);
    mintext   =uicontrol('Style', 'text', 'String', 'Min Cluster Size',...
        'Position', [200 95 85  20], 'parent', f);
    Minh      =uicontrol('Style', 'edit', 'String',...
        4,'Position', [290 95 30 20], 'parent', f);

    Maxtext   =uicontrol('Style', 'text', 'String', 'Max Cluster Size',...
        'Position', [200 65 85 20], 'parent', f);
    Maxh      =uicontrol('Style', 'edit', 'String',...
        100,'Position', [290 65 30 20], 'parent', f);

    Ecctext   =uicontrol('Style', 'text', 'String', 'Eccentricity',...
        'Position', [200 35 85 20], 'parent', f);
    Ecch      =uicontrol('Style', 'edit', 'String',...
        0.97,'Position', [290 35 30 20], 'parent', f);

    uicontrol('Style', 'pushbutton', 'String', 'Imcontrast',...
        'Position', [95 5 75 20], 'parent', f,'Callback',@ContrastSynapses );

    uicontrol('Style', 'pushbutton', 'String', 'Reset Image',...
        'Position', [10 5 75 20], 'parent', f,'Callback',@ResetImage);

    %remove image overlay
    Removetxtim = uicontrol('Style', 'checkbox', 'String', 'remove Overlay','Value',1,'Position', [350 10 120 25], 'parent', f,'Callback',@Segment_synapses);

    % Done button
    Runh = uicontrol('Style', 'pushbutton','String','Run',...
        'Position', [350 70 50 50], 'Callback',@Segment_synapses, 'parent', f);
    Doneh = uicontrol('Style', 'togglebutton','String','Done',...
        'Position', [425 70 50 50], 'Callback',@Segment_synapses, 'parent', f);
    %     Batchh = uicontrol('Style', 'togglebutton','String','Send to BatchFile',...
    %         'Position', [350 40 120 25], 'Callback',@Segment_synapses, 'parent', f);

        function Segment_synapses(f,~)

            W = round(get(slideW,'Value'));
            set( Wh, 'String',W)
            kSD = (get(slidekSD,'Value'));
            set( kSDh, 'String',kSD)
            BW1 = f_Adaptive_Threshold_Morpho(gdata.Synapses,W,kSD,str2num(get(Minh,'String')),str2num(get(Maxh,'String')),str2num(get(Ecch,'String')));
            figure(10)
            imshow(BW1,[])
            text(10,10,'Adaptive threshold It','color','w');
            set(10,'Position',[841   219   833   754])
            J = ind2rgb(gdata.Synapses,gray(2^16)); BW1RGB = ind2rgb(BW1,[0 0 0 ; 1 0 0]);
            %remove image overlay
            if (get(Removetxtim,'Value') == get(Removetxtim,'Max'))
                figure(11)
                imshow(BW1RGB+J,'InitialMagnification',100)
                text(10,10,'Overlay It','color','w');
            else
                % Checkbox is not checked-take approriate action
                figure(11)
                imshow(J,'InitialMagnification',100)
                text(10,10,'Adaptive threshold It','color','w');
            end
            set(11,'Position',[1 219 833 754])
            %Done button
            try
                button_state = get(Doneh,'Value');
                if button_state == get(Doneh,'Max')

                    gdata.Synapses_BW = BW1;
                    guidata(hexo,gdata)
                    savename = strcat('Synapses_BW.tif');
                    imwrite(BW1,savename,'tif','Compression','none')
                    close(findobj('tag', 'Segemntaion parameters'))
                    close(11)
                    close(10)
                end
            catch ME
            end
        end

        function ContrastSynapses(f,eventdata)
            figure(10)
            imshow(gdata.Synapses,[])
            imch= imcontrast(10);
            waitfor (imch)
            Synapses = getimage(10);
            close(10)
            gdata.Synapses = Synapses;
            guidata(hexo,gdata);
            imwrite(gdata.Synapses,'Synapses.tif','tif', 'Compression','none')
            Segment_synapses
        end

        function ResetImage(f,eventdata)

            info = imfinfo(QD.SynName);
            for i=1:numel(info)
                QD.ISyn(i).data = imread(QD.SynName,i);
            end

            %addition des images
            I_sum=zeros(size(QD.ISyn(1).data,1),size(QD.ISyn(1).data,2));
            for l=1:numel(QD.ISyn)
                I2= double(QD.ISyn(l).data);
                I_sum = I_sum + I2;
            end
            QD.Synapses_sum = uint16((I_sum / max(I_sum(:)))*2^16);
            guidata (findobj('tag','MainMenu'),QD);
            imwrite(QD.Synapses_sum,strcat([QD.SynName(1:end-4) '_sum.tif']),'tif', 'Compression','none')
            Segment_synapses

        end

    end

    function GetBleach(hexo,~)
    gdata = guidata(hexo);
    %         BW_Mask= gdata.BW_Mask;
    BW_Mask =  imread('Homer-DsRed_001_Projection_Mask.tif');
    BW_Mask = double(BW_Mask ./ max(BW_Mask(:)));

    Synapses = gdata.Synapses;
    Green_Channel = gdata.Green_Channel;
    Red_Channel = gdata.Red_Channel;
    I_Mean =[];

    figure('name','Overlay')
    imshow(BW_Mask,[])

    for i=1:numel(Green_Channel)
        I(i).data = double(Green_Channel(i).data) .*  BW_Mask ;
        I_Mean = [ I_Mean; mean(I(i).data(I(i).data > 0))];
    end
    % Delta F
    F0 = mean(I_Mean(1:2));
    DF  = ((I_Mean-F0)./ F0)*100;

    figure(1112)
    plot(DF)


    end

    function GetDeltaF(hexo,~)
    %%
    gdata = guidata(hexo);
    Keepers=gdata.Keepers;
    DeltaF = gdata.DeltaF;
    
    IdedDeltaF = [];
    
    figure(1111)
    imshow(gdata.Red_Channel,[])
    
    if size(DeltaF,1)>0
        
        hold on
        for i=1:max(DeltaF(:,4))
            Ided = DeltaF(DeltaF(:,4)==i,:); pos=Ided(1,1:2);
            plot(pos(1),pos(2),'sr','Markersize',6)
            text(pos(1),pos(2),num2str(i),'color','y')
            
            try
                Ided2 = Ided;
                % Facteur de taille en x.
                xscale= 0.5;
                yscale = 50;
                Pos = pos(1):xscale:((pos(1)+size(Ided2,1)-1));
                if Ided(1,6) == 1 % Synaptic
                    line(Pos(1:size(Ided2,1)),(pos(2)-yscale*Ided2(:,5)),'color','r')
                elseif Ided(1,6) == 2 %Dendritic
                    line(Pos(1:size(Ided2,1)),(pos(2)-yscale*Ided2(:,5)),'color','g')
                elseif Ided(1,6) == 0 % NaN
                    line(Pos(1:size(Ided2,1)),(pos(2)-yscale*Ided2(:,5)),'color','y')
                end
                text(10,10,'Synaptic Events','color','r')
                text(10,25,'Dendritic Events','color','y')
            catch ME
                disp(['Missed trace ' num2str(i)])
            end
        end
        
        % Scale Bar
        line(10,yscale*(0:0.01:1)+475,'color','y')
        text(12,475,'DeltaF (100%)','color','y')
        line(xscale*(0:0.1:100)+10,475+yscale,'color','y')
        text(100*xscale+12,yscale+475,'Time (10s)','color','y')
        text(10,10,cd,'color','k','fontsize',4);
    end
    

    
    saveas(1111,'Events on Homer','pdf')
    
    guidata(hexo,gdata)

    end

    function GetDecayAmplitude(hexo,~)

    gdata = guidata(hexo);
    IdedDeltaF = gdata.DeltaF;
     % Keeper = [x y frame particule localisation Amplitude Amp1 tau1 Amp2 tau2]
    Keepers = gdata.Keepers;
    for i=1:max(Keepers(:,4))
        Amplitude=Keepers(i,6);
        b=Keepers(i,7:10);
        peak= Keepers(i,3);

        deltaF = IdedDeltaF(IdedDeltaF(:,4)==i,[3 5] );
        if  findobj('name','Fitted Event')>0
            close( findobj('name','Fitted Event'))
        end
         figure9=figure(9);set(9,'name','Fitted Event');
         axes1 = axes('Parent',figure9);
         
        plot(deltaF(:,1),deltaF(:,2),'-k')
%         try
            fun = @(c,x) c(1)*exp(-x/c(2))+c(3)*exp(-x/c(4));
% %             fun = @(c,x) c(1)*exp(-x/c(2));
%             b0 = [max(deltaF3) 500 0.75*max(deltaF3) 10];
            X = 1:size(deltaF(find(deltaF(:,1)==peak):end,1),1);
            fitted = fun(b,X);
            hold on
            plot(deltaF(:,1),[zeros(1,size(deltaF(:,1),1)-size(fitted,2)) fitted],'r')
            
            nume = round(max([b(2) b(4)]));
            deno = round(min([b(2) b(4)]));
            
            Ratio = nume / deno;
            l= ['Amplitude=' num2str(Amplitude) '\nAmp1=' num2str(b(1)) ' tau1=' num2str(b(2)) '\nAmp2=' num2str(b(3)) ' tau2=' num2str(b(4)) '\nRatio=' num2str(Ratio)];
            legend(sprintf(l)); 
            p=cd;
            text(deltaF(1,1),-0.2,[p 'Keeper ' num2str(Keepers(i,4))],'FontSize',4)
            hold off
            ylim(axes1,[-0.25 2]);
            xlim(axes1,[0 1000]);
%             axis([min(deltaF(:,1)) max(deltaF(:,1)) -0.25 2])
           
           
% 
%             legend(['Amplitude' num2str(Amplitude) ' Max = ' num2str(b(1)) ' ' num2str(b(3)) ' Tau% = ' num2str(b(2)) ' ' num2str(b(4))])
%         catch ME
%             disp('argh!!!!')
%         end
            
            % save images       
            a=dir('DeltaF');
            if size(a,1)==0
                mkdir('DeltaF')
            end
            cd([p '\DeltaF'])
            saveas(findobj('name','Fitted Event'),['DeltaF_Event_' num2str(Keepers(i,4))],'pdf')
            cd(p)

        pause(0.1)
    end

    gdata.Keepers = Keepers;
    guidata(hexo,gdata)

%     save([gdata.StreamName(1:end-4) '_Keepers.txt'],'Keepers','-ascii')
    end

    function Final_Analysis(hexo,~)
    gdata = guidata(hexo);
    K = dir('*_Keepers.txt');
    if size(K,1) > 0
        Keepers = load(K(1).name);
    end
    
    Spots2 = Keepers;
    set(findobj('Tag','Keepbutton'),'enable','on')
    set(findobj('Tag','KeepbuttonS'),'enable','on')
    set(findobj('Tag','KeepbuttonD'),'enable','on')
    set(findobj('Tag','Deletebutton'),'enable','on')
    
    set(findobj('Tag','Type1button'),'enable','on')
    set(findobj('Tag','Type2button'),'enable','on')
    set(findobj('Tag','Type3button'),'enable','on')
    
    Keepers = [];
    IdedDeltaF = [];
    k=1;
    i=1;
    refitValue = 0;
    refitPeakValue =0;
    resetend=0;
    while i <= size(Spots2,1)
%     for i = 1:size(Spots2,1)
        Spots2(i,:)
         if refitValue == 0 || resetend==1;
            % Get DeltaF first approx.
            framebefore = 5;
            frameafter = 10;
            if Spots2(i,3)+frameafter > numel(gdata.I)
                DeltaFrame = (Spots2(i,3)-framebefore:numel(gdata.I))';
            elseif Spots2(i,3)-framebefore <= 0
                DeltaFrame = (1:Spots2(i,3)+frameafter)';
            else
                DeltaFrame = (Spots2(i,3)-framebefore:Spots2(i,3)+frameafter)';
            end
            %   Ided(:,StartFrame:EndFrame) = [x y frame]
            Ided =  [ones(size(DeltaFrame,1),1).* Spots2(i,1)  ones(size(DeltaFrame,1),1).*Spots2(i,2) DeltaFrame];
            % DeltaF
            [DeltaF]  = f_DeltaF_SquareD(Ided, gdata.I,2, 3);
            Ided = [Ided  DeltaF];

            % trouve l'évènement d'insertion le plus grand
            variation = DeltaF-circshift(DeltaF,1);
            peak = find(variation==max(variation));

            %Get Full DeltaF
            framebefore = 15;
            if Ided(peak,3)-framebefore <= 0
                DeltaFrame = (1:numel(gdata.I))';
            else
                DeltaFrame = (Ided(peak,3)-framebefore:numel(gdata.I))';
            end
            %   Ided(:,StartFrame:EndFrame) = [x y frame]
            Ided2 =  [ones(size(DeltaFrame,1),1).* Ided(1,1)  ones(size(DeltaFrame,1),1).*Ided(1,2) DeltaFrame ];
            [DeltaF2]  = f_DeltaF_SquareD(Ided2, gdata.I,2,3);
            Ided2 = [Ided2  DeltaF2];
            
            if resetend==0
                Ided2 = [Ided2 bwlabel(Ided2(1:end,4)>0)];
                eventnum = Ided2(Ided2(:,4) == max(Ided2(:,4)),5);
                zeroframe =  Ided2(Ided2(:,5) == eventnum(1)+1,3); 
                if size(zeroframe,1) > 0
                    Ided2 = Ided2(Ided2(:,3)<=zeroframe(end),:);
                end
            end
            if refitValue > 0 
                Ided2 = Ided2(Ided2(:,3) <= refitValue,:);
            end
         else
            Ided2 = Ided2(Ided2(:,3) <= refitValue,:);
         end
        
        %Fit for decay DeltaF
        variation = Ided2(:,4)-circshift(Ided2(:,4),1);
        if  refitPeakValue==0
            peak = find(variation==max(variation));
            DeltaF3 = Ided2(peak:end,4);
        else
            peak = find(Ided2(:,3)==refitPeakValue);
            DeltaF3 = Ided2(peak:end,4);
%             Ided2 = Ided2(peak-framebefore:end,:);
        end
        rr =0;
        exp2 = 1;
        try
            fun = @(c,x) c(1)*exp(-x/c(2))+c(3)*exp(-x/c(4));
            b0 = [max(DeltaF3) 1 max(DeltaF3) 1];
            X = 1:size(DeltaF3,1);
            [b,r] = nlinfit(X',DeltaF3,fun,b0);
            fitted = fun(b,X); 
            ssresid = sum(r.^2);
            sstotal = sum(DeltaF3);
            ssreg =sstotal- ssresid;
            rr = ssreg/sstotal;
        catch ME
            exp2 = 0;
            disp('Double Exponential failed')
        end
        try
        if exp2==0
            fun = @(c,x) c(1)*exp(-x/c(2));
            b0 = [max(DeltaF3) 1];
            X = 1:size(DeltaF3,1);
            [b,r] = nlinfit(X',DeltaF3,fun,b0);
            b(3)=0;b(4)=0;
            fitted = fun(b,X);
            ssresid = sum(r.^2); %http://office.microsoft.com/en-ca/excel-help/linest-function-HP010069838.aspx
            sstotal = sum(DeltaF3);
            ssreg =sstotal- ssresid;
            rr = ssreg/sstotal;
        end
        catch ME
            disp('Single Exponential failed too')
            fitted = zeros(1,size(Ided2(:,3),1)-size([zeros(1,peak-1)],2));
            b(1)=0; b(2)=0; b(3)=0; b(4)=0;
        end
        
        if size(DeltaF3,1) < 10
            Amplitude = max(DeltaF3(1:size(DeltaF3,1)));
        else
            Amplitude = max(DeltaF3(1:10));
        end
%             DeltaF = IdedDeltaF(IdedDeltaF(:,5)==i,:);
            % Keeper = [x y frame particule localisation Amplitude Amp1 tau1 Amp2 tau2]
%             Keepers = [Keepers ;DeltaF(peak,1:3) DeltaF(peak,5:6)  Amplitude b(1) b(2) b(3) b(4) ];
% % 
        if size(findobj('name','DeltaF Event'),1)>0
            close(findobj('name','DeltaF Event'))
        end
        hDeltaF = figure('Name','DeltaF Event');
        hold on      
        set(hDeltaF,'Position',[1168 92 512 420])
        uicontrol('Style','ToggleButton','Position',[5 5 50 50],'String','Refitend','Value',0,'tag','RefitEndTag')
        uicontrol('Style','ToggleButton','Position',[5 105 50 50],'String','Refitpeak','Value',0,'tag','RefitPeakTag')
        uicontrol('Style','ToggleButton','Position',[5 55 50 50],'String','ResetEnd','Value',0,'tag','ResetEndTag')
        uicontrol('Style','ToggleButton','Position',[5 155 50 50],'String','Back','Value',0,'tag','BackTag')
% %        
        plot(Ided2(:,3),Ided2(:,4),'-k')
        plot(Ided2(:,3),[zeros(1,peak-1) fitted],'r')
        l= ['Amplitude=' num2str(Amplitude) '\nAmp1=' num2str(b(1)) ' tau1=' num2str(b(2)) '\nAmp2=' num2str(b(3)) ' tau2=' num2str(b(4))...
            '\nRatio=' num2str(b(2)/b(4)) '\nR2=' num2str(rr)];
        hlegend = legend(sprintf(l));
        if b(2) > 1000 || b(4) > 1000  
            set(hlegend, 'Color', [1 0 0])
        elseif b(2) < 0 || b(4) < 0
            set(hlegend, 'Color', [1 1 0])
        else
            set(hlegend, 'Color', 'None')
        end
        axis([Ided2(1,3) Ided2(end,3) -0.5 5])
        hold off
        
        % Suggest Decay type
        Ratio=round(b(2)/b(4));
        if Ratio ==1 || Ratio == Inf
            set(findobj('Tag','Type1button'),'Value',1)
        else
            set(findobj('Tag','Type2button'),'Value',1)
        end
        
        
        %% display evnets in movie
        squareD=10;
        Roi= [];
        p=1;
        for time=min(Ided2(:,3)):max(Ided2(:,3))
            Roi(p).data = f_cut_square_on_Image(Ided2(Ided2(:,3)==time,1:2), gdata.I(time).data, squareD);
            Roi(p).data = imresize(Roi(p).data,2);
            p=p+1;
        end
        f_StackDisplay('Event',Roi,'data') 

        % Montage SEP
        Montage = f_MakeMontage(Ided, gdata.I,20);
        Montage = imresize(Montage,2);
        if size(findobj('name','Montage Event'),1)>0
            close(findobj('name','Montage Event'))
        end
        hmontage = figure('Name','Montage Event');
%         set(hmontage,'Position',[200 768 1476 168])
        screensize = get( 0, 'Screensize' );
        set(hmontage,'Position',[0.0179*screensize(3) 0.7314*screensize(4) 0.8786*screensize(3)  0.16*screensize(4)]); 
        imshow(Montage,[],'initialMagnification',200)
        text(10,10,[num2str(i) '/' num2str(size(Spots2,1))],'color','w')

        % Overlay Synapses
        if size(findobj('name','Synapses Montage'),1)>0
            close(findobj('name','Synapses Montage'))
        end
        hmontage2 = figure('Name','Synapses Montage');
        set(hmontage2,'Position',[0.0202   318   488   410])
        set(hmontage,'Position',[0.0179*screensize(3)  0.3029*screensize(4) 0.2905*screensize(3)  0.3905*screensize(4)]); 
        RoiHomer = cut_square_on_Image(Ided, gdata.Synapses, 40);
        % Make subtraction from two first images
        id = Ided(Ided(:,4) == max(Ided(:,4)),1:3);
        RoiSEP = cut_square_on_Image(id, gdata.I(id(3)).data, 40);
        RoiSEP0 = cut_square_on_Image(Ided, gdata.I(Ided(1,3)).data, 40);
        RoiSEP1 = cut_square_on_Image(Ided, gdata.I(Ided(2,3)).data, 40);
        RoiSEP0 = imadd(RoiSEP0,RoiSEP1)/2;
        RoiSEP=RoiSEP-RoiSEP0;
        Homer_RGB=ind2rgb(RoiHomer, makecolormaps(RoiHomer,'Red'));
        SEP_RGB=ind2rgb(RoiSEP,makecolormaps(RoiSEP, 'Green'));
        imshow(Homer_RGB+SEP_RGB,'InitialMagnification',400)

        % Set events on player
        set(findobj('tag','Movieeditslider'),'String',num2str(Spots2(i,3)))
        set(findobj('tag','Movieslider'),'Value',Spots2(i,3))
        figure(findobj('name','Movie Player'))
        hold on
        plot(Spots2(i,1),Spots2(i,2),'sr','Markersize',6,'tag','eventPlot')
        hold off

        while sum([get(findobj('Tag','Keepbutton'),'Value') get(findobj('Tag','Deletebutton'),'Value') get(findobj('Tag','KeepbuttonS'),'Value') get(findobj('Tag','KeepbuttonD'),'Value')...
            get(findobj('tag','RefitEndTag'),'Value') get(findobj('tag','RefitPeakTag'),'Value')  get(findobj('tag','ResetEndTag'),'Value') get(findobj('tag','BackTag'),'Value')]) == 0
            waitforbuttonpress;
            pause(0.2)
        end
        
        if get(findobj('tag','RefitEndTag'),'Value')==1
            [x]=ginput(1);
           refitValue = floor(x(1));
           set(findobj('tag','RefitEndTag','Value',0))
        end
        
        if get(findobj('tag','RefitPeakTag'),'Value')==1
            [x]=ginput(1);
            refitPeakValue = floor(x(1));
            set(findobj('tag','RefitPeakTag','Value',0))
        end
        
        if get(findobj('tag','ResetEndTag'),'Value')==1
            refitValue =0;
            refitPeakValue = 0;
            resetend=1;
            set(findobj('tag','ResetEndTag','Value',0))
        end
        
        if get(findobj('tag','BackTag'),'Value')==1
            resetend=0;
            refitValue =0;
            refitPeakValue =0;
            if i>1
                i=i-1;
            end
            set(findobj('tag','BackTag','Value',0))
        end
        
        if get(findobj('Tag','Keepbutton'),'Value')==1 || get(findobj('Tag','KeepbuttonS'),'Value')==1 || get(findobj('Tag','KeepbuttonD'),'Value')==1
            % Get Localisation
            if get(findobj('Tag','Keepbutton'),'Value')
                localization = 0; %'NaN'
            elseif get(findobj('Tag','KeepbuttonS'),'Value')
                localization = 1; % Synaptic
            elseif get(findobj('Tag','KeepbuttonD'),'Value')
                localization = 2; %Dendritique
            end
            %Get Deacy type
            if get(findobj('Tag','Type1button'),'Value')
                Decaytype = 1; %'NaN'
            elseif get(findobj('Tag','Type2button'),'Value')
                Decaytype = 2; % Synaptic
            elseif get(findobj('Tag','Type3button'),'Value')
                Decaytype = 3; %Dendritique
            end
            
            % Keeper = [x y frame particule localisation Decaytype Amplitude Amp1 tau1 Amp2 tau2]
            Keepers = [Keepers ;Ided2(peak,1:3) k localization Decaytype Amplitude b(1) b(2) b(3) b(4)];
            IdedDeltaF = [IdedDeltaF ; Ided2(:,1:3) k*ones(size(Ided2,1),1) Ided2(:,4) localization*ones(size(Ided2,1),1)];
            k=k+1;
            delete(findobj('tag','eventPlot'))
            figure(findobj('name','Movie Player'))
            hold on
            plot(Spots2(i,1),Spots2(i,2),'sy','Markersize',6)
            text(Spots2(i,1)+5,Spots2(i,2)+5,['E ' num2str(i)  'F ' num2str(Spots2(i,3))],'color','r','fontsize',6)
            hold off
            % save images
            p=cd;
            a=dir('Montages');
            if size(a,1)==0
                mkdir('Montages')
            end
            cd([p '\Montages'])
            saveas(hDeltaF,['Deltaf_Fit ' num2str(i) '.pdf'],'pdf')
            imwrite(uint16((Montage*65335)/max(Montage(:))),['Montage ' num2str(i) '.tif'],'Compression','none')
            imwrite(Homer_RGB+SEP_RGB,['Overlay ' num2str(i) '.tif'],'Compression','none')
            cd(p)
            i=i+1;
            refitValue = 0;
            refitPeakValue=0;
            resetend=0;
        end

        if get(findobj('Tag','Deletebutton'),'Value')==1
            delete(findobj('tag','eventPlot'))
            i=i+1;
            refitValue = 0;
            refitPeakValue=0;
        end
        set(findobj('Tag','Keepbutton'),'Value',0)
        set(findobj('Tag','KeepbuttonS'),'Value',0)
        set(findobj('Tag','KeepbuttonD'),'Value',0)
        set(findobj('Tag','Deletebutton'),'Value',0)
    end

%           Keepers = [x y frame amplitude]
    if isempty(Keepers)==1
        Keepers = [];
        save([gdata.StreamName(1:end-4) '_Keepers_Final.txt'],'Keepers','-ascii')
        Freq = 0;
        save([gdata.StreamName(1:end-4) '_Freq_Final.txt'], 'Freq', '-ascii', '-double')
        IdedDeltaF = [];
        save([gdata.StreamName(1:end-4) '_DeltaF_Final.txt'], 'IdedDeltaF', '-ascii', '-double')
    else
%         Keepers = removeRecEvents2(Keepers,10,10);
        save([gdata.StreamName(1:end-4) '_Keepers_Final.txt'],'Keepers','-ascii')
        %Calcule fréquence
        numevents = size(Keepers,1);
        lenghtD =  sum(gdata.BW_Skel(:)) * 0.16; % in um
        time = (numel(gdata.I)) * 0.1; %in sec
        Freq = (numevents/time)
%         FreqNorm =  (numevents/time)/lenghtD;
        save([gdata.StreamName(1:end-4) '_Freq_Final.txt'], 'Freq', '-ascii', '-double')
%         save([gdata.StreamName(1:end-4) '_FreqNorm.txt'], 'FreqNorm', '-ascii', '-double')
        save([gdata.StreamName(1:end-4) '_DeltaF_Final.txt'], 'IdedDeltaF', '-ascii', '-double')
    end

    gdata.DeltaF = IdedDeltaF;
    gdata.Keepers = Keepers;
    gdata.Freq = Freq;
    guidata(hexo,gdata)

    try close(4); close(5); close(6);
    catch ME

    end
    
    GetDeltaF(hexo)
    
    
    
    
    end
