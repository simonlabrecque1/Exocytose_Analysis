% function Exocytose_PeakDetector
% cd('Y:\Simon\Netrin\Simon_Submission_3\Matlab\Data_Set')
[StreamName,directory_name] = uigetfile(('*.tif'),'Select the Movie');
cd(directory_name)

% Read Stream movie
info = imfinfo(StreamName);
h = waitbar(0,'Please wait reading Stream...');
I= [];
for i=1:numel(info)
    I(i).data = imread(StreamName,i);
    waitbar(i/numel(info))
end
close(h)
I_max = f_StackProjections(I,'Maximum');
imwrite(I_max,[StreamName(1:end-4) '_MaxProjection.tif'],'tif','Compression','none')

 
% Get BCKG (defined as the fisrt 25 frames of the movie)
I_Bckg=[];
info = imfinfo(StreamName);
for i=1:25
    I_Bckg(i).data = imread(StreamName,i);
end
I_Bckg = f_StackProjections(I_Bckg,'Mean');
% figure
% imshow(I_Bckg,[])

h = waitbar(0,'Please wait subtract SEP...');
for i=1:numel(I)
   I(i).data = I(i).data-I_Bckg;
   waitbar(i/numel(info))
end
close(h)

I_maxEvents = f_StackProjections(I,'Maximum');
imwrite(I_maxEvents,[StreamName(1:end-4) '_MaxEvents.tif'],'tif','Compression','none')

% % Get Spots
infoSEP = imfinfo(StreamName);
Spots = zeros(1,4);
h = waitbar(0,'Please wait computing Spots...');
for i=1:numel(infoSEP)-1
    set(h,'Name',['Detect Spots. Frame ' num2str(i) 'of ' num2str(numel(infoSEP)-1)] )
    SEP1 = I(i).data;
    SEP2 = I(i+1).data;
    BS   = SEP2- SEP1;
    [Spot] = f_ExoSpots(BS);
    Spot = [Spot zeros(size(Spot,1),1)+i];
    Spots = [Spots ; Spot];
    waitbar(i/numel(infoSEP))
end
close(h)  
Spots = [Spots(:,1:2) Spots(:,4)  Spots(:,3)];   
Spots = sortrows(Spots,-4);
save([StreamName(1:end-4) '_Spots.mat'], 'Spots')

%% DF/F comptutation
% At this point we identified appearing spots on the movie. Next, we
% compute the DeltaF over F value for each spots. We will use this curve to
% discriminate noise from real exocytose events.

% Go back to original data
info = imfinfo(StreamName);
h = waitbar(0,'Please wait reading Stream...');
for i=1:numel(info)
    I(i).data = imread(StreamName,i);
    waitbar(i/numel(info))
end
close(h)

h = waitbar(0,'Please wait computing DeltaFs...');
Spots_Area=[];
j=1;
for i = 1:size(Spots,1)
    if Spots(i,3)+10 > numel(I)
       DeltaFrame = (Spots(i,3)-2:numel(I))';
    elseif Spots(i,3)-2 <= 0
        DeltaFrame = (1:Spots(i,3)+10)';
    else
        DeltaFrame = (Spots(i,3)-2:Spots(i,3)+10)';
    end
    %   Ided(:,StartFrame:EndFrame) = [x y frame sum(DeltaF)]
    Ided =  [ones(size(DeltaFrame,1),1).* Spots(i,1)  ones(size(DeltaFrame,1),1).*Spots(i,2) DeltaFrame];
    % DeltaF
    [DeltaF]  = f_DeltaF_SquareD(Ided, I,2, 3);

    Spots(i,1:4)=[Spots(i,1:3) sum(DeltaF)];
    waitbar(i/size(Spots,1))
%     figure(1)
%     Montage = f_MakeMontage(Ided, I,20);
%     Montage = imresize(Montage,2);
% %     set(1,'Position',[44 68 1230 168])
%     imshow(Montage,[],'initialMagnification',500)
%     figure(2)
%     plot(DeltaF)
%     legend(num2str(sum(DeltaF)))
%     pause
end
% Spots= Spots_Area;
save([StreamName(1:end-4) '_Spots_wDeltaF.mat'], 'Spots')
close(h)
%%


Exocytose_Analysis

