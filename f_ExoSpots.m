function [all_Events] = f_ExoSpots(BS)

% Detection based in part on the Spot (fluorophore transition event) detection
% algorithm, in 'Single-molecule-based super-resolution images in the
% presence of multiple fluorophores, Paul D. Simonson, Eli Rothenberg,
% Paul R. Selvin, Nano Lett.  11, 11, 5090-5096'.

    %1. Change all pixels values that are above 0
    %     BS((BS > 0)) =0;
    %     figure(2); imshow(BS,[]);colorbar; text(10,10,'Background Subtracted Image', 'color','w')

    %2. Dilate the image. Each pixel intensity value is replace with the
    %Intensity of the smallest value of 8 neighbors. Do 2X for low SNR.
    BS_dilate =imdilate(BS,ones(3));
    BS_dilate =imdilate(BS_dilate,ones(3));
    %     figure(3); imshow(BS_dilate,[])

    %3. Find the Local Minima; positions of all pixels not affected by the dilatation.
    local_minima = BS_dilate-BS;
    local_minima = local_minima==0;
    local_minima = double(local_minima) .* double(BS);
    Even = sortrows(double(local_minima(:)));
    BW2 = ismember(local_minima,Even(end-100:end)); % consider the 100 brightest events per frame
    %     BW2 = ismember(local_minima,Even(1:end)); % all events could be
    %     big...
    %     figure(4); imshow(local_minima,[]); colorbar ; text(10,10,'Local Minima', 'color','w')
    %     local_minima = local_minima > 0;
    %     %local_minima = local_minima > Minima_Thres; % doit autmatiser le threshold sur minima
    %     figure(5); imshow(local_minima,[]); colorbar; text(10,10,'Local Minima', 'color','w')

    %4.Find the Spot Average Pixel Intensity. Determine average pixel intensity of all
    %  pixels within radius 'R'.

    %5.Find the Average Background Intensity. Determine average pixel intensity of all
    %  pixels between radius 'R + dR'.

    % Get dF for deltaF
    % %%
    [y x] = find(BW2>0);
    all_Events = [];
    if size(x,1) ~= size(BS,1) .* size(BS,2) % do not compute if all pixels are positive (prevents black frame)
        for i=1:size(x,1)
            Roi = f_cut_square_on_Image([x(i) y(i)], BS, 3);
            dFlocal  = mean(double(Roi(:)));
            event = [ x(i) y(i) dFlocal];
            all_Events = [all_Events ; event ];
        end
    end
end
