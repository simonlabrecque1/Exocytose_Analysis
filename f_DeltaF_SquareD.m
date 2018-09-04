function [DeltaF]  = f_DeltaF_SquareD(Ided, SEP,squareD, numF0)
    %Define Roi intensity on Ided.
    mRois = [];
    % stdRois = [];
    for time=min(Ided(:,3)):max(Ided(:,3))
        Roi = f_cut_square_on_Image(Ided(Ided(:,3)==time,1:2), SEP(time).data, squareD);
        mRoi = mean(double(Roi(:)));
        mRois = [mRois ; mRoi];
        %     stdRoi = std(double(Roi(:)));
        %     stdRois = [stdRois ; stdRoi];
    end
    % Compute DeltaF
    F0 = mean(mRois(1:numF0));
    DeltaF = (mRois-F0)./ F0;
end