function Qafilt = field_filter(Qa, threshold, window, sfilt)

[Z, R] = meshgrid(1:size(Qa,2), 1:size(Qa,1));

TF  = del_outliers(Qa, threshold, window);
% figure; imshow(flip(TF)); colorbar;

PT = Qa(~TF);
ZT = Z(~TF);
RT = R(~TF);

PT = scatteredInterpolant(ZT(:),RT(:),PT(:),'linear','nearest');
P2 = PT(Z,R);

Qafilt = imgaussfilt(P2,[sfilt*size(P2,1)/100, sfilt*size(P2,2)/100]);

Qafilt(Qafilt<1e-10) = 1e-10;
end

function TF = del_outliers(F, threshold, window)
sumTFi = 1; TF = 0*F;
while ~(sumTFi==0)
    M1 = movmedian(F, window,1,'omitnan');
    M2 = movmedian(F, window,2,'omitnan');
    TFi = (F > threshold * M1 | F > threshold * M2) & F > prctile(F(:),75);
    sumTFi = sum(TFi(:));
%     figure; imshow(flip(TFi)); colorbar;
    TF = TFi | TF;
    F(TF) = nan;
end
disp(sum(TF(:)));
end


% function TF = del_outliers2(F, threshold, window)
% sumTFi = 1; TF = 0*F;
% while ~(sumTFi==0)
%     M = movmean(movmean(F, window,1,'omitnan'),2,'omitnan') * threshold;
%     TFi = F > M & F > prctile(F(:),80);
%     sumTFi = sum(TFi(:));
%     disp(sumTFi/numel(TFi));
% %     figure; imshow(flip(TFi)); colorbar;
%     TF = TFi | TF;
%     F(TF) = nan;
% end

