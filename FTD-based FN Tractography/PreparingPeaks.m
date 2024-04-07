function [peaks,WeightedPeaks] = PreparingPeaks(Data,Maxima,Thr_Maxima,Norm)

if nargin == 3
    Norm = 0;
end

% Reshape Peaks
DataSize = size(Data);
Data(isnan(Data)) = 0;
peaks = reshape(Data, DataSize(1), DataSize(2), DataSize(3), 3, DataSize(4)/3);
peaks = permute(peaks, [1,2,3,5,4]);
peaks = peaks(:,:,:,1:3,:);

% Filtering Small-Value Maxima Peaks
PeaksNormSet = repmat(sqrt(sum(peaks.^2,5)), 1, 1, 1, 1, size(peaks,4));
PeaksNormSet_VectorValue = repmat(squeeze(PeaksNormSet(:,:,:,1,:)), 1, 1, 1, 1, size(peaks,4));

peaks(PeaksNormSet_VectorValue<Maxima) = 0;

% Filtering Small-Value Lower than Thr * MaximaPeaks
for i = 2:size(PeaksNormSet, 4)
    MinPeaksNorm = zeros(size(PeaksNormSet)); PeaksNorm = zeros(size(PeaksNormSet));
    MinPeaksNorm(:,:,:,i,:) = Thr_Maxima * PeaksNormSet(:,:,:,1,:);
    PeaksNorm(:,:,:,i,:) = PeaksNormSet(:,:,:,i,:);
    Dvalue = PeaksNorm - MinPeaksNorm;
    Dvalue(Dvalue>=0)=1; Dvalue(Dvalue<0)=0;
    peaks = peaks.*Dvalue;
    PeaksNormSet = PeaksNormSet.*Dvalue;
end

% peaks Normalization
if Norm
    WeightedPeaks = peaks./PeaksNormSet_VectorValue;
    WeightedPeaks(isnan(WeightedPeaks)) = 0;
else
    WeightedPeaks = peaks;
end
peaks = peaks./PeaksNormSet;
peaks(isnan(peaks)) = 0;

end
