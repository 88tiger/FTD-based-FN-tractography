%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create in HangZhou, on Oct.16, 2023
% AONR (Automatic Optic Nerves Reconstruction)
% Author: Qiming Hu
% Mail: 2112103115@zjut.edu.cn
% Version: MATLAB R2023a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add environment
clear;clc;
addpath("NIFTI");

%% Parameters
step_size = 0.3;
FiberCount = 10000; % number of fibers for each subdivision
FiberCup = 0; % Is fibercup data?
multi = 1;

%% Load Data
name = 'your folder path';
path = ['NIFTI/',name];
DWIFile = ['NIFTI/',name,'/data_eddy.nii.gz'];
PeaksFile = ['NIFTI/',name,'/peaks.nii.gz'];
StreamFiles = ['NIFTI/',name,'/fiber1000_iFOD1_FA05_angle30.tck'];
tckmapFile = ['NIFTI/',name,'/fiber1000_iFOD1_FA05_angle30.nii.gz'];
IACFile = ['NIFTI/',name,'/ROI_2.nii.gz'];
CPAFile = ['NIFTI/',name,'/ROI_1.nii.gz'];
ftrackname = ['NIFTI/',name,'/FN'];

%% Read file
DWIdata = niftiread(DWIFile); DWIinfo = niftiinfo(DWIFile);
DWIaffine = DWIinfo.Transform.T; DWIaffine = DWIaffine'; DWIaffineInv = DWIaffine^-1;
peaks = niftiread(PeaksFile);
DWIdim = size(DWIdata); DWIdim = DWIdim(1:3);
tckmap = niftiread(tckmapFile); 
RGVPregion = tckmap;
RGVPregion((RGVPregion>=1)) = 1;RGVPregion((RGVPregion<1)) = 0;
RGVPregion = uint8(RGVPregion);
mask = RGVPregion;
IACpoint = niftiread(IACFile);
CPApoint = niftiread(CPAFile);
endpoint = IACpoint.*mask;
endind = find(endpoint==1);
[endx, endy, endz] = ind2sub(DWIdim, endind);
endposition = [endx,endy,endz];

%% Start point and IAC (end) point

[CPApointx,CPApointy,CPApointz] = ind2sub(size(CPApoint),find(CPApoint==1));
startposition = [CPApointx,CPApointy,CPApointz];

% Find IAC (start) points
[IACpointx,IACpointy,IACpointz] = ind2sub(size(IACpoint),find(IACpoint==1));
IACposition = [IACpointx,IACpointy,IACpointz];


%% Preprocess peaks
[Peaks, WeightedPeaks] = PreparingPeaks(peaks, 0.005, 0.6, 0);%peaks
Peaks(:,:,:,:,1) = -Peaks(:,:,:,:,1);
WeightedPeaks(:,:,:,:,1) = -WeightedPeaks(:,:,:,:,1);

savefilename = strcat(ftrackname,'_order4.tck');
flowfilename = strcat(ftrackname,'_flow');
selectfilename = strcat(ftrackname,'_select');

[CPApositionx,CPApositiony,CPApositionz] = ind2sub(size(CPApoint),find(CPApoint==1));
CPAposition = [CPApositionx,CPApositiony,CPApositionz];

%% FN streamlines for directions guidence
% startposition = leftstartposition;
tckimg = read_mrtrix_tracks(StreamFiles); oristreamlines = tckimg.data; oristreamlines = oristreamlines';
oristreamlines = transform(oristreamlines, DWIaffineInv);

%conpute CPA ROI mean
CPApositionmean = mean(CPAposition);
% Change start point if streamlines(1) not in the start ROI
for i=1:size(oristreamlines,1)
    Tract = oristreamlines{i};
    intTract = round(Tract);
    intTractyrange = [max(intTract(1:30,2)):-1:max(intTract(1:30,2))-20];
    
    if sqrt(sum((CPApositionmean - intTract(1,:)).^2)) < sqrt(sum((CPApositionmean - intTract(end,:)).^2))
        oristreamlines{i} = Tract(end:-1:1,:);
    end
    
    [Lia, Locb] = ismember(round(oristreamlines{i}), IACposition, 'rows');
    value = min(Locb(Locb > 0));
    if isempty(value)
        oristreamlines{i} = [];
        continue
    end
    p = max(find(Locb == value));
    oristreamlines{i}(1:p,:) = [];
   
end
clear CPApositionmean Lia Locb value p;

%% Get Streamlines peaks
for i=1:size(oristreamlines,1)
    oristreamlines{i} = [oristreamlines{i};0 0 0];
end

oristreamlines = cell2mat(oristreamlines);
oristreamlines = oristreamlines+[0.5,0.5,0.5];
BundlePeaks = oristreamlines(2:end,:)-oristreamlines(1:end-1,:);
BundleROIs = round(oristreamlines(2:end,:));
BundlePeaksnorm = sqrt(sum(BundlePeaks.^2,2));
BundlePeaks(BundlePeaksnorm>1,:) = [];
BundleROIs(BundlePeaksnorm>1,:) = [];
Bundlepeaks = BundlePeaks./repmat(sqrt(sum(BundlePeaks.^2,2)),1,3);


%% Intersect regions between tckmap and mask (segmented by Li et al method)
tckmapind = find(tckmap>0);
maskind = find(mask==1);
submaskind = intersect(maskind, tckmapind);
submask = zeros(size(mask));
submask(submaskind) = 1;
submask=single(submask);

%% Correct directions by streamlines
DirsROI = zeros(0,3); WeightedDirsROI = zeros(0,3); ROIpositions = zeros(0,3);
[submaskx,submasky,submaskz] = ind2sub(size(submask),find(RGVPregion==1));
submaskposition = [submaskx,submasky,submaskz];

BundleROIsind = sub2ind(size(mask),BundleROIs(:,1),BundleROIs(:,2),BundleROIs(:,3));
[C,IA,IC] = unique(BundleROIsind);
Peaks_reshape = reshape(Peaks, size(Peaks,1)*size(Peaks,2)*size(Peaks,3),3,3);
WeightedPeaks_reshape = reshape(WeightedPeaks, size(Peaks,1)*size(Peaks,2)*size(Peaks,3),3,3);
Peaks_ROIs = Peaks_reshape(BundleROIsind,:,:);
BundlePeaks_reshape = reshape(BundlePeaks,size(BundlePeaks,1),1,3);
repBundlePeaks_reshape = repmat(BundlePeaks_reshape,1,3,1);
cosval_streampeak = repBundlePeaks_reshape.*Peaks_ROIs;
cosval_streampeak = sum(cosval_streampeak,3);
abscosval_streampeak = abs(cosval_streampeak);
[cosvalmax,cosvalind] = max(abscosval_streampeak,[],2);

a = cosvalind; a(cosvalind==1) = 1; a(cosvalind~=1) = 0;
b = cosvalind; b(cosvalind==2) = 1; b(cosvalind~=2) = 0; 
c = cosvalind; c(cosvalind==3) = 1; c(cosvalind~=3) = 0;
abc = [a,b,c];
cosvalori = sum(cosval_streampeak.*abc,2);

virtual_cosvalind(cosvalind==1) = 100000000;
virtual_cosvalind(cosvalind==2) = 10000;
virtual_cosvalind(cosvalind==3) = 1;

final_cosvalind = zeros(size(C));
check_inverse = zeros(size(C));
for i=1:size(C,1)
    if size(find(BundleROIsind==C(i)),1)>1 % max(tckmap(:))/2.
        final_cosvalind(i) = sum(virtual_cosvalind(find(BundleROIsind==C(i))));
        check_inverse(i) = sum(cosvalori(find(BundleROIsind==C(i))));
    end
end
check_inverse(check_inverse>0) = 1; check_inverse(check_inverse<0) = -1;
rep_check_interse = repmat(check_inverse,1,3);
firnumber = fix(final_cosvalind/100000000.);
secnumber = fix((final_cosvalind-firnumber*100000000)/10000);
thinumber = fix(((final_cosvalind-firnumber*100000000-secnumber*10000)));

[~,maxind] = max([firnumber,secnumber,thinumber],[],2);
maxindnum = [1:size(maxind,1)];maxindnum = maxindnum';


[Cx,Cy,Cz] = ind2sub(size(mask),C);
Cposition = [Cx,Cy,Cz];
Cposition(:,1) = 126-Cposition(:,1);
C_ind = sub2ind(size(mask),Cposition(:,1)',Cposition(:,2)',Cposition(:,3)');

final_ind = size(mask,1)*size(mask,2)*size(mask,3)*(maxind-1)+C;

Peaks_reshape_re = reshape(Peaks_reshape,size(Peaks_reshape,1)*3,3);
WeightedPeaks_reshape_re = reshape(WeightedPeaks_reshape,size(WeightedPeaks_reshape,1)*3,3);
DirsROI = Peaks_reshape_re(final_ind,:); DirsROI = DirsROI.*rep_check_interse;
WeightedDirsROI = WeightedPeaks_reshape_re(final_ind,:); WeightedDirsROI = WeightedDirsROI.*rep_check_interse;
repC = repmat(C,1,3);
repC = reshape(repC,size(repC,1)*3,1);

[ROIpositionsx,ROIpositionsy,ROIpositionsz] = ind2sub(size(mask),C);
ROIpositions = [ROIpositionsx,ROIpositionsy,ROIpositionsz];
[~,~,interROIpositionind] = intersect(submaskposition,ROIpositions,'rows');
ROIpositions = ROIpositions(interROIpositionind,:);
DirsROI = DirsROI(interROIpositionind,:);
WeightedDirsROI = WeightedDirsROI(interROIpositionind,:);

ROIpositions = ROIpositions(any(DirsROI~=0,2),:);
WeightedDirsROI = WeightedDirsROI(any(DirsROI~=0,2),:);
DirsROI = DirsROI(any(DirsROI~=0,2),:);

%% selected directions threshold
dirind = zeros(size(DirsROI,1),1);
for number=1:size(DirsROI,1)
    pos = ROIpositions(number,:);
    posind = sub2ind(size(mask),pos(1),pos(2),pos(3));
    bpeaks = BundlePeaks(find(IC==find(C==posind)),:);
    posval = mean(DirsROI(number,:)*bpeaks');
    if posval<0.05
        dirind(number,1) = 1;
    end
end
ROIpositions(find(dirind==1),:) = [];
WeightedDirsROI(find(dirind==1),:) = [];
DirsROI(find(dirind==1),:) = [];

%% Get seeds
seeds = GenSeeds(IACposition,FiberCount);

%% Calculate A

%     [A1, fvalIntra1] = GetATernaryCubic(ROIpositions, DirsROI, WeightedDirsROI);
    [A1, fvalIntra1] = GetATernaryForth(ROIpositions, DirsROI, WeightedDirsROI);
%     [A1, fvalIntra1] = GetATernaryFifth(ROIpositions, DirsROI, WeightedDirsROI);
%     [A1, fvalIntra1] = GetATernarySeventh(ROIpositions, DirsROI, WeightedDirsROI);
%     [A1, fvalIntra1] = GetATernaryTenth(ROIpositions, DirsROI, WeightedDirsROI);
    Aall = cell(1,1); 
    Aall{1} = A1; 

    flowpeaks = zeros(size(peaks));
    [ROIx,ROIy,ROIz] = ind2sub(size(mask),maskind);
    ROI = [ROIx,ROIy,ROIz];
    Dirs = GetFlowPeaks(A1,ROI,FiberCup);
    for i=1:size(Dirs,1)
        flowpeaks(ROI(i,1),ROI(i,2),ROI(i,3),1:3) = Dirs(i,:);
    end


selectedpeaks = zeros(size(peaks));
for i=1:size(DirsROI,1)
    DirsROI(i,2) = -DirsROI(i,2);
    selectedpeaks(ROIpositions(i,1),ROIpositions(i,2),ROIpositions(i,3),1:3) = DirsROI(i,:);  
end
flowpeaks = zeros(size(peaks));
[ROIx,ROIy,ROIz] = ind2sub(size(mask),maskind);
ROI = [ROIx,ROIy,ROIz];
Dirs = GetFlowPeaks(A1,ROI,FiberCup);
for i=1:size(Dirs,1)
    flowpeaks(ROI(i,1),ROI(i,2),ROI(i,3),1:3) = Dirs(i,:);
end
info = niftiinfo(PeaksFile);
flowpeaks(:,:,:,2) = -flowpeaks(:,:,:,2);
info.Description = 'Modified using MATLAB R2023a';
info.Datatype = 'double';
niftiwrite(flowpeaks,flowfilename,info);
info.Datatype = 'double';
niftiwrite(selectedpeaks,selectfilename,info);

%% Tracking
Tractsall = cell(1,1); 
wholemask = ones(size(mask));
for k=1:size(Aall,2)
    Tracts = cell(1, FiberCount);
    A = Aall{k};
    A(3,:) = A(3,:)*DWIinfo.PixelDimensions(1)/DWIinfo.PixelDimensions(3); 
    Num = 0;
    for FiberNum = 1:size(seeds,1)
        start_point = seeds(FiberNum,:);
        if ~isempty(A) && any(A(:)~=0)
            FTrace = ftrack(A,start_point,step_size,mask,RGVPregion,CPApoint);
        end
        % Except empty Trace
        if ~isempty(FTrace)
            Tract = FTrace;
        else
            Tract = [];
        end

        intTract = ceil(Tract);

        if size(Tract, 1) <= 20
            continue
        elseif ~isempty(intersect(CPAposition,squeeze(intTract(:,:)),'rows'))
            Tract = Tract-[0.5,0.5,0.5]; % Mrtrix is start from 0, Matlab is start from 1.
            Tracts{FiberNum} = Tract; Num = Num + 1;
            fprintf('%5d/%d/%d\n', Num, FiberNum, FiberCount);
        else
            fprintf('%5d/%d/%d\n', Num, FiberNum, FiberCount);
        end
    end
    Tracts(cellfun('length', Tracts)==0) = [];
    Tracts = Tracts';
    Tractsall{k} = Tracts;
end

Tracts = transform(Tracts, DWIaffine);
savetckimg = tckimg;
savetckimg.command_history = 'MATLAB 2023a';
savetckimg.data = Tracts';
write_mrtrix_tracks(savetckimg,savefilename);
clear Tracts;

%% Save peaks
    flowpeaksinfo = niftiinfo(PeaksFile);
    flowpeaksinfo.Description = 'Modified using MATLAB R2023a';
    flowpeaksinfo.Datatype = 'double';
    filename = 'flowpeaks.nii'; % must be 'XXX' rather than "XXX"
    niftiwrite(flowpeaks, filename, flowpeaksinfo)
