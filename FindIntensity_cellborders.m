function [Polarized,BW2,Ineg,results] = FindIntensity_cellborders(A,f_dapi,channel,cellbw,cell_border)

no_clusters = 4;
Polarized = true;
dapi = A(:,:,f_dapi);
dapi_bin = imbinarize(dapi);
dapi_inside = imerode(dapi_bin , strel('disk',5));     % how much thickeness should be considered


I = A(:,:,channel); % original image to be analysed
background = imopen(I,strel('disk',40)); % similar to rolling ball average. Change the 40 if you see artifacts
I2 = imsubtract(I,background); % back ground subtraction
thresholdval = multithresh(I2,no_clusters); % how many clusters of intensity you expect.

Ipos = (double(I2>thresholdval(no_clusters)).*double(cell_border)); % Image with High intensity signal
Ineg = (double(I2<thresholdval(no_clusters)).*double(cell_border));  % Image with low intenstiy signal


% mean intensity of entire cell
MeanIntCell = (double(I2).*cellbw);
MeanIntCell(MeanIntCell==0) = NaN;
MeanIntCell = nanmean(MeanIntCell(:));

% mean intensity of bright area
Int_in_bright_area= double(I2);
Int_in_bright_area(Ipos==0) = 0;         % only consider area that is bright in cell border area
cellbwt = cellbw;
cellbwt(Int_in_bright_area==0) = 0;
BW2 = bwareaopen(cellbwt, 45);          % how big amount of pixels should be removed
if sum(BW2(:)) == 0
    Polarized = false;
    results = [];
    return
end
cc = bwconncomp(BW2,4);         % count the number of objects in the bwimage
no_objects = cc.NumObjects;

for ii = 1 : no_objects
    len_objects(ii) = length(cc.PixelIdxList{ii});
end

sorted_objects = sort(len_objects);
if length(sorted_objects) > 1
    if sorted_objects(end-1) > sorted_objects(end) - 0.2* sorted_objects(end)
        Polarized = false;
        results = [];
        return
    end
end

%% Here we check if there are two clusters. if secondary clusters is less than 15 percent of original cluster intensity, ignore it. 

T = bwconncomp(BW2,18);         % count the number of objects in the bwimage
tobj = T.NumObjects;
T1= regionprops(BW2,(double(I2)),'MeanIntensity');
T2 = regionprops(BW2,'all');
[x,primary] = max([T2.Area]);

if tobj == 2
    if primary ==1
        secondary = 2;
    else 
        secondary = 1;
    end
    
tmp3 = pdist2(T2(secondary).Centroid,T2(primary).Centroid); % find distance from primary to secondary clusters
tmp4 = regionprops(cellbw,'MajorAxisLength'); % find major axis which is reference

if tmp3 < tmp4.MajorAxisLength/2   % if clusters are less than 1 radius away, check them 
    
    tmp = T1(primary).MeanIntensity - (T1(primary).MeanIntensity*0.10);  % if secondary cluster intensity is <90% of primary cluster, ignore it 
    if T1(secondary).MeanIntensity < tmp
        BW2(T2(secondary).PixelIdxList) = 0;
    end
    
    tmp2 = (T2(primary).Area*0.25); % if secondary cluster is 25% lower than primary. 
    if T2(secondary).Area < tmp2
        BW2(T2(secondary).PixelIdxList) =0;
    end
    
end
end

%%

cell_inside = imerode(cellbw , strel('disk',10));     % how much thickeness should be considered

if length(sorted_objects) > 2
    if sum(sum(bwconvhull(BW2))) / sum(cellbw(:)) > 0.6
        Polarized = false;
        results = [];
        return
    end
end
subplot(2,2,1)
imagesc(I2)
subplot(2,2,2)
imagesc(imfuse(BW2,Ineg))



h_int_details = regionprops(BW2,'all');
Int_in_bright_area = (double(I2).*BW2);
Int_in_bright_area(Int_in_bright_area==0) = NaN;

no_pixels_bright = sum(~isnan(Int_in_bright_area(:)));
Int_in_bright_area = nanmean(Int_in_bright_area(:));
% mean intensity of low intensity area
Int_in_dark_area = double(I2);
Int_in_dark_area(Ineg==0) = NaN;
no_pixels_dark = sum(~isnan(Int_in_dark_area(:)));
Int_in_dark_area = nanmean(Int_in_dark_area(:));

% start with conditions. set thresholds here...
if Int_in_bright_area < 1.25*Int_in_dark_area || Int_in_bright_area <  MeanIntCell*1.25
    Polarized = false;
    results = [];
    return;
end

if isnan(Int_in_bright_area)
    Polarized = false;
    results = [];
    return;
end

if no_pixels_bright < 50 || no_pixels_bright / no_pixels_dark > 0.4
    Polarized = false;
    display('ratio screwed up')
    results = [];
    return
end

% lets see if there are lot of dots around
cellcentre = regionprops(cellbw,'Centroid');
r=7; x_c=cellcentre(1).Centroid(1); y_c=cellcentre(1).Centroid(2);
imageSizeX = size(A(:,:,1),2);
imageSizeY = size(A(:,:,1),1);
[columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
circlePixels = (rowsInImage - y_c).^2 ...
    + (columnsInImage - x_c).^2 <= r.^2;

if sum(sum((bwconvhull(BW2)).* circlePixels)) ~= 0
    Polarized = false;
    results = [];
    display('circle is problem')
end

forcal = bwareafilt(BW2,1);
angleofregion = regionprops(forcal,'Centroid');

subplot(2,2,3)
imagesc(double(I2).*cell_border)

subplot(2,2,4)
t = linspace(1, 2*pi);
r = 8; x = r*cos(t); y = r*sin(t);
if Polarized
    patch(x, y, 'g')
else
    patch(x, y, 'r')
end
xticks([]); yticks([])


temp= (Ineg.*double(I));
temp(temp==0) = NaN;
ma_cell_temp = regionprops(bwconvhull(Ineg),'MajorAxisLength');
ma_cell_temp = ma_cell_temp.MajorAxisLength(1);
cell_intensityMean  = mean(temp(:),'omitnan')/ma_cell_temp;

temp = (BW2.*double(I));
temp(temp==0) = NaN;
ma_reg_temp = regionprops(bwconvhull(BW2),'MajorAxisLength');
ma_reg_temp = ma_reg_temp.MajorAxisLength(1);
region_intensityMean  = mean(temp(:),'omitnan')/ma_reg_temp;

ratio = region_intensityMean / cell_intensityMean;
ratio = ratio;
ch = string(channel);
results = array2table([Polarized,angleofregion.Centroid(1),angleofregion.Centroid(2),cell_intensityMean,region_intensityMean,ratio]);
% results. = cell2table(results,'VariableNames',...
%     [strcat('Polar',ch),strcat('CentroidX',ch),strcat('CentroidY',ch),strcat('CellIntensity',ch),strcat('RegionIntensity',ch),strcat('Ratio',ch)]);
results.Properties.VariableNames = [strcat('Polar',ch),strcat('CentroidX',ch),strcat('CentroidY',ch),strcat('CellIntensity',ch),strcat('RegionIntensity',ch),strcat('Ratio',ch)];
%results.Properties.VariableNames = {'Polar','CentroidX','CentroidY','CellIntensity','RegionIntensity','Ratio'};

if exist('results','var')==0
    results = [];
end

end


