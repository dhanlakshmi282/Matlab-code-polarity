%Dhanalaxmi code to detect polarized cells from multiple stainings.
%Lokesh Pimpale 2022

clear all ;close all;
fileList = dir('*.tif');
Polarized = nan(1,length(fileList));

for i = 1 : length(fileList)
    %% DO CHECK THIS
    f_dapi = 1;   %3
    f_refch1 = 2;
    f_refch2 = 3;  %1
    f_tub = 4;     %2
    dic = 5;

    %%
    figure(1)
    Polarized(i) = true;
    htemp = length(imfinfo(fileList(i).name));
    A = cell(htemp,1);    % create a empty array
    for q = 1:htemp
        A{q} = imread(fileList(i).name,q);
    end
    A = cat(3,A{:});

    mip = max(cat(3,A(:,:,1:end)),[],3);                           % Maximum Intensity Projection
    %mip = A(:,:,dic);
    %backgroundt = imopen(mip,strel('disk',20));

    %background = imopen(mip,strel('disk',30));       % similar to rolling ball average. Change the 40 if you see artifacts
    %I2 = imsubtract(mip,backgroundt);                     % back ground subtraction
    I2 = imadjust(mip);
    thresholdval = multithresh(I2,3);                       % how many clusters of intensity you expect.
    %1 will be cell, 2 is overall cell inside and 3 is high intensity pixels

    cellbw1 = (I2>thresholdval(1));                          % this is cell mask
    cellbw1 = imclearborder(cellbw1);
    cellbw2 = imbinarize(I2);
    cellbw2 = imclearborder(cellbw2);
    if sum(cellbw2(:)) > sum(cellbw1(:))
        cellbw = cellbw2;
    else
        cellbw = cellbw1;
    end
    cellbw = imdilate(cellbw,strel('disk',2));
    cellbw = bwareafilt(cellbw,1);                            % consider the biggest object which is the cell
    cellbw = imfill(cellbw,'holes');
    cellbw = imerode(cellbw,strel('disk',8));            % default is 5
    cellbw = bwconvhull(cellbw);
    cellinfo = regionprops(cellbw,'all');                   % this is cell information, centroid, area bladibla

    cell_inside = imerode(cellbw , strel('disk',15));  % how much thickeness should be considered. was 15
    cell_border = ~cell_inside.*cellbw;                    % find cell border

    I2b = I2; % Create a backup for I2 just in case...
    x = cellbw.*double(A(:,:,f_dapi));
    %x = cellbw.*double(mip);
    lineplot1 = (mean(imgaussfilt(x,0.25)));
    lineplot2 = (mean(imgaussfilt(x',0.25)));
    for i = 1 : htemp-1
        tmp = findpeaks(mean(imgaussfilt(A(:,:,i),2)));
        x_temp(i) = length(tmp);
    end

    %% here we will move to floroscent channel - 3 and 4
    clear a b c d e f g h j k l m
    [a,b,c,d] = FindIntensity_cellborders(A,f_dapi,f_refch2,cellbw,cell_border);
    [e,f,g,h] = FindIntensity_cellborders(A,f_dapi,f_tub,cellbw,cell_border);
    [j,k,l,m] = FindIntensity_cellborders(A,f_dapi,f_refch1,cellbw,cell_border);

    if isempty(d)
        ch = string(f_refch2);
        d = array2table(zeros(1,6));
        d.Properties.VariableNames = [strcat('Polar',ch), strcat('CentroidX',ch), strcat('CentroidY',ch), ...
            strcat('CellIntensity',ch), strcat('RegionIntensity',ch), strcat('Ratio',ch)];
    end

    if isempty(h)
        ch = string(f_tub);
        h = array2table(zeros(1,6));
        h.Properties.VariableNames = [strcat('Polar',ch), strcat('CentroidX',ch), strcat('CentroidY',ch), ...
            strcat('CellIntensity',ch), strcat('RegionIntensity',ch), strcat('Ratio',ch)];
    end

    if isempty(m)
        ch = string(f_refch1);
        m = array2table(zeros(1,6));
        m.Properties.VariableNames = [strcat('Polar',ch), strcat('CentroidX',ch), strcat('CentroidY',ch), ...
            strcat('CellIntensity',ch), strcat('RegionIntensity',ch), strcat('Ratio',ch)];
    end

    clear centroid t
    %figure();imagesc(cellbw)
    if sum(cellbw(:)) >1
        centroid.cell = cellinfo.Centroid;
        t = regionprops(k);
        if ~isempty(t)
            centroid.ch2 =  regionprops(k).Centroid;
        else
            centroid.ch2 = [];
        end
        t = regionprops(b);
        if ~isempty(t)
            centroid.ch3 =  regionprops(b).Centroid;
        else
            centroid.ch3 = [];
        end
        t = regionprops(f);
        if ~isempty(t)
            centroid.ch4 =  regionprops(f).Centroid;
        else
            centroid.ch4 = [];
        end
        %pause
    end
    %%



    %% check for directionality

    clear directionCh3lr directionCh4lr Table.AngleBetweenSpots Table.AngleBetweenSpots2
    directionCh2lr = 0;
    directionCh2up = 0;
    directionCh3lr = 0;
    directionCh4lr = 0;
    directionCh3ud = 0;
    directionCh4ud = 0;
    angle_between_spots(i) = 0;
    angle_between_spots2(i) = 0;

    if sum(cellbw(:)) >1

        if ~isempty(centroid.cell)

            if ~isempty(centroid.ch2)
                if (centroid.cell(1)-centroid.ch2(1)) >0
                    directionCh2lr = -1;                             % its on the left of centre
                else
                    directionCh2lr = 1;                               % its on the right of centre
                end

                if (centroid.cell(2)-centroid.ch2(2)) >0
                    directionCh2ud = -1;                           % its below the centre
                else
                    directionCh2ud = 1;                             % its above the centre
                end
            end

            if ~isempty(centroid.ch3)
                if (centroid.cell(1)-centroid.ch3(1)) >0
                    directionCh3lr = -1;                             % its on the left of centre
                else
                    directionCh3lr = 1;                               % its on the right of centre
                end

                if (centroid.cell(2)-centroid.ch3(2)) >0
                    directionCh3ud = -1;                           % its below the centre
                else
                    directionCh3ud = 1;                             % its above the centre
                end
            end

            if  isfield(centroid,'ch4')
                if ~isempty(centroid.ch4)
                    if (centroid.cell(1)-centroid.ch4(1)) >0
                        directionCh4lr = -1;
                    else
                        directionCh4lr = 1;
                    end
                    if (centroid.cell(2)-centroid.ch4(2)) >0
                        directionCh4ud = -1;                           % its on the bttom
                    else
                        directionCh4ud = 1;
                    end
                end
            end
        end

        if  isfield(centroid,'ch4') && isfield(centroid,'ch2')
            if ~isempty(centroid.ch4) && ~isempty(centroid.ch2)
                cellcent = scatter(cellinfo(1).Centroid(1),cellinfo(1).Centroid(2));
                P0 = centroid.cell;
                P1 = centroid.ch2;
                P2 = centroid.ch4;
                angle_between_spots2(i) = atan2d(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0));
            end
        end

        if  isfield(centroid,'ch4') && isfield(centroid,'ch3')
            if ~isempty(centroid.ch4) && ~isempty(centroid.ch3)
                cellcent = scatter(cellinfo(1).Centroid(1),cellinfo(1).Centroid(2));
                P0 = centroid.cell;
                P1 = centroid.ch3;
                P2 = centroid.ch4;
                angle_between_spots(i) = atan2d(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0));
            end
        end

    end


    if exist('my_table','var') == 0
        % Make N by 2 matrix of fieldname + value type
        variable_names_types = ...
            [[strcat('Polar',string(f_tub)),'double'];[strcat('CentroidX',string(f_tub)), "double"];[strcat('CentroidY',string(f_tub)), "double"];[strcat('CellIntensity',string(f_tub)), "double"];[strcat('RegionIntensity',string(f_tub)), "double"];[strcat('Ratio',string(f_tub)), "double"]; ...
            [strcat('Polar',string(f_refch2)),'double'];[strcat('CentroidX',string(f_refch2)), "double"];[strcat('CentroidY',string(f_refch2)), "double"];[strcat('CellIntensity',string(f_refch2)), "double"];[strcat('RegionIntensity',string(f_refch2)), "double"];[strcat('Ratio',string(f_refch2)), "double"];...
            [strcat('Polar',string(f_refch1)),'double'];[strcat('CentroidX',string(f_refch1)), "double"];[strcat('CentroidY',string(f_refch1)), "double"];[strcat('CellIntensity',string(f_refch1)), "double"];[strcat('RegionIntensity',string(f_refch1)), "double"];[strcat('Ratio',string(f_refch1)), "double"];...
            ["AngleBetweenSpots",'double']; ["AngleBetweenSpots2",'double']];
        % Make table using fieldnames & value types from above
        my_table = table('Size',[1,size(variable_names_types,1)],...
            'VariableNames', variable_names_types(:,1),...
            'VariableTypes', variable_names_types(:,2));
    end

    %%
    if isempty(d) + isempty(h) < 2 || isempty(m) + isempty(h) < 2

        Table = [m,d,h];
        if exist("angle_between_spots",'var')
            Table.AngleBetweenSpots = angle_between_spots(i) ;
        else
            Table.AngleBetweenSpots = nan;
        end
        if exist("angle_between_spots2",'var')
            Table.AngleBetweenSpots2 = angle_between_spots2(i) ;
        else
            Table.AngleBetweenSpots2 = nan;
        end
        my_table = vertcat(my_table, Table);
    end

    t = linspace(1, 2*pi);
    r = 1; x = r*cos(t); y = r*sin(t);

    subplot(2,3,1)
    tmp1 = imfuse(b,c);
    thickness = 10;
    imagesc(tmp1)
    if a
        rectangle('Position',[thickness/2 thickness/2 size(tmp1,2)-thickness size(tmp1,1)-thickness],'EdgeColor','g', 'LineWidth',round(thickness));
    else
        rectangle('Position',[thickness/2 thickness/2 size(tmp1,2)-thickness size(tmp1,1)-thickness],'EdgeColor','r', 'LineWidth',round(thickness));
    end
    subplot(2,3,2)
    tmp1 = imfuse(f,g);
    imagesc(tmp1)
    if e
        rectangle('Position',[thickness/2 thickness/2 size(tmp1,2)-thickness size(tmp1,1)-thickness],'EdgeColor','g', 'LineWidth',round(thickness));
    else
        rectangle('Position',[thickness/2 thickness/2 size(tmp1,2)-thickness size(tmp1,1)-thickness],'EdgeColor','r', 'LineWidth',round(thickness));
    end
    subplot(2,3,3)
    tmp1 = imfuse(k,l);
    imagesc(tmp1)
    if e
        rectangle('Position',[thickness/2 thickness/2 size(tmp1,2)-thickness size(tmp1,1)-thickness],'EdgeColor','g', 'LineWidth',round(thickness));
    else
        rectangle('Position',[thickness/2 thickness/2 size(tmp1,2)-thickness size(tmp1,1)-thickness],'EdgeColor','r', 'LineWidth',round(thickness));
    end
    subplot(2,3,4)
    imagesc(A(:,:,f_refch2));
    subplot(2,3,5)
    imagesc(A(:,:,f_tub));
    subplot(2,3,6)
    imagesc(A(:,:,f_refch1));
    %pause
end
my_table(1,:) = [];

str1 = strcat('Polar',string(f_refch2));
str2 = strcat('Polar',string(f_tub));
str3 = strcat('Polar',string(f_refch1));



tmpC = array2table(eval(strcat('my_table.',str1)).*eval(strcat('my_table.',str2)));
tmpC.Properties.VariableNames = {'Polarity1'};
tmpD = array2table(eval(strcat('my_table.',str3)).*eval(strcat('my_table.',str2)));
tmpD.Properties.VariableNames = {'Polarity2'};

names = cell2table({fileList.name}');
names.Properties.VariableNames = {'File'};
my_table = [tmpC my_table];
my_table = [tmpD my_table];

my_table = [names my_table];
close (figure(1))
subplot(1,3,1)
hist(eval(strcat('my_table.',str1)));
xlim([-0.5 1.5])
xticks([0 1])
pbaspect([0.5 1 0.5])
title('Polarity ch3')
subplot(1,3,2)
hist(eval(strcat('my_table.',str2)));
xlim([-0.5 1.5])
xticks([0 1])
pbaspect([0.5 1 0.5])
title('Polarity ch4')
subplot(1,3,3)
hist(eval(strcat('my_table.',str3)));
xlim([-0.5 1.5])
xticks([0 1])
pbaspect([0.5 1 0.5])
title('Polarity ch2')
saveas(gcf,'stats.pdf')


%% plot the percentage plots
figure(2)
tp1 = nnz(eval(strcat('my_table.',str1))==1)/numel(eval(strcat('my_table.',str1)))*100;
tp2 = nnz(eval(strcat('my_table.',str2))==1)/numel(eval(strcat('my_table.',str2)))*100;
tp3 = nnz(eval(strcat('my_table.',str3))==1)/numel(eval(strcat('my_table.',str3)))*100;

hB = bar([tp1 tp2 tp3]);
hAx=gca;            % get a variable for the current axes handle
tstr = {str1,str2,str3};
hAx.XTickLabel=tstr; % label the ticks
hT=[];              % placeholder for text object handles
for i=1:length(hB)  % iterate over number of bar objects
    hT=[hT text(hB(i).XData+hB(i).XOffset,hB(i).YData,num2str(hB(i).YData.','%.3f'), ...
        'VerticalAlignment','bottom','horizontalalign','center')];
end
saveas(gcf,'percentages.pdf')


%%
writetable(my_table,'results.xlsx')

Mangle1 = nanmean(my_table.AngleBetweenSpots(my_table.Polarity1==1));
Sangle1 = nanstd(my_table.AngleBetweenSpots(my_table.Polarity1==1));
Mangle2 = nanmean(my_table.AngleBetweenSpots(my_table.Polarity2==1));
Sangle2 = nanstd(my_table.AngleBetweenSpots(my_table.Polarity2==1));
figure(45)
er = bar(1,Mangle1); hold on
er = errorbar(1,[Mangle1],[Sangle1]);
er.Color = [0 0 0];
er.LineStyle = 'none';
ylabel('Average angle of Polarity - degrees')







































%% MERGED CODE INTO one file

function [Polarized,BW2,Ineg,results] = FindIntensity_cellborders(A,f_dapi,channel,cellbw,cell_border)

%% how many clusters do we want to distinguish the pixels in ?
no_clusters = 5;
Polarized = true;
dapi = A(:,:,f_dapi);
dapi_bin = imbinarize(dapi);
%dapi_inside = imerode(dapi_bin , strel('disk',5));     % how much thickeness should be considered


I = A(:,:,channel); % original image to be analysed
background = imopen(I,strel('disk',40)); % similar to rolling ball average. Change the 40 if you see artifacts
I2 = imsubtract(I,background); % back ground subtraction
thresholdval = multithresh(I2,no_clusters); % how many clusters of intensity you expect.

Ipos = (double(I2>thresholdval(no_clusters)).*double(cell_border)); % Image with High intensity signal
Ineg = (double(I2<thresholdval(no_clusters)).*double(cell_border));  % Image with low intenstiy signal


%% mean intensity of entire cell
MeanIntCell = (double(I2).*cellbw);
MeanIntCell(MeanIntCell==0) = NaN;
MeanIntCell = mean(MeanIntCell(:),'omitnan');

%% mean intensity of bright area
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
[~,primary] = max([T2.Area]);

if tobj == 2
    if primary ==1
        secondary = 2;
    else
        secondary = 1;
    end

    tmp3 = pdist2(T2(secondary).Centroid,T2(primary).Centroid); % find distance from primary to secondary clusters
    tmp4 = regionprops(cellbw,'MajorAxisLength'); % find major axis which is reference

    if tmp3 < tmp4.MajorAxisLength/2   % if clusters are less than 1 radius away, check them

        tmp = T1(primary).MeanIntensity - (T1(primary).MeanIntensity*0.05);  % if secondary cluster intensity is <90% of primary cluster, ignore it  default 0.10
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

%cell_inside = imerode(cellbw , strel('disk',10));     % how much thickeness should be considered

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



%h_int_details = regionprops(BW2,'all');
Int_in_bright_area = (double(I2).*BW2);
Int_in_bright_area(Int_in_bright_area==0) = NaN;

no_pixels_bright = sum(~isnan(Int_in_bright_area(:)));
Int_in_bright_area = mean(Int_in_bright_area(:),'omitnan');

%% mean intensity of low intensity area
Int_in_dark_area = double(I2);
Int_in_dark_area(Ineg==0) = NaN;
no_pixels_dark = sum(~isnan(Int_in_dark_area(:)));
Int_in_dark_area = mean(Int_in_dark_area(:),'omitnan');

%% start with conditions. set thresholds here...
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
    disp('ratio screwed up')
    results = [];
    return
end

%% lets see if there are lot of dots around
% find out the centroid of the cell
temp = regionprops(cellbw,'MajorAxisLength');
cellcentre = regionprops(cellbw,'Centroid');
r = round(temp.MajorAxisLength/2/10); % this is the radius of the inner cicle - now dependent on cell size
%r=5; % this is the radius of the inner cicle -
x_c=cellcentre(1).Centroid(1); y_c=cellcentre(1).Centroid(2);
imageSizeX = size(A(:,:,1),2);
imageSizeY = size(A(:,:,1),1);
[columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
circlePixels = (rowsInImage - y_c).^2 ...
    + (columnsInImage - x_c).^2 <= r.^2;

if sum(sum((bwconvhull(BW2)).* circlePixels)) ~= 0
    Polarized = false;
    results = [];
    disp('circle is problem')
end


allmasks = bwconvhull(BW2);
allmask = sum(allmasks(:));
allcell = sum(cellbw(:));
if allmask/allcell < 0.5
    disp('here');
end
%imagesc(convhull(BW2))

forcal = bwareafilt(BW2,1);
angleofregion = regionprops(forcal,'Centroid');

subplot(2,2,3)
imagesc(double(I2).*cell_border)

subplot(2,2,4)
t = linspace(1, 2*pi);
r = 3; x = r*cos(t); y = r*sin(t);
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
%ratio = ratio;
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

