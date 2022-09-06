
% plots receptive fields for a selection of footprints



[parentdir,~,~] = fileparts(pwd);
load(fullfile(parentdir, 'supporting_data/position_tracking.mat'));
load(fullfile(parentdir, 'supporting_data/extract_curated.mat'));





nbucket = 100;

ptile = 80;
speed = abs(diff(sunwrap));
speed(end+1)=speed(end);
speed = movmean(speed, 30);
p = prctile(speed(:), ptile);
isrunning = speed > p;
srun = s(isrunning);
fiT = iT(isrunning, :);

dir = diff(sunwrap);
dir(end+1) = dir(end);
dir = movmean(dir,30);
dir = -1+2*double(dir>0);



smin = prctile(srun, 1);
smax = prctile(srun, 99);
sbucket = (srun-smin)/(smax-smin);
sbucket(sbucket < 0) = 0;
sbucket(sbucket > .999) = .999;
sbucket = floor(sbucket*nbucket);



pthresh = 0.01;
ethresh = pthresh*max(fiT, [], 1);
events = double(fiT > repmat(ethresh, size(fiT, 1), 1));

placeD = zeros(nbucket, size(iT, 2));
eventD = zeros(nbucket, size(iT, 2));
resD = zeros(nbucket, 1);
for i = 1:numel(sbucket)
    placeD(sbucket(i)+1, :) = placeD(sbucket(i)+1, :) + fiT(i, :);
    eventD(sbucket(i)+1, :) = eventD(sbucket(i)+1, :) + events(i, :);
    resD(sbucket(i)+1) = resD(sbucket(i)+1) + 1;
end

normD = placeD./repmat(resD, 1, size(placeD, 2));
normD = normD./repmat(max(fiT, [], 1), nbucket, 1);

placeD = placeD./repmat(resD, 1, size(placeD, 2));
placeD = placeD./repmat(sum(placeD, 1), nbucket, 1);

eventD = eventD./repmat(resD, 1, size(placeD, 2));

e = zeros(size(placeD,2),1);
for i = 1:numel(e)
    p = placeD(:,i);
    e(i) = -sum(p(p>0).*log2(p(p>0)));
end

[B,I] = sort(e);



r=5;

xdim = ceil(max(x)-min(x));
ydim = ceil(max(y)-min(y));
xbin = floor(x(isrunning)-min(x))+1;
ybin = floor(y(isrunning)-min(y))+1;
rxbin = ceil(xbin/r);
rybin = ceil(ybin/r);
amap1 = zeros(ydim, xdim, numel(e),'single');
amap2 = zeros(ydim, xdim, numel(e),'single');
a1filt = amap1;
a2filt = amap2;
tmap1 = amap1;
tmap2 = amap2;
t1filt = amap1;
t2filt = amap2;
resmap = zeros(ydim, xdim, 'single');

rmap1 = zeros(ceil(ydim/r), ceil(xdim/r), numel(e),'single');
rmap2 = zeros(ceil(ydim/r), ceil(xdim/r), numel(e),'single');
rresmap = zeros(ceil(ydim/r), ceil(xdim/r),'single');

for j = 1:size(fiT, 1)
    resmap(ybin(j), xbin(j)) = resmap(ybin(j), xbin(j)) + 1;
    rresmap(rybin(j), rxbin(j)) = rresmap(rybin(j), rxbin(j)) + 1;
end

for i = 1:numel(e)
    for j = 1:size(fiT, 1)
        if dir(j) < 0
            amap1(ybin(j), xbin(j), i) = amap1(ybin(j), xbin(j), i) + events(j, i);
            tmap1(ybin(j), xbin(j), i) = tmap1(ybin(j), xbin(j), i) + fiT(j, i);
            rmap1(rybin(j), rxbin(j), i) = rmap1(rybin(j), rxbin(j), i) + events(j,i);
        else
            amap2(ybin(j), xbin(j), i) = amap2(ybin(j), xbin(j), i)+ events(j, i);
            tmap2(ybin(j), xbin(j), i) = tmap2(ybin(j), xbin(j), i)+ fiT(j, i);
            rmap2(rybin(j), rxbin(j), i) = rmap2(rybin(j), rxbin(j), i) + events(j,i);
        end
    end
end





amap1 = amap1./repmat(resmap, 1, 1, numel(e));
amap2 = amap2./repmat(resmap, 1, 1, numel(e));

tmap1 = tmap1./repmat(resmap, 1, 1, numel(e));
tmap2 = tmap2./repmat(resmap, 1, 1, numel(e));


amap1(isnan(amap1)) = 0;
amap2(isnan(amap2)) = 0;

tmap1(isnan(tmap1)) = 0;
tmap2(isnan(tmap2)) = 0;

sigma=10;

nannorm = imgaussfilt(double(resmap>0), sigma, 'FilterSize', 101);

fresmap = imgaussfilt(resmap, sigma, 'FilterSize', 101, 'FilterDomain', 'spatial');

for i = 1:numel(e)
    a1filt(:,:,i) = imgaussfilt(amap1(:,:,i), sigma, 'FilterSize', 101, 'FilterDomain', 'spatial');
    a2filt(:,:,i) = imgaussfilt(amap2(:,:,i), sigma, 'FilterSize', 101, 'FilterDomain', 'spatial');
    t1filt(:,:,i) = imgaussfilt(tmap1(:,:,i), sigma, 'FilterSize', 101, 'FilterDomain', 'spatial');
    t2filt(:,:,i) = imgaussfilt(tmap2(:,:,i), sigma, 'FilterSize', 101, 'FilterDomain', 'spatial');
end

fudge=.001;
a1filt = a1filt./(repmat(nannorm, 1,1,numel(e))+fudge);
a2filt = a2filt./(repmat(nannorm, 1,1,numel(e))+fudge);
t1filt = t1filt./(repmat(nannorm, 1,1,numel(e))+fudge);
t2filt = t2filt./(repmat(nannorm, 1,1,numel(e))+fudge);

a=zeros(1,1,size(fiT,2));
a(1,1,:) = max(fiT,[],1);
t1filt = t1filt./repmat(a,size(t1filt,1),size(t1filt,2),1);
t2filt = t2filt./repmat(a,size(t2filt,1),size(t2filt,2),1);




atot = a1filt+a2filt;
amax = squeeze(max(atot, [], [1,2]));


[m,aI] = sort(amax);
aI = flipud(aI);




%{

for i = 1:numel(I)
    nI = aI(i);
    f = figure('visible', 'off');
    tiledlayout(1,2);
    nexttile;
    imagesc(t1filt(:,:,nI)+t2filt(:,:,nI));
    colorbar;
    axis image;
    title('average activity');
    
    nexttile;
    imagesc(a1filt(:,:,nI)+a2filt(:,:,nI));
    colorbar;
    axis image;
    title('average event rate');
    
    sgtitle(strcat('Neuron',' ',num2str(nI)));
    
    saveas(f,strcat('E:/place_cells/peak/',num2str(i),'_neuron',num2str(nI),'.png'));
    close(f);
    disp(i);
end
%}


% selected place cells
ifig =[363,1467,1805,1136,693,2001,1712,2331,2286,782,351,1190];
iplot=ifig;

plotOut = a1filt+a2filt;
cmax = max(plotOut(:,:,iplot),[],'all');
figure;
t = tiledlayout('flow');
t.TileSpacing = 'compact';
t.Padding = 'compact';
for i = 1:numel(iplot)
    nexttile;
    imagesc(a1filt(:,:,iplot(i))+a2filt(:,:,iplot(i)));
    colorbar;
    colormap parula;
    title(['Neuron',' ',num2str(iplot(i))]);
    axis image;
    axis off;
end
shg;






%cmap = hsv(3);
cmap = [0,0,0; .3, .3, .3; 0,1,1; 1,.3,.3];
%cmap = [0,0,0; .3, .3, .3; 0,1,1];
plotOut = double(sum(S,3)>.1);
plotOut = plotOut + sum(double(S(:,:,ifig)>.1),3);
figure;
imagesc(plotOut);
colormap(cmap);

%imwrite(plotOut+1, cmap, 'place_cell_locations_new.png');





