

[parentdir,~,~] = fileparts(pwd);
load(fullfile(parentdir, 'supporting_data/extract_curated.mat'));
load(fullfile(parentdir, 'supporting_data/position_tracking.mat'));


fiT = iT;



centroids = zeros(size(S,3), 2);
for i = 1:size(S,3)
    [row, col] = find(S(:,:,i)==1);
    centroids(i,1) = row(1);
    centroids(i,2) = col(1);
end

sortD = centroids(:,1);

xlim = [140,700];
ylim = [1,608];

mapscale = 5;


amap = NaN*ones(round((ylim(2)-ylim(1))/mapscale), round((xlim(2)-xlim(1))/mapscale));





Nclass = 8;


Nsamp=40;

holdOut = .5;

tLinear = templateLinear('Lambda', 6.75e-8, 'Learner', 'svm', 'Regularization', 'ridge');


ptile = 50;


figure;
imagesc(sum(S(ylim(1):ylim(2), xlim(1):xlim(2), :), 3));

aplot = figure;

rsize = 50;
nsize = 64;

speed = abs(diff(sunwrap));
speed(end+1)=speed(end);
speed = movmean(speed, 30);
fiT = movmax(fiT,30);
p = prctile(speed(:), ptile);
isrunning = speed>p;

dir = diff(sunwrap);
dir(end+1) = dir(end);
dir = movmean(dir,30);
dir = -1+2*double(dir>0);

dir = dir(isrunning);
fiT = fiT(isrunning, :);
fq = octant(isrunning);
%fq = quadrant(isrunning);

holdoutbuffer = round(size(fiT, 1)*.01);
%holdoutbuffer = 0;




N=round(size(fiT,1)*holdOut);




for row = 1:size(amap, 1)
    for col = 1:size(amap, 2)

        coord = [row*mapscale + ylim(1), col*mapscale + xlim(1)];
        dists = sqrt(sum((centroids-repmat(coord, size(centroids, 1), 1)).^2, 2));
        if sum(dists<=rsize) >= nsize

            [B, Idist] = sort(dists);

            ncluster = Idist(1:nsize);

            fT = fiT(:, ncluster);
            

            htot=zeros(Nclass,Nclass,Nsamp);
            avgcount=zeros(1,Nclass, Nsamp);
            parfor j = 1:Nsamp

                rs = round(rand*size(fT,1));
                fitfT = circshift(fT, rs);
                fitq = circshift(fq, rs+0);
                %rs = round(rand*size(fT,1));
                %fitsclass = circshift(fitsclass, rs);
                mdl=fitcecoc(fitfT(1:N-holdoutbuffer,:),fitq(1:N-holdoutbuffer),'Learners',tLinear,'Coding','onevsall','Verbose',0);
                labels=predict(mdl,fitfT(N+1:end-holdoutbuffer,:));

                bins=linspace(.5,Nclass+.5,Nclass+1);
                [h,xe,ye]=histcounts2(labels,fitq(N+1:end-holdoutbuffer),bins,bins);
                norm=repmat(sum(h,1),Nclass,1);
                avgcount(:,:,j) = double(sum(h,1)>0);
                norm(norm==0)=1;
                htot(:,:,j) = h./norm;

            end

            avgh = sum(htot,3)./repmat(sum(avgcount, 3),Nclass,1);
            amap(row, col) = median(diag(avgh));


            %save('accuracy_map_octant_64neuron_highres.mat', 'amap');
            
            figure(aplot);
            cla;
            im = imagesc(amap);
            colorbar;
            set(im, 'AlphaData', ~isnan(amap));

        end
    end
end

im = imagesc(amap, [0.125,max(amap(~isnan(amap)))]);
colorbar;
set(im, 'AlphaData', ~isnan(amap));

