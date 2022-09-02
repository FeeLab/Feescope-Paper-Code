
% constructs a series of decoders for octant to generate
% the decoding matrix in fig 3.
% plots the average decoding matrix and updates as more decoders are
% trained



[parentdir,~,~] = fileparts(pwd);
load(fullfile(parentdir, 'supporting_data/extract_curated.mat'));
load(fullfile(parentdir, 'supporting_data/position_tracking.mat'));


Nclass = 8;

Nsamp=200;
figure;

holdOut = .5;

tLinear = templateLinear('Lambda', 6.75e-8, 'Learner', 'svm', 'Regularization', 'ridge');

fiT = iT;

%{
for i = 1:size(fiT, 2)
    fiT(:,i) = circshift(fiT(:,i), randi(size(fiT, 1)));
end
%}



ptile = [50];

plabels = {'1', '2', '3', '4','5', '6','7', '8'};

ssmooth = [30];
tsmooth = [30];



speed = abs(diff(sunwrap));
speed(end+1)=speed(end);
speed = movmean(speed, ssmooth);
fT = movmax(fiT,tsmooth);
p = prctile(speed(:), ptile);
isrunning = speed>p;

dir = diff(sunwrap);
dir(end+1) = dir(end);
dir = movmean(dir,30);
dir = -1+2*double(dir>0);

dir = dir(isrunning);
fT = fT(isrunning, :);
fq = octant(isrunning);

holdoutbuffer = round(size(fT, 1)*.01);


N=round(size(fT,1)*holdOut);

htot=zeros(Nclass,Nclass,Nsamp);
avgcount=zeros(1,Nclass);
nexttile;
for j = 1:Nsamp
    
    %shuffle direction
    %{
    fq = quadrant(isrunning);
    dir = circshift(dir, randi(numel(dir)));
    fq(fq==4) = 7 + double(dir(fq==4)>0);
    fq(fq==3) = 5 + double(dir(fq==3)>0);
    fq(fq==2) = 3 + double(dir(fq==2)>0);
    fq(fq==1) = 1 + double(dir(fq==1)>0);
    %}
                
    rs = round(rand*size(fT,1));
    fitfT = circshift(fT, rs);
    fitq = circshift(fq, rs+0);
    %rs = round(rand*size(fT,1));
    %fitsclass = circshift(fitsclass, rs);
    
    mdl=fitcecoc(fitfT(1:N-holdoutbuffer,:),fitq(1:N-holdoutbuffer),'Learners',tLinear,'Coding','onevsall','Verbose',0,'Options', statset('UseParallel', true));
    labels=predict(mdl,fitfT(N+1:end-holdoutbuffer,:));

    bins=linspace(.5,Nclass+.5,Nclass+1);
    [h,xe,ye]=histcounts2(labels,fitq(N+1:end-holdoutbuffer),bins,bins);
    norm=repmat(sum(h,1),Nclass,1);
    avgcount = avgcount+double(sum(h,1)>0);
    norm(norm==0)=1;
    htot(:,:,j) = h./norm;
    disp(j);
    cla
    
    
    plotOut = flipud(sum(htot,3)./repmat(avgcount, Nclass, 1));
    
    maxVal = max((max(plotOut(:))-1/Nclass), -(min(plotOut(:)-1/Nclass)));
    imagesc(plotOut, [-maxVal, maxVal]+1/Nclass);
    
    xticks(1:Nclass);
    yticks(1:Nclass);
    set(gca, 'XTick', 1:Nclass, 'XTickLabel', plabels);
    set(gca, 'YTick', 1:Nclass, 'YTickLabel', fliplr(plabels));
    colormap(redbluecmap(256));
    cbar = colorbar;
    cbar.Limits = [0, maxVal+1/Nclass];
    
    
    shg
end



xlabel('actual state');
ylabel('predicted state');
title('Octant prediction, 50% holdout');



%fitcecoc(fT, fq, 'Learners', templateSVM('KernelFunction', 'polynomial'), 'OptimizeHyperparameters', 'all', 'HyperparameterOptimizationOptions', struct('UseParallel', true, 'Repartition', true, 'Holdout', 0.5));
