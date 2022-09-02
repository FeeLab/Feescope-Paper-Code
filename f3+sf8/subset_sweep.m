
% computes the average decoding matrix for Nsubsamp independent randomly
% chosen subpopulations of footprints, and plots the distribution of median
% decoding accuracies for each sampled decoder size (subN)

[parentdir,~,~] = fileparts(pwd);
load(fullfile(parentdir, 'supporting_data/extract_curated.mat'));
load(fullfile(parentdir, 'supporting_data/position_tracking.mat'));

Nclass = 8;

Nsamp=100;

holdOut = .5;

tLinear = templateLinear('Lambda', 6.75e-8, 'Learner', 'svm', 'Regularization', 'ridge');

fiT = iT;

%{
for i = 1:size(fiT, 2)
    fiT(:,i) = circshift(fiT(:,i), randi(size(fiT, 1)));
end
%}


ptile = 50;

Nsubsamp = 64;
subN = 2.^(1:11);
Nsubsizes = numel(subN);

accuracies_raw = zeros(Nsubsizes, Nsubsamp, Nclass);
accuracies_shuffled = zeros(Nsubsizes, Nsubsamp, Nclass);


aplot = figure;
for l = 1:Nsubsizes
    
    
    parfor k = 1:Nsubsamp
        

        speed = abs(diff(sunwrap));
        speed(end+1)=speed(end);
        speed = movmean(speed, 30);
        fT = movmax(iT,30);
        p = prctile(speed(:), ptile);
        isrunning = speed>p;
        
        dir = diff(sunwrap);
        dir(end+1) = dir(end);
        dir = movmean(dir,30);
        dir = -1+2*double(dir>0);
        
        dir = dir(isrunning);
        fT = fT(isrunning, :);
        %fq = octant(isrunning);
        fq = quadrant(isrunning);
        
        
        fT = fT(:, randsample(size(fT,2), subN(l)));
        
        holdoutbuffer = round(size(fT, 1)*.01);
        %holdoutbuffer = 0;
        
        fq(fq==4) = 7 + double(dir(fq==4)>0);
        fq(fq==3) = 5 + double(dir(fq==3)>0);
        fq(fq==2) = 3 + double(dir(fq==2)>0);
        fq(fq==1) = 1 + double(dir(fq==1)>0);
        
        
        %RANDOM SHUFFLE
        %fq = circshift(fq, randi(numel(fq)));

        %{
        for i = 1:size(fT, 2)
            fT(:, i) = circshift(fT(:,i), randi(size(fT, 1)));
        end
        %}
        
        %fT=movmax(fiT,1+1*(l-1));
        N=round(size(fT,1)*holdOut);

        htot=zeros(Nclass,Nclass,Nsamp);
        avgcount=zeros(1,Nclass);
        for j = 1:Nsamp
            
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
            avgcount = avgcount+double(sum(h,1)>0);
            norm(norm==0)=1;
            htot(:,:,j) = h./norm;
            
        end
        
        avgh = sum(htot,3)./repmat(avgcount,Nclass,1);
        accuracies_raw(l,k,:) = diag(avgh);
        

        %random shuffle
        for i = 1:size(fT, 2)
            fT(:, i) = circshift(fT(:,i), randi(size(fT, 1)));
        end
        
        htot=zeros(Nclass,Nclass,Nsamp);
        avgcount=zeros(1,Nclass);
        for j = 1:Nsamp
            
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
            avgcount = avgcount+double(sum(h,1)>0);
            norm(norm==0)=1;
            htot(:,:,j) = h./norm;
            
        end
        avgh = sum(htot,3)./repmat(avgcount,Nclass,1);
        accuracies_shuffled(l,k,:) = diag(avgh);
        disp(k);
    end
    

    %{
    figure(aplot);
    cla;
    plotOut = median(accuracies(1:l, :, :), 3);
    boxplot(plotOut', 'Labels', subN(1:l), 'Symbol', '.');
    %}

    
    figure(aplot);
    cla;
    plotOut = median(accuracies_raw(1:l, :, :), 3);
    boxplot(plotOut', 'Labels', subN(1:l), 'Symbol', '.', 'boxstyle','filled','widths',.2,'Positions',(1:l)-.1,'color','b','symbol','b.','Whisker',3);
    hold on;
    plotOut = median(accuracies_shuffled(1:l, :, :), 3);
    boxplot(plotOut', 'Labels', subN(1:l), 'Symbol', '.', 'boxstyle','filled','widths',.2,'Positions',(1:l)+.1,'color','r','symbol','b.','Whisker',3);
    ylim([0,0.5]);








end

%save('accuracies_quadrant_holdout_shuffled.mat', 'subN', 'accuracies');

%{
figure;
for i = 1:numel(subN)
    swarmchart(categorical(subN(i)*ones(Nsubsamp,1)), median(accuracies(i,:,:), 3));
    hold on;
end
%}
