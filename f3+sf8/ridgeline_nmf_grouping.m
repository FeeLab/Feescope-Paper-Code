
% factors extracted 1p calcium activity into Nc components and uses the
% factorization to order the traces for the ridgeline plot in fig 3.
% for the ordering, each footprint is assigned to the factor that it has
% the largest component in

% outputs:
% Wtot: timeseries weights of each NMF factor
% Htot: footprint weights of each NMF factor
% groupings: factor number assigned to each footprint
% orderT: plot order for ridgeline plot


[parentdir,~,~] = fileparts(pwd);
load(fullfile(parentdir, 'supporting_data/extract_curated.mat'));

Nc = 8;


[Wtot,Htot] = nnmf(iT, Nc, 'Algorithm', 'mult', 'Replicates', 4, 'Options', statset('Display', 'final', 'TolX', 1e-7, 'TolFun', 1e-7, 'UseParallel', true));
[B,I] = sort(Htot');


groupings = zeros(size(iT, 2),1);
for i = 1:numel(groupings)
    [r,c] = find(I==i);
    [gm, gi] = max(r);
    groupings(i) = gi;
end

orderT = [];
for i = 1:Nc
    k = find(groupings==i);
    add = k(randperm(numel(k)));
    orderT = [orderT; add];
end
