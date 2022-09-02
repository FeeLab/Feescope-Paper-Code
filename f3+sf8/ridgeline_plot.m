
% generates ridgeline plot of kiloscope data in fig 3 using ordering from
% NMF extraction output (orderT)

% note the plot will take a few minutes to render on most personal machines



[parentdir,~,~] = fileparts(pwd);
load(fullfile(parentdir, 'supporting_data/extract_curated.mat'));
load(fullfile(parentdir, 'supporting_data/NMF_ordering.mat'));


set(0, 'DefaultFigureRenderer', 'OpenGL');

dr=10;
delta = 0.01;
x = (0:(size(iT, 1)/dr+2))/30;

figure;
for i = 1:numel(orderT)
    y = iT(:,orderT(i));
    y = decimate(double(y), dr);
    y = [0; y; 0];
    y = y+i*delta;
    ax = gca;
    h = fill(ax,x,y,'b');
    %h = plot(T(:,i)+j*delta, 'color', 'b', 'LineWidth', 0.01);
    %h.Color(4) = 0.01;
    set(h, 'facealpha', 0.5);
    set(h, 'edgealpha', 0);
    hold on;
end
ax = gca;

xlim([0, (size(iT, 1)/dr+2)/30]);
ylim([0,delta*size(iT,2)*1.01]);

%ax.YTickLabel = 0:200:2200;
shg