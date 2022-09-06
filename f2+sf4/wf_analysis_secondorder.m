
% fits a second-order quantized gaussian model to the rois of the image
% stack specified in wf_roi.csv

images=dir('wf/beads/*.png');

T = struct2table(images);
sortedT = sortrows(T, 'date');
images = table2struct(sortedT);

imStack=uint8([]);
imFilt=[];
for i=1:numel(images)
    im=imread(strcat(images(i).folder,'/',images(i).name));
    imStack(:,:,i)=im;
end

T = readtable('wf_roi.csv');

dpixel=5.35;

center=[259,420];

roi=[];
%/(2*pi*x(4)*x(6)*sqrt(1-x(2)^2))
%fun = @(x,xdata) x(1) + x(7)/(x(4)*x(6)*sqrt(1-x(2)^2)) * exp( (-1/(2*(1-(x(2)^2)))) * ( (xdata(:,1)-x(3)).^2/x(4)^2 + (xdata(:,2)-x(5)).^2/x(6)^2 - 2*x(2)*(xdata(:,1)-x(3)).*(xdata(:,2)-x(5))/(x(4)*x(6))));

a = @(x) cos(x(2))^2/(2*x(4)^2) + sin(x(2))^2/(2*x(6)^2);
b = @(x) -sin(2*x(2))/(4*x(4)^2) + sin(2*x(2))/(4*x(6)^2);
c = @(x) sin(x(2))^2/(2*x(4)^2) + cos(x(2))^2/(2*x(6)^2);

gaussExp = @(x,xdata) exp( -( a(x)*(xdata(:,1)-x(3)).^2 + 2*b(x)*(xdata(:,1)-x(3)).*(xdata(:,2)-x(5)) + c(x)*(xdata(:,2)-x(5)).^2 ) );

f0 = @(x,xdata) gaussExp(x,xdata);

fx2 = @(x,xdata) 2 * gaussExp(x,xdata) .* (2 * a(x)^2 * (xdata(:,1)-x(3)).^2 + a(x)*(-1 + 4*b(x)*(xdata(:,1)-x(3)).*(xdata(:,2)-x(5))) + 2*b(x)^2*(xdata(:,1)-x(3)).^2 );

fy2 = @(x,xdata) 2 * gaussExp(x,xdata) .* (2 * c(x)^2 * (xdata(:,2)-x(5)).^2 + c(x)*(-1 + 4*b(x)*(xdata(:,2)-x(5)).*(xdata(:,1)-x(3))) + 2*b(x)^2*(xdata(:,2)-x(5)).^2 );

fun = @(x,xdata) x(1) + x(7) * (f0(x,xdata) + fx2(x,xdata)*dpixel^2/24 + fy2(x,xdata)*dpixel^2/24);


options = optimoptions('lsqcurvefit','Display','off','MaxFunctionEvaluations',100000,'MaxIter',10000,'FunctionTolerance',1e-8);

widths=zeros(size(T,1),2);
r=zeros(size(T,1),1);


amps=zeros(size(T,1),1);
p=zeros(size(T,1),1);


minWidth=0;
parfor i=1:size(T,1)
    roi=imStack(T.BY(i)+1:T.BY(i)+T.Height(i),T.BX(i)+1:T.BX(i)+T.Width(i),T.Slice(i));
    
    imVec=double(roi(:));
    
    hh=(T.Height(i)-1)/2;
    hw=(T.Width(i)-1)/2;
    
    xgrid = repmat(dpixel*(-hw:hw),T.Height(i),1);
    ygrid = repmat(dpixel*(-hh:hh)',1,T.Width(i));
    grid=[xgrid(:),ygrid(:)];
    
    [m, xi] = max(imVec);
    xi = grid(xi,1);
    [m, yi] = max(imVec);
    yi = grid(yi,2);
    
    x0=[median(imVec),0,xi,2,yi,2,max(imVec)];
    lb=[0,-Inf,-hw*dpixel,minWidth,-hh*dpixel,minWidth,0];
    ub=[Inf,Inf,hw*dpixel,Inf,hh*dpixel,Inf,Inf];
    
    f=lsqcurvefit(fun,x0,grid,imVec,lb,ub,options);

    widths(i,:) = 2.355*sort([f(4),f(6)]);
    r(i)=dpixel*sqrt( (center(1)-(T.BY(i)+hh+1)).^2 + (center(2)-(T.BX(i)+hw+1)).^2 );
    amps(i)=f(7);
    p(i)=f(2);
    
    disp(i);
    
    roiWrite=reshape(fun(f,grid),T.Height(i),T.Width(i));
    roiWrite = (roiWrite-min(roiWrite(:)))/(max(roiWrite(:))-min(roiWrite(:)));
    %imwrite(roiWrite,strcat('models/',num2str(i),'.png'));
end


scatter(r,widths(:,1));hold on;scatter(r,widths(:,2));
ylim([0,15]);
xlim([-100,1600]);

%groups=[50,150,250,350,450];
groups=[0,250,500,750,1000,1250,1500];
[m,I]=min(abs(repmat(r,1,numel(groups))-repmat(groups,numel(r),1)),[],2);
grouping=groups(I)';

figure;
scatter(r,T.Slice);
hold on;
roc=50000;
plot(0:1700,sqrt(roc^2-(0:1700).^2)-roc+50);


figure;
boxplot(widths(:,1),grouping,'boxstyle','filled','widths',.2,'Positions',(1:numel(groups))-.1,'color','b','symbol','b.','Whisker',5);
hold on;
boxplot(widths(:,2),grouping,'boxstyle','filled','widths',.2,'Positions',(1:numel(groups))+.1,'color','r','symbol','r.','Whisker',5);
ylim([0,10]);

rVals=groups;
rBins=cell(numel(groups),1);
for i=1:numel(r)
    rBins{find(grouping(i)==groups)}=[rBins{find(grouping(i)==groups)};widths(i,:)];
end

rMean=[];
rErr=[];
for i=1:numel(rVals)
    rMean(i,:)=mean(rBins{i});
    rErr(i,:)=std(rBins{i});
end

figure;
errorbar(rVals,rMean(:,1),2*rErr(:,1));
hold on;
errorbar(rVals,rMean(:,2),2*rErr(:,2));
xlim([-50,1550]);
ylim([0,10]);
