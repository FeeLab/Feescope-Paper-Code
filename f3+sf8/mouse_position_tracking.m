
% processes behavioral tracking video and outputs position of mouse in
% maze as a function of time.
% correlates timestamps of 1p imaging and behavioral streams, interpolates
% position to yield simultaneous datapoints for decoding analysis

% outputs:
% s: estimated 1d position of mouse body in maze
% sunwrap: unwrapped 1d position of mouse body in maze
% x: x coordinate of mouse body in maze
% y: y coordinate of mouse body in maze
% quadrant: current quadrant occupied by mouse (fig 3)
% octant: current octant occupied by mouse (sfig 8)
% iT: extracted neural activity, truncated to overlap with behavioral data


[parentdir,~,~] = fileparts(pwd);

pos=readmatrix(fullfile(parentdir,'supporting_data/','home_pos-speed-in_2021-07-26T13_50_50.csv'),'Delimiter',' ');

posT=readmatrix(fullfile(parentdir,'supporting_data/','home_pos-speed-in_2021-07-26T13_50_50.csv'),'Delimiter',' ','OutputType','string');
caT=readmatrix(fullfile(parentdir,'supporting_data/','invivo_2021-07-26T13_50_50.csv'),'Delimiter',' ','OutputType','string');

load(fullfile(parentdir, 'supporting_data/extract_output.mat'));
T = exOut.temporal_weights;


caT = str2double(extractBetween(caT,18,23))+60*str2double(extractBetween(caT,15,16))+3600*str2double(extractBetween(caT,12,13));
posT = str2double(extractBetween(posT(:,1),18,23))+60*str2double(extractBetween(posT(:,1),15,16))+3600*str2double(extractBetween(posT(:,1),12,13));


pos = pos(:,2:3);


xmin = 101;
xmax = 378;
ymin = 229;
ymax = 560;

rAngle = -3.21+90; %rotate video by degrees

red = zeros(1,1,3);
%red(1) = 174;
%red(2) = 81;
%red(3) = 67;

red(1) = 0;
red(2) = 3;
red(3) = 4;
redlim = 1600;



v = VideoReader(fullfile(parentdir,'supporting_data/','home_arena_2021-07-26T13_50_50.avi'));

v.CurrentTime = 64; %start at 1'04 to remove period of experimenter futzing

i=1;
vPos = zeros(int16(v.Duration-v.CurrentTime)*v.FrameRate,2);
while hasFrame(v)
    frame = readFrame(v);
    frame = imrotate(frame, rAngle, 'bilinear');
    frame = frame(ymin:ymax,xmin:xmax,:);
    track = double(frame) - repmat(red, size(frame,1), size(frame,2), 1);
    track = track(:,:,1).^2 + track(:,:,2).^2 + track(:,:,3).^2;
    BW = track < redlim;
    s = regionprops('table', BW, 'Area', 'Centroid');
    s = table2array(s);
    [M,I] = max(s(:,1));
    vPos(i,:) = s(I,2:3);
    disp(i);
    i = i+1;
end


isgood = ones(size(T,2),1);

posT = posT(end-size(vPos,1)+1:end);
a = caT>posT(1);
a = diff(a);
iT = T(find(a)+1:end,isgood==1);
icaT = caT(find(a)+1:end);
iPos = interp1(posT, vPos, icaT);



%%
x = iPos(:,1);
y = iPos(:,2);

xmin = 44;
xmax = 228;
ymin = 34;
ymax = 270;
xl = xmax-xmin;
yl = ymax-ymin;

d = [abs(x-xmin), abs(y-ymin), abs(x-xmax), abs(y-ymax)];

[m, quadrant] = min(d, [], 2, 'linear');
[r,quadrant] = ind2sub(size(d), quadrant);

s = zeros(numel(quadrant), 1);


s(quadrant==1) = ymax - y(quadrant==1);
s(quadrant==2) = yl + x(quadrant==2) - xmin;
s(quadrant==3) = yl + xl + y(quadrant==3) - ymin;
s(quadrant==4) = 2*yl + xl + xmax - x(quadrant==4);

%{
jumps = [0; diff(s)];

phase = 0;
for i = 2:numel(s)
    if jumps(i) > xl + yl
        phase = phase - (2*xl + 2*yl);
    elseif jumps(i) < -xl - yl
        phase = phase + (2*xl + 2*yl);
    end
    s(i) = s(i) + phase;
end
%}

lf = .508/290; %50.8 cm per 290 pixel
s = s*lf;

ymid = mean([ymin, ymax]);
xmid = mean([xmin, xmax]);

octant = quadrant;
octant(quadrant==1) = 1 + double(y(quadrant==1) < ymid);
octant(quadrant==2) = 3 + double(x(quadrant==2) > xmid);
octant(quadrant==3) = 5 + double(y(quadrant==3) > ymid);
octant(quadrant==4) = 7 + double(x(quadrant==4) < xmid);


sunwrap = unwrap(s*100, max(s)*100/2)/100;

%save('position_tracking_output.mat', 's', 'sunwrap', 'quadrant', 'octant', 'x', 'y', 'iT');
