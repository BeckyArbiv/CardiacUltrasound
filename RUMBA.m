% Code last updated on 04/27/2021 by:
% Eric Wei, Michelle Dantzler, Becky Arbiv

%% Initialize and Define Variables
clear;
close all
X = resampleDicom('03.dcm');

% Parameters

% 26 for x for dataset 7 frame 1
% 33 for x for dataset 6 frame 5
% 27 for x for dataset 5 frame 1
% 30 for x for dataset 4 frame 1
% 37 for x for dataset 3 frame 1
% 30 for x for dataset 2 frame 1
% 28 for x for dataset 1 frame 1
Xslice_num = 37; % Change as needed

% 25 for y for dataset 7 frame 1
% 30 for y for dataset 6 frame 5
% 27 for y for dataset 5 frame 1
% 22 for y for dataset 4 frame 1
% 29 for y for dataset 3 frame 1
% 28 for y for dataset 2 frame 1
% 25 for y for dataset 1 frame 1
Yslice_num = 29; % Change as needed

% 15 for z for dataset 6 frame 1
% 10 for z for dataset 6 frame 1
% 10 for z for dataset 5 frame 1
% 20 for z for dataset 4 frame 1
% 13 for z for dataset 3 frame 1
% 15 for z for dataset 2 frame 1
% 15 for z for dataset 1 frame 1
Zslice_num = 13; % Change as needed
avgslice = 5; % Change as needed
frameselect = 1; % Single Frame Value for Figures
dataselector = 3; % Single Dataset Value for Figures

% Other Parameters
NumVolumes = size(X.data, 4);

%% Edge tracking lets gooo
xdata1 = imrotate3(X.data(:, :, :, 1), 0, [1 0 0]);
[xavgslice,yavgslice,zavgslice] = avg_slices(xdata1, avgslice);
% SLICE SELECTION COMMENT
%[i] = slice_selection(frameselect, dataselector, xavgslice, yavgslice, zavgslice);  % Comment out when not selecting slices
[Yrow, Ycol, Ydistance,Yz_input,Yx_input,Yindex_x, Yindex_y] = edge_tracking(yavgslice,Yslice_num, dataselector); %Y SLICE
[row, col, distance,z_input,y_input,index_x, index_y] = edge_tracking(xavgslice, Xslice_num, dataselector); %X SLICE
[Zx_input,Zy_input] = short_axis(zavgslice,Zslice_num,dataselector);% Z SLICE

[Longaxis, Magnitude] = longaxis_1(Zx_input, Zy_input, Yz_input, Yx_input, z_input, y_input, Xslice_num, Yslice_num, avgslice)
for t = 2:X.NumVolumes
    xdata1 = imrotate3(X.data(:, :, :, t), 0, [1 0 0]);
    [xavgslice,yavgslice,zavgslice] = avg_slices(xdata1, avgslice);
    [Yindex_x, Yindex_y] = track_over_frames(yavgslice,Yindex_x, Yindex_y, Yslice_num,t, dataselector) %Y SLICE
    [index_x, index_y] = track_over_frames(xavgslice,index_x, index_y, Xslice_num,t, dataselector) %X SLICE
    [longaxis, magnitude] = longaxis_calc(Zx_input, Zy_input, Yz_input, Yx_input, z_input, y_input, index_y, index_x, Yindex_x, Yindex_y, t)
    Longaxis = [Longaxis; longaxis];
    Magnitude = [Magnitude; magnitude];
end
[k] = printLongAxis(Magnitude, NumVolumes);

function [x_input,y_input] = short_axis(zavgslice,slice_num, dataselector)
figure(107); clf
hold on
% find optimal slice number
J = wiener2(squeeze(zavgslice(slice_num,:,:)),[7 7]);
B = ordfilt2(J,10,true(8));
filt_x = imgaussfilt(B, 1);
[~,threshold] = edge(filt_x,'Roberts');
fudgeFactor = 0.2;
BG = edge(filt_x,'Roberts',threshold * fudgeFactor);

se90 = strel('line',3,90);
se0 = strel('line',3,0);
BWsdilX = imdilate(BG,[se90 se0]);

imagesc(BWsdilX)
text(5,5,sprintf('Slice: %d',slice_num),'color','white');
text(5,20,sprintf('Frame: %d',1),'color','white');
text(5,35,sprintf('Dataset: %d',dataselector),'color','white');
xlabel('x')
ylabel('y')
click = ginput(1)
x_input = round(click(:,1));
y_input = round(click(:,2));
plot(x_input, y_input,'or','MarkerSize',7, 'MarkerFaceColor', 'r')
hold off
end
function [row, col, distance,x_input,y_input,index_x, index_y] = edge_tracking(direc_avgslice,slice_num, dataselector)
%X SLICE
figure(1); clf
hold on
J = wiener2(squeeze(direc_avgslice(slice_num,:,:)),[7 7]);
B = ordfilt2(J,10,true(8));
filt_x = imgaussfilt(B, 1);
[~,threshold] = edge(filt_x,'Roberts');
fudgeFactor = 0.2;
BG = edge(filt_x,'Roberts',threshold * fudgeFactor);

se90 = strel('line',3,90);
se0 = strel('line',3,0);
BWsdilX = imdilate(BG,[se90 se0]);

imagesc(BWsdilX)
text(5,5,sprintf('Slice: %d',slice_num),'color','white');
text(5,20,sprintf('Frame: %d',1),'color','white');
text(5,35,sprintf('Dataset: %d',dataselector),'color','white');
xlabel('z')
ylabel('y axis (x slice) and x axis (y slice)')
click = ginput(3)
x_input = round(click(:,1));
y_input = round(click(:,2));
plot(x_input, y_input,'or','MarkerSize',7, 'MarkerFaceColor', 'r')
pause(0.5);
hold off

% define ROI
roiSize = 5;
roiEdge = (roiSize-1)/2;
roiCent = roiEdge+1;
dim=size(BWsdilX,1);
roi = zeros(roiSize,roiSize,length(click));
for i = 1:length(click)
    roi(:,:,i) = flip(int8(BWsdilX([y_input(i)-roiEdge:y_input(i)+roiEdge],[x_input(i)-roiEdge:x_input(i)+roiEdge])));
end

if any(all(all(roi == 1)))
    error("Please select an edge")
elseif any(all(all(roi == 0)))
    error("Please select an edge")
end
row = zeros(1,3);
col = zeros(1,3);
for q = 1:3
    for p = 1:roiEdge
        if any(roi(roiCent, roiCent,q) ~= BWsdilX([y_input(q)-p:y_input(q)+p],[x_input(q)-p:x_input(q)+p]),'all')
            [edge_row, edge_col] = find((roi(roiCent, roiCent,q) ~= BWsdilX([y_input(q)-p:y_input(q)+p],[x_input(q)-p:x_input(q)+p])))
            if p == 1
                roiCent = 2;
            end
            distance = sqrt( (edge_row-roiCent).^2 + (edge_col-roiCent).^2 );
            closest_point = find(distance== min(distance))
            row(q) = edge_row(closest_point(1));
            col(q) = edge_col(closest_point(1));
            break
            %elseif any(roi(roiCent, roiCent) == BWsdil([y_input-p:y_input+p],[x_input-p:x_input+p]),'all')
        end
    end
end
index_y = y_input' + (-row+roiCent);
index_x = x_input' + (col-roiCent);
end
% Average Slices
function [xavgslice,yavgslice,zavgslice] = avg_slices(xdata1, avgslice)

for k = 1:(size(xdata1,1)/avgslice)
    for i = 1:avgslice
        xdata1avg(i,:,:) = xdata1((k-1)*avgslice+i, :, :);
        
    end
    xavgslice(k,:,:) = mean(xdata1avg);
end

xdata1y = permute(xdata1, [2 1 3]);
for k = 1:(size(xdata1y,1)/avgslice)
    for i = 1:avgslice
        ydata1avg(i,:,:) = xdata1y((k-1)*avgslice+i, :, :);
        
    end
    yavgslice(k,:,:) = mean(ydata1avg);
end

xdata1z = permute(xdata1, [3 2 1]);
for k = 1:(size(xdata1z,1)/avgslice)
    for i = 1:avgslice
        zdata1avg(i,:,:) = xdata1z((k-1)*avgslice+i, :, :);
        
    end
    zavgslice(k,:,:) = mean(zdata1avg);
end
end

function [rumba_x, rumba_y] = track_over_frames(xavgslice,index_x, index_y,sliceNum,t, dataselector)
figure(t); clf
hold on
slice = sliceNum; % find optimal slice number
J = wiener2(squeeze(xavgslice(slice,:,:)),[7 7]);
B = ordfilt2(J,10,true(8));
filt_x = imgaussfilt(B, 1);
[~,threshold] = edge(filt_x,'Roberts');
fudgeFactor = 0.2;
BG = edge(filt_x,'Roberts',threshold * fudgeFactor);

se90 = strel('line',3,90);
se0 = strel('line',3,0);
BWsdilX = imdilate(BG,[se90 se0]);

imagesc(BWsdilX)
text(5,5,sprintf('Slice: %d',slice),'color','white');
text(5,20,sprintf('Frame: %d',t),'color','white');
text(5,35,sprintf('Dataset: %d',dataselector),'color','white');
xlabel('z')
ylabel('y')


roiSize = 3;
roiEdge = (roiSize-1)/2;
roiCent = roiEdge+1;

roi = zeros(roiSize,roiSize,length(index_x));
for i = 1:length(index_x)
    roi(:,:,i) = flip(int8(BWsdilX([index_y(i)-roiEdge:index_y(i)+roiEdge],[index_x(i)-roiEdge:index_x(i)+roiEdge])));
end
y_input = index_y;
x_input = index_x;

rumba = 0;
row = zeros(1,3);
col = zeros(1,3);
%plot([x_input-roiEdge:x_input+roiEdge],[y_input-roiEdge:y_input+roiEdge],'*')
for q = 2:3
    rumba = 0;
    while rumba == 0
        if any(roi(:,:,q) ~= roi(roiCent, roiCent,q),'all')
            %if  any(roi(roiCent, roiCent,q) ~= BWsdilX([y_input(q)-roiEdge:y_input(q)+roiEdge],[x_input(q)-roiEdge:x_input(q)+roiEdge]),'all')
            %[edge_row, edge_col] = find((roi(roiCent, roiCent,q) ~= BWsdilX([y_input(q)-roiEdge:y_input(q)+roiEdge],[x_input(q)-roiEdge:x_input(q)+roiEdge])))
            [edge_row, edge_col] = find((roi(roiCent, roiCent,q) ~= roi(:,:,q)))
            distance = sqrt( (edge_row-roiCent).^2 + (edge_col-roiCent).^2 );
            closest_point = find(distance== min(distance))
            row(q) = edge_row(closest_point(1));
            col(q) = edge_col(closest_point(1));
            rumba = 1;
        else
            roiSize = roiSize + 2;
            roiEdge = (roiSize-1)/2;
            roiCent = roiEdge+1;
            roi = zeros(roiSize,roiSize,length(index_x));
            roi(:,:,q) = flip(int8(BWsdilX([y_input(q)-roiEdge:y_input(q)+roiEdge],[x_input(q)-roiEdge:x_input(q)+roiEdge])));
        end
    end
end

rumba_y = y_input + (-row+roiCent)
rumba_x = x_input + (col-roiCent)
rumba_y(1) = y_input(1);
rumba_x(1) = x_input(1);
plot(rumba_x,rumba_y,'g*')
plot(x_input, y_input, 'r*')
hold off
print -dpng Dataset
end

function [LongAxis, MagLongAxis] = longaxis_calc(Zx_input, Zy_input, Yz_input, Yx_input, z_input, y_input, index_y, index_x, Yindex_x, Yindex_y, t)
Xslice_num = 33;
Yslice_num = 30;
Zslice_num = 10;
avgslice = 5;
A = [(Zx_input(1)+Yx_input(1))/2 (y_input(1) + Zy_input(1))/2 (z_input(1) + Yz_input(1))/2];
AC = [Xslice_num*avgslice index_y(2) index_x(2)];
PC = [Xslice_num*avgslice index_y(3) index_x(3)];
AA = [Yindex_y(2) Yslice_num*avgslice Yindex_x(2)];
PA = [Yindex_y(3) Yslice_num*avgslice Yindex_x(3)];

APC = PC-A;
AAC = AC-A;
APA = PA-A;
AAA = AA-A;

LongAxis = (APC + AAC + APA + AAA)/4
MagLongAxis = norm(LongAxis)

end

function [i] = slice_selection(frameselect, dataselector, xavgslice, yavgslice, zavgslice)
% X slice selector/viewer
figure(104); clf

se90 = strel('line',1,90);
se0 = strel('line',1,0);
for i = 1:size(xavgslice,1)
    J = wiener2(squeeze(xavgslice(i,:,:)),[7 7]);
    B = ordfilt2(J,10,true(8));
    [~,threshold] = edge(B,'Canny');
    fudgeFactor = .8;
    BWs = edge(B,'Canny',threshold * fudgeFactor);
    BWsdil = imdilate(BWs,[se90 se0]);
    colormap gray
    imagesc(B);
    xlabel('z')
    ylabel('y')
    text(5,5,sprintf('X-slice: %d',i),'color','white');
    text(5,20,sprintf('Frame: %d',frameselect),'color','white');
    text(5,35,sprintf('Dataset: %d',dataselector),'color','white');
    frame = getframe(figure(104));
    gif{i}= frame2im(frame);
    pause(0.25);
end
% Write first frame separately, all others will be appended
delay = 1.25; % The delay between successive frames (in seconds)
[A, map] = rgb2ind(gif{1},256);
imwrite(A, map, 'Xsliceselector.gif', 'gif', 'LoopCount', Inf, 'DelayTime', delay);
for i = 2:size(xavgslice,1)
    [A, map] = rgb2ind(gif{i},256);
    imwrite(A, map, 'Xsliceselector.gif', 'gif', 'WriteMode', 'append', 'DelayTime', delay);
end
% 33 for x for dataset 6 frame 5
% 27 for x for dataset 5 frame 1

% Y slice selector/viewer
figure(105)

for i = 1:size(yavgslice,1)
    J = wiener2(squeeze(yavgslice(i,:,:)),[7 7]);
    B = ordfilt2(J,10,true(8));
    [~,threshold] = edge(B,'Canny');
    fudgeFactor = .8;
    BWs = edge(B,'Canny',threshold * fudgeFactor);
    BWsdil = imdilate(BWs,[se90 se0]);
    colormap gray
    imagesc(B);
    xlabel('z')
    ylabel('x')
    text(5,5,sprintf('Y-slice: %d',i),'color','white');
    text(5,20,sprintf('Frame: %d',frameselect),'color','white');
    text(5,35,sprintf('Dataset: %d',dataselector),'color','white');
    frame = getframe(figure(105));
    gif{i}= frame2im(frame);
    pause(0.25);
end
% Write first frame separately, all others will be appended
delay = 1.25; % The delay between successive frames (in seconds)
[A, map] = rgb2ind(gif{1},256);
imwrite(A, map, 'Ysliceselector.gif', 'gif', 'LoopCount', Inf, 'DelayTime', delay);
for i = 2:size(yavgslice,1)
    [A, map] = rgb2ind(gif{i},256);
    imwrite(A, map, 'Ysliceselector.gif', 'gif', 'WriteMode', 'append', 'DelayTime', delay);
end
% 30 for y for dataset 6 frame 5
% 27 for y for dataset 5 frame 1

% Z slice selector/viewer
figure(106)

for i = 1:size(zavgslice,1)
    J = wiener2(squeeze(zavgslice(i,:,:)),[7 7]);
    B = ordfilt2(J,10,true(8));
    [~,threshold] = edge(B,'Canny');
    fudgeFactor = .8;
    BWs = edge(B,'Canny',threshold * fudgeFactor);
    BWsdil = imdilate(BWs,[se90 se0]);
    colormap gray
    imagesc(B);
    xlabel('x')
    ylabel('y')
    colormap gray
    text(5,5,sprintf('Z-slice: %d',i),'color','white');
    text(5,20,sprintf('Frame: %d',frameselect),'color','white');
    text(5,35,sprintf('Dataset: %d',dataselector),'color','white');
    frame = getframe(figure(106));
    gif{i}= frame2im(frame);
    pause(0.25);
end
% Write first frame separately, all others will be appended
delay = 1.25; % The delay between successive frames (in seconds)
[A, map] = rgb2ind(gif{1},256);
imwrite(A, map, 'Zsliceselector.gif', 'gif', 'LoopCount', Inf, 'DelayTime', delay);
for i = 2:size(zavgslice,1)
    [A, map] = rgb2ind(gif{i},256);
    imwrite(A, map, 'Zsliceselector.gif', 'gif', 'WriteMode', 'append', 'DelayTime', delay);
end
% 10 for z for dataset 6 frame 1

end

%Long axis 1
function [Longaxis, Magnitude] = longaxis_1(Zx_input, Zy_input, Yz_input, Yx_input, z_input, y_input, Xslice_num, Yslice_num, avgslice)
A = [(Zx_input(1)+Yx_input(1))/2 (y_input(1) + Zy_input(1))/2 (z_input(1) + Yz_input(1))/2];
AC = [Xslice_num*avgslice y_input(2) z_input(2)];
PC = [Xslice_num*avgslice y_input(3) z_input(3)];
AA = [Yx_input(2) Yslice_num*avgslice Yz_input(2)];
PA = [Yx_input(3) Yslice_num*avgslice Yz_input(3)];

APC = PC-A;
AAC = AC-A;
APA = PA-A;
AAA = AA-A;

LongAxis = (APC + AAC + APA + AAA)/4;
MagLongAxis = norm(LongAxis);
Longaxis = [LongAxis];
Magnitude = [MagLongAxis];
end

function [k] = printLongAxis(Magnitude, NumVolumes)

% Fractional shortening
MaxMag = max(Magnitude);
for k = 1:NumVolumes
    FractMag(k) = Magnitude(k)./MaxMag;
end

figure(100); clf
plot(1:NumVolumes, Magnitude, 'b', 'LineWidth', 2)
title('Magnitude of Long-Axis')
xlabel('Frame #')
ylabel('Long-Axis Magnitude (pixels)')
axis([0 NumVolumes min(Magnitude) max(Magnitude)])
print -dpng DataSet3_Magnitude

figure(101); clf
plot(1:NumVolumes, FractMag(1:NumVolumes), 'b', 'LineWidth', 2)
title('Fractional Shortening of Long-Axis')
xlabel('Frame #')
ylabel('Percent Long-Axis')
axis([0 NumVolumes min(FractMag) 1])
print -dpng DataSet3_Fraction

end