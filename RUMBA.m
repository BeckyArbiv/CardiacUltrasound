% Code last updated on 04/19/2021 by Eric Wei

%% Initialize and Define Variables
clear;
close all
X = resampleDicom('06.dcm');

Xslice_num = 33;
Yslice_num = 30;
Zslice_num = 10;
avgslice = 5;
frameselect = 1; % Change as needed
dataselector = 6; % Change as needed

%% Edge tracking lets gooo
xdata1 = imrotate3(X.data(:, :, :, 1), 0, [1 0 0]);
[xavgslice,yavgslice,zavgslice] = avg_slices(xdata1);
[i] = slice_selection(frameselect, dataselector, xavgslice, yavgslice, zavgslice);
[row, col, distance,z_input,y_input,index_x, index_y] = edge_tracking(xavgslice, Xslice_num); %X SLICE
[Yrow, Ycol, Ydistance,Yz_input,Yx_input,Yindex_x, Yindex_y] = edge_tracking(yavgslice,Yslice_num); %Y SLICE
[Zx_input,Zy_input] = short_axis(zavgslice,Zslice_num);

[Longaxis, Magnitude] = longaxis_1(Zx_input, Zy_input, Yz_input, Yx_input, z_input, y_input, Xslice_num, Yslice_num, avgslice)
for t = 2:X.NumVolumes
xdata1 = imrotate3(X.data(:, :, :, t), 0, [1 0 0]);
[xavgslice,yavgslice,zavgslice] = avg_slices(xdata1);
[index_x, index_y] = track_over_frames(xavgslice,index_x, index_y,t) %X SLICE
[Yindex_x, Yindex_y] = track_over_frames(yavgslice,Yindex_x, Yindex_y,t) %Y SLICE
[longaxis, magnitude] = longaxis_calc(Zx_input, Zy_input, Yz_input, Yx_input, z_input, y_input, index_y, index_x, Yindex_x, Yindex_y, t)
Longaxis = [Longaxis; longaxis];
Magnitude = [Magnitude; magnitude];
end
[k] = printLongAxis(Magnitude)

function [x_input,y_input] = short_axis(zavgslice,slice_num)
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
xlabel('x')
ylabel('y')
click = ginput(1)
x_input = round(click(:,1));
y_input = round(click(:,2));
plot(x_input, y_input,'or','MarkerSize',7, 'MarkerFaceColor', 'r')
hold off
end
function [row, col, distance,x_input,y_input,index_x, index_y] = edge_tracking(xavgslice,slice_num)
%X SLICE
figure(7); clf
hold on
 % find optimal slice number
J = wiener2(squeeze(xavgslice(slice_num,:,:)),[7 7]);
B = ordfilt2(J,10,true(8));
filt_x = imgaussfilt(B, 1);
[~,threshold] = edge(filt_x,'Roberts');
fudgeFactor = 0.2;
BG = edge(filt_x,'Roberts',threshold * fudgeFactor);

se90 = strel('line',3,90);
se0 = strel('line',3,0);
BWsdilX = imdilate(BG,[se90 se0]);

imagesc(BWsdilX)
% xlabel('z')
% ylabel('y')
click = ginput(3)
x_input = round(click(:,1));
y_input = round(click(:,2));
plot(x_input, y_input,'or','MarkerSize',7, 'MarkerFaceColor', 'r')
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

% figure(7)
% plot([x_input-roiEdge:x_input+roiEdge],[y_input-roiEdge:y_input+roiEdge],'*')


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
function [xavgslice,yavgslice,zavgslice] = avg_slices(xdata1)
avgslice = 5;

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

function [rumba_x, rumba_y] = track_over_frames(xavgslice,index_x, index_y,g)
figure(g); clf
hold on
xslice = 33; % find optimal slice number
J = wiener2(squeeze(xavgslice(xslice,:,:)),[7 7]);
B = ordfilt2(J,10,true(8));
filt_x = imgaussfilt(B, 1);
[~,threshold] = edge(filt_x,'Roberts');
fudgeFactor = 0.2;
BG = edge(filt_x,'Roberts',threshold * fudgeFactor);

se90 = strel('line',3,90);
se0 = strel('line',3,0);
BWsdilX = imdilate(BG,[se90 se0]);

imagesc(BWsdilX)
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

se90 = strel('line',3,90);
se0 = strel('line',3,0);
for i = 1:size(xavgslice,1)
    J = wiener2(squeeze(xavgslice(i,:,:)),[7 7]);
    B = ordfilt2(J,10,true(8));
    [~,threshold] = edge(B,'Canny');
    fudgeFactor = .8;
    BWs = edge(B,'Canny',threshold * fudgeFactor);
    BWsdil = imdilate(BWs,[se90 se0]);
    colormap gray
    imagesc(BWsdil);
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
    imagesc(BWsdil);
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
for i = 2:size(xavgslice,1)
    [A, map] = rgb2ind(gif{i},256);
    imwrite(A, map, 'Ysliceselector.gif', 'gif', 'WriteMode', 'append', 'DelayTime', delay);
end
% 30 for y for dataset 6

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
for i = 2:size(xavgslice,1)
    [A, map] = rgb2ind(gif{i},256);
    imwrite(A, map, 'Zsliceselector.gif', 'gif', 'WriteMode', 'append', 'DelayTime', delay);
end
% 10 for z for dataset 6

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

function [k] = printLongAxis(Magnitude)
figure(100); clf
plot(1:X.NumVolumes, Magnitude)
print -dpng Dataset6_FullyAutomaticLongAxis

% Fractional shortening
MaxMag = max(Magnitude);
for k = 1:X.NumVolumes 
    FractMag(k) = Magnitude(k)./MaxMag; 
end
figure(101); clf
plot(1:X.NumVolumes, FractMag, 'b', 'LineWidth', 2)
title('Fractional Shortening of Long-Axis')
xlabel('Frame #')
ylabel('Percent Long-Axis')
print -dpng DataSet6_Fraction

end