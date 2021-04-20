% Code last updated on 04/19/2021 by Eric Wei

%% Initialize and Define Variables
clear;
x = resampleDicom('06.dcm')

figure(1); clf
xdata1 = imrotate3(x.data(:, :, :, 5), 0, [1 0 0]);
sliceViewer(xdata1, 'SliceDirection', 'X');

figure(2); clf
sliceViewer(xdata1, 'SliceDirection', 'Y');

figure(3); clf
sliceViewer(xdata1, 'SliceDirection', 'Z');

%% Average Slices

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

figure(4); clf

se90 = strel('line',3,90);
se0 = strel('line',3,0);
for i = 1:size(xavgslice,1)
    J = wiener2(squeeze(xavgslice(i,:,:)),[7 7]);
    B = ordfilt2(J,10,true(8));
    % filt_x = imgaussfilt(B, 1);
    % [~,threshold] = edge(filt_x,'sobel');
    % fudgeFactor = 0.2;
    % BG = edge(filt_x,'sobel',threshold * fudgeFactor);
    % imagesc(BG)
    [~,threshold] = edge(B,'Canny');
    fudgeFactor = .8;
    BWs = edge(B,'Canny',threshold * fudgeFactor);
    BWsdil = imdilate(BWs,[se90 se0]);
    colormap gray
    imagesc(BWsdil);
    xlabel('z')
    ylabel('y')
    text(5,5,sprintf('X-slice: %d',i),'color','white');
    pause(0.25);
end
% 33 for x for dataset 6 frame 5

figure(5)

for i = 1:size(yavgslice,1)
    imagesc(squeeze(yavgslice(i,:,:)));
    %  colormap gray
    text(5,5,sprintf('Y-slice: %d',i),'color','white');
    pause(0.25);
end
% 30 for y for dataset 6

% figure(6)
% 
%     for i = 1:size(zavgslice,1)
%  imagesc(squeeze(zavgslice(i,:,:)));
%  colormap gray
%  text(5,5,sprintf('Z-slice: %d',i),'color','white');
%  pause(0.25);
%     end
% 10 for z for dataset 6

figure(7); clf
J = wiener2(squeeze(xavgslice(33,:,:)),[7 7]);
B = ordfilt2(J,10,true(8));
% filt_x = imgaussfilt(B, 1);
% [~,threshold] = edge(filt_x,'sobel');
% fudgeFactor = 0.2;
% BG = edge(filt_x,'sobel',threshold * fudgeFactor);
% imagesc(BG)
imagesc(B)
colormap gray
xlabel('z')
ylabel('y')
[z,y]  = ginput(3)

figure(8); clf
J = wiener2(squeeze(yavgslice(30,:,:)),[7 7]);
B = ordfilt2(J,10,true(8));
imagesc(B)
colormap gray
xlabel('z')
ylabel('x')
[z,x]  = ginput(2)

figure(9); clf
J = wiener2(squeeze(zavgslice(10,:,:)),[7 7]);
B = ordfilt2(J,10,true(8));
imagesc(B)
colormap gray
xlabel('x')
ylabel('y')
[x,y]  = ginput(1)

A = [33*5; 143.328; 44.8318];
AC = [33*5; 113.554; 122.4724];
PC = [33*5; 157.8120; 135.8917];
AA = [149.1283; 30*5; 135.8917];
PA = [197.0729; 30*5; 115.7627];
 
APC = PC-A;
AAC = AC-A;
APA = PA-A;
AAA = AA-A;
 
LongAxis = (APC + AAC + APA + AAA)/4
MagLongAxis = norm(LongAxis)

% A = [28.1731; 150; 172.8532];
% PC = [134.9954; 150; 156.1926];
% AC = [150.0898; 150; 196.8315];
% PA = [177; 124.9932; 130.0587];
% AA = [177; 160.1820; 128.2815];
% 
% APC = PC-A;
% AAC = AC-A;
% APA = PA-A;
% AAA = AA-A;
% 
% LongAxis = (APC + AAC + APA + AAA)/4
% MagLongAxis = norm(LongAxis)

figure(10); clf
line([0 LongAxis(1)], [0 LongAxis(2)], [0 LongAxis(3)]);
view(3)