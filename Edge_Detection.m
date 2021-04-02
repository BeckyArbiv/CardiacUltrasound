% Code last updated on 03/28/2021 by Eric Wei

%% Initialize and Define Variables
clear;
file = ('/Users/beckyarbiv/Documents/BME/BME 543/06.dcm');
x = resampleDicom(file)

figure(1); clf;
xdata1 = imrotate3(x.data(:, :, :, 5), 0, [1 0 0]);
sliceViewer(xdata1, 'SliceDirection', 'X');
title('Unfiltered X')
%[x,y]  = ginput(3)


figure(5);clf;
imshow(xdata1(:,:,100))
title('Raw Image')

figure(8);clf;
J = wiener2(xdata1(:,:,100),[7 7]);
%J = imadjust(J);
B = ordfilt2(J,10,true(8));
imshow(B)
title('2D Wiener Noise-Removal Filter')

figure(7);clf;
J = wiener2(xdata1(:,:,100),[7 7]);
J = imadjust(J);
B = ordfilt2(J,14,true(8));
imshow(B)
title('2D Wiener Noise-Removal Filter With Contrast Enhanced')



%% Binary Gradient Mask
figure(4); clf;
filt_x = imgaussfilt(J, 1);
[~,threshold] = edge(filt_x,'sobel');
fudgeFactor = 0.5;
BWs = edge(filt_x,'sobel',threshold * fudgeFactor);
imshow(BWs)
title('Binary Gradient Mask')

%% Bilated Gradient Mask
se90 = strel('line',3,90);
se0 = strel('line',3,0);

BWsdil = imdilate(BWs,[se90 se0]);
figure(7);clf;
imshow(BWsdil)
title('Sobel Dilated Gradient Mask')

%% Prewitt
figure(8); clf;
[~,threshold] = edge(filt_x,'Prewitt');
fudgeFactor = 1;
BWs = edge(filt_x,'Prewitt',threshold * fudgeFactor);
BWsdil = imdilate(BWs,[se90 se0]);
imshow(BWsdil)
title('Prewitt Dilated Gradient Mask')

%% Roberts
figure(9); clf;
[~,threshold] = edge(filt_x,'Roberts');
fudgeFactor = 1;
BWs = edge(filt_x,'Roberts',threshold * fudgeFactor);
BWsdil = imdilate(BWs,[se90 se0]);
imshow(BWsdil)
title('Roberts Dilated Gradient Mask')

%% Canny
figure(10); clf;
[~,threshold] = edge(filt_x,'Canny');
fudgeFactor = 5;
BWs = edge(filt_x,'Canny',threshold * fudgeFactor);
BWsdil = imdilate(BWs,[se90 se0]);
imshow(BWsdil)
title('Canny Dilated Gradient Mask')

%sliceViewer(BWs, 'SliceDirection', 'X');
%title('Gaussian Filter X')

% figure(2); clf
% sliceViewer(xdata1, 'SliceDirection', 'Y');
% 
% figure(3); clf
% sliceViewer(xdata1, 'SliceDirection', 'Z');
%[x,y]  = ginput(2)

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