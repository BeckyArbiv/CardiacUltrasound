% Code last updated on 03/28/2021 by Eric Wei

%% Initialize and Define Variables
clear;
% x = resampleDicom('cylinder.dcm')
file = ('/Users/beckyarbiv/Documents/BME/BME 543/06.dcm');
x = resampleDicom(file)
% Frame 1
figure(1); clf
xdata1 = imrotate3(x.data(:, :, :, 1), 0, [0 0 1]);
%a = [];
% for i = 1:length(xdata1)
%     b = xdata1(i,:,:);
%     b = reshape(b, [276,208]);
%     image(b);pause(.1);clf;
%     a =  [a; b];
%     
% end

for k = 1:(size(xdata1,1)/5)
    for i = 1:5
        xdata1avg(i,:,:) = xdata1((k-1)*5+i, :, :);
        
    end
    xavgslice(k,:,:) = mean(xdata1avg);
end

for k = 1:(size(xdata1,2)/5)
    for i = 1:5
        ydata1avg(:,i,:) = xdata1(:,(k-1)*5+i, :);
        
    end
    yavgslice(:,k,:) = mean(reshape(ydata1avg, [5,299,208]));
end
yavgslice = reshape(yavgslice, [55,299,208]);
s = sliceViewer(xdata1, 'SliceDirection', 'X');
hAx = getAxesHandle(s);
s.SliceNumber = 105;
I1 = getframe(hAx);
[indI1,cm1] = rgb2ind(I1.cdata,256);

% Frame 2
figure(2); clf
xdata2 = imrotate3(x.data(:, :, :, 2), 10, [0 0 1]);
s = sliceViewer(xdata2, 'SliceDirection', 'X');
hAx = getAxesHandle(s);
s.SliceNumber = 105;
I2 = getframe(hAx);
[indI2,cm2] = rgb2ind(I2.cdata,256);

% Frame 3
figure(3); clf
xdata3 = imrotate3(x.data(:, :, :, 3), 20, [0 0 1]);
s = sliceViewer(xdata3, 'SliceDirection', 'X');
hAx = getAxesHandle(s);
s.SliceNumber = 105;
I3 = getframe(hAx);
[indI3,cm3] = rgb2ind(I3.cdata,256);

% Frame 4
figure(4); clf
xdata4 = imrotate3(x.data(:, :, :, 4), 30, [0 0 1]);
s = sliceViewer(xdata4, 'SliceDirection', 'X');
hAx = getAxesHandle(s);
s.SliceNumber = 105;
I4 = getframe(hAx);
[indI4,cm4] = rgb2ind(I4.cdata,256);

% Frame 5
figure(5); clf
xdata5 = imrotate3(x.data(:, :, :, 5), 40, [0 0 1]);
s = sliceViewer(xdata5, 'SliceDirection', 'X');
hAx = getAxesHandle(s);
s.SliceNumber = 105;
I5 = getframe(hAx);
[indI5,cm5] = rgb2ind(I5.cdata,256);

% Frame 6
figure(6); clf
xdata6 = imrotate3(x.data(:, :, :, 6), 45, [0 0 1]);
s = sliceViewer(xdata6, 'SliceDirection', 'X');
hAx = getAxesHandle(s);
s.SliceNumber = 105;
I6 = getframe(hAx);
[indI6,cm6] = rgb2ind(I6.cdata,256);

% Frame 7
figure(7); clf
xdata7 = imrotate3(x.data(:, :, :, 7), 50, [0 0 1]);
s = sliceViewer(xdata7, 'SliceDirection', 'X');
hAx = getAxesHandle(s);
s.SliceNumber = 105;
I7 = getframe(hAx);
[indI7,cm7] = rgb2ind(I7.cdata,256);

% Frame 8
figure(8); clf
xdata8 = imrotate3(x.data(:, :, :, 8), 60, [0 0 1]);
s = sliceViewer(xdata8, 'SliceDirection', 'X');
hAx = getAxesHandle(s);
s.SliceNumber = 105;
I8 = getframe(hAx);
[indI8,cm8] = rgb2ind(I8.cdata,256);

% Frame 9
figure(9); clf
xdata9 = imrotate3(x.data(:, :, :, 9), 70, [0 0 1]);
s = sliceViewer(xdata9, 'SliceDirection', 'X');
hAx = getAxesHandle(s);
s.SliceNumber = 105;
I9 = getframe(hAx);
[indI9,cm9] = rgb2ind(I9.cdata,256);

% Frame 10
figure(10); clf
xdata10 = imrotate3(x.data(:, :, :, 10), 80, [0 0 1]);
s = sliceViewer(xdata10, 'SliceDirection', 'X');
hAx = getAxesHandle(s);
s.SliceNumber = 105;
I10 = getframe(hAx);
[indI10,cm10] = rgb2ind(I10.cdata,256);

% Frame 11
figure(11); clf
xdata11 = imrotate3(x.data(:, :, :,11), 90, [0 0 1]);
s = sliceViewer(xdata11, 'SliceDirection', 'X');
hAx = getAxesHandle(s);
s.SliceNumber = 105;
I11 = getframe(hAx);
[indI11,cm11] = rgb2ind(I11.cdata,256);


figure(12); clf
subplot(3, 4, 1)
imshow(indI1,cm1)
title('Frame 1');
subplot(3, 4, 2)
imshow(indI2,cm2)
title('Frame 2');
subplot(3, 4, 3)
imshow(indI3,cm3)
title('Frame 3');
subplot(3, 4, 4)
imshow(indI4,cm4)
title('Frame 4');
subplot(3, 4, 5)
imshow(indI5,cm5)
title('Frame 5');
subplot(3, 4, 6)
imshow(indI6,cm6)
title('Frame 6');
subplot(3, 4, 7)
imshow(indI7,cm7)
title('Frame 7');
subplot(3, 4, 8)
imshow(indI8,cm8)
title('Frame 8');
subplot(3, 4, 9)
imshow(indI9,cm9)
title('Frame 9');
subplot(3, 4, 10)
imshow(indI10,cm10)
title('Frame 10');
subplot(3, 4, 11)
imshow(indI11,cm11)
title('Frame 11');
sgtitle('Dataset 2: Frames 1-11 YZ Slice')
print -dpng Tutorial4

figure(13); clf
grid on
volumeViewer(xdata8);
axis on