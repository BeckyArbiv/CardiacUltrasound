% Code last updated on 03/28/2021 by Eric Wei

%% Initialize and Define Variables
clear;
file = ('/Users/beckyarbiv/Documents/BME/BME 543/06.dcm');
x = resampleDicom(file)

xdata1 = imrotate3(x.data(:, :, :, 1), 0, [0 0 1]);
s = sliceViewer(xdata1, 'SliceDirection', 'Y');
hAx = getAxesHandle(s);
s.SliceNumber = 155;
I1 = getframe(hAx);
[indI1,cm1] = rgb2ind(I1.cdata,256);
figure(1);clf;
imshow(indI1',cm1)

J = squeeze(x.data(:,:,:,1)); %Data at NFrame = 1

% % for frames 35-60
% a = 1;
% for i = 35:60
%     A(:,:,a) = J(:,:,i);
%     a = a+1;
% end
% imgA = mean(A,3)/(max(max(mean(A,3))));
% figure(1); clf; 
% imshow(imgA');
% % for frames 61-86
% b = 1;
% for i = 61:86
%     B(:,:,b) = J(:,:,i);
%     b = b+1;
% end
% imgB = mean(B,3)/(max(max(mean(B,3))));
% figure(2); clf; 
% imshow(imgB');
% % for frames 87-112
% c = 1;
% for i = 87:112
%     C(:,:,c) = J(:,:,i);
%     c = c+1;
% end
% imgC = mean(C,3)/(max(max(mean(C,3))));
% figure(3); clf; 
% imshow(imgC');
% 
% % for frames 113-138
% d = 1;
% for i = 87:112
%     D(:,:,d) = J(:,:,i);
%     d = d+1;
% end
% imgD = mean(D,3)/(max(max(mean(D,3))));
% figure(4); clf; 
% imshow(imgD');

% for frames 139:165
e = 1;
for i = 139:165
    E(:,:,e) = J(:,:,i);
    e = e+1;
end
imgE = mean(E,3)/(max(max(mean(E,3))));

figure(1); clf; 
hold on
imshow(imgE');
title('Thick Slice')
hold off


J = wiener2(imgE',[7 7]);
B = ordfilt2(J,10,true(8));
figure(2);clf;
imshow(B)
title('2D Wiener Noise-Removal Filter')



J = imadjust(J);
B = ordfilt2(J,14,true(8));
figure(3);clf;
imshow(B)
title('2D Wiener Noise-Removal Filter With Contrast Enhanced')


J = squeeze(x.data(:,:,:,1));
figure(4); clf;
imshow(J(:,:,165)')
title('Unfiltered X')
%[x,y]  = ginput(3)

% ====================== %
% figure(5);clf;
% imshow(xdata1(:,:,100))
% title('Raw Image')
% 
% figure(8);clf;
% J = wiener2(xdata1(:,:,100),[7 7]);
% %J = imadjust(J);
% B = ordfilt2(J,10,true(8));
% imshow(B)
% title('2D Wiener Noise-Removal Filter')
% 
% figure(7);clf;
% J = wiener2(xdata1(:,:,100),[7 7]);
% J = imadjust(J);
% B = ordfilt2(J,14,true(8));
% imshow(B)
% title('2D Wiener Noise-Removal Filter With Contrast Enhanced')



%% Binary Gradient Mask
figure(5); clf;
filt_x = imgaussfilt(B, 1);
[~,threshold] = edge(filt_x,'sobel');
fudgeFactor = 0.2;
BG = edge(filt_x,'sobel',threshold * fudgeFactor);
imshow(BG)
title('Binary Gradient Mask')

% %% Bilated Gradient Mask
se90 = strel('line',3,90);
se0 = strel('line',3,0);
% 
% BWsdil = imdilate(BWs,[se90 se0]);
% figure(7);clf;
% imshow(BWsdil)
% title('Sobel Dilated Gradient Mask')
% 
% %% Prewitt
% figure(8); clf;
% [~,threshold] = edge(filt_x,'Prewitt');
% fudgeFactor = 1;
% BWs = edge(filt_x,'Prewitt',threshold * fudgeFactor);
% BWsdil = imdilate(BWs,[se90 se0]);
% imshow(BWsdil)
% title('Prewitt Dilated Gradient Mask')
% 
%% Roberts
figure(9); clf;
[~,threshold] = edge(B,'Roberts');
fudgeFactor = .3;
R = edge(B,'Roberts',threshold * fudgeFactor);
BWsdil = imdilate(R,[se90 se0]);
imshow(BWsdil)
title('Roberts Dilated Gradient Mask')

%% Canny
figure(6); clf;
[~,threshold] = edge(B,'Canny');
fudgeFactor = .8;
BWs = edge(B,'Canny',threshold * fudgeFactor);
BWsdil = imdilate(BWs,[se90 se0]);
imshow(BWsdil)
title('Canny Dilated Gradient Mask')

figure(10);clf;
[~,threshold] = edge(R,'Canny');
fudgeFactor = 1.3;
BWs = edge(B,'Canny',threshold * fudgeFactor);
BWsdil = imdilate(BWs,[se90 se0]);
imshow(BWsdil)
title('Roberts -> Canny Dilated Gradient Mask')

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