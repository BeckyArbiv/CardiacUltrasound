% Code last updated on 03/28/2021 by Eric Wei

%% Initialize and Define Variables
clear;
x = resampleDicom('06.dcm')

figure(1); clf
xdata1 = imrotate3(x.data(:, :, :, 5), 0, [1 0 0]);
sliceViewer(xdata1, 'SliceDirection', 'X');
%[x,y]  = ginput(3)

figure(2); clf
sliceViewer(xdata1, 'SliceDirection', 'Y');

figure(3); clf
sliceViewer(xdata1, 'SliceDirection', 'Z');
[x,y]  = ginput(2)

A = [28.1731; 150; 172.8532];
PC = [134.9954; 150; 156.1926];
AC = [150.0898; 150; 196.8315];
PA = [177; 124.9932; 130.0587];
AA = [177; 160.1820; 128.2815];

APC = PC-A;
AAC = AC-A;
APA = PA-A;
AAA = AA-A;

LongAxis = (APC + AAC + APA + AAA)/4
MagLongAxis = norm(LongAxis)