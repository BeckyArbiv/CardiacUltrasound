% Code last updated on 04/19/2021 by Eric Wei

%% Initialize and Define Variables
clear;
close all
%x = resampleDicom('06.dcm')
file = ('/Users/beckyarbiv/Documents/BME/BME 543/06.dcm');
X = resampleDicom(file);
MagLongAxis = zeros(1, X.NumVolumes);
% % % for t = 1:X.NumVolumes
% % % % figure(1); clf
 %xdata1 = imrotate3(X.data(:, :, :, 1), 0, [1 0 0]);
% % % % sliceViewer(xdata1, 'SliceDirection', 'X');
% % % % 
% % % % figure(2); clf
% % % % sliceViewer(xdata1, 'SliceDirection', 'Y');
% % % % 
% % % % figure(3); clf
% % % % sliceViewer(xdata1, 'SliceDirection', 'Z');
% % % 

% figure(4); clf

% se90 = strel('line',3,90);
% se0 = strel('line',3,0);
% for i = 1:size(xavgslice,1)
%     J = wiener2(squeeze(xavgslice(i,:,:)),[7 7]);
%     B = ordfilt2(J,10,true(8));
%     % filt_x = imgaussfilt(B, 1);
%     % [~,threshold] = edge(filt_x,'sobel');
%     % fudgeFactor = 0.2;
%     % BG = edge(filt_x,'sobel',threshold * fudgeFactor);
%     % imagesc(BG)
%     [~,threshold] = edge(B,'Canny');
%     fudgeFactor = .8;
%     BWs = edge(B,'Canny',threshold * fudgeFactor);
%     BWsdil = imdilate(BWs,[se90 se0]);
%     colormap gray
%     imagesc(BWsdil);
%     xlabel('z')
%     ylabel('y')
%     text(5,5,sprintf('X-slice: %d',i),'color','white');
%     pause(0.25);
% end
% 33 for x for dataset 6 frame 5

% figure(5)
% 
% for i = 1:size(yavgslice,1)
%     imagesc(squeeze(yavgslice(i,:,:)));
%     %  colormap gray
%     text(5,5,sprintf('Y-slice: %d',i),'color','white');
%     pause(0.25);
% end
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
%% Edge tracking lets gooooo
xdata1 = imrotate3(X.data(:, :, :, 1), 0, [1 0 0]);
[xavgslice,yavgslice,zavgslice] = avg_slices(xdata1);
[row, col, distance,x_input,y_input,index_x, index_y] = edge_tracking(xavgslice);
for t = 2:X.NumVolumes
xdata1 = imrotate3(X.data(:, :, :, t), 0, [1 0 0]);
[xavgslice,yavgslice,zavgslice] = avg_slices(xdata1);
[index_x, index_y] = track_over_frames(xavgslice,x_input,y_input,index_x, index_y,t)
     
end

function [row, col, distance,x_input,y_input,index_x, index_y] = edge_tracking(xavgslice)
figure(7); clf
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
BWsdil = imdilate(BG,[se90 se0]);

imagesc(BWsdil)
xlabel('z')
ylabel('y')
click = ginput(1)
x_input = round(click(1));
y_input = round(click(2));
plot(x_input, y_input,'or')
% define ROI
roiSize = 5;
roiEdge = (roiSize-1)/2;
roiCent = roiEdge+1;
dim=size(BWsdil,1);
roi = int8(BWsdil([y_input-roiEdge:y_input+roiEdge],[x_input-roiEdge:x_input+roiEdge]));
xkernal = [click(1)-1 click(1)+1 click(1)+1 click(1)-1 click(1)+1 click(1) click(1) click(1)+1 click(1)-1];
ykernal = [click(2)-1 click(2)-1 click(2)+1 click(2)+1 click(2)-1 click(2)+1 click(2)-1 click(2) click(2)];
plot([x_input-roiEdge:x_input+roiEdge],[y_input-roiEdge:y_input+roiEdge],'*')
hold off

if all(roi == 1)
    error("Please select an edge")
elseif all(roi == 0)
    error("Please select an edge")
end

for p = 1:roiEdge
    if any(roi(roiCent, roiCent) ~= BWsdil([y_input-p:y_input+p],[x_input-p:x_input+p]),'all')
        [edge_row, edge_col] = find((roi(roiCent, roiCent) ~= BWsdil([y_input-p:y_input+p],[x_input-p:x_input+p])))
        if p == 1
            roiCent = 2;
        end
        distance = sqrt( (edge_row-roiCent).^2 + (edge_col-roiCent).^2 );
        closest_point = find(distance== min(distance))
        row = edge_row(closest_point);
        col = edge_col(closest_point);
        break
    %elseif any(roi(roiCent, roiCent) == BWsdil([y_input-p:y_input+p],[x_input-p:x_input+p]),'all')
    end
end
index_y = y_input + (-row+roiCent);
index_x = y_input + (col-roiCent);
end
%% Average Slices
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

function [rumba_x, rumba_y] = track_over_frames(xavgslice,x_input,y_input,index_x, index_y,g)
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
BWsdil = imdilate(BG,[se90 se0]);

imagesc(BWsdil)
xlabel('z')
ylabel('y')

plot(x_input,y_input,'g*')


roiSize = 3;
roiEdge = (roiSize-1)/2;
roiCent = roiEdge+1;

roi = flip(int8(BWsdil([y_input-roiEdge:y_input+roiEdge],[x_input-roiEdge:x_input+roiEdge])),1);
rumba = 0;
plot([x_input-roiEdge:x_input+roiEdge],[y_input-roiEdge:y_input+roiEdge],'*')
hold off
while rumba == 0
    if  any(roi(roiCent, roiCent) ~= BWsdil([y_input-roiEdge:y_input+roiEdge],[x_input-roiEdge:x_input+roiEdge]),'all')
        [edge_row, edge_col] = find((roi(roiCent, roiCent) ~= BWsdil([y_input-roiEdge:y_input+roiEdge],[x_input-roiEdge:x_input+roiEdge])))
        distance = sqrt( (edge_row-roiCent).^2 + (edge_col-roiCent).^2 );
        closest_point = find(distance== min(distance))
        row = edge_row(closest_point);
        col = edge_col(closest_point);
        rumba = 1;
    else
        roiSize = roiSize + 2;
        roiEdge = (roiSize-1)/2;
        roiCent = roiEdge+1;
        roi = int8(BWsdil([y_input-roiEdge:y_input+roiEdge],[x_input-roiEdge:x_input+roiEdge]));
    end
end

rumba_y = y_input + (-row+roiCent);
rumba_x = x_input + (col-roiCent);
end
%%
% % figure(7); clf
% % xslice = 33; % find optimal slice number
% % J = wiener2(squeeze(xavgslice(xslice,:,:)),[7 7]);
% % B = ordfilt2(J,10,true(8));
% % filt_x = imgaussfilt(B, 1);
% % [~,threshold] = edge(filt_x,'sobel');
% % fudgeFactor = 0.2;
% % BG = edge(filt_x,'sobel',threshold * fudgeFactor);
% % imagesc(BG)
% % % imagesc(B)
% % % colormap gray
% % xlabel('z')
% % ylabel('y')
% % [z1,y1]  = ginput(3)
% % 
% % figure(8); clf
% % yslice = 30; % find optimal slice number
% % J = wiener2(squeeze(yavgslice(yslice,:,:)),[7 7]);
% % B = ordfilt2(J,10,true(8));
% % imagesc(B)
% % colormap gray
% % xlabel('z')
% % ylabel('x')
% % [z2,x1]  = ginput(3)
% % 
% % figure(9); clf
% % zslice = 10; % find optimal slice number
% % J = wiener2(squeeze(zavgslice(zslice,:,:)),[7 7]);
% % B = ordfilt2(J,10,true(8));
% % imagesc(B)
% % colormap gray
% % xlabel('x')
% % ylabel('y')
% % [x2,y2]  = ginput(1)
% % 
% % A = [(x1(1)+x2(1))/2 (y1(1) + y2(1))/2 (z1(1) + z2(1))/2];
% % AC = [xslice*avgslice y1(2) z1(2)];
% % PC = [xslice*avgslice y1(3) z1(3)];
% % AA = [x1(2) yslice*5 z2(2)];
% % PA = [x1(3) yslice*5 z2(3)];
% % 
% % APC = PC-A;
% % AAC = AC-A;
% % APA = PA-A;
% % AAA = AA-A;
% % 
% % LongAxis = (APC + AAC + APA + AAA)/4
% % MagLongAxis(t) = norm(LongAxis);
% % 
% % % A = [28.1731; 150; 172.8532];
% % % PC = [134.9954; 150; 156.1926];
% % % AC = [150.0898; 150; 196.8315];
% % % PA = [177; 124.9932; 130.0587];
% % % AA = [177; 160.1820; 128.2815];
% % % 
% % % APC = PC-A;
% % % AAC = AC-A;
% % % APA = PA-A;
% % % AAA = AA-A;
% % % 
% % % LongAxis = (APC + AAC + APA + AAA)/4
% % % MagLongAxis = norm(LongAxis)
% % 
% % % figure(10);
% % % hold on
% % % line([A(1) LongAxis(1)], [A(2) LongAxis(2)], [A(3) LongAxis(3)],'Color',[.6 .9 1/t]);
% % % legend(['Frame: ' t])
% % % view(3)
% % 
% % figure(11);clf
% % hold on
% % line([A(3) (LongAxis(3)+A(3))], [A(2) (LongAxis(2)+A(2))]);axis([0 206 0 277]);
% % end
