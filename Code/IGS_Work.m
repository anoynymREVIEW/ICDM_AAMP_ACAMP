% To show the problem of sigma = 0 in original STOMP algorithm and how the
% proposed algorithm "AAMP" is capable to overcome this issue

% To run this code, you need to download the dataset from link : 
% https://drive.google.com/drive/folders/10WHOK5qEaUiZREW5Tf9dReRJoAua4dMy

% Author information omitted for ICDM review.
% For details of the code, see:
% "Efficient Matrix Profile Algorithms for Normalized and Non-Normalized Distances", submitted to ICDM 2020.

clear
close all
clc


targetFilePath = '/home/mondal/Documents/All_Safrans/igs-repro2/0194.dat';

getTable = readtable(targetFilePath);

noOfOutliers = 10;
winLen = 15;

Varib2 = getTable.Var2;
Varib3 = getTable.Var3;
Varib4 = getTable.Var4;


[~,fileNamOnly,~] = fileparts(targetFilePath);



% The first column
% [pro_mul_1, pro_idx_1] = STOMP(Varib2, winLen);
[pro_mul_1, pro_idx_1] = AAMP(Varib2', winLen);

pro_mul_1(isinf(pro_mul_1(:,1)), 1) = 0; % if there are some inf values then replace them by zeros so that they couldn't occur in outliers 
[sortVal1, sortIndx1] = sort(pro_mul_1, 'descend');



% The second column
%  [pro_mul_2, pro_idx_2] = STOMP(Varib3, winLen);
[pro_mul_2, pro_idx_2] = AAMP(Varib3', winLen);

pro_mul_2(isinf(pro_mul_2(:,1)), 1) = 0; % if there are some inf values then replace them by zeros so that they couldn't occur in outliers 
[sortVal2, sortIndx2] = sort(pro_mul_2, 'descend');



% The third column
% [pro_mul_3, pro_idx_3] = STOMP(Varib4, winLen);
[pro_mul_3, pro_idx_3] = AAMP(Varib4', winLen);

pro_mul_3(isinf(pro_mul_3(:,1)), 1) = 0; % if there are some inf values then replace them by zeros so that they couldn't occur in outliers 
[sortVal3, sortIndx3] = sort(pro_mul_3, 'descend');



plotTheGraph(Varib2, sortIndx1, noOfOutliers, pro_mul_1, winLen, 'Longitude', fileNamOnly);
% plotTheGraph(Varib3, sortIndx2, noOfOutliers, pro_mul_2, winLen, 'Latitude', fileNamOnly);
plotTheGraph(Varib4, sortIndx3, noOfOutliers, pro_mul_3, winLen, 'Height', fileNamOnly);




function plotTheGraph(Varib2, sortIndx1, noOfOutliers, pro_mul_1, winLen, str, fileNamOnly)

hFig1 = figure();
% h(1) = subplot(1,1,1);
plt1 = plot(1:length(Varib2), Varib2, 'k-', 'LineWidth',1);
hold on;
for iOut = 1:1:noOfOutliers
    outlierIndx = (sortIndx1(iOut)); % here we get the paticular sub-sequence and then we repaint it accordingly
    plot(outlierIndx:(outlierIndx+winLen-1), Varib2(outlierIndx:(outlierIndx+winLen-1)),'Color','r', 'LineWidth',2);  
end
hold off;
% hleg1 = legend(plt1,{str});
% set(hleg1,'Location','Best')
% set(hleg1,'FontSize',12)

fileNamDesig = strcat(fileNamOnly,str );
fileNamDesig = strcat(fileNamDesig,'Outliers');
imageNam = strcat(fileNamDesig,'.png');
saveas(hFig1,imageNam);
close(hFig1);

hFig1_1 = figure();
% h(1) = subplot(1,1,1);


plt1_1 = plot(pro_mul_1, 'r-', 'LineWidth',1);
matrixStr = strcat('Matrix Profile-',str);


% hleg1_1 = legend(plt1_1,{matrixStr});
% set(hleg1_1,'Location','Best')
% set(hleg1_1,'FontSize',12)

fileNamDesig = strcat(fileNamOnly,str );
fileNamDesig = strcat(fileNamDesig,'MP');
imageNam = strcat(fileNamDesig,'.png');
saveas(hFig1_1,imageNam);
close(hFig1_1);
end