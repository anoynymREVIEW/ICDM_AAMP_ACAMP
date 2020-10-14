
clear
clc
% close all


subSeqLen = 100;  % The length of the sub-sequence
kNN_Uwant = 10;   % How many kNN do you want, just tell me

trainData = '/home/mondal/Music/Data_Backup/Miscellaneous/Workspace_MatLab/UCR_DTW_Tests/UCR_TS_Archive_2015/ArrowHead/ArrowHead_TRAIN';
testData = '/home/mondal/Music/Data_Backup/Miscellaneous/Workspace_MatLab/UCR_DTW_Tests/UCR_TS_Archive_2015/ArrowHead/ArrowHead_TEST';

TRAIN = load(trainData); % Only these two lines need to be changed to test a different dataset.?%
TEST  = load(testData); % Only these two lines need to be changed to test a different dataset.?%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TRAIN_class_labels = TRAIN(:,1);    % Pull out the class labels.
TRAIN(:,1) = [];                    % Remove class labels from training set.
TEST_class_labels = TEST(:,1);      % Pull out the class labels.
TEST(:,1) = [];                     % Remove class labels from testing set.
correct = 0; %


firstClass = 0;
secondClass = 1;
keepFirstClassData = zeros(1,size(TEST,2));
keepFirstClassDataIndx = zeros(1,1);

keepSecondClassData = zeros(1,size(TEST,2));
keepSecondClassDataIndx = zeros(1,1);

classCnt1 = 1;
classCnt2 = 1;
for ii = 1:1:size(TEST,1)
    if(( TEST_class_labels(ii)  == firstClass ) && (classCnt1 < 11) )
        keepFirstClassData(classCnt1, :) = TEST(ii, :);
        keepFirstClassDataIndx(classCnt1,1) = ii;
        
        classCnt1 = classCnt1 +1;
    end
    
    if( ( TEST_class_labels(ii)  == secondClass ) && (classCnt2 < 11) )
        keepSecondClassData(classCnt2, :) = TEST(ii, :);
        keepSecondClassDataIndx(classCnt2,1) = ii;
        
        classCnt2 = classCnt2 +1;
    end
end

[concatFirstClassData, concatFirstFileInfo] = concatenateData(keepFirstClassData);
[concatSecondClassData, concatSecondFileInfo] = concatenateData(keepSecondClassData);

[STOMP_pro_mul, STOMP_pro_idx] = STOMP_Youcef(concatFirstClassData, subSeqLen); % PAA ; A = concatFirstClassData

[STOMP_pro_mul_indep, STOMP_pro_idx_indep] = ...
                            STOMP_Independent_Join_Simple(concatSecondClassData, concatFirstClassData, subSeqLen); % PAB; A = concatFirstClassData, B = concatSecondClassData

[AAMP_pro_mul, AAMP_pro_idx] = AAMP_Tan(concatFirstClassData', subSeqLen); % PAA ; A = concatFirstClassData
[AAMP_pro_mul_indep, AAMP_pro_idx_indep] = ...
                            AAMP_1NN_IndependentJoin(concatSecondClassData, concatFirstClassData, subSeqLen); % PAB; A = concatFirstClassData, B = concatSecondClassData


diffDist_STOMP = ( STOMP_pro_mul_indep(1,:) - STOMP_pro_mul(:,1)');
shapelet_STOMP = [559, 828];

diffDist_AAMP = ( AAMP_pro_mul_indep(1,:) - AAMP_pro_mul(1,:) );
shapelet_AAMP = [534, 2038];

colorString = {'r','b','k','m','c','g','y',[.5 .6 .7],[.8 .2 .6], [0.72 0.52 0.04], [0.6 0.8 0.2], [0.619 1 0.4 ],...
    [0.9 0.9 0.9], [0.4 0.4 0.4], [0.8 0.8 0.8], [0.6 0.6 0.6], [0.2 0.2 0.2], [0 0.75 0.75], [0 0.5 0],...
    [0.4 0.58 0.9], [0.75 0 0.75], [0.8 0.7 0.6], [0.6 0.5 0.4 ], [0.8 0.6 1 ], [0 1 1], [1 0.6 0.8]};

% Original data 
figure();
hold on;
plot(1:length(concatFirstClassData(:)), concatFirstClassData(:),'color', colorString{3}, 'LineWidth',1);
plot(1:length(concatSecondClassData(:)), concatSecondClassData(:),'color', colorString{4}, 'LineWidth',1);
hleg1 = legend('Original time series A', ...
                'Original time series B');        
set(hleg1,'Location','NorthWest')
set(hleg1,'FontSize',16)

hold off;

% Plot for MP of STOMP and MP of STOMP_independent
figure();
hold on;
plot(1:length(STOMP_pro_mul(:)), STOMP_pro_mul(:),'color', colorString{1}, 'LineWidth',1);
plot(1:length(STOMP_pro_mul_indep(:)), STOMP_pro_mul_indep(:),'color', colorString{2}, 'LineWidth',1);
hleg1 = legend('Result of STOMP by using time series join J-AA', ...
                'Result of STOMP by using time series join J-AB');        
set(hleg1,'Location','NorthWest')
set(hleg1,'FontSize',16)
hold off;

% Plot for difference between MP of STOMP and MP of STOMP_independent
figure();
hold on;
plot(1:length(diffDist_STOMP(:)), diffDist_STOMP(:),'color', colorString{5}, 'LineWidth',2);
hleg1 = legend('Result difference between matrix profile J-AA and J-AB by STOMP');        
set(hleg1,'Location','NorthWest')
set(hleg1,'FontSize',16)
hold off;

% Plot for MP of AAMP and MP of AAMP_independent
figure();
hold on;
plot(1:length(AAMP_pro_mul(:)), AAMP_pro_mul(:),'color', colorString{1}, 'LineWidth',1);
plot(1:length(AAMP_pro_mul_indep(:)), AAMP_pro_mul_indep(:),'color', colorString{2}, 'LineWidth',1);
hleg1 = legend('Result of AAMP by using time series join J-AA', ...
                'Result of AAMP by using time series join J-AB');        
set(hleg1,'Location','NorthWest')
set(hleg1,'FontSize',16)
hold off;

% Plot for difference between MP of AAMP and MP of AAMP_independent
figure();
hold on;
plot(1:length(diffDist_AAMP(:)), diffDist_AAMP(:),'color', colorString{5}, 'LineWidth',2);
hleg1 = legend('Result difference between matrix profile J-AA and J-AB by AAMP');       
set(hleg1,'Location','NorthWest')
set(hleg1,'FontSize',16)
hold off;


figure();
hold on;
plot(1:length(concatFirstClassData(:)), concatFirstClassData(:),'color', colorString{3}, 'LineWidth',1);
% hleg1 = legend('Original time series A');        
% set(hleg1,'Location','NorthWest')
% set(hleg1,'FontSize',16)

for iOut = 1:1:length(shapelet_STOMP)
    outlierIndx = (shapelet_STOMP(iOut)); % here we get the paticular sub-sequence and then we repaint it accordingly
    plot(outlierIndx:(outlierIndx+subSeqLen-1), concatFirstClassData(outlierIndx:(outlierIndx+subSeqLen-1)),'Color','r', 'LineWidth',2);  
end
hold off;

figure();
hold on;
plot(1:length(concatFirstClassData(:)), concatFirstClassData(:),'color', colorString{3}, 'LineWidth',1);
for iOut = 1:1:length(shapelet_AAMP)
    outlierIndx = (shapelet_AAMP(iOut)); % here we get the paticular sub-sequence and then we repaint it accordingly
    plot(outlierIndx:(outlierIndx+subSeqLen-1), concatFirstClassData(outlierIndx:(outlierIndx+subSeqLen-1)),'Color','b', 'LineWidth',2);  
end
hold off;


function [keepAllTargetTogether, keepDataFileInfo] = concatenateData(subFoldersTarget)

% Get all the target sequences together and merged
keepAllTargetTogether = zeros(1,1);
fullPtCnt = 1;
keepDataFileInfo = cell(1,1);

getGoodFileCnt = 1;
for lTarget = 1:1:size(subFoldersTarget,1)
    C1Target = subFoldersTarget(lTarget,:);
    getLengthTarget = length(C1Target);
    
    C1TargetArr = zeros(getLengthTarget,1);
    
    C1TargetArr(:,1) = C1Target(1, :);
    
    clearvars C1Target
    
    if (getGoodFileCnt == 1)
        keepAllTargetTogether(1:(getLengthTarget),1) = C1TargetArr(1:end);
        
        keepDataFileInfo{getGoodFileCnt,1}.FileNum = lTarget;
        keepDataFileInfo{getGoodFileCnt,1}.DataStart = 1;
        keepDataFileInfo{getGoodFileCnt,1}.DataEnd = (getLengthTarget);
        fullPtCnt = fullPtCnt + (getLengthTarget-1);
    else
        keepAllTargetBackup = keepAllTargetTogether;
        keepAllTargetTogether = zeros ((size(keepAllTargetBackup,1)+size(C1TargetArr,1)),1);
        keepAllTargetTogether(1:size(keepAllTargetBackup,1),:) = keepAllTargetBackup(:,:);
        
        clearvars keepAllTargetBackup
        keepAllTargetTogether(fullPtCnt:(fullPtCnt+(getLengthTarget-1)),1) = C1TargetArr(:,1);
        
        keepDataFileInfo{getGoodFileCnt,1}.FileNum = lTarget;
        keepDataFileInfo{getGoodFileCnt,1}.DataStart = fullPtCnt;
        keepDataFileInfo{getGoodFileCnt,1}.DataEnd = (fullPtCnt+(getLengthTarget-1));
        
        if(( keepDataFileInfo{getGoodFileCnt,1}.DataEnd - ...
                keepDataFileInfo{getGoodFileCnt,1}.DataStart) > getLengthTarget)
            error('need to check it, there could be some problem here');
        end
        fullPtCnt = fullPtCnt + (getLengthTarget-1);
    end
    getGoodFileCnt = getGoodFileCnt + 1;
end

end