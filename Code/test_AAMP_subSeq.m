function test_AAMP_subSeq()
clear
clc
close all;


% DetectMotifsDiscordsProtien()
% DetectMotifsDiscordsSheep()

% PlotDataSubSeq();
 PlotDataSubSeq_All_Mean();

% zNormalizationProblemProteinData();
% zNormalizationProblemSheepData();

% matchDiffProblemProteinData()

% subSequenceExpSheepData();
% subSequenceExpProteinData();

% subSequenceExpSeismicData();
% subSequenceExpRandomWalkData();

end

function PlotDataSubSeq_All_Mean()

textFilePathSubSeq_Protien_1 = 'AAMP_ALL_Time/AAMP_All_Protein_SubSeq_1.txt';
textFilePathSubSeq_Sheep_1 = 'AAMP_ALL_Time/AAMP_All_Sheep_SubSeq_1.txt';

textFilePathSubSeq_Protien_2 = 'AAMP_ALL_Time/AAMP_All_Protein_SubSeq_2.txt';
textFilePathSubSeq_Sheep_2 = 'AAMP_ALL_Time/AAMP_All_Sheep_SubSeq_2.txt';

textFilePathSubSeq_Protien_3 = 'AAMP_ALL_Time/AAMP_All_Protein_SubSeq_3.txt';
textFilePathSubSeq_Sheep_3 = 'AAMP_ALL_Time/AAMP_All_Sheep_SubSeq_3.txt';

keepDataProtienSubSeq_1 = ReadFile_1(textFilePathSubSeq_Protien_1);
keepDataSheepSubSeq_1 = ReadFile_1(textFilePathSubSeq_Sheep_1);

keepDataProtienSubSeq_2 = ReadFile_1(textFilePathSubSeq_Protien_2);
keepDataSheepSubSeq_2 = ReadFile_1(textFilePathSubSeq_Sheep_2);

keepDataProtienSubSeq_3 = ReadFile_1(textFilePathSubSeq_Protien_3);
keepDataSheepSubSeq_3 = ReadFile_1(textFilePathSubSeq_Sheep_3);


keepDataProtienSubSeqMean = keepDataProtienSubSeq_1;
keepDataProtienSubSeqMean(:, 2:end) =  (keepDataProtienSubSeq_1(:, 2:end) + keepDataProtienSubSeq_2(:, 2:end) + keepDataProtienSubSeq_3(:, 2:end))./3;

keepDataSheepSubSeqMean = keepDataSheepSubSeq_1;
keepDataSheepSubSeqMean(:, 2:end) =  (keepDataSheepSubSeq_1(:, 2:end) + keepDataSheepSubSeq_2(:, 2:end) + keepDataSheepSubSeq_3(:, 2:end))./3;

ShowOnlyFigure(keepDataProtienSubSeqMean); % for protien data
ShowOnlyFigure(keepDataSheepSubSeqMean); % for sheep data


end


function PlotDataSubSeq()

textFilePathSubSeq_Protien = 'AAMP_SubSeq_Protien.txt';
textFilePathSubSeq_Sheep = 'AAMP_SubSeq_Sheep.txt';
textFilePathSubSeq_Seismic = 'AAMP_SubSeq_Seismic.txt';
textFilePathSubSeq_RandomWalk = 'AAMP_SubSeq_RandomWalk.txt';


keepDataProtienSubSeq = ReadFile_1(textFilePathSubSeq_Protien);
keepDataSheepSubSeq = ReadFile_1(textFilePathSubSeq_Sheep);
keepDataSeismicSubSeq = ReadFile_1(textFilePathSubSeq_Seismic);
keepDataRandomWalkSubSeq = ReadFile_2(textFilePathSubSeq_RandomWalk);

ShowOnlyFigure(keepDataProtienSubSeq); % for protien data
ShowOnlyFigure(keepDataSheepSubSeq); % for sheep data
ShowOnlyFigure(keepDataSeismicSubSeq); % for seismic data
ShowOnlyFigure(keepDataRandomWalkSubSeq); % for randomwalk data

end

function ShowOnlyFigure(keepData)
figure();
subplot(1,1,1);

hold on;
plot(keepData(:,1),keepData(:,2),'r-*');
hold on;
plot(keepData(:,1),keepData(:,3),'b-d');
hold on;
plot(keepData(:,1),keepData(:,4),'m-s');
hold on;
plot(keepData(:,1),keepData(:,5),'k-o');
hold on;
plot(keepData(:,1),keepData(:,6),'g-o');
hold on;
plot(keepData(:,1),keepData(:,7),'c-x');
hold on;


hleg1 = legend('AAMP','STOMP', 'SCRIMP', 'SCRIMP++', 'ACAMP-Optimized', 'ACAMP');
set(hleg1,'Location','NorthEast')
set(hleg1,'FontSize',10)
set(gca,'XTick', keepData(:,1));
grid on;
set(gca,'FontSize',10);
xl = xlabel('Sub-sequence Length');
yl = ylabel('Time Needed (sec.)');
set(xl,'FontSize',12,'FontWeight','bold','FontName','Courier');
set(yl,'FontSize',14,'FontWeight','bold','FontName','Courier');
hold off;
disp('see me');
end

function matchDiffProblemProteinData()
clear
subSeqLen = 256;  % The length of the sub-sequence
load ('ProteinData.mat');

[keepAllTargetTogether, ~] = ConcatenateAllSeries(realData(1:150,:)); % concatenate all the time series
[STOMP_pro_mul, STOMP_pro_idx] = STOMP_Youcef(keepAllTargetTogether, subSeqLen);
[AAMP_pro_mul, AAMP_pro_idx] = AAMP_Tan(keepAllTargetTogether', subSeqLen);

[sortedCloseVals, sortedCloseIndx] = sort(AAMP_pro_mul, 'descend');

onlyPickImpIndexes = zeros(30,1);
goodIndxCnt = 1;
for iTake = 1:1:length(sortedCloseVals)
    getIndx = sortedCloseIndx(iTake);
    
    if(iTake == 1)
        onlyPickImpIndexes(1,1) = getIndx;
    else
        
        smallDist = Inf;
        
        for ik = 1:1:(goodIndxCnt)
            getDist =  abs(onlyPickImpIndexes(ik,1)-getIndx);
            if (getDist < smallDist )
                smallDist = getDist;
                % closeIndex = onlyPickImpIndexes(ik,1);
            end
        end
        
        % after the for loop, I will get the closest index and it's
        % distance
        if(smallDist > round(subSeqLen/2))
            goodIndxCnt = goodIndxCnt +1;
            onlyPickImpIndexes(goodIndxCnt, 1) = getIndx;
        end
    end
    
    if(goodIndxCnt > 30)
        break;
    end
end

JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, STOMP_pro_idx, keepAllTargetTogether, 1, 10, subSeqLen);
JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, AAMP_pro_idx, keepAllTargetTogether, 1, 10, subSeqLen);

JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, STOMP_pro_idx, keepAllTargetTogether, 11, 20, subSeqLen);
JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, AAMP_pro_idx, keepAllTargetTogether, 11, 20, subSeqLen);

JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, STOMP_pro_idx, keepAllTargetTogether, 21, 30, subSeqLen);
JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, AAMP_pro_idx, keepAllTargetTogether, 21, 30, subSeqLen);

end

function zNormalizationProblemProteinData()
clear
subSeqLen = 256;  % The length of the sub-sequence
load ('ProteinData.mat');

[keepAllTargetTogether, ~] = ConcatenateAllSeries(realData(1:150,:)); % concatenate all the time series
[STOMP_pro_mul, STOMP_pro_idx] = STOMP_Youcef(keepAllTargetTogether, subSeqLen);
[SCRIMP_pro_mul, SCRIMP_pro_idx] = SCRIMP_Codde_Test(keepAllTargetTogether, subSeqLen, 2);
[ACAMP_mul, ACAMP_pro_idx] = ACAMP(keepAllTargetTogether, subSeqLen); % ACAMP algo
[AAMP_pro_mul, AAMP_pro_idx] = AAMP_Tan(keepAllTargetTogether', subSeqLen);

keepAllDiffMedians = zeros(length(STOMP_pro_mul),1);
for iDetect = 1:1:length(STOMP_pro_mul)
    matchIndx = STOMP_pro_idx(iDetect);
    realSubSeq = keepAllTargetTogether(iDetect:iDetect+subSeqLen-1, 1);
    matchSubSeq = keepAllTargetTogether(matchIndx:matchIndx+subSeqLen-1, 1);
    
    realSubSeqMedian = median(realSubSeq);
    matchSubSeqMedian = median(matchSubSeq);
    diffMedian = abs(matchSubSeqMedian-realSubSeqMedian);
    
    keepAllDiffMedians(iDetect,1) = diffMedian;
end

[sortedMedianVals, sortedMedianIndx] = sort(keepAllDiffMedians, 'descend');

onlyPickImpIndexes = zeros(30,1);
goodIndxCnt = 1;
for iTake = 1:1:5000
    getIndx = sortedMedianIndx(iTake);
    
    if(iTake == 1)
        onlyPickImpIndexes(1,1) = getIndx;
    else
        
        % calculate the entry with all other previous elements in the array
        % and compute the distances with them so that we can find the
        % closest index. Now see whether the closest index is very close
        % with the new one or not
        smallDist = Inf;
        
        for ik = 1:1:(goodIndxCnt)
            getDist =  abs(onlyPickImpIndexes(ik,1)-getIndx);
            if (getDist < smallDist )
                smallDist = getDist;
                % closeIndex = onlyPickImpIndexes(ik,1);
            end
        end
        
        % after the for loop, I will get the closest index and it's
        % distance
        if(smallDist > round(subSeqLen/2))
            goodIndxCnt = goodIndxCnt +1;
            onlyPickImpIndexes(goodIndxCnt, 1) = getIndx;
        end
    end
    
    if(goodIndxCnt > 30)
        break;
    end
end
indxArr = [17, 19, 23, 27];
JustPlotFor_Z_NOrmalizationSelected(onlyPickImpIndexes, STOMP_pro_idx, AAMP_pro_idx, keepAllTargetTogether, indxArr, subSeqLen);

% JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, STOMP_pro_idx, keepAllTargetTogether, 1, 10, subSeqLen);
% JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, AAMP_pro_idx, keepAllTargetTogether, 1, 10, subSeqLen);
% JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, SCRIMP_pro_idx, keepAllTargetTogether, 1, 10, subSeqLen);
% JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, ACAMP_pro_idx, keepAllTargetTogether, 1, 10, subSeqLen);
% 
% JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, STOMP_pro_idx, keepAllTargetTogether, 11, 20, subSeqLen);
% JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, AAMP_pro_idx, keepAllTargetTogether, 11, 20, subSeqLen);
% JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, SCRIMP_pro_idx, keepAllTargetTogether, 11, 20, subSeqLen);
% JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, ACAMP_pro_idx, keepAllTargetTogether, 11, 20, subSeqLen);
% 
% JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, STOMP_pro_idx, keepAllTargetTogether, 21, 30, subSeqLen);
% JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, AAMP_pro_idx, keepAllTargetTogether, 21, 30, subSeqLen);
% JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, SCRIMP_pro_idx, keepAllTargetTogether, 21, 30, subSeqLen);
% JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, ACAMP_pro_idx, keepAllTargetTogether, 21, 30, subSeqLen);

end

function JustPlotFor_Z_NOrmalizationSelected(onlyPickImpIndexes, pro_idx, AAMP_pro_idx, keepAllTargetTogether, indxArr, subSeqLen)
% for AAMP plot
% take top 10
hFig = figure;
subplot(3,4,1);

plotCntOrig = 1;
plotCntMatch_stomp = 5;
plotCntMatch_aamp = 9;

for iTake = 1:1:length(indxArr)
    getTake = indxArr(iTake);
    
    getIndx = onlyPickImpIndexes(getTake);
    stomp_matchIndx = pro_idx(getIndx);
    aamp_matchIndx = AAMP_pro_idx(getIndx);
    
    realSubSeq = keepAllTargetTogether(getIndx:getIndx+subSeqLen-1, 1);
    stomp_matchSubSeq = keepAllTargetTogether(stomp_matchIndx:stomp_matchIndx+subSeqLen-1, 1);
    aamp_matchSubSeq = keepAllTargetTogether(aamp_matchIndx:aamp_matchIndx+subSeqLen-1, 1);
        
    hold on;
    subplot(3,4,plotCntOrig);
    plot(realSubSeq, 'color', 'r', 'LineWidth',2);
    hleg1 = legend('Query Signal');        
    set(hleg1,'Location','best')
    set(hleg1,'FontSize',12)
    
    subplot(3,4,plotCntMatch_stomp);
    plot(stomp_matchSubSeq, 'color', 'k', 'LineWidth',2);
    hleg1 = legend('STOMP Match');        
    set(hleg1,'Location','best')
    set(hleg1,'FontSize',12)
    
    subplot(3,4,plotCntMatch_aamp);
    plot(aamp_matchSubSeq, 'color', 'b', 'LineWidth',2);
    hleg1 = legend('AAMP Match');        
    set(hleg1,'Location','best')
    set(hleg1,'FontSize',12)
    
    plotCntOrig = plotCntOrig +1;
    plotCntMatch_stomp = plotCntMatch_stomp +1;
    plotCntMatch_aamp = plotCntMatch_aamp +1;
end
hold off;
end

function JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, pro_idx, keepAllTargetTogether, st, endLp, subSeqLen)
% for AAMP plot
% take top 10
hFig = figure;
subplot(2,10,1);

plotCntOrig = 1;
plotCntMatch = 11;

for iTake = st:1:endLp
    getIndx = onlyPickImpIndexes(iTake);
    matchIndx = pro_idx(getIndx);
    realSubSeq = keepAllTargetTogether(getIndx:getIndx+subSeqLen-1, 1);
    matchSubSeq = keepAllTargetTogether(matchIndx:matchIndx+subSeqLen-1, 1);
    
    hold on;
    subplot(2,10,plotCntOrig);
    plot(realSubSeq);
    
    subplot(2,10,plotCntMatch);
    plot(matchSubSeq);
    
    plotCntOrig = plotCntOrig +1;
    plotCntMatch = plotCntMatch +1;
end
hold off;
end


function zNormalizationProblemSheepData()
clear

% Loading the sheep data
load('SheepDataFull.mat')

subSeqLen = 256;  % The length of the sub-sequence


[keepAllTargetTogether, ~] = ConcatenateAllSeries(keepAllData(1:20,:)); % concatenate all the time series
[STOMP_pro_mul, STOMP_pro_idx] = STOMP_Youcef(keepAllTargetTogether, subSeqLen);
[AAMP_pro_mul, AAMP_pro_idx] = AAMP_Tan(keepAllTargetTogether', subSeqLen);

keepAllDiffMedians = zeros(length(STOMP_pro_mul),1);
for iDetect = 1:1:length(STOMP_pro_mul)
    matchIndx = STOMP_pro_idx(iDetect);
    realSubSeq = keepAllTargetTogether(iDetect:iDetect+subSeqLen-1, 1);
    matchSubSeq = keepAllTargetTogether(matchIndx:matchIndx+subSeqLen-1, 1);
    
    realSubSeqMedian = median(realSubSeq);
    matchSubSeqMedian = median(matchSubSeq);
    diffMedian = abs(matchSubSeqMedian-realSubSeqMedian);
    
    keepAllDiffMedians(iDetect,1) = diffMedian;
end

[sortedMedianVals, sortedMedianIndx] = sort(keepAllDiffMedians, 'descend');

onlyPickImpIndexes = zeros(30,1);
goodIndxCnt = 1;
for iTake = 1:1:5000
    getIndx = sortedMedianIndx(iTake);
    
    if(iTake == 1)
        onlyPickImpIndexes(1,1) = getIndx;
    else
        
        % calculate the entry with all other previous elements in the array
        % and compute the distances with them so that we can find the
        % closest index. Now see whether the closest index is very close
        % with the new one or not
        smallDist = Inf;
        
        for ik = 1:1:(goodIndxCnt)
            getDist =  abs(onlyPickImpIndexes(ik,1)-getIndx);
            if (getDist < smallDist )
                smallDist = getDist;
                % closeIndex = onlyPickImpIndexes(ik,1);
            end
        end
        
        % after the for loop, I will get the closest index and it's
        % distance
        if(smallDist > round(subSeqLen/2))
            onlyPickImpIndexes(goodIndxCnt, 1) = getIndx;
            goodIndxCnt = goodIndxCnt +1;
        end
    end
    
    if(goodIndxCnt > 30)
        break;
    end
end

indxArr = [1, 8, 12, 15];
JustPlotFor_Z_NOrmalizationSelected(onlyPickImpIndexes, STOMP_pro_idx, AAMP_pro_idx, keepAllTargetTogether, indxArr, subSeqLen);

% JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, STOMP_pro_idx, keepAllTargetTogether, 1, 10, subSeqLen);
% JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, AAMP_pro_idx, keepAllTargetTogether, 1, 10, subSeqLen);
% 
% JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, STOMP_pro_idx, keepAllTargetTogether, 11, 20, subSeqLen);
% JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, AAMP_pro_idx, keepAllTargetTogether, 11, 20, subSeqLen);
% 
% JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, STOMP_pro_idx, keepAllTargetTogether, 21, 30, subSeqLen);
% JustPlotFor_Z_NOrmalization(onlyPickImpIndexes, AAMP_pro_idx, keepAllTargetTogether, 21, 30, subSeqLen);
end


function DetectMotifsDiscordsProtien()
load ('ProteinData.mat');
[keepAllTargetTogether, ~] = ConcatenateAllSeries(realData(1:100,:)); % concatenate all the time series
apply_AAMP(keepAllTargetTogether, 500);   % working for the 2nd column here

end

function DetectMotifsDiscordsSheep()
load('SheepDataFull.mat')
[keepAllTargetTogether, ~] = ConcatenateAllSeries(realData(1:100,:)); % concatenate all the time series
apply_AAMP(keepAllTargetTogether, 500);   % working for the 2nd column here

end


function subSequenceExpProteinData()
clear

subSeqLen = 2000;  % The length of the sub-sequence
load ('ProteinData.mat');

[keepAllTargetTogether, ~] = ConcatenateAllSeries(realData(1:100,:)); % concatenate all the time series

while (subSeqLen < ( round((length(keepAllTargetTogether)*52)/100) ))
    fprintf('Protien Data : The sub-seqeunce length is : %d \n', subSeqLen);
    
    apply_AAMP(keepAllTargetTogether, subSeqLen);   % working for the 2nd column here
    subSeqLen = subSeqLen + 3000;
end
end


function subSequenceExpSheepData()
clear

% Loading the sheep data
load('SheepDataFull.mat')


[keepAllTargetTogether, ~] = ConcatenateAllSeries(keepAllData(1:100,:)); % concatenate all the time series
subSeqLen = 2000;  % The length of the sub-sequence

while (subSeqLen < ( round((length(keepAllTargetTogether)*52)/100) ))
    fprintf('Sheep Data : The sub-seqeunce length is : %d \n', subSeqLen);
    apply_AAMP(keepAllTargetTogether, subSeqLen);   % working for the 2nd column here
    
    subSeqLen = subSeqLen + 3000;
end
end


function subSequenceExpSeismicData()
clear
% Loading the sheep data
load('seismic_50000.mat');
subSeqLen = 30;  % The length of the sub-sequence

keepAllData = seeMe1;
clear seeMe1;

queryIndex = 307;
queryComplete = keepAllData(queryIndex, 1:end-1);

noOfZerosInQuery = sum(queryComplete(:)==0);
if(noOfZerosInQuery > (length(queryComplete) /2))
    error ('The number of zeros in the subsequence is more than the half of \n query size. Please choose another query');
end

while (subSeqLen < ( round((length(queryComplete)*90)/100) ))
    fprintf('Seismic Data : The sub-seqeunce length is : %d \n', subSeqLen);
    
    apply_AAMP(keepAllData(1:100,:), subSeqLen);   % working for the 2nd column here
    subSeqLen = subSeqLen + 20;
end

end


function subSequenceExpRandomWalkData()
clear
% Loading the sheep data
load('randomWalk_50000.mat');
subSeqLen = 30;  % The length of the sub-sequence

keepAllData = seeMe1;
clear seeMe1;

queryIndex = 307;
queryComplete = keepAllData(queryIndex, 1:end-1);

noOfZerosInQuery = sum(queryComplete(:)==0);
if(noOfZerosInQuery > (length(queryComplete) /2))
    error ('The number of zeros in the subsequence is more than the half of \n query size. Please choose another query');
end

while (subSeqLen < ( round((length(queryComplete)*90)/100) ))
    fprintf('Random Walk Data : The sub-seqeunce length is : %d \n', subSeqLen);
    
    apply_AAMP(keepAllData(1:100,:), subSeqLen);   % working for the 2nd column here
    subSeqLen = subSeqLen + 20;
end

end


function [keepAllTargetTogether, keepDataFileInfo] = ConcatenateAllSeries(subFoldersTarget )
% Get all the target sequences together and merged
keepAllTargetTogether = zeros(1,1);
fullPtCnt = 1;
keepDataFileInfo = cell(1,1);

getGoodFileCnt = 1;
for lTarget = 1:1:size(subFoldersTarget,1) % get the target files
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
        fullPtCnt = fullPtCnt + (getLengthTarget);
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
        fullPtCnt = fullPtCnt + (getLengthTarget);
    end
    getGoodFileCnt = getGoodFileCnt + 1;
end

return
end


function keepData = ReadFile_1(textFilePath)
keepData = zeros(1,2);
fid = fopen(textFilePath);
tline = fgetl(fid);
lnCnt = 1;
while ischar(tline)
    if(~isempty(tline))
        newStr = split(tline,' ');
        %     subSeqlen = (newStr{6});
        if(strcmp(newStr{6}, 'length'))
            subSeqlen = str2double(newStr{9});
        end
        keepData(lnCnt, 1) = subSeqlen;
        tlineAAMP = fgetl(fid);
        newStrAAMP = split(tlineAAMP,' ');
        
        tlineSTOMP = fgetl(fid);
        newStrSTOMP = split(tlineSTOMP,' ');
        
        tlineSCRIMP = fgetl(fid);
        newStrSCRIMP = split(tlineSCRIMP,' ');
        
        tlineSCRIMPPlusPlus = fgetl(fid);
        newStrSCRIMP_PlusPLus = split(tlineSCRIMPPlusPlus,' ');
        
        tlineACAMP = fgetl(fid);
        newStrACAMP = split(tlineACAMP,' ');
        
        tlineACAMP_Colez = fgetl(fid);
        newStrACAMP_Colez = split(tlineACAMP_Colez,' ');
        
        keepData(lnCnt, 2) = str2double(newStrAAMP{7});
        keepData(lnCnt, 3) = str2double(newStrSTOMP{7});
        keepData(lnCnt, 4) = str2double(newStrSCRIMP{7});
        keepData(lnCnt, 5) = str2double(newStrSCRIMP_PlusPLus{9});
        
        keepData(lnCnt, 6) = str2double(newStrACAMP{7});
        keepData(lnCnt, 7) = str2double(newStrACAMP_Colez{7});
        
        lnCnt = lnCnt +1;
    end
    tline = fgetl(fid);
end
fclose(fid);
return;
end


function keepData = ReadFile_2(textFilePath)
keepData = zeros(1,2);
fid = fopen(textFilePath);
tline = fgetl(fid);
lnCnt = 1;
while ischar(tline)
    if(~isempty(tline))
        newStr = split(tline,' ');
        %     subSeqlen = (newStr{6});
        if(strcmp(newStr{7}, 'length'))
            subSeqlen = str2double(newStr{10});
        end
        keepData(lnCnt, 1) = subSeqlen;
        tlineAAMP = fgetl(fid);
        newStrAAMP = split(tlineAAMP,' ');
        
        tlineSTOMP = fgetl(fid);
        newStrSTOMP = split(tlineSTOMP,' ');
        
        tlineSCRIMP = fgetl(fid);
        newStrSCRIMP = split(tlineSCRIMP,' ');
        
        tlineSCRIMPPlusPlus = fgetl(fid);
        newStrSCRIMP_PlusPLus = split(tlineSCRIMPPlusPlus,' ');
        
        tlineACAMP = fgetl(fid);
        newStrACAMP = split(tlineACAMP,' ');
        
        tlineACAMP_Colez = fgetl(fid);
        newStrACAMP_Colez = split(tlineACAMP_Colez,' ');
        
        keepData(lnCnt, 2) = str2double(newStrAAMP{7});
        keepData(lnCnt, 3) = str2double(newStrSTOMP{7});
        keepData(lnCnt, 4) = str2double(newStrSCRIMP{7});
        keepData(lnCnt, 5) = str2double(newStrSCRIMP_PlusPLus{9});
        
        keepData(lnCnt, 6) = str2double(newStrACAMP{7});
        keepData(lnCnt, 7) = str2double(newStrACAMP_Colez{7});
        
        lnCnt = lnCnt +1;
    end
    tline = fgetl(fid);
end
fclose(fid);
return;
end
