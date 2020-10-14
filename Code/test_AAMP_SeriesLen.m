function test_AAMP_SeriesLen()
% clear
% clc
% close all;


% PlotDataSubSeq();
PlotDataSubSeq_All_Mean();
 
% subSequenceExpSheepData();
% subSequenceExpProteinData();


% subSequenceExpSeismicData();
% subSequenceExpRandomWalkData();

end

function PlotDataSubSeq_All_Mean()

textFilePathSeriesLen_Protien_1 = 'AAMP_ALL_Time/AAMP_All_Protein_SeriesLen_1.txt';
textFilePathSeriesLen_Sheep_1 = 'AAMP_ALL_Time/AAMP_All_Sheep_SeriesLen_1.txt';

textFilePathSeriesLen_Protien_2 = 'AAMP_ALL_Time/AAMP_All_Protein_SeriesLen_2.txt';
textFilePathSeriesLen_Sheep_2 = 'AAMP_ALL_Time/AAMP_All_Sheep_SeriesLen_2.txt';

textFilePathSeriesLen_Protien_3 = 'AAMP_ALL_Time/AAMP_All_Protein_SeriesLen_3.txt';
textFilePathSeriesLen_Sheep_3 = 'AAMP_ALL_Time/AAMP_All_Sheep_SeriesLen_3.txt';

keepDataProtienSeriesLen_1 = ReadFile_1(textFilePathSeriesLen_Protien_1);
keepDataSheepSeriesLen_1 = ReadFile_2(textFilePathSeriesLen_Sheep_1);

keepDataProtienSeriesLen_2 = ReadFile_1(textFilePathSeriesLen_Protien_2);
keepDataSheepSeriesLen_2 = ReadFile_2(textFilePathSeriesLen_Sheep_2);

keepDataProtienSeriesLen_3 = ReadFile_1(textFilePathSeriesLen_Protien_3);
keepDataSheepSeriesLen_3 = ReadFile_2(textFilePathSeriesLen_Sheep_3);

keepDataProtienSeriesLenMean = keepDataProtienSeriesLen_1;
keepDataProtienSeriesLenMean(:, 2:end) =  (keepDataProtienSeriesLen_1(:, 2:end) ...
                    + keepDataProtienSeriesLen_2(:, 2:end) + keepDataProtienSeriesLen_3(:, 2:end))./3;

keepDataSheepSeriesLenMean = keepDataSheepSeriesLen_1;
keepDataSheepSeriesLenMean(:, 2:end) =  (keepDataSheepSeriesLen_1(:, 2:end) ...
                + keepDataSheepSeriesLen_2(:, 2:end) + keepDataSheepSeriesLen_3(:, 2:end))./3;


ShowOnlyFigure(keepDataProtienSeriesLenMean); % for protien data
ShowOnlyFigure(keepDataSheepSeriesLenMean); % for sheep data

end

function PlotDataSubSeq()

textFilePathSeriesLen_Protien = 'AAMP_SeriesLen_Protien.txt';
textFilePathSeriesLen_Sheep = 'AAMP_SeriesLen_Sheep.txt';
% textFilePathSeriesLen_Seismic = 'AAMP_SeriesLen_Seismic.txt';
% textFilePathSeriesLen_RandomWalk = 'AAMP_SeriesLen_RandomWalk.txt';


keepDataProtienSeriesLen = ReadFile_1(textFilePathSeriesLen_Protien);
keepDataSheepSeriesLen = ReadFile_2(textFilePathSeriesLen_Sheep);
% keepDataSeismicSeriesLen = ReadFile_3(textFilePathSeriesLen_Seismic);
% keepDataRandomWalkSeriesLen = ReadFile_3(textFilePathSeriesLen_RandomWalk);

ShowOnlyFigure(keepDataProtienSeriesLen); % for protien data
ShowOnlyFigure(keepDataSheepSeriesLen); % for sheep data
% ShowOnlyFigure(keepDataSeismicSeriesLen); % for seismic data
% ShowOnlyFigure(keepDataRandomWalkSeriesLen); % for randomwalk data

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
xl = xlabel('Time Series Length');
yl = ylabel('Time Needed (sec.)');
set(xl,'FontSize',12,'FontWeight','bold','FontName','Courier');
set(yl,'FontSize',14,'FontWeight','bold','FontName','Courier');
hold off;
disp('see me');
end

function subSequenceExpProteinData()
clear

subSeqLen = 50;  % The length of the sub-sequence
load ('ProteinData.mat');

[keepAllTargetTogether, ~] = ConcatenateAllSeries(realData(1:100,:)); % concatenate all the time series

takenRws = 20;
while (takenRws < 150)
    fprintf('Protien Data : The numer of taken rows are  : %d \n', takenRws);
    
    apply_AAMP(keepAllTargetTogether, subSeqLen);   % working for the 2nd column here
    takenRws = takenRws + 20;
end
end


function subSequenceExpSheepData()
clear

% Loading the sheep data
load('SheepDataFull.mat')


[keepAllTargetTogether, ~] = ConcatenateAllSeries(keepAllData(1:100,:)); % concatenate all the time series

subSeqLen = 50;  % The length of the sub-sequence
takenRws = 20;
while (takenRws < 150)
    fprintf('Sheep Data : The numer of taken rows are %d \n', takenRws);
    apply_AAMP(keepAllTargetTogether, subSeqLen);   % working for the 2nd column here
    
    takenRws = takenRws + 20;
end
end


function subSequenceExpSeismicData()
clear

% Loading the sheep data
load('seismic_50000.mat');

keepAllData = seeMe1;
clear seeMe1;

queryIndex = 307;
queryComplete = keepAllData(queryIndex, 1:end-1);

noOfZerosInQuery = sum(queryComplete(:)==0);
if(noOfZerosInQuery > (length(queryComplete) /2))
    error ('The number of zeros in the subsequence is more than the half of \n query size. Please choose another query');
end

subSeqLen = 50;  % The length of the sub-sequence
takenRws = 20;
while (takenRws < 150)
    fprintf('Seismic Data :  The numer of taken rows are : %d \n', takenRws);
    
    apply_AAMP(keepAllData(1:takenRws,:), subSeqLen);   % working for the 2nd column here
    takenRws = takenRws + 20;
end

end


function subSequenceExpRandomWalkData()
clear

load('randomWalk_50000.mat');
keepAllData = seeMe1;
clear seeMe1;

queryIndex = 307;
queryComplete = keepAllData(queryIndex, 1:end-1);

noOfZerosInQuery = sum(queryComplete(:)==0);
if(noOfZerosInQuery > (length(queryComplete) /2))
    error ('The number of zeros in the subsequence is more than the half of \n query size. Please choose another query');
end

subSeqLen = 50;  % The length of the sub-sequence
takenRws = 20;
while (takenRws < 150)
    fprintf('Random Walk Data : The numer of taken rows are : %d \n', takenRws);
    
    apply_AAMP(keepAllData(1:takenRws,:), subSeqLen);   % working for the 2nd column here
    takenRws = takenRws + 20;
end

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
        if(strcmp(newStr{8}, 'rows'))
            seriesLen = str2double(newStr{12});
        end
        keepData(lnCnt, 1) = seriesLen;
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
        if(strcmp(newStr{8}, 'rows'))
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

function keepData = ReadFile_3(textFilePath)
keepData = zeros(1,2);
fid = fopen(textFilePath);
tline = fgetl(fid);
lnCnt = 1;
while ischar(tline)
    if(~isempty(tline))
        newStr = split(tline,' ');
        %     subSeqlen = (newStr{6});
        if(strcmp(newStr{9}, 'rows'))
            subSeqlen = str2double(newStr{12});
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