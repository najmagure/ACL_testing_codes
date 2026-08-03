clear
clc
close all

%% SIMPLE AUTOMATIC DIC–ACUMEN STRAIN PLOT
%
% This version keeps the structure of the original StrainPlot.m script.
%
% Automatic processing:
%   1. DIC: find the final step-and-hold event lasting 20–30 seconds.
%   2. DIC: remove everything before the loading step and after unloading.
%   3. Acumen: detect and remove the first 10 preconditioning cycles.
%   4. Acumen: find the main loading step after preconditioning.
%   5. Acumen: retain no more than 30 seconds.
%   6. Align both cropped signals at t = 0 and calculate strain.
%
% Fixed acquisition frequencies:
dicFs = 25;
acFs = 100;

% Fixed test-protocol settings:
numberOfPrecycles = 10;
minimumHoldTime = 20;
maximumHoldTime = 30;

%% Acumen resampling frequency
% This is the only analysis-settings dialog.

answer = inputdlg( ...
    {'Acumen resampling frequency (Hz)'}, ...
    'Analysis Settings', ...
    [1 45], ...
    {'25'});

if isempty(answer)
    error('Analysis cancelled.');
end

targetFs = str2double(answer{1});

if ~isfinite(targetFs) || targetFs <= 0
    error('The resampling frequency must be a positive number.');
end

%% Select and read DIC file

[fileDIC,pathDIC] = uigetfile( ...
    {'*.csv;*.xlsx','DIC File (*.csv, *.xlsx)'}, ...
    'Select DIC File');

if isequal(fileDIC,0)
    error('No DIC file selected.');
end

dicFile = fullfile(pathDIC,fileDIC);
raw = readcell(dicFile);

% DIC coordinate data begin in row 3, column 5.
DIC = cell2mat(raw(3:end,5:end));

%% Insertion centroid calculation

xIns = mean(DIC([1 4 7],:),1,'omitnan');
yIns = mean(DIC([2 5 8],:),1,'omitnan');
zIns = mean(DIC([3 6 9],:),1,'omitnan');

%% Origin centroid calculation

xOrg = mean(DIC([19 22 25],:),1,'omitnan');
yOrg = mean(DIC([20 23 26],:),1,'omitnan');
zOrg = mean(DIC([21 24 27],:),1,'omitnan');

%% DIC distance calculation

distance = sqrt( ...
    (xIns-xOrg).^2 + ...
    (yIns-yOrg).^2 + ...
    (zIns-zOrg).^2 );

distance = distance(:);
distance = fillmissing(distance,'linear','EndValues','nearest');

%% Automatically remove DIC ramp tests and unwanted sections

% Smooth only for detecting sharp loading and unloading steps.
distanceSmooth = smoothdata(distance,'movmean',5);

% Change in distance from one frame to the next.
dicChange = gradient(distanceSmooth);

% Detect substantial positive loading steps and negative unloading steps.
stepProminence = 0.10 * range(distanceSmooth);
minimumStepSpacing = round(5*dicFs);

[positiveSteps,positiveFrames] = findpeaks( ...
    dicChange, ...
    'MinPeakProminence',stepProminence, ...
    'MinPeakDistance',minimumStepSpacing);

[~,negativeFrames] = findpeaks( ...
    -dicChange, ...
    'MinPeakProminence',stepProminence, ...
    'MinPeakDistance',minimumStepSpacing);

% A stress-relaxation test is defined here as:
%   positive step -> 20 to 30 second hold -> negative step.
validStarts = [];
validEnds = [];
validPeakValues = [];

for i = 1:length(positiveFrames)

    earliestEnd = positiveFrames(i) + round(minimumHoldTime*dicFs);
    latestEnd = positiveFrames(i) + round(maximumHoldTime*dicFs);

    possibleEnds = negativeFrames( ...
        negativeFrames >= earliestEnd & ...
        negativeFrames <= latestEnd);

    if ~isempty(possibleEnds)
        validStarts(end+1,1) = positiveFrames(i); %#ok<SAGROW>
        validEnds(end+1,1) = possibleEnds(1); %#ok<SAGROW>
        validPeakValues(end+1,1) = positiveSteps(i); %#ok<SAGROW>
    end
end

if isempty(validStarts)
    error(['No DIC step-and-hold event lasting 20–30 seconds was found. ' ...
        'Check the DIC data or the fixed timing settings.']);
end

% Use the final valid step-and-hold event. This removes earlier tests and
% excludes the later ramp test.
dicStepFrame = validStarts(end);
dicEndFrame = validEnds(end)-1;
dicStepValue = validPeakValues(end);

% Move backward from the derivative peak to the beginning of the loading
% edge so the sharp step begins at t = 0.
dicStartFrame = dicStepFrame;
startThreshold = 0.05*dicStepValue;

while dicStartFrame > 1 && ...
        dicChange(dicStartFrame) > startThreshold
    dicStartFrame = dicStartFrame-1;
end

% Reference length is the median distance during the 0.5 seconds directly
% before the detected loading step.
referenceWindow = round(0.5*dicFs);
referenceStart = max(1,dicStartFrame-referenceWindow);

L0 = median(distance(referenceStart:dicStartFrame),'omitnan');

% Crop DIC to only the stress-relaxation step and hold.
distance = distance(dicStartFrame:dicEndFrame);

dicDisplacement = distance-L0;
dicStrain = dicDisplacement./L0;
dicTime = (0:length(distance)-1)'/dicFs;

fprintf('\nDIC stress-relaxation frames: %d to %d\n', ...
    dicStartFrame,dicEndFrame);
fprintf('DIC duration: %.2f s\n',dicTime(end));
fprintf('Reference length L0: %.4f mm\n',L0);

%% Select and read Acumen file

[fileAC,pathAC] = uigetfile( ...
    {'*.txt;*.csv;*.xlsx','Acumen File (*.txt, *.csv, *.xlsx)'}, ...
    'Select Acumen File');

if isequal(fileAC,0)
    error('No Acumen file selected.');
end

ACFile = fullfile(pathAC,fileAC);
[~,~,acExtension] = fileparts(ACFile);

if strcmpi(acExtension,'.txt')

    % Acumen text export:
    % row 1 = names, row 2 = units, row 3 onward = numeric data.
    opts = delimitedTextImportOptions('NumVariables',3);
    opts.DataLines = [3 Inf];
    opts.Delimiter = '\t';
    opts.VariableNames = ...
        {'AxialForce','AxialDisplacement','Time'};
    opts.VariableTypes = {'double','double','double'};
    opts.ExtraColumnsRule = 'ignore';

    ACTable = readtable(ACFile,opts);

else
    ACTable = readtable(ACFile);
end

% Locate columns without depending on spaces in their names.
columnNames = lower(regexprep( ...
    ACTable.Properties.VariableNames,'[^a-zA-Z0-9]',''));

forceColumn = find(contains(columnNames,'axialforce'),1);
dispColumn = find(contains(columnNames,'axialdisplacement'),1);

if isempty(forceColumn) || isempty(dispColumn)
    error('Axial Force or Axial Displacement was not found.');
end

acForce = ACTable{:,forceColumn};
acDisp = ACTable{:,dispColumn};

acForce = acForce(:);
acDisp = acDisp(:);

validRows = isfinite(acForce) & isfinite(acDisp);
acForce = acForce(validRows);
acDisp = acDisp(validRows);

%% Automatically remove Acumen preconditioning cycles

% Smooth only for event detection.
acDispSmooth = smoothdata(acDisp,'movmean',15);

% Locate the first 10 preconditioning peaks.
[~,precycleFrames] = findpeaks( ...
    acDispSmooth, ...
    'MinPeakProminence',0.20, ...
    'MinPeakDistance',round(1.0*acFs));

if length(precycleFrames) < numberOfPrecycles
    error('Fewer than 10 Acumen preconditioning peaks were detected.');
end

lastPrecycleFrame = precycleFrames(numberOfPrecycles);

% Find the strongest positive loading step after the preconditioning cycles.
acChange = gradient(acDispSmooth);
searchStart = lastPrecycleFrame + round(1.0*acFs);

[maximumStep,relativeFrame] = max(acChange(searchStart:end));
acStepFrame = searchStart + relativeFrame-1;

% Move backward to the beginning of the loading edge.
acStartFrame = acStepFrame;
acStartThreshold = 0.05*maximumStep;

while acStartFrame > searchStart && ...
        acChange(acStartFrame) > acStartThreshold
    acStartFrame = acStartFrame-1;
end

% Keep no more than 30 seconds after the loading step.
acEndFrame = min( ...
    length(acDisp), ...
    acStartFrame + round(maximumHoldTime*acFs));

% Crop and zero the Acumen displacement.
acDisp = acDisp(acStartFrame:acEndFrame);
acForce = acForce(acStartFrame:acEndFrame);
acDisp = acDisp-acDisp(1);

fprintf('\nAcumen stress-relaxation samples: %d to %d\n', ...
    acStartFrame,acEndFrame);

%% Resample Acumen

[p,q] = rat(targetFs/acFs);

dispResampled = resample(acDisp,p,q);
forceResampled = resample(acForce,p,q);

timeResampled = (0:length(dispResampled)-1)'/targetFs;
acumenStrain = dispResampled./L0;

%% Match the final durations

commonDuration = min([ ...
    dicTime(end), ...
    timeResampled(end), ...
    maximumHoldTime]);

dicKeep = dicTime <= commonDuration;
acKeep = timeResampled <= commonDuration;

dicTime = dicTime(dicKeep);
dicStrain = dicStrain(dicKeep);

timeResampled = timeResampled(acKeep);
acumenStrain = acumenStrain(acKeep);
forceResampled = forceResampled(acKeep);

%% Plot

figure

plot(dicTime,dicStrain, ...
    'LineWidth',2)

hold on

plot(timeResampled,acumenStrain, ...
    'LineWidth',2)

xlabel('Time (s)')
ylabel('Strain')
title('DIC Strain vs Acumen Strain')
legend('DIC','Acumen','Location','best')
grid on
xlim([0 commonDuration])
