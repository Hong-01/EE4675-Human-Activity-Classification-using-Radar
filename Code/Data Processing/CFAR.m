% Author: Yanqi Hong, Alexandros Theocharous
% CFAR - Constant False Alarm Rate detection algorithm
%
% Syntax: [RT,dections] = CFAR(Data, guardBand, trainingBand, False_alarm, copy_edge, flag_norm)
%
% Inputs:
%    Data - Radar data matrix
%    guardBand - Size of the guard band [range, doppler]
%    trainingBand - Size of the training band [range, doppler]
%    False_alarm - Probability of false alarm
%    copy_edge - Number of rows/columns to copy from the edge to handle edge effects
%    flag_norm - Flag to normalize the output in the detection region, 1 to normalize, 0 to not normalize
%
% Outputs:
%    RT - Range-time matrix with CFAR detections
%    detections - Detection matrix
% Example:
%    Data = rand(100, 100); % Example radar data
%    guardBand = [30, 6]; % Guard band size
%    trainingBand = [100, 20]; % Training band size
%    False_alarm = 0.3; % Probability of false alarm
%    copy_edge = 50; % Number of rows/columns to copy from the edge
%    flag_norm = 1; % Normalize the output
%    [RT,detections] = CFAR(Data, guardBand, trainingBand, False_alarm, copy_edge,1); % Run CFAR detection

function [RT,detections] = CFAR(Data, guardBand, trainingBand, False_alarm, copy_edge,flag_norm)

    Data_range_MTI = Data;
    
    % Copy rows/columns from the edge to handle edge effects
    Data_range_MTI_edge = Data_range_MTI;
    Data_range_MTI_edge = [Data_range_MTI(end-copy_edge+1:end, :); Data_range_MTI_edge];
    Data_range_MTI_edge = [Data_range_MTI_edge; Data_range_MTI(1:copy_edge, :)];
    Data_range_MTI_edge = [Data_range_MTI_edge(:, end-copy_edge+1:end), Data_range_MTI_edge];
    Data_range_MTI_edge = [Data_range_MTI_edge, Data_range_MTI_edge(:, 1:copy_edge)];
    
    radar_data = Data_range_MTI_edge;
    
    % Create CFAR detector object
    cfar2D = phased.CFARDetector2D('GuardBandSize', guardBand, 'TrainingBandSize', trainingBand, ...
        'ProbabilityFalseAlarm', False_alarm, 'ThresholdFactor', 'Auto');
    
    JJ = size(radar_data, 2);
    II = size(radar_data, 1);
    Grid = transpose(1:max(size(Data)));

    totalGuardTraining = max(guardBand) + max(trainingBand);

    validRangeStart = totalGuardTraining + 1;
    validRangeEnd = II - totalGuardTraining;

    validTimeStart = totalGuardTraining + 1;
    validTimeEnd = JJ - totalGuardTraining;

    % Ensure the indices are within the valid range
    [~, rangeIndx] = min(abs(Grid - [0 (II-1)]));
    rangeIndx(1) = max(rangeIndx(1), validRangeStart);
    rangeIndx(2) = max(rangeIndx(2), validRangeEnd);

    [~, timeIndx] = min(abs(Grid - [0 (JJ-1)]));
    timeIndx(1) = max(timeIndx(1), validTimeStart);
    timeIndx(2) = max(timeIndx(2), validTimeEnd);

    % Generate indices for CUT cells within the valid range
    [columnInds, rowInds] = meshgrid(timeIndx(1):timeIndx(2), rangeIndx(1):rangeIndx(2));

    % Run CFAR detector
    CUTIdx = [rowInds(:) columnInds(:)]';
    detections = cfar2D(Data_range_MTI_edge, CUTIdx);
    detections = reshape(detections, [sqrt(length(detections)), sqrt(length(detections))]);
    diff1 = (size(detections, 1) - size(Data, 1)) / 2;
    diff2 = (size(detections, 2) - size(Data, 2)) / 2;
    detections = detections(diff1+1:end-diff1, diff2+1:end-diff2);
    detections = flip(detections);

    RT = flip(detections .* flip(Data_range_MTI));

    if flag_norm == 1   % Normalize the output in the detection region
        RT_nonzero = RT;
        RT_nonzero(RT_nonzero == 0) = Inf;
        min_value=min(RT_nonzero(:));
        RT = ((RT - min_value).*flip(detections))/ (max(RT(:)) - min_value) ;
    end
end
