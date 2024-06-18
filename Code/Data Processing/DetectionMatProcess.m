% Author: Yanqi Hong, Alexandros Theocharous
% DetectionMatProcess - Process the raw data using a detection matrix.
%
% Syntax:
%   ProcessedData = DetectionMatProcess(RawData, DetectionMat, Norm_flag)
%
% Inputs:
%   - RawData: The raw data to be processed.
%   - DetectionMat: The detection matrix to be applied to the raw data.
%   - Norm_flag: A flag indicating whether to normalize the output in the detection region.
%
% Output:
%   - ProcessedData: The processed data after applying the detection matrix.
%
% Description:
%   This function processes the raw data by applying a detection matrix. The detection matrix is multiplied element-wise with the raw data, and the result is returned as the processed data. If the Norm_flag is set to 1, the output in the detection region is normalized.
%
% Example:
%   raw_data = [1 2 3; 4 5 6; 7 8 9];
%   detection_mat = [1 0 1; 0 1 0; 1 0 1];
%   norm_flag = 1;
%   processed_data = DetectionMatProcess(raw_data, detection_mat, norm_flag);
%
%   The output processed_data will be:
%   processed_data = [0 0 0; 0 5 0; 0 0 0];
%
%   If norm_flag is set to 0, the output will not be normalized.
function ProcessedData = DetectionMatProcess(RawData,DetectionMat, Norm_flag)

    ProcessedData = flip(DetectionMat .* flip(RawData));

    if Norm_flag == 1   % Normalize the output in the detection region
        Data_nonzero = ProcessedData;
        Data_nonzero(Data_nonzero == 0) = Inf;
        min_value=min(Data_nonzero(:));
        ProcessedData = ((ProcessedData - min_value).*flip(DetectionMat))/ (max(ProcessedData(:)) - min_value) ;
    end
end