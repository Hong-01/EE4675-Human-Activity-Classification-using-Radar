% Author: Alexandros Theocharous, Yanqi Hong 

% This script processes radar data and performs object detection using CFAR (Constant False Alarm Rate) algorithm.
% It loads radar data from a specified folder, applies range-time, doppler-time, and doppler-range processing techniques,
% and generates RGB images representing the processed data. The script also provides options to save the detection results
% and the generated images (256*256).

clear; clc; close all;
%%-------------------------------------------------------------------------%
% process_mode choice, NOTE: the cfar detection process is time consuming and it is recommended to save the detection result
% 1: Load data and process cfar detection
% 2: Load data and process with saved cfar detection matrix result
process_mode = 2;

% save image choice
% 0: do not save image
% 1: save image
save_image = 0;

% save cfar detection result choice, only valid when process_mode = 1
% 0: do not save cfar detection result
% 1: save cfar detection result
save_detection = 0;

%Path 
folder_path = 'E:/DATA/TUD/Master/TUD_Master_Y1/Q4/EE4675 Object Classification with Radar/Project/Dataset_848';
detection_DR_path = 'E:/DATA/TUD/Master/TUD_Master_Y1/Q4/EE4675 Object Classification with Radar/Project/detection_DR.mat';   %path of detection_DR (DR: Doppler-Range)
detection_DT_path = 'E:/DATA/TUD/Master/TUD_Master_Y1/Q4/EE4675 Object Classification with Radar/Project/detection_DT.mat';   %path of detection_DT (DT: Doppler-Time)
detection_RT_path = 'E:/DATA/TUD/Master/TUD_Master_Y1/Q4/EE4675 Object Classification with Radar/Project/detection_RT.mat';   %path of detection_RT (RT: Range-Time)
image_save_path = 'E:/DATA/TUD/Master/TUD_Master_Y1/Q4/EE4675 Object Classification with Radar/Project/image_data/image_data';   %path of image save

% Load data
inputFolder = folder_path;
inputFiles = dir(fullfile(inputFolder, '*.dat'));

% Initialize variables
if process_mode == 1
  detection_DR = zeros(length(inputFiles), 256, 256);
  detection_DT = zeros(length(inputFiles), 256, 256);
  detection_RT = zeros(length(inputFiles), 256, 256);
end

if process_mode==2
  detection_DR = load(detection_DR_path).detection_DR;
  detection_DT = load(detection_DT_path).detection_DT;
  detection_RT = load(detection_RT_path).detection_RT;
end

%%-------------------------------------------------------------------------%
%  Process data
for z =1:length(inputFiles)    % z: number of files in the folder
  clear Data_range
  % progress bar, print out the file number being processed
  if mod(z, 50) == 0
    fprintf('Processing file %d of %d\n', z, length(inputFiles));
  end

    % Input file
    inputFile = fullfile(inputFolder, inputFiles(z).name);
    radarData = dlmread(inputFile);

    % Define variables 
    clearvars fileID dataArray ans;
    fc = radarData(1); % Center frequency
    Tsweep = radarData(2); % Sweep time in ms
    Tsweep=Tsweep/1000; %then in sec
    NTS = radarData(3); % Number of time samples per sweep
    Bw = radarData(4); % FMCW Bandwidth. For FSK, it is frequency step;
    % For CW, it is 0.
    Data = radarData(5:end); % raw data in I+j*Q format
    fs=NTS/Tsweep; % sampling frequency ADC
    record_length=length(Data)/NTS*Tsweep; % length of recording in s
    nc=record_length/Tsweep; % number of chirps
    
    %-------
    %Range-Time Domian
    %-------

    %% Reshape data into chirps and plot Range-Time
    Data_time=reshape(Data, [NTS nc]);
    win = ones(NTS,size(Data_time,2));
    %Part taken from Ancortek code for FFT and IIR filtering
    tmp = fftshift(fft(Data_time.*win),1);
    Data_range(1:NTS/2,:) = tmp(NTS/2+1:NTS,:);
    ns = oddnumber(size(Data_range,2))-1;
    Data_range_MTI = zeros(size(Data_range,1),ns);
    [b,a] = butter(4, 0.0075, 'high');
    [h, f1] = freqz(b, a, ns);
    for k=1:size(Data_range,1)
      Data_range_MTI(k,1:ns) = filter(b,a,Data_range(k,1:ns));
    end
    freq =(0:ns-1)*fs/(2*ns); 
    range_axis=(freq*3e8*Tsweep)/(2*Bw);
    Data_range_MTI=Data_range_MTI(2:size(Data_range_MTI,1),:);
    Data_range=Data_range(2:size(Data_range,1),:);
    Data_RT=20*log10(abs(Data_range_MTI));  % log compression
    Data_RT = imresize(Data_RT,[256,256]);  % resize to 256x256
    if process_mode == 1
      % CFAR processing
      [Data_RT_norm,detections]=CFAR(Data_RT,[30,6],[100,20],0.3,200,1);
      %save detections
      detection_RT(z, :, :) = detections;
    end

    if process_mode == 2
      Data_RT_norm = DetectionMatProcess(Data_RT,squeeze(detection_RT(z,:,:)), 1);
    end

    %-------
    %Doppler-Time Domain
    %-------

    % Spectrogram processing for 2nd FFT to get Doppler
    % This selects the range bins where we want to calculate the spectrogram
    bin_indl = 10;
    bin_indu = 30;
    
    MD.PRF=1/Tsweep;
    MD.TimeWindowLength = 200;
    MD.OverlapFactor = 0.95;
    MD.OverlapLength = round(MD.TimeWindowLength*MD.OverlapFactor);
    MD.Pad_Factor = 4;
    MD.FFTPoints = MD.Pad_Factor*MD.TimeWindowLength;
    MD.DopplerBin=MD.PRF/(MD.FFTPoints);
    MD.DopplerAxis=-MD.PRF/2:MD.DopplerBin:MD.PRF/2-MD.DopplerBin;
    MD.WholeDuration=size(Data_range_MTI,2)/MD.PRF;
    MD.NumSegments=floor((size(Data_range_MTI,2)-MD.TimeWindowLength)/floor(MD.TimeWindowLength*(1-MD.OverlapFactor)));
      
    Data_spec_MTI2=0;
    Data_spec2=0;
    for RBin=bin_indl:1:bin_indu
        Data_MTI_temp = fftshift(spectrogram(Data_range_MTI(RBin,:),MD.TimeWindowLength,MD.OverlapLength,MD.FFTPoints),1);
        Data_spec_MTI2=Data_spec_MTI2+abs(Data_MTI_temp);                                
        Data_temp = fftshift(spectrogram(Data_range(RBin,:),MD.TimeWindowLength,MD.OverlapLength,MD.FFTPoints),1);
        Data_spec2=Data_spec2+abs(Data_temp);
    end
    MD.TimeAxis=linspace(0,MD.WholeDuration,size(Data_spec_MTI2,2));
    Data_spec_MTI2=flipud(Data_spec_MTI2);
    Data_DT=20*log10(Data_spec_MTI2);
    Data_DT = imresize(Data_DT,[256,256]);

    if process_mode == 1
      % CFAR processing
      [Data_DT_norm,detections]=CFAR(Data_DT,[10,2],[50,10],0.31,100,1);
      %save detections
      detection_DT(z, :, :) = detections;
    end

    if process_mode == 2
      Data_DT_norm = DetectionMatProcess(Data_DT,squeeze(detection_DT(z,:,:)), 1);
    end

    %-------
    %Doppler-Range Domain
    %-------

    % Reshape data into range-Doppler matrix
    Data_DR = transpose(fftshift(fft(Data_range_MTI,[],2),2));

    Data_DR=20*log10(abs(Data_DR));
    Data_DR = imresize(Data_DR,[256,256]);

    if process_mode == 1
      % CFAR processing
      [Data_DR_norm,detections]=CFAR(Data_DR,[30,8],[50,15],0.31,100,1);
      %save detections
      detection_DR(z, :, :) = detections;
    end

    if process_mode == 2
      Data_DR_norm=DetectionMatProcess(Data_DR,squeeze(detection_DR(z,:,:)), 1);
    end
 
    % Move the data to the top of the image
    Data_DR_norm = padarray(Data_DR_norm, [70, 0], 'pre');
    Data_DR_norm = Data_DR_norm(1:end-70, :);

    %-------
    %RGB Image and Save Image
    %-------

    [height, width] = size(Data_RT_norm);
    rgb_image = zeros(height, width, 3); 
    
    rgb_image(:,:,1) = Data_DT_norm; % Red channel

    rgb_image(:,:,2) = Data_DR_norm; % Green channel

    rgb_image(:,:,3) = Data_RT_norm; % Blue channel

    % Extract label
    [number1,number2,number3] = Label_extract4(inputFile);  %number2: type of movement
    number2=cell2mat(number2);
    if length(number2)~=2
        number2 = char('0',number2);
    end

    movement = {'walk', 'sit', 'stand',  'pick', 'drink','fall'};
    movement = movement(str2num(number2(2)));

    % save image
    if save_image == 1
      % save image
      imwrite(flip(rgb_image), sprintf('%s/%s/image_%s.png', image_save_path, cell2mat(movement), num2str(z)));
    end
    % imwrite(flip(rgb_image), sprintf('C:/Users/alext/Desktop/radar/Dataset_848/image_data/%s/image_%s.png', cell2mat(number2), num2str(z)));
    % imwrite(flip(rgb_image), sprintf('E:/DATA/TUD/Master/TUD_Master_Y1/Q4/EE4675 Object Classification with Radar/Project/image_data_phase/image_data/%s/image_%s.png',  cell2mat(movement), num2str(z)));
end

% save detections
if process_mode == 1 && save_detection == 1
  % save detections
  save(detection_DR_path, 'detection_DR');
  save(detection_DT_path, 'detection_DT');
  save(detection_RT_path, 'detection_RT');
end



