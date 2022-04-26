%% spike detection main function (32-channels)

clear all;
% close all;
fclose('all');

%% general setting
mode = 2;                                               % 1 for software, 2 for hardware
test_type = 1;                                          % 1 for SINICA data, 2 for synthetic data (https://leicester.figshare.com/articles/dataset/Simulated_dataset/11897595)
result_dir = 'D:\kilosort2\SINICA\32channels_design\';  % output file directory (%'D:\kilosort2\SINICA\32channels_SWdesign\';%'D:\kilosort2\Simulator\results\')
thS_factor = 8;                                         % https://ieeexplore.ieee.org/document/5740375 (mu * thS_factor)
thH_factor = 5;                                         % https://ieeexplore.ieee.org/document/6070974 (mu + thH_factor*std)
sigma = 3;                                              % NEO parameter (neo_value = data[n]^2 - data[n-sigma]*data[n+sigma])

%% SINICA or synthetic
if test_type == 1                                       % for SINICA data
    trial = 40;                                         % if trail = 40, B3D71S40.h5 is used
    tstart = 0;                                         % start time
    tcount = 10;                                        % duration
    ch = 32;                                            % total channel number
    plot_ch_begin = 1;                                  % analyze from which channel
    plot_ch_end = 32;                                   % analyze end at which channel
else                                                    % for synthetic data
    test_data = 'C_Difficult1_noise01';
    result_dir = strcat(result_dir, test_data, '\');
    tstart = 0;
    tcount = 1;
    plot_ch_begin = 1;
    plot_ch_end = 1;
    ch = 1;
end

%% files management (CAREFUL)
if ~isfolder(result_dir)                                % if the result_dir does not exist, create one.
    disp('Creating new directory...');
    mkdir(result_dir);
else
    disp('Deleting existing txt files...');             % if the result_dir exists, delete all the .txt files in it.
    delete(sprintf('%s', strcat(result_dir, '*.txt')));    
end

%% read raw data (this part should be modified under different file hierarchy)
disp('Loading raw data...');
if test_type == 1
    [dataRAW, samplingInterval] = rmdata(strcat('D:\kilosort2\SINICA\data\B3D71S', num2str(trial), '.h5'), ch, tstart, tcount);
else    
    load(strcat(result_dir(1:23), 'Simulator\', test_data, '.mat'));
    if ~exist('tcount')
        tcount = size(data, 2)*samplingInterval/1000;
        dataRAW = data;
    else
        dataRAW = data(1, 1:(tstart+tcount)/samplingInterval*1000);
    end
    clear data OVERLAP_DATA spike_class spike_times startData
end

fs = 1000/samplingInterval;
time = tstart:1/fs:tstart+tcount-1/fs;

%% Processing raw data
tic
step = 'raw';
% % Calculating NEO for RAW data
% for i = plot_ch_begin:1:plot_ch_end
%     disp(strcat('Processing data, calculating NEO1 for ch', num2str(i), '...'));    
%     exe_neo(result_dir, step, i, sigma, dataRAW(i, :), thS_factor, thH_factor);
% end

%% substract mean inside single channel
% dataRAW = dataRAW - mean(dataRAW, 1);

%% substract median cross channels
if test_type == 1
    disp(strcat('Processing data, subtracting median...'));    
    if mode == 2
        for i = 1:1:size(time, 2)                
            dataRAW(:, i) = dataRAW(:, i) - median(dataRAW(:, i));  % model hw behavior
        end
    else
        dataRAW = dataRAW - median(dataRAW, 1);                     % model sw behavior
    end
%     % Calculating NEO for data after substracting median value
%     step = 'median';
%     for i = plot_ch_begin:1:plot_ch_end
%         disp(strcat('Processing data, calculating NEO2 for ch', num2str(i), '...'));    
%         exe_neo(result_dir, step, i, sigma, dataRAW(i, :), thS_factor, thH_factor);
%     end
else
    disp('Running single-channel data, no need for substracting median...');    
    % there is only one channel in synthetic data, no need to do cross-channel preprocessing
end

%% high pass filter with cuttoff frequency = 150 Hz
b1 = [0.9697, -2.9092, 2.9092, -0.9697];
a1 = [1, -2.9385, 2.8790, -0.9404];

disp('Processing data, filtering...');
dataRAW = filter(b1, a1, dataRAW);

% Calculating NEO for data after substracting median and high pass filter
step = 'filter';
for i = plot_ch_begin:1:plot_ch_end
    disp(strcat('Processing data, calculating NEO3 for ch', num2str(i), '...'));        
    exe_neo(result_dir, step, i, sigma, dataRAW(i, :), thS_factor, thH_factor);    
end

%% time-reverse high pass filter with cuttoff frequency = 150 Hz
% disp('Processing data, reverse filtering...');
% dataRAW = filter(b1, a1, flipud(dataRAW));
% dataRAW = flipud(dataRAW);

% % Calculating NEO for data after substracting median, high pass filter and reverse high pass filter
% step = 're-filter';
% for i = plot_ch_begin:1:plot_ch_end
%     disp(strcat('Processing data, calculating NEO4 for ch', num2str(i), '...'));
%     exe_neo(result_dir, step, i, sigma, dataRAW(i, :), thS_factor, thH_factor);
% end

%% normalized by STD
% disp('Processing data, normalizing by STD...');
% dataRAW = dataRAW ./ std(dataRAW, 1, 1);

% % Calculating NEO for data after substracting median, high pass filter, reverse high pass filter and normalization
% step = 'normal';
% for i = plot_ch_begin:1:plot_ch_end
%     disp(strcat('Processing data, calculating NEO5 for ch', num2str(i), '...'));
%     exe_neo(result_dir, step, i, sigma, dataRAW(i, :), thS_factor, thH_factor);
% end

%% spike alignment
disp('Spike alignment...');
train = zeros(ch, size(time, 2));
loc = zeros(ch, size(time, 2));
Strain = zeros(ch, size(time, 2));
Sloc = zeros(ch, size(time, 2));

for i = plot_ch_begin:1:plot_ch_end
    T = readtable(strcat(result_dir, 'filter_', num2str(i), '.txt'));
    [train(i, :), loc(i, :)] = alignment(T.RAW, T.NEO, T.THH);
    [Strain(i, :), Sloc(i, :)] = alignment(T.RAW, T.SNEO, T.STHH);
end

loc_num = sum(loc, 2)';
train_num = sum(train, 2)';
Sloc_num = sum(Sloc, 2)';
Strain_num = sum(Strain, 2)';

disp('Done');
toc

%% This part is just ploting the detection results.

% plot_ch = 2;
% T = readtable(strcat(result_dir, 'filter_', num2str(plot_ch), '.txt'));
% % figure(1);
% % plot(time, T.RAW, time, T.NEO, time, T.THS, time, T.THH); legend('RAW', 'NEO', 'THS', 'THH');
% %
% figure(plot_ch);
% plot(time, T.RAW);
% for i = 1:1:size(Strain, 2)
%     if Strain(plot_ch, i) == 1
%         xline(i/fs, '--r');
% %         disp(strcat('plot', num2str(i), 'at', num2str(Strain(i)));
%     end
% end
% 
% for i = 1:1:size(Sloc, 2)
%     if Sloc(plot_ch, i) == 1
%         xline(i/fs, '--b');
%     end
% end
% title('window = 64');
        