clc; clear; close all;

%% ================== PARAMETERS ==================
ROLL_NUMBER   = 244; 
fs            = 2.4e6;             
duration      = 5;           
samples_per_frame = 256*1024;
output_bin    = 'selected_station.bin';
meta_file     = 'record_meta.mat';

fm_start  = 88e6;
fm_stop   = 108e6;
freq_step = 2.0e6; 
sweep_freqs = fm_start : freq_step : fm_stop;

%% ================== SDR SETUP ==================
try
    rx = comm.SDRRTLReceiver( ...
        'SampleRate', fs, ...
        'SamplesPerFrame', samples_per_frame, ...
        'OutputDataType', 'single');
    disp("SDR Initialized successfully.");
catch
    error("SDR not found. Please check connection.");
end

%% ================== WIDEBAND PSD SCAN ==================
disp("Step 1: Sweeping FM Band (88-108 MHz)...");
wide_psd = [];
wide_freq = [];

for f = sweep_freqs
    rx.CenterFrequency = f;
    for i=1:5, rx(); end % Flush buffer
    data = rx();
    
    [Pxx, F] = pwelch(data, 1024, 512, 1024, fs, 'centered');
    wide_psd = [wide_psd; 10*log10(Pxx + 1e-12)];
    wide_freq = [wide_freq; (F + f)];
end

[wide_freq, sortIdx] = sort(wide_freq);
wide_psd = wide_psd(sortIdx);

%% ================== IDENTIFY STATIONS ==================
stations = [91.1e6; 91.9e6; 92.7e6; 93.5e6; 94.3e6; 98.3e6; 100.1e6; 101.3e6; 102.9e6; 104.0e6];
power_vals = zeros(length(stations), 1);

for k = 1:length(stations)
    [~, idx_f] = min(abs(wide_freq - stations(k)));
    bin_range = floor(100e3 / (fs/1024)); 
    power_vals(k) = max(wide_psd(max(1,idx_f-bin_range) : min(length(wide_psd),idx_f+bin_range)));
end

[sorted_power, idx] = sort(power_vals, 'descend');
top5_indices = idx(1:5);

% Roll Number logic
station_target_idx = mod(ROLL_NUMBER, 5) + 1;
selected_fc = stations(top5_indices(station_target_idx));

fprintf('\nSelected (Index %d): %.2f MHz\n', station_target_idx, selected_fc/1e6);

%% ================== RECORD SELECTED STATION ==================
rx.CenterFrequency = selected_fc;
num_frames = ceil((fs * duration) / samples_per_frame);

psd_matrix = [];
time_axis = [];
recorded_iq_data = []; 
frame_time = samples_per_frame / fs;

% Capture start timestamp
capture_timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

fprintf('\nRecording for %d seconds starting at %s...\n', duration, capture_timestamp);
fid = fopen(output_bin,'wb');

for k = 1:num_frames
    data = rx();
    if isempty(data), continue; end
    
    % FIX: Remove invalid SDR samples
    data(~isfinite(data)) = 0;
    
    data(abs(data) > 10) = 0;
    
    recorded_iq_data = [recorded_iq_data; data];
    
    % Save binary
    interleaved = zeros(2*length(data),1,'single');
    interleaved(1:2:end) = real(data);
    interleaved(2:2:end) = imag(data);
    fwrite(fid, interleaved, 'float32');
    
    % Spectrogram data
    [Pxx_slice, F_slice] = pwelch(data, 2048, 1024, 2048, fs, 'centered');
    psd_matrix(:,k) = 10*log10(Pxx_slice + 1e-12);
    time_axis(k) = (k-1) * frame_time;
    
    if mod(k,10) == 0, fprintf('Progress: %.0f%%\n', (k/num_frames)*100); end
end

fclose(fid);
release(rx);

%% ================== SAVE META & DATA ==================
save(meta_file, 'wide_freq', 'wide_psd', 'psd_matrix', 'selected_fc', ...
                'time_axis', 'recorded_iq_data', 'capture_timestamp');

fprintf('\nDONE! Files saved to %s\n', meta_file);

%% ================== VISUALIZATION ==================

% --- Figure 1: Wideband FM Scan ---
figure('Name','Figure 1 - Wideband FM Scan','NumberTitle','off','Position',[50 600 900 380]);
subplot(1,2,1);
plot(wide_freq/1e6, wide_psd, 'Color', [0.2 0.5 0.9], 'LineWidth', 0.8); hold on;
plot(selected_fc/1e6, max(wide_psd), 'rv', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
for k = 1:length(stations)
    text(stations(k)/1e6, power_vals(k)+1, sprintf('%.1f', stations(k)/1e6), ...
        'FontSize', 7, 'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5]);
end
xlabel('Frequency (MHz)'); ylabel('Power (dB)');
title(sprintf('FM Band Scan — Selected: %.1f MHz', selected_fc/1e6));
grid on; box on;

subplot(1,2,2);
b_bar = bar(1:5, sorted_power(1:5), 0.6, 'FaceColor', 'flat');
b_bar.CData(station_target_idx,:) = [0.9 0.2 0.2]; % highlight selected
for k = 1:5
    set(gca, 'XTick', 1:5, 'XTickLabel', ...
        arrayfun(@(x) sprintf('%.1f', x/1e6), stations(top5_indices(1:5)), 'UniformOutput', false));
    text(k, sorted_power(k)+0.5, sprintf('%.1f dB', sorted_power(k)), ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
end
xlabel('Station (MHz)'); ylabel('Peak Power (dB)');
title(sprintf('Top 5 Stations  |  Index %d Selected (Roll %d)', station_target_idx, ROLL_NUMBER));
grid on; box on;
sgtitle(sprintf('Task 1: Wideband FM Scan  |  %s', capture_timestamp), 'FontWeight', 'bold', 'FontSize', 12);

% --- Figure 2: Captured IQ & Baseline Demod ---
% Quick baseline demod for visualization only
iq_viz = double(recorded_iq_data(:));
iq_viz(~isfinite(iq_viz)) = 0;
iq_viz = iq_viz - mean(iq_viz);
demod_viz = angle(iq_viz(2:end) .* conj(iq_viz(1:end-1)));
demod_viz(~isfinite(demod_viz)) = 0;
audio_viz = decimate(decimate(double(demod_viz), 10), 5);
audio_viz = audio_viz / (max(abs(audio_viz)) + 1e-10);

figure('Name','Figure 2 - Captured IQ & Baseline Demod','NumberTitle','off','Position',[50 150 900 380]);
subplot(1,2,1);
t_iq = (0:min(4999, length(iq_viz)-1)) / fs * 1000;
plot(t_iq, abs(iq_viz(1:length(t_iq))), 'Color', [0.2 0.5 0.9], 'LineWidth', 0.8);
xlabel('Time (ms)'); ylabel('Magnitude');
title(sprintf('Raw IQ Magnitude — %.2f MHz', selected_fc/1e6)); grid on; box on;

subplot(1,2,2);
t_raw = (0:min(4999, length(audio_viz)-1)) / 48000 * 1000;
plot(t_raw, audio_viz(1:length(t_raw)), 'Color', [0.9 0.4 0.2], 'LineWidth', 0.8);
xlabel('Time (ms)'); ylabel('Amplitude');
title('Baseline Demod (No Filter)'); grid on; box on;
sgtitle(sprintf('Task 2: IQ Capture @ %.2f MHz  |  %d samples  |  %.2f s', ...
    selected_fc/1e6, length(iq_viz), length(iq_viz)/fs), 'FontWeight', 'bold', 'FontSize', 12);

% --- Figure 3: Capture Spectrogram ---
figure('Name','Figure 3 - Capture Spectrogram','NumberTitle','off','Position',[100 100 1000 420]);
subplot(1,2,1);
imagesc(time_axis, F_slice/1e6, psd_matrix);
axis xy; colormap('hot'); colorbar;
xlabel('Time (s)'); ylabel('Offset from centre (MHz)');
title(sprintf('Spectrogram @ %.2f MHz', selected_fc/1e6)); box on;

subplot(1,2,2);
mean_psd = mean(psd_matrix, 2);
plot(F_slice/1e6, mean_psd, 'Color', [0.9 0.3 0.1], 'LineWidth', 1.5);
xlabel('Offset from centre (MHz)'); ylabel('Power (dB)');
title('Time-averaged PSD'); grid on; box on;
xline(0, 'k--', 'Centre', 'FontSize', 9);
sgtitle(sprintf('Task 2: Capture Spectrogram  |  %s', capture_timestamp), 'FontWeight', 'bold', 'FontSize', 12);