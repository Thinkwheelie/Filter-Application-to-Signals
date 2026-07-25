clc; clear; close all;

%% ================================================================
%% PROJECT: DSP for SDR FM Receivers
%% ================================================================
ROLL_NUMBER   = 244;
station_index = mod(ROLL_NUMBER, 5) + 1;
fprintf('Roll Number: %d  →  Station Index: %d\n', ROLL_NUMBER, station_index);

%% ================================================================
%% TASK 1 - BAND SCANNING (commented out - use with live RTL-SDR)
%% PDF Deliverable 1.1: PSD plot with top 5 peaks labeled
%% PDF Deliverable 1.2: Table of fc and estimated bandwidth
%% PDF Deliverable 1.3: Roll Number → Station Index calculation
%%   Index = mod(244, 5) + 1 = 4 + 1 = 5  → Station Index = 5
%% ================================================================
% fm_start = 88e6; fm_end = 108e6; scan_bw = 2.4e6; step = 1e6;
% freqs = fm_start:step:fm_end;
% psd_peaks = zeros(size(freqs));
% for i = 1:length(freqs)
%     rx = comm.SDRRTLReceiver('CenterFrequency',freqs(i),'SampleRate',scan_bw,'SamplesPerFrame',65536);
%     data = rx(); release(rx);
%     psd_peaks(i) = mean(abs(fft(data)).^2);
% end
% [sorted_psd,sort_idx] = sort(psd_peaks,'descend');
% top5_freqs = freqs(sort_idx(1:5));
% top5_psd   = sorted_psd(1:5);
% figure;
% plot(freqs/1e6, 10*log10(psd_peaks),'b'); hold on;
% plot(top5_freqs/1e6, 10*log10(top5_psd),'rv','MarkerSize',10,'MarkerFaceColor','r');
% for k=1:5
%     text(top5_freqs(k)/1e6, 10*log10(top5_psd(k))+1, sprintf('%.1f MHz',top5_freqs(k)/1e6),'Color','r','FontSize',9);
% end
% xlabel('Frequency (MHz)'); ylabel('PSD (dB)'); title('FM Band Scan - Top 5 Stations'); grid on;
% selected_fc = top5_freqs(station_index);
% fprintf('Selected: %.2f MHz\n', selected_fc/1e6);
%
% --- PDF Table 1.2: Top 5 Discovered Stations ---
% Rank | fc (MHz) | Est. BW (kHz) | Notes
%  1   |   xx.x   |     200       | Strongest station
%  2   |   xx.x   |     200       |
%  3   |   xx.x   |     200       |
%  4   |   xx.x   |     200       |
%  5   |   xx.x   |     200       | ← Assigned (Roll 244, Index 5)

%% ================================================================
%% TASK 2 - LOAD IQ + TIMESTAMP
%% ================================================================
fs        = 2.4e6;
audio_fs  = 48000;

meta_file = 'record_meta.mat';
if exist(meta_file, 'file')
    load(meta_file);
    if ~exist('capture_timestamp', 'var'), capture_timestamp = 'Unknown'; end
else
    error('Metadata file %s not found.', meta_file);
end
fprintf('Data timestamp: %s\n', capture_timestamp);

iq = double(recorded_iq_data(:));   % cast to double for filtfilt compatibility
iq(~isfinite(iq)) = 0;              % remove any NaN/Inf from SDR dropouts
iq = iq - mean(iq);
fprintf('Loaded %d samples (%.2f s)\n', length(iq), length(iq)/fs);

% Baseline demod - no filtering at all
demod_raw = angle(iq(2:end) .* conj(iq(1:end-1)));
demod_raw(~isfinite(demod_raw)) = 0;   % sanitize before filtfilt inside decimate
audio_raw = decimate(decimate(double(demod_raw), 10), 5);
audio_raw = audio_raw / (max(abs(audio_raw)) + 1e-10);
audiowrite('audio_raw.wav', audio_raw, audio_fs);
fprintf('Saved: audio_raw.wav\n');

% Figure 1
figure('Name','Figure 1 - Task 2: Data Acquisition','NumberTitle','off','Position',[50 600 900 380]);
subplot(1,2,1);
t_iq = (0:4999)/fs*1000;
plot(t_iq, abs(iq(1:5000)),'Color',[0.2 0.5 0.9],'LineWidth',0.8);
xlabel('Time (ms)'); ylabel('Magnitude');
title('Raw IQ Signal - Magnitude envelope'); grid on; box on;
subplot(1,2,2);
t_raw = (0:min(4999,length(audio_raw)-1))/audio_fs*1000;
plot(t_raw, audio_raw(1:length(t_raw)),'Color',[0.9 0.4 0.2],'LineWidth',0.8);
xlabel('Time (ms)'); ylabel('Amplitude');
title('Baseline FM Demod - No filtering applied'); grid on; box on;
sgtitle(sprintf('Task 2: Raw IQ Data (timestamp: %s) & Baseline FM Demodulation', capture_timestamp),'FontWeight','bold','FontSize',12);

%% ================================================================
%% TASK 3 - FIR vs IIR FILTER COMPARISON
%% PDF Task 3: FIR (Hamming window) vs IIR (Chebyshev Type II)
%% Deliverables: magnitude+phase overlay, group delay, 3-panel spectrogram, SNR bar
%% ================================================================
Wn      = 100e3 / (fs/2);         % normalised cutoff = 100 kHz
Rs      = 60;                      % Chebyshev Type II stopband attenuation (dB)
b_fir   = fir1(101, Wn, hamming(102));           % FIR - Hamming window, order 101
[b_iir, a_iir] = cheby2(6, Rs, Wn);             % IIR - Chebyshev Type II, order 6

iq_fir = filter(b_fir, 1,      iq);
iq_iir = filter(b_iir, a_iir,  iq);

[H_fir, w] = freqz(b_fir, 1,      2048, fs);
[H_iir, ~] = freqz(b_iir, a_iir,  2048, fs);
gd_fir     = grpdelay(b_fir, 1,      2048);
gd_iir     = grpdelay(b_iir, a_iir,  2048);

% RAW IQ → FIR filter → demod → audio  (independent pipeline)
d_fir = angle(iq_fir(2:end) .* conj(iq_fir(1:end-1)));
d_iir = angle(iq_iir(2:end) .* conj(iq_iir(1:end-1)));
d_fir(~isfinite(d_fir)) = 0;
d_iir(~isfinite(d_iir)) = 0;
a_fir_q = decimate(decimate(double(d_fir), 10), 5);
a_iir_q = decimate(decimate(double(d_iir), 10), 5);
a_fir_q = a_fir_q / (max(abs(a_fir_q)) + 1e-10);
a_iir_q = a_iir_q / (max(abs(a_iir_q)) + 1e-10);
snr_fir = snr(a_fir_q);
snr_iir = snr(a_iir_q);
fprintf('FIR SNR: %.1f dB  |  IIR (Cheby2) SNR: %.1f dB\n', snr_fir, snr_iir);
audiowrite('audio_FIR.wav', a_fir_q, audio_fs);
audiowrite('audio_IIR.wav', a_iir_q, audio_fs);
fprintf('Saved: audio_FIR.wav, audio_IIR.wav\n');

% --- Figure 2: Filter Magnitude, Phase, Group Delay & SNR ---
% PDF Deliverable 3.1: Overlay magnitude+phase, highlight phase linearity difference
figure('Name','Figure 2 - Task 3: FIR vs IIR Filter Analysis','NumberTitle','off','Position',[50 150 1000 600]);

subplot(2,2,1);
plot(w/1e3, 20*log10(abs(H_fir)+1e-10),'b','LineWidth',1.8); hold on;
plot(w/1e3, 20*log10(abs(H_iir)+1e-10),'r--','LineWidth',1.8);
xlabel('Frequency (kHz)'); ylabel('Magnitude (dB)');
title('(a) Magnitude Response');
legend('FIR - Hamming window, order 101','IIR - Chebyshev Type II, order 6','Location','southwest');
xlim([0 300]); ylim([-80 5]); grid on; box on;
xline(100,'k--','100 kHz cutoff','LabelVerticalAlignment','bottom','FontSize',8);

subplot(2,2,2);
plot(w/1e3, unwrap(angle(H_fir)),'b','LineWidth',1.8); hold on;
plot(w/1e3, unwrap(angle(H_iir)),'r--','LineWidth',1.8);
xlabel('Frequency (kHz)'); ylabel('Phase (radians)');
title('(b) Phase Response — FIR linear, Cheby2 nonlinear');
legend('FIR (linear phase)','Cheby2 (nonlinear phase)','Location','southwest');
xlim([0 300]); grid on; box on;

subplot(2,2,3);
% PDF Deliverable 3.1: Highlight group delay difference
plot(w/1e3, gd_fir,'b','LineWidth',1.8); hold on;
plot(w/1e3, gd_iir,'r--','LineWidth',1.8);
xlabel('Frequency (kHz)'); ylabel('Group delay (samples)');
title('(c) Group Delay — FIR constant, Cheby2 varies near cutoff');
legend('FIR (flat/constant)','Cheby2 (varies)','Location','northeast');
xlim([0 300]); grid on; box on;

subplot(2,2,4);
b_bar = bar([snr_fir, snr_iir], 0.5, 'FaceColor','flat');
b_bar.CData(1,:) = [0.2 0.4 0.8]; b_bar.CData(2,:) = [0.8 0.2 0.2];
set(gca,'XTickLabel',{'FIR (Hamming)','IIR (Cheby2)'},'FontSize',10);
ylabel('SNR (dB)'); title('(d) Resulting Audio SNR Comparison');
text(1, snr_fir+0.3, sprintf('%.1f dB',snr_fir),'HorizontalAlignment','center','FontWeight','bold');
text(2, snr_iir+0.3, sprintf('%.1f dB',snr_iir),'HorizontalAlignment','center','FontWeight','bold');
grid on; box on;
sgtitle('Task 3: FIR (Hamming) vs IIR (Chebyshev Type II) — cutoff = 100 kHz', ...
    'FontWeight','bold','FontSize',12);

% --- Figure 3: 3-panel channel isolation spectrograms ---
% PDF Deliverable 3.2: (a) unfiltered, (b) FIR isolated, (c) IIR isolated
n_spec = min(50000, length(iq));
figure('Name','Figure 3 - Task 3: Channel Isolation Spectrograms','NumberTitle','off','Position',[100 100 900 750]);

subplot(3,1,1);
spectrogram(real(iq(1:n_spec)), 256, 200, 256, fs, 'yaxis');
title('(a) Centered but Unfiltered — adjacent channels visible on both sides');
cb = colorbar; cb.Label.String = 'Power (dB)';

subplot(3,1,2);
spectrogram(real(iq_fir(1:n_spec)), 256, 200, 256, fs, 'yaxis');
title('(b) Isolated via FIR Low-Pass Filter (Hamming) — linear phase preserved');
cb = colorbar; cb.Label.String = 'Power (dB)';

subplot(3,1,3);
spectrogram(real(iq_iir(1:n_spec)), 256, 200, 256, fs, 'yaxis');
title('(c) Isolated via IIR Low-Pass Filter (Chebyshev Type II) — phase nonlinearity introduced');
cb = colorbar; cb.Label.String = 'Power (dB)';

sgtitle('Task 3: Baseband Channel Isolation — FIR vs Chebyshev Type II IIR', ...
    'FontWeight','bold','FontSize',12);

%% ================================================================
%% MAIN PIPELINE - best quality audio via PLL
%% ================================================================
b_main  = fir1(1023, 80e3/(fs/2), kaiser(1024,10));
iq_filt = filter(b_main, 1, iq);

sig_pwr      = mean(abs(iq_filt).^2);
noise_pwr    = mean(abs(iq - iq_filt).^2);
snr_predemod = 10*log10(sig_pwr / (noise_pwr + 1e-12));
fprintf('Pre-demod SNR: %.1f dB\n', snr_predemod);

Kp = 0.02; Ki = 0.0005;
phase = 0; integrator = 0;
N = length(iq_filt);
demod_pll = zeros(N,1);
for n = 1:N
    err        = angle(iq_filt(n) * exp(-1j*phase));
    if ~isfinite(err), err = 0; end   % PLL safety gate
    integrator = integrator + Ki * err;
    freq_inc   = Kp * err + integrator;
    phase      = phase + freq_inc;
    demod_pll(n) = freq_inc * 0.4;
end
demod_pll(~isfinite(demod_pll)) = 0;  % sanitize before decimation

audio_240k = decimate(demod_pll, 10);
audio_48k  = decimate(audio_240k, 5);

tau      = 75e-6;
d        = exp(-1 / (audio_fs * tau));
audio_de = filter((1-d), [1,-d], audio_48k);

[b_hp, a_hp] = butter(3,    40/(audio_fs/2), 'high');
[b_lp, a_lp] = butter(6, 13000/(audio_fs/2), 'low');
audio_bp = filter(b_lp, a_lp, filter(b_hp, a_hp, audio_de));

%% ================================================================
%% TASK 4 - WIENER FILTER   H(f) = Pss / (Pss + Pnn)
%% Independent pipeline: RAW IQ → baseline demod → Wiener directly
%% (No FIR/IIR filtering applied — Wiener operates on unfiltered audio)
%% ================================================================
% audio_raw = baseline demod of raw IQ (already computed in Task 2)
% Wiener is applied directly to this — separate from FIR/IIR paths
wiener_input = audio_raw;   % RAW IQ → demod → Wiener (no bandpass pre-filter)
N_audio   = length(wiener_input);
frame_len = 2048;
hop       = 512;
win_w     = sqrt(hann(frame_len,'periodic'));

% Find 20 quietest frames for noise PSD estimate
n_total = floor((N_audio - frame_len) / hop);
fpwr    = zeros(n_total, 1);
for k = 1:n_total
    seg = wiener_input((k-1)*hop+1 : (k-1)*hop+frame_len);
    fpwr(k) = mean(seg.^2);
end
[~, quiet_idx] = mink(fpwr, 20);
noise_acc = zeros(frame_len, 1);
for k = quiet_idx.'
    seg = wiener_input((k-1)*hop+1:(k-1)*hop+frame_len) .* win_w;
    noise_acc = noise_acc + abs(fft(seg)).^2;
end
Pnn = noise_acc / numel(quiet_idx);

% Apply Wiener: H = Pss / (Pss + Pnn)
audio_wiener = zeros(N_audio, 1);
norm_w       = zeros(N_audio, 1);
for k = 1:hop:N_audio-frame_len
    seg = wiener_input(k:k+frame_len-1) .* win_w;
    F   = fft(seg);
    Pss = abs(F).^2;
    H   = Pss ./ (Pss + Pnn + 1e-12);
    out = real(ifft(H .* F)) .* win_w;
    audio_wiener(k:k+frame_len-1) = audio_wiener(k:k+frame_len-1) + out;
    norm_w(k:k+frame_len-1)       = norm_w(k:k+frame_len-1) + win_w.^2;
end
valid = norm_w > 0.1;
audio_wiener(valid) = audio_wiener(valid) ./ norm_w(valid);
audio_wiener(~isfinite(audio_wiener)) = 0;
audio_wiener = audio_wiener / (max(abs(audio_wiener)) + 1e-10);
audiowrite('audio_wiener.wav', audio_wiener, audio_fs);
fprintf('Saved: audio_wiener.wav\n');

snr_raw    = snr(wiener_input);   % raw demod SNR (Wiener input)
snr_before = snr_raw;
snr_after  = snr(audio_wiener);
fprintf('\n--- SNR Summary ---\n');
fprintf('Pre-demod SNR (IQ):    %.1f dB\n', snr_predemod);
fprintf('FIR audio SNR:         %.1f dB\n', snr_fir);
fprintf('IIR audio SNR:         %.1f dB\n', snr_iir);
fprintf('Raw demod SNR:         %.1f dB\n', snr_raw);
fprintf('Wiener output SNR:     %.1f dB\n', snr_after);
fprintf('Wiener improvement:    %.1f dB\n', snr_after - snr_before);

% --- Figure 4: Wiener PSD comparison + 3-way SNR bar ---
figure('Name','Figure 4 - Task 4: Wiener Filter Analysis','NumberTitle','off','Position',[100 100 1100 480]);
[pxx_before, f_psd] = pwelch(wiener_input, [], [], [], audio_fs);
[pxx_after,  ~    ] = pwelch(audio_wiener, [], [], [], audio_fs);

subplot(1,3,1);
plot(f_psd/1e3, 10*log10(pxx_before),'b','LineWidth',1.5); hold on;
plot(f_psd/1e3, 10*log10(pxx_after), 'r--','LineWidth',1.5);
xlabel('Frequency (kHz)'); ylabel('PSD (dB/Hz)');
title('(a) PSD: Raw Demod vs After Wiener');
legend('Raw demod (input)','After Wiener','Location','northeast');
grid on; box on;

subplot(1,3,2);
freq_bins = (0:frame_len-1) * audio_fs / frame_len / 1e3;
plot(freq_bins(1:frame_len/2), 10*log10(Pnn(1:frame_len/2)+1e-12),'k','LineWidth',1.5);
xlabel('Frequency (kHz)'); ylabel('Power (dB)');
title('(b) Estimated Noise PSD P_{nn}(f)');
xlim([0 audio_fs/2/1e3]); grid on; box on;

subplot(1,3,3);
% 3-way SNR bar: FIR | IIR | Wiener — all from same raw IQ source
all_snr = [snr_fir, snr_iir, snr_after];
b4 = bar(all_snr, 0.55, 'FaceColor','flat');
b4.CData(1,:) = [0.2 0.4 0.8];   % FIR - blue
b4.CData(2,:) = [0.8 0.2 0.2];   % IIR - red
b4.CData(3,:) = [0.1 0.65 0.3];  % Wiener - green
set(gca,'XTickLabel',{'FIR','IIR (Cheby2)','Wiener'},'FontSize',10);
ylabel('SNR (dB)');
title('(c) SNR Comparison: FIR vs IIR vs Wiener');
for k = 1:3
    text(k, all_snr(k)+0.3, sprintf('%.1f dB', all_snr(k)), ...
        'HorizontalAlignment','center','FontWeight','bold','FontSize',9);
end
grid on; box on;

sgtitle('Task 4: Wiener Filter — H(f) = P_{ss}(f) / [P_{ss}(f) + P_{nn}(f)]  |  Noise Floor Suppression', ...
    'FontWeight','bold','FontSize',12);

%% ================================================================
%% TASK 5 - STFT ANALYSIS
%% 3 pipelines: RAW demod | FIR demod | Wiener on raw demod
%% ================================================================
figure('Name','Figure 5 - Task 5: Time-Frequency Analysis','NumberTitle','off','Position',[100 100 1100 750]);
N_show  = min(length(audio_raw), length(audio_wiener));
t_audio = (0:N_show-1) / audio_fs;
n_show  = min(length(audio_raw), N_show);

subplot(3,2,1);
spectrogram(audio_raw(1:n_show),512,400,512,audio_fs,'yaxis');
title('(a) STFT — Raw demod (no filter, Wiener input)'); cb=colorbar; cb.Label.String='dB';

subplot(3,2,3);
spectrogram(a_fir_q(1:min(n_show,length(a_fir_q))),512,400,512,audio_fs,'yaxis');
title('(b) STFT — After FIR filter (Hamming, 100 kHz)'); cb=colorbar; cb.Label.String='dB';

subplot(3,2,5);
spectrogram(audio_wiener(1:n_show),512,400,512,audio_fs,'yaxis');
title('(c) STFT — After Wiener denoising (applied to raw demod)'); cb=colorbar; cb.Label.String='dB';

subplot(3,2,2);
plot(t_audio, audio_raw(1:N_show),'Color',[0.5 0.5 0.5],'LineWidth',0.6);
xlabel('Time (s)'); ylabel('Amplitude'); title('Waveform — Raw demod (unfiltered)'); grid on; box on;

subplot(3,2,4);
n_fir = min(N_show, length(a_fir_q));
plot((0:n_fir-1)/audio_fs, a_fir_q(1:n_fir),'b','LineWidth',0.6);
xlabel('Time (s)'); ylabel('Amplitude'); title('Waveform — After FIR filter'); grid on; box on;

subplot(3,2,6);
plot(t_audio, audio_wiener(1:N_show),'r','LineWidth',0.6);
xlabel('Time (s)'); ylabel('Amplitude'); title('Waveform — After Wiener filter'); grid on; box on;

sgtitle('Task 5: STFT Time-Frequency Analysis — Raw vs FIR vs Wiener (all from same RAW IQ)','FontWeight','bold','FontSize',12);

%% ================================================================
%% SEPARATE FIGURES FOR FIR / IIR / WIENER (for presentation)
%% ================================================================

% --- Figure 7: FIR audio alone ---
figure('Name','Figure 7 - FIR Filter Result','NumberTitle','off','Position',[100 100 1000 600]);

subplot(2,2,1);
plot(w/1e3, 20*log10(abs(H_fir)+1e-10),'b','LineWidth',2);
xlabel('Frequency (kHz)'); ylabel('Magnitude (dB)');
title('FIR Magnitude Response (Hamming, order 101)');
xline(100,'r--','Cutoff','FontSize',9);
xlim([0 300]); ylim([-80 5]); grid on; box on;

subplot(2,2,2);
plot(w/1e3, unwrap(angle(H_fir)),'b','LineWidth',2);
xlabel('Frequency (kHz)'); ylabel('Phase (radians)');
title('FIR Phase Response - perfectly linear');
xlim([0 300]); grid on; box on;

subplot(2,2,3);
n_wave = min(5000, length(a_fir_q));
t_w = (0:n_wave-1)/audio_fs*1000;
plot(t_w, a_fir_q(1:n_wave),'b','LineWidth',0.8);
xlabel('Time (ms)'); ylabel('Amplitude');
title('FIR demodulated audio waveform'); grid on; box on;

subplot(2,2,4);
[pxx_fir, f_fir] = pwelch(a_fir_q,[],[],[],audio_fs);
plot(f_fir/1e3, 10*log10(pxx_fir),'b','LineWidth',1.5);
xlabel('Frequency (kHz)'); ylabel('PSD (dB/Hz)');
title(sprintf('FIR audio spectrum  |  SNR = %.1f dB', snr_fir));
grid on; box on;
sgtitle('FIR Low-Pass Filter - Audio Result','FontWeight','bold','FontSize',13);

% --- Figure 8: IIR audio alone ---
figure('Name','Figure 8 - IIR Filter Result','NumberTitle','off','Position',[100 100 1000 600]);

subplot(2,2,1);
plot(w/1e3, 20*log10(abs(H_iir)+1e-10),'r','LineWidth',2);
xlabel('Frequency (kHz)'); ylabel('Magnitude (dB)');
title('IIR Magnitude Response (Chebyshev Type II, order 6, Rs=60 dB)');
xline(100,'b--','Cutoff','FontSize',9);
xlim([0 300]); ylim([-80 5]); grid on; box on;

subplot(2,2,2);
plot(w/1e3, unwrap(angle(H_iir)),'r','LineWidth',2);
xlabel('Frequency (kHz)'); ylabel('Phase (radians)');
title('IIR Phase Response — Cheby2 nonlinear near cutoff');
xlim([0 300]); grid on; box on;

subplot(2,2,3);
n_wave = min(5000, length(a_iir_q));
plot(t_w, a_iir_q(1:n_wave),'r','LineWidth',0.8);
xlabel('Time (ms)'); ylabel('Amplitude');
title('IIR demodulated audio waveform'); grid on; box on;

subplot(2,2,4);
[pxx_iir, f_iir] = pwelch(a_iir_q,[],[],[],audio_fs);
plot(f_iir/1e3, 10*log10(pxx_iir),'r','LineWidth',1.5);
xlabel('Frequency (kHz)'); ylabel('PSD (dB/Hz)');
title(sprintf('IIR audio spectrum  |  SNR = %.1f dB', snr_iir));
grid on; box on;
sgtitle('IIR Low-Pass Filter (Chebyshev Type II) - Audio Result','FontWeight','bold','FontSize',13);

% --- Figure 9: Wiener filter alone ---
figure('Name','Figure 9 - Wiener Filter Result','NumberTitle','off','Position',[100 100 1000 600]);

subplot(2,2,1);
plot(f_psd/1e3, 10*log10(pxx_before),'b','LineWidth',1.5); hold on;
plot(f_psd/1e3, 10*log10(pxx_after),'r','LineWidth',1.5);
xlabel('Frequency (kHz)'); ylabel('PSD (dB/Hz)');
title('PSD: Before vs After Wiener');
legend('Input (audio\_bp)','After Wiener','Location','northeast');
grid on; box on;

subplot(2,2,2);
freq_bins = (0:frame_len-1) * audio_fs / frame_len / 1e3;
plot(freq_bins(1:frame_len/2), 10*log10(Pnn(1:frame_len/2)+1e-12),'k','LineWidth',1.5);
xlabel('Frequency (kHz)'); ylabel('Power (dB)');
title('Estimated noise PSD (Pnn) from quietest frames');
xlim([0 audio_fs/2/1e3]); grid on; box on;

subplot(2,2,3);
n_wave = min(5000, N_audio);
t_w2 = (0:n_wave-1)/audio_fs*1000;
plot(t_w2, audio_bp(1:n_wave),'b','LineWidth',0.8); hold on;
plot(t_w2, audio_wiener(1:n_wave),'r','LineWidth',0.8);
xlabel('Time (ms)'); ylabel('Amplitude');
title('Waveform before (blue) vs after (red) Wiener');
legend('Before','After'); grid on; box on;

subplot(2,2,4);
b3 = bar([snr_before, snr_after], 0.6, 'FaceColor','flat');
b3.CData(1,:) = [0.2 0.4 0.8]; b3.CData(2,:) = [0.1 0.65 0.3];
set(gca,'XTickLabel',{'Before Wiener','After Wiener'},'FontSize',11);
ylabel('SNR (dB)'); title('SNR improvement from Wiener filter');
text(1,snr_before+0.3,sprintf('%.1f dB',snr_before),'HorizontalAlignment','center','FontWeight','bold');
text(2,snr_after+0.3, sprintf('%.1f dB',snr_after), 'HorizontalAlignment','center','FontWeight','bold');
grid on; box on;
sgtitle('Wiener Filter - H(f) = Pss/(Pss+Pnn)  |  Noise Suppression Result','FontWeight','bold','FontSize',13);

%% ================================================================
%% FIGURE 6 - REAL-TIME SPECTROGRAM
%% ================================================================
audio_final = tanh(2 * audio_wiener) / tanh(2);
audio_final = audio_final / (max(abs(audio_final)) + 1e-10);
audiowrite('final_audio.wav', audio_final, audio_fs);

fft_len=1024; hop_rt=256;
win_rt=sqrt(hann(fft_len,'periodic'));
num_cols=300;
f_kHz=linspace(0,audio_fs/2000,fft_len/2+1);
wfall=-80*ones(fft_len/2+1,num_cols);

fig=figure('Name','Figure 6 - Real-Time Spectrogram','NumberTitle','off','Color','k','Position',[50 50 1000 450]);
ax=axes(fig);
img_h=imagesc(ax,wfall); axis(ax,'xy');
colormap(ax,hot); clim(ax,[-70 -10]);
ylabel(ax,'Frequency (kHz)','Color','w');
xlabel(ax,'Time →','Color','w');
title(ax,sprintf('Real-Time FM Audio Spectrogram  |  SNR: %.1f dB  →  %.1f dB  →  %.1f dB (after Wiener)', snr_predemod, snr_before, snr_after),'Color','w','FontSize',10);
set(ax,'Color','k','XColor','w','YColor','w','XTick',[],'TickDir','out');
tick_idx=round(linspace(1,fft_len/2+1,7));
set(ax,'YTick',tick_idx,'YTickLabel',arrayfun(@(x)sprintf('%.0f',x),f_kHz(tick_idx),'UniformOutput',false));
cb=colorbar(ax); cb.Color='w'; cb.Label.String='Power (dB)'; cb.Label.Color='w';

fprintf('\nPlaying final audio...\n');
player=audioplayer(audio_final,audio_fs);
play(player); ptr=1;
while isplaying(player) && ptr+fft_len<=length(audio_final)
    seg=audio_final(ptr:ptr+fft_len-1).*win_rt;
    spec_db=20*log10(abs(fft(seg,fft_len))/fft_len+1e-10);
    wfall=[wfall(:,2:end),spec_db(1:fft_len/2+1)];
    set(img_h,'CData',wfall);
    drawnow limitrate;
    ptr=ptr+hop_rt;
end

fprintf('\n✓ All files saved:\n');
fprintf('  selected_station.mat - raw IQ with timestamp\n');
fprintf('  audio_raw.wav        - baseline, no filtering\n');
fprintf('  audio_FIR.wav        - FIR filtered\n');
fprintf('  audio_IIR.wav        - IIR filtered\n');
fprintf('  audio_wiener.wav     - Wiener denoised\n');
fprintf('  final_audio.wav      - final soft-clipped output\n');