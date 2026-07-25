
clear; clc; close all;

filename = '.wav'; 
[iq_data, fs_sdr] = audioread(filename);
complex_signal = iq_data(:,1) + 1i*iq_data(:,2); 

bw_fm = 200e3; 
lpf_iq = designfilt('lowpassfir', 'PassbandFrequency', bw_fm/2, ...
                    'StopbandFrequency', (bw_fm/2) + 20e3, 'SampleRate', fs_sdr);
complex_signal_filtered = filter(lpf_iq, complex_signal);

demod_raw = diff(unwrap(angle(complex_signal_filtered)));

audio_fs = 48e3;
audio_baseband = resample(demod_raw, audio_fs, fs_sdr);

tau = 75e-6; 
alpha = 1 / (2 * pi * audio_fs * tau);
audio_clear = filter(alpha, [1, -(1-alpha)], audio_baseband);

audio_clear = audio_clear / max(abs(audio_clear));
soundsc(audio_clear, audio_fs);