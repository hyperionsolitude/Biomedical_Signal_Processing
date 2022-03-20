clc;
clear;
close all;
ecg_1=load('ecg_1.mat');
ecg_2=load('ecg_2.mat');
ecg_3=load('ecg_3.mat');
%% Q1.a
Fs=1000;% 1000 Samples per second(1 second)
%% Q1.b
ecg_1_len=length(ecg_1.ecg_lfn);% getting length of first matrix
ecg_2_len=length(ecg_2.ecg_hfn);% getting length of second matrix
ecg_3_len=length(ecg_3.ecg_noisy);% getting length of third matrix
record_time_1=ecg_1_len/Fs;%23.4840 second
record_time_2=ecg_2_len/Fs;%8.5680 second
record_time_3=ecg_3_len/Fs;%10.6960 second
%% Q1.c
t1=0:1:ecg_1_len-1;
t2=0:1:ecg_2_len-1;
t3=0:1:ecg_3_len-1;
% We have 3000 Samples for this interval [2-5)*1000
figure;
plot(t1((2001:5000)),ecg_1.ecg_lfn(2001:5000));
xlabel('Time(ms)');
ylabel('Amplitude');
title('Ecg_1 Signal at Time interval [2,5)');
figure;
plot(t2((2001:5000)),ecg_2.ecg_hfn(2001:5000));
xlabel('Time(ms)');
ylabel('Amplitude');
title('Ecg_2 Signal at Time interval [2,5)');
figure;
plot(t3((2001:5000)),ecg_3.ecg_noisy(2001:5000));
xlabel('Time(ms)');
ylabel('Amplitude');
title('Ecg_3 Signal at Time interval [2,5)');
%% Q1.d
ecg1_freq1=abs(fftshift(fft(ecg_1.ecg_lfn(2001:5000))));ecg1_freq1=ecg1_freq1/length(ecg1_freq1);
ecg2_freq1=abs(fftshift(fft(ecg_2.ecg_hfn(2001:5000))));ecg2_freq1=ecg2_freq1/length(ecg2_freq1);
ecg3_freq1=abs(fftshift(fft(ecg_3.ecg_noisy(2001:5000))));ecg3_freq1=ecg3_freq1/length(ecg3_freq1);
ecg1_phase1=unwrap(angle(fft(ecg_1.ecg_lfn(2001:5000))/length(ecg1_freq1)));
ecg2_phase1=unwrap(angle(fft(ecg_2.ecg_hfn(2001:5000))/length(ecg2_freq1)));
ecg3_phase1=unwrap(angle(fft(ecg_3.ecg_noisy(2001:5000))/length(ecg3_freq1)));

f1=linspace(-Fs/2,Fs/2,length(ecg1_freq1));
f2=linspace(-Fs/2,Fs/2,length(ecg2_freq1));
f3=linspace(-Fs/2,Fs/2,length(ecg3_freq1));
figure;
subplot(211);
plot(f1,ecg1_freq1);
xlabel('Frequency(Hz)');
ylabel('Magnitude');
title('Ecg_1 Signal Frequency Spectrum at Time interval [2,5)');
subplot(212);
plot(f1,ecg1_phase1);
xlabel('Frequency(Hz)');
ylabel('Phase angle');
title('Ecg_1 Signal Phase Spectrum at Time interval [2,5)');
figure;
subplot(211);
plot(f2,ecg2_freq1);
xlabel('Frequency(Hz)');
ylabel('Magnitude');
title('Ecg_2 Signal Frequency Spectrum at Time interval [2,5)');
subplot(212);
plot(f2,ecg2_phase1);
xlabel('Frequency(Hz)');
ylabel('Phase angle');
title('Ecg_2 Signal Phase Spectrum at Time interval [2,5)');
figure;
subplot(211);
plot(f3,ecg3_freq1);
xlabel('Frequency(Hz)');
ylabel('Magnitude');
title('Ecg_3 Signal Frequency Spectrum at Time interval [2,5)');
subplot(212);
plot(f3,ecg3_phase1);
xlabel('Frequency(Hz)');
ylabel('Phase angle');
title('Ecg_3 Signal Phase Spectrum at Time interval [2,5)');
%There are zero angle phases

%% Q1.e
figure;
plot(t1((5001:8000)),ecg_1.ecg_lfn(5001:8000));
xlabel('Time(ms)');
ylabel('Amplitude');
title('Ecg_1 Signal at Time interval [5,8)');
figure;
plot(t2((5001:8000)),ecg_2.ecg_hfn(5001:8000));
xlabel('Time(ms)');
ylabel('Amplitude');
title('Ecg_2 Signal at Time interval [5,8)');
figure;
plot(t3((5001:8000)),ecg_3.ecg_noisy(5001:8000));
xlabel('Time(ms)');
ylabel('Amplitude');
title('Ecg_3 Signal at Time interval [5,8)');
ecg1_freq2=abs(fftshift(fft(ecg_1.ecg_lfn(5001:8000))));ecg1_freq2=ecg1_freq2/length(ecg1_freq2);
ecg2_freq2=abs(fftshift(fft(ecg_2.ecg_hfn(5001:8000))));ecg2_freq2=ecg2_freq2/length(ecg2_freq2);
ecg3_freq2=abs(fftshift(fft(ecg_3.ecg_noisy(5001:8000))));ecg3_freq2=ecg3_freq2/length(ecg3_freq2);

ecg1_phase2=unwrap(angle(fft(ecg_1.ecg_lfn(5001:8000))/length(ecg1_freq2)));
ecg2_phase2=unwrap(angle(fft(ecg_2.ecg_hfn(5001:8000))/length(ecg2_freq2)));
ecg3_phase2=unwrap(angle(fft(ecg_3.ecg_noisy(5001:8000))/length(ecg3_freq2)));

f12=linspace(-Fs/2,Fs/2,length(ecg1_freq2));
f22=linspace(-Fs/2,Fs/2,length(ecg2_freq2));
f32=linspace(-Fs/2,Fs/2,length(ecg3_freq2));
figure;
subplot(211);
plot(f12,ecg1_freq2);
xlabel('Frequency(Hz)');
ylabel('Magnitude');
title('Ecg_1 Signal Frequency Spectrum at Time interval [5,8)');
subplot(212);
plot(f12,ecg1_phase2);
xlabel('Frequency(Hz)');
ylabel('Phase angle');
title('Ecg_1 Signal Phase Spectrum at Time interval [5,8)');
figure;
subplot(211);
plot(f22,ecg2_freq2);
xlabel('Frequency(Hz)');
ylabel('Magnitude');
title('Ecg_2 Signal Frequency Spectrum at Time interval [5,8)');
subplot(212);
plot(f22,ecg2_phase2);
xlabel('Frequency(Hz)');
ylabel('Phase angle');
title('Ecg_2 Signal Phase Spectrum at Time interval [5,8)');
figure;
subplot(211);
plot(f32,ecg3_freq2);
xlabel('Frequency(Hz)');
ylabel('Magnitude');
title('Ecg_3 Signal Frequency Spectrum at Time interval [5,8)');
subplot(212);
plot(f32,ecg3_phase2);
xlabel('Frequency(Hz)');
ylabel('Phase angle');
title('Ecg_3 Signal Phase Spectrum at Time interval [5,8)');
%There are high band noises because hearth beats frequency response should
%be between 0.05 and 150Hz but we have companents further than 150 Hz So It
%should be High frequency noises. But there is also powerline noises(50-60 Hz) so this can
%be counted as bandlimited frequency noises and also there are Baseline Wander which is
%around 0.5Hz this can be counted as Low frequency noise. Thus we have all
%kinds of noise in our data.
%There are nonlinear angle phases in our data.

%% Q1.f
% There is magnitude differences in frequency domain This can be caused by noise(electrodes,powerline etc.) or simply condition of the patient has been changed 
%% Q2.a
rect_window=rectwin(256);
figure;
spectrogram(ecg_1.ecg_lfn,rect_window,128);
title('Spectrogram of ecg_1 with rectangular window %50 overlap');
figure;
spectrogram(ecg_2.ecg_hfn,rect_window,128);
title('Spectrogram of ecg_2 with rectangular window %50 overlap');
figure;
spectrogram(ecg_3.ecg_noisy,rect_window,128);
title('Spectrogram of ecg_3 with rectangular window %50 overlap');
%% Q2.b
hamm_window=hamming(256);
figure;
spectrogram(ecg_1.ecg_lfn,hamm_window,128);
title('Spectrogram of ecg_1 with hamming window %50 overlap');
figure;
spectrogram(ecg_2.ecg_hfn,hamm_window,128);
title('Spectrogram of ecg_2 with hamming window %50 overlap');
figure;
spectrogram(ecg_3.ecg_noisy,hamm_window,128);
title('Spectrogram of ecg_3 with hamming window %50 overlap');