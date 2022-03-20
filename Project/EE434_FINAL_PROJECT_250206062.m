%% UTKU ACAR 250206062
%  Final Project
close all
clear 
clc
fs=1000;
fc=150;
ord=9;
Wp=2*fc/fs;
ecg = load('ecg_data_org.mat');
ecgorg=cell2mat(struct2cell(ecg));
ecgorg=ecgorg.ecg_data_org;
ecgorg=ecgorg(1:2001);
ecg_data=ecgorg;
t=1:1:length(ecg_data);
t=t/fs;
f=linspace(-fs/2,fs/2,length(ecg_data));
ecg_data=awgn(ecg_data,20,'measured');
%ecg_data=ecgorg;

%% Time domain representation of ECG Signal

figure;
plot(t,ecgorg);
title('Time domain Representatiom of ECG Signal');
xlabel('Time(sec)');
ylabel('Amplitude');

%% Magnitude and Phase Response of ECG signal

fft_ecg=fft(ecg_data)/length(ecg_data);
figure;
subplot(211);
plot(f,abs(fftshift(fft_ecg)));
title('FFT of Signal');
xlabel('Frequency(Hz)');
ylabel('Magnitude');
subplot(212);
plot(f,unwrap(angle(fft_ecg)));
title('Phase of Signal');
xlabel('Frequency(Hz)');
ylabel('Degree');

%% Filtering

gaussianmask=[1,2,1]/4;
Gaussian_Filtered_ECG=filtfilt(gaussianmask,1,ecg_data);
fft_ecg_gauss_filtered=fft(Gaussian_Filtered_ECG)/length(Gaussian_Filtered_ECG);

MA_kernel_size=4;
MA_kernel=ones(1,MA_kernel_size)/MA_kernel_size;
MA_Filtered_ECG=filtfilt(MA_kernel,1,ecg_data);
fft_ecg_MA_filtered=fft(MA_Filtered_ECG)/length(MA_Filtered_ECG);

[b1,a1]=butter(ord,Wp,"low");
Butter_Filtered_ECG=filtfilt(b1,a1,ecg_data);
fft_ecg_Butter_filtered=fft(Butter_Filtered_ECG)/length(Butter_Filtered_ECG);

[b2,a2]=fir1(ord,Wp,"low",rectwin(ord+1));
Fir_Filtered_ECG=filtfilt(b2,a2,ecg_data);
fft_ecg_Fir_filtered=fft(Fir_Filtered_ECG)/length(Fir_Filtered_ECG);

%% Filtering Results

figure;
subplot(211);
plot(f,abs(fftshift(fft_ecg)));
title('FFT of Signal');
xlabel('Frequency(Hz)');
ylabel('Magnitude');
hold on
plot(f,abs(fftshift(fft_ecg_gauss_filtered)));
plot(f,abs(fftshift(fft_ecg_MA_filtered)));
plot(f,abs(fftshift(fft_ecg_Butter_filtered)));
plot(f,abs(fftshift(fft_ecg_Fir_filtered)));
legend("Original","Gauss Filtered","MA Filtered","Butter Filtered","FIR Filtered");
hold off
subplot(212);
plot(f,unwrap(angle(fft_ecg)));
title('Phase of Signal');
xlabel('Frequency(Hz)');
ylabel('Degree');
hold on
plot(f,unwrap(angle(fft_ecg_gauss_filtered)));
plot(f,unwrap(angle(fft_ecg_MA_filtered)));
plot(f,unwrap(angle(fft_ecg_Butter_filtered)));
plot(f,unwrap(angle(fft_ecg_Fir_filtered)));
legend("Original","Gauss Filtered","MA Filtered","Butter Filtered","FIR Filtered");
hold off
figure;
plot(t,ecgorg);
title('Time domain Representation of Signal');
xlabel('Time(sec)');
ylabel('Amplitude');
hold on
plot(t,ecg_data);
plot(t,Gaussian_Filtered_ECG);
plot(t,MA_Filtered_ECG);
plot(t,Butter_Filtered_ECG);
plot(t,Fir_Filtered_ECG);
legend("Original","Noisy","Gauss Filtered","MA Filtered","Butter Filtered","FIR Filtered");
hold off

%% Spectrogram(Short Time Fourier Transform)

% Rectangular Window

figure;
stft(ecgorg,fs,'Window',rectwin(256),'OverlapLength',128);
title('Spectrogram of Original ECG with Rectangular Window %50 Overlap');
figure;
stft(ecg_data,fs,'Window',rectwin(256),'OverlapLength',128);
title('Spectrogram of Noisy ECG with Rectangular Window %50 Overlap');
figure;
stft(Gaussian_Filtered_ECG,fs,'Window',rectwin(256),'OverlapLength',128);
title('Spectrogram of Gaussian Filtered with Rectangular Window %50 Overlap');
figure;
stft(MA_Filtered_ECG,fs,'Window',rectwin(256),'OverlapLength',128);
title('Spectrogram of MA Filtered with Rectangular Window %50 Overlap');
figure;
stft(Butter_Filtered_ECG,fs,'Window',rectwin(256),'OverlapLength',128);
title('Spectrogram of Butter Filtered with Rectangular Window %50 Overlap');
figure;
stft(ecgorg,fs,'Window',rectwin(256),'OverlapLength',128);
title('Spectrogram of FIR filtered with Rectangular Window %50 Overlap');

% Hamming Window

figure;
stft(ecgorg,fs,'Window',hamming(256),'OverlapLength',128);
title('Spectrogram of Original ECG with Hamming Window %50 Overlap');
figure;
stft(ecg_data,fs,'Window',hamming(256),'OverlapLength',128);
title('Spectrogram of Noisy ECG with Hamming Window %50 Overlap');
figure;
stft(Gaussian_Filtered_ECG,fs,'Window',hamming(256),'OverlapLength',128);
title('Spectrogram of Gaussian Filtered with Hamming Window %50 Overlap');
figure;
stft(MA_Filtered_ECG,fs,'Window',hamming(256),'OverlapLength',128);
title('Spectrogram of MA Filtered with Hamming Window %50 Overlap');
figure;
stft(Butter_Filtered_ECG,fs,'Window',hamming(256),'OverlapLength',128);
title('Spectrogram of Butter Filtered with Hamming Window %50 Overlap');
figure;
stft(ecgorg,fs,'Window',hamming(256),'OverlapLength',128);
title('Spectrogram of FIR filtered with Hamming Window %50 Overlap');

%% Numeric Error Calculations(Mean Square Error)

RMSE_Noisy=sqrt(mean((ecgorg-ecg_data).^2));
RMSE_Gaussian=sqrt(mean((ecgorg-Gaussian_Filtered_ECG).^2));
RMSE_MA=sqrt(mean((ecgorg-MA_Filtered_ECG).^2));
RMSE_Butter=sqrt(mean((ecgorg-Butter_Filtered_ECG).^2));
RMSE_FIR=sqrt(mean((ecgorg-Fir_Filtered_ECG).^2));

%% Periodogram 

figure;
subplot(211)
periodogram(ecgorg);
title('PSD by Periodogram with Rectangular window');
subplot(212)
periodogram(ecgorg,hamming(length(ecgorg)));
title('PSD by Periodogram with Hamming window');

%% Welch Method

figure;
subplot(211)
pwelch(ecgorg,rectwin(256));
title('PSD by Welch Method with Rectangular window with l=256');
subplot(212)
pwelch(ecgorg,hamming(256));
title('PSD by Welch Method with Hamming window with l=256');

%% Second Step Filtering

[b3,a3]=cheby1(ord,30,Wp);
[b4,a4]=cheby2(ord,50,Wp);
[b5,a5]=fir1(ord,Wp,"low",kaiser(ord+1,1));
[b6,a6]=fir1(ord,Wp,"low",hamming(ord+1));

Butter_Cheby1_Filtered_ECG=filtfilt(b3,a3,Butter_Filtered_ECG);
Butter_Cheby2_Filtered_ECG=filtfilt(b4,a4,Butter_Filtered_ECG);
Fir_Kaiser_Filtered_ECG=filtfilt(b5,a5,Fir_Filtered_ECG);
Fir_Hamming_Filtered_ECG=filtfilt(b6,a6,Fir_Filtered_ECG);
MA_Kaiser_Filtered_ECG=filtfilt(b5,a5,MA_Filtered_ECG);
MA_Hamming_Filtered_ECG=filtfilt(b6,a6,MA_Filtered_ECG);
Gaussian_Kaiser_Filtered_ECG=filtfilt(b5,a5,Gaussian_Filtered_ECG);
Gaussian_Hamming_Filtered_ECG=filtfilt(b6,a6,Gaussian_Filtered_ECG);

RMSE_Butter_Cheby1=sqrt(mean((ecgorg-Butter_Cheby1_Filtered_ECG).^2));
RMSE_Butter_Cheby2=sqrt(mean((ecgorg-Butter_Cheby2_Filtered_ECG).^2));
RMSE_FIR_Kaiser=sqrt(mean((ecgorg-Fir_Kaiser_Filtered_ECG).^2));
RMSE_FIR_Hamming=sqrt(mean((ecgorg-Fir_Hamming_Filtered_ECG).^2));
RMSE_MA_Kaiser=sqrt(mean((ecgorg-MA_Kaiser_Filtered_ECG).^2));
RMSE_MA_Hamming=sqrt(mean((ecgorg-MA_Hamming_Filtered_ECG).^2));
RMSE_Gaussian_Kaiser=sqrt(mean((ecgorg-Gaussian_Kaiser_Filtered_ECG).^2));
RMSE_Gaussian_Hamming=sqrt(mean((ecgorg-Gaussian_Hamming_Filtered_ECG).^2));

%% Second Filtering Results
fft_butter_cheby1=fft(Butter_Cheby1_Filtered_ECG)/length(Butter_Cheby1_Filtered_ECG);
fft_butter_cheby2=fft(Butter_Cheby1_Filtered_ECG)/length(Butter_Cheby1_Filtered_ECG);
fft_fir_kaiser=fft(Butter_Cheby1_Filtered_ECG)/length(Butter_Cheby1_Filtered_ECG);
fft_fir_hamming=fft(Butter_Cheby1_Filtered_ECG)/length(Butter_Cheby1_Filtered_ECG);
fft_ma_kaiser=fft(Butter_Cheby1_Filtered_ECG)/length(Butter_Cheby1_Filtered_ECG);
fft_ma_hamming=fft(Butter_Cheby1_Filtered_ECG)/length(Butter_Cheby1_Filtered_ECG);
fft_gaussian_kaiser=fft(Butter_Cheby1_Filtered_ECG)/length(Butter_Cheby1_Filtered_ECG);
fft_gaussian_hamming=fft(Butter_Cheby1_Filtered_ECG)/length(Butter_Cheby1_Filtered_ECG);
figure;
subplot(211);
plot(f,abs(fftshift(fft_ecg)));
title('FFT of Signal');
xlabel('Frequency(Hz)');
ylabel('Magnitude');
hold on
plot(f,abs(fftshift(fft_ecg_Butter_filtered)));
plot(f,abs(fftshift(fft_butter_cheby1)));
plot(f,abs(fftshift(fft_butter_cheby2)));
legend("Original","Butter","Butter and Cheby1","Butter and Cheby2");
hold off
subplot(212);
plot(f,unwrap(angle(fft_ecg)));
title('Phase of Signal');
xlabel('Frequency(Hz)');
ylabel('Degree');
hold on
plot(f,unwrap(angle(fft_ecg_Butter_filtered)));
plot(f,unwrap(angle(fft_butter_cheby1)));
plot(f,unwrap(angle(fft_butter_cheby2)));
legend("Original","Butter","Butter and Cheby1","Butter and Cheby2");
hold off

figure;
subplot(211);
plot(f,abs(fftshift(fft_ecg)));
title('FFT of Signal');
xlabel('Frequency(Hz)');
ylabel('Magnitude');
hold on
plot(f,abs(fftshift(fft_ecg_Fir_filtered)));
plot(f,abs(fftshift(fft_fir_kaiser)));
plot(f,abs(fftshift(fft_fir_hamming)));
legend("Original","FIR","FIR and Kaiser","FIR and Hamming");
hold off
subplot(212);
plot(f,unwrap(angle(fft_ecg)));
title('Phase of Signal');
xlabel('Frequency(Hz)');
ylabel('Degree');
hold on
plot(f,unwrap(angle(fft_ecg_Fir_filtered)));
plot(f,unwrap(angle(fft_fir_kaiser)));
plot(f,unwrap(angle(fft_fir_hamming)));
legend("Original","FIR","FIR and Kaiser","FIR and Hamming");
hold off

figure;
subplot(211);
plot(f,abs(fftshift(fft_ecg)));
title('FFT of Signal');
xlabel('Frequency(Hz)');
ylabel('Magnitude');
hold on
plot(f,abs(fftshift(fft_ecg_MA_filtered)));
plot(f,abs(fftshift(fft_ma_kaiser)));
plot(f,abs(fftshift(fft_ma_hamming)));
legend("Original","MA","MA and Kaiser","MA and Hamming");
hold off
subplot(212);
plot(f,unwrap(angle(fft_ecg)));
title('Phase of Signal');
xlabel('Frequency(Hz)');
ylabel('Degree');
hold on
plot(f,abs(fftshift(fft_ecg_MA_filtered)));
plot(f,unwrap(angle(fft_ma_kaiser)));
plot(f,unwrap(angle(fft_ma_hamming)));
legend("Original","MA","MA and Kaiser","MA and Hamming");
hold off

figure;
subplot(211);
plot(f,abs(fftshift(fft_ecg)));
title('FFT of Signal');
xlabel('Frequency(Hz)');
ylabel('Magnitude');
hold on
plot(f,abs(fftshift(fft_ecg_gauss_filtered)));
plot(f,abs(fftshift(fft_gaussian_kaiser)));
plot(f,abs(fftshift(fft_gaussian_hamming)));
legend("Original","Gaussian","Gaussian and Kaiser","Gaussian and Hamming");
hold off
subplot(212);
plot(f,unwrap(angle(fft_ecg)));
title('Phase of Signal');
xlabel('Frequency(Hz)');
ylabel('Degree');
hold on
plot(f,abs(fftshift(fft_ecg_gauss_filtered)));
plot(f,unwrap(angle(fft_gaussian_kaiser)));
plot(f,unwrap(angle(fft_gaussian_hamming)));
legend("Original","Gaussian","Gaussian and Kaiser","Gaussian and Hamming");
hold off

figure;
plot(t,ecgorg);
title('Time domain Representation of Signal');
xlabel('Time(sec)');
ylabel('Amplitude');
hold on
plot(t,Butter_Filtered_ECG);
plot(t,Butter_Cheby1_Filtered_ECG);
plot(t,Butter_Cheby2_Filtered_ECG);
legend("Original","Butter","Butter and Cheby1","Butter and Cheby2");
hold off

figure;
plot(t,ecgorg);
title('Time domain Representation of Signal');
xlabel('Time(sec)');
ylabel('Amplitude');
hold on
plot(t,Fir_Filtered_ECG);
plot(t,Fir_Kaiser_Filtered_ECG);
plot(t,Fir_Hamming_Filtered_ECG);
legend("Original","FIR","FIR and Kaiser","FIR and Hamming");
hold off

figure;
plot(t,ecgorg);
title('Time domain Representation of Signal');
xlabel('Time(sec)');
ylabel('Amplitude');
hold on
plot(t,MA_Filtered_ECG);
plot(t,MA_Kaiser_Filtered_ECG);
plot(t,MA_Hamming_Filtered_ECG);
legend("Original","MA","MA and Kaiser","MA and Hamming");
hold off

figure;
plot(t,ecgorg);
title('Time domain Representation of Signal');
xlabel('Time(sec)');
ylabel('Amplitude');
hold on
plot(t,Gaussian_Filtered_ECG);
plot(t,Gaussian_Kaiser_Filtered_ECG);
plot(t,Gaussian_Hamming_Filtered_ECG);
legend("Original","Gaussian","Gaussian and Kaiser","Gaussian and Hamming");
hold off
