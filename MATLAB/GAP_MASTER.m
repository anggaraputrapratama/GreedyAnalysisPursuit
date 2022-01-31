tic
clc;
clear all;
close all;
load('512 Sampel/eeg.mat');  %raw data;
load('512 Sampel/BT50.mat');%trans_matrix
A = BT50;
N = length(eeg);
D = zeros(size(A,2));
% D(1,:)=eeg;
f = @(x) A* x;%function handle
x=fft(eeg); %Trasnformasi Sparsity (input data)
y=A*x;  


%------Reconstruction------%
xr= GAP(y,f,D,N); %Reconstruction GAP
xr_i = ifft(xr); %ifft
xreal = real(xr_i); %real GAP

%-----Plot------%
subplot(3,1,1); plot(eeg);
title('Original signal')
subplot(3,1,2); plot(xreal);
title('Greedy Analysis Pursuit (Reconstructed)')

subplot(3,1,3);
stem(real(eeg), 'blue');hold on;stem(xreal,'filled');hold on;
legend('Original','GAP');

%------------MEAN ABSOLUTE PERCENTAGE ERROR (MAPE)--------------%
abseeg = max(abs(eeg));
MAE = sum(abs(eeg - xreal))./ numel(xreal);
MAPE = (MAE/abseeg);
fprintf('MAE: %d\n', MAE);
fprintf('MAPE: %d\n', MAPE);

%------------MEAN SQUARE ERROR--------------%

MSE = sum(((eeg - xreal).^2)/length(xreal));
% RootSquare = sqrt (MSE);
MSE1 = immse(eeg, xreal);
fprintf('MSE: %d\n', MSE);

S = sum(eeg);
noi = sum(xreal);
Noise = S-noi;
SNR = mag2db(S/Noise);
SNR1 = db(S/Noise);
fprintf('SNR: %d\n', SNR);
fprintf('SNR: %d\n', SNR1);
toc
