close all; clear all; clc


%% Create 2D Signal

N = 32; M = 64;
[X,Y] = meshgrid(1:N,1:M); 
Z = 1/2*(1+sin(2*pi/N*3*X) + sin(2*pi/M*10*Y)); 



figure; sgtitle('Input Signal X: 2D Sinusoidal Image')
subplot(1, 2, 1); imagesc(Z); 
colorbar; title('X')
subplot(1, 2, 2); imagesc(abs(fftshift(fft2(Z)))); 
colorbar; title('FFT(X)')


%% Call 2D Wigner Distribution 

% Call Function
[W1] = WD_2D(Z);

% Plot the output
figure; 
subplot(2, 2, 1); imagesc(Z.^2); 
colorbar; title('X^2')

subplot(2, 2, 2); imagesc(squeeze(sum(W1, [3, 4]))); 
colorbar; title('Primal Domain Marginal')

% Normalize by number of samples to account for energy
subplot(2, 2, 3); imagesc(abs(fftshift(fft2(Z)).^2)./(prod(size(Z))));
colorbar; title('|FFT(X)^2|')

subplot(2, 2, 4); imagesc(squeeze(sum(W1, [1, 2]))); 
colorbar; title('Dual Domain Marginal')






