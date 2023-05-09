function [W1] = WD_2D(Ex)

%% Summary
% Inspired by 1D Implementation: 
% https://www.mathworks.com/matlabcentral/fileexchange/15637-calculate-wigner-distribution
% 
% Input is EX: A 2D signal/image that has an even number of samples in x and y, 
%              Otherwise the code throws an error
%
% Output: W1 is the 4D wigner distribution of the input image, 
%         The frequency domain sampling rate is twice that of the original
%         sample, which gives smears the spectral energy into different frequency bins, 
%         as evidenced by the dual domain marginal
% 
% Requires the Matlab Version >= R2022a due to the use of tensorprod
% function:% https://www.mathworks.com/help/matlab/ref/tensorprod.html
%%

if rem(prod(size(Ex)-1), 2) == 0
    error('One or both of the image dimensions have an odd number of samples')
end



sx = 1; sy = 1;



%WD_2D_MYWIGNER Summary of this function goes here
%   Detailed explanation goes here
%   Ex is a 2D function/image

[M, N] = size(Ex);

% Generate Frequency Vector for each dimension
fx = ifftshift(((0:N-1)'-N/2)*2*pi/(N-1));
fy = ifftshift(((0:M-1)'-M/2)*2*pi/(M-1));

% Generate Frequency Grids, n indicates Grids
[fxn, fyn] = meshgrid(fx, fy);


% Generate Shifts
SX = (0:1/sx:N-1/sx)-N/2;
SY = (0:1/sy:M-1/sy)-M/2;


N1 = length(SX); M1 = length(SY);

% Generate Shift Grid
[SXn, SYn] = meshgrid(SX, SY);

% Create Copies of the Fourier transform
T1 = repmat(fft2(Ex), [1, 1, sx*M, sy*N]); 

% All Possible +ve X shifts
EX1 = ifft2(T1.*...
    exp(1i.*tensorprod(+fxn, SXn/2)).*exp(1i.*tensorprod(+fyn,  SYn/2)));  % Only need half shifts

% All Possible -ve X shifts
EX2 = ifft2(T1.*...
    exp(1i.*tensorprod(-fxn, SXn/2)).*exp(1i.*tensorprod(-fyn,  SYn/2)));  %%  % Only need half shifts


% Combine the negative and positive shifts, and transform back
T45 = fft(fft(fftshift(fftshift(EX1.*conj(EX2), 3), 4), [], 3), [], 4)./(M1*N1);
W1 = real(fftshift(fftshift(T45 , 3), 4));
W1 = W1(:,:, [M1/2-M/2+1:M1/2+M/2], [N1/2-(N/2-1):N1/2+N/2]);

end