function [F,P2]=computeFFT(X,Y)
%COMPUTEFFT Summary of this function goes here
%   Detailed explanation goes here

% doEx=1;
% if doEx
%    X=linspace(-10*pi,10*pi,1E3);
%    k=.2;
%    Y=exp(1i*k*X);
%     
% end


dX=X(2)-X(1);               % X spacing
Fs=1/dX;                    % Sampling Frequency
dF=Fs/length(Y);            % Delta frequency
F=-Fs/2+dF:dF:Fs/2;      % Frequency vector

F=2*F;
% F=F;

Y = fftshift(fft(Y));       % Double sided FFT
P2 = abs(Y/length(X));      % magnitudes

% figure
% plot(F*pi,P2);

% if psi=exp(1i*p*x/hb) and H=-hb^2/2m dx^2 then
% H*psi=-p^2/(hb^2)*(-hb^2/2m)*psi=p^2/2m
end

