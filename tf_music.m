%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [1] B. Boashash, A. Aissa-El-Bey, Multisensor Time-Frequency Signal Processing:
%     A tutorial review with illustrations in selected application areas, Digital
%     Signal Processing, In Press.
% [2] B. Boashash, A. Aissa-El-Bey, M. F. Al-Sa'd, Multisensor time-frequency
%     signal processing software Matlab package: An analysis tool for multichannel
%     non-stationary data , SoftwareX, In Press.
%
% Last Modification: 25-04-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 DOA estimation using TF MUSIC algorithm
%
% Syntax : P = tf_music(Ds, n, m, lamda, d, theta)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <INPUTs>
% Ds       : Selected MTFD Matrix (m x m).
% n        : Number of sources.
% m        : Number of sensors.
% lamda    : Wavelength of the received signal.
% d        : ULA elements spacing.
% theta    : Solution space for the 'MUSIC'.
%
% <OUTPUTs>
% P        : Estimated Spectrum for 'MUSIC'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function P = tf_music(Ds, n, m, lamda, d, theta)
theta   = theta/180*pi;
theta_N = length(theta);
[vv,~]  = svd(Ds);        % Find the eigenvalues and eigenvectors of STFD
UN_TF   = vv(:,n+1:m);    % Estimate/Selection of noise subspace

%% TF MUSIC Main
P = zeros(1,theta_N);
for ii = 1:theta_N
    a_theta = exp(-1j*2*pi*(d/lamda)*cos(theta(ii))*(0:m-1));
    PP_TF   = conj(a_theta)*(UN_TF*UN_TF')*a_theta.';
    P(ii)   = abs(1/PP_TF);
end
P = P/max(P);
end