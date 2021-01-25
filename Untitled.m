close all;
clear all;
addpath('D:\tfsa_5-5\windows\win64_bin');
%crossing componentsi8
n=0:127;
s1=exp(2*pi*1i*(0.05*n+0.45*n.^3/(128*128*3)));
%s2=exp(2*pi*1i*(0.2*n+0.1*n.^2/(2*128)+0.15*n.^3/(128*128*3)));
s3=exp(2*pi*1i*(0.45*n-0.45*n.^3/(128*128*3)));%+exp(2*pi*1i*(0.45*n+0*0.1*n.^2/(2*128)-0*0.5*n.^3/(128*128*3)));
%s4=exp(2*pi*1i*(0.45*n+0.1*n.^2/(2*128)-0.1*n.^3/(128*128*3)));
%crossing components LFM sources

s = [(s1.')  (s3.') ];

% s1=exp(2*pi*1i*(0.05*n+0*0.2*n.^2/(2*128)));
s2=exp(2*pi*1i*(0.1*n+0.2*n.^2/(2*128)));
s3=exp(2*pi*1i*(0.3*n-0.2*n.^2/(2*128)));
% s4=exp(2*pi*1i*(0.45*n-0*0.2*n.^2/(2*128)));

%s = [(s2.')  (s3.') ];
% s1=exp(2*pi*1i*(0.05*n+0*0.2*n.^2/(2*128)));
% s2=exp(2*pi*1i*(0.15*n+0*0.3*n.^2/(2*128)));
% s3=exp(2*pi*1i*(0.3*n-0*0.3*n.^2/(2*128)));
% s4=exp(2*pi*1i*(0.45*n-0*0.2*n.^2/(2*128)));
%s = [(s1.') (s2.') (s3.') (s4.')];

%+exp(2*pi*1i*(0.45*n+0*0.1*n.^2/(2*128)-0*0.5*n.^3/(128*128*3)));
N_sensors=5;
perc=0.4;

%s = [(s1') (s2') (s3') (s4')];
%s = [(s1.') (s2.') (s3.') ];
n_sources=2;
N_C=2;
s_orig=s;

% set mixing matrix A
%theta = [15,30,50]*pi/180;   % sensor separation angles in radians
LL=50;

theta = [-5,5]*pi/180;   % sensor separation angles in radians

A = exp(1j*pi*[0:N_sensors-1].'*sin(theta));  % mixing matrix A


X = A*s.';
SNR=0;
sigma = 10^(-SNR/20);
w = sigma*(randn(N_sensors,length(n)) + 1j*(randn(N_sensors,length(n))))/sqrt(2); % noise

X=X+w;
D   = mtfd(X, 'ckd',1, 0.05, 0.05, length(X));
%%% Averaged Auto-TFD
D_avg = zeros(length(X), length(X));
for mm = 1:N_sensors, D_avg = D{mm,mm} + D_avg; end
D_avg = D_avg./N_sensors;
%%% Selection of high-energy (t,f) points
thr = perc*max(max(D_avg));
Tr = abs(D_avg) >= thr;
[F_trace, ~] = find(Tr);
n_p = length(F_trace);
D_s = zeros(N_sensors, N_sensors);
for m1 = 1:N_sensors
    for m2 = 1:N_sensors
        D_s(m1,m2) = (1/n_p).*sum(sum(D{m1,m2}.*Tr));
    end
end
theta1=-90:1:90;

%%% DOA Estimation
P_tf_music_ckd = tf_music(D_s, n_sources, N_sensors, 2,1, theta1);

plot(P_tf_music_ckd);
hold on;
[D_avg,D,~]= SADTFD(X,2,20,84);
D_avg(D_avg<0)=0;
thr = 2*perc*max(max(D_avg));
Tr = abs(D_avg) >= thr;
[F_trace, ~] = find(Tr);
n_p = length(F_trace);
D_s = zeros(N_sensors, N_sensors);
for m1 = 1:N_sensors
    for m2 = 1:N_sensors
        D_s(m1,m2) = (1/n_p).*sum(sum(D{m1,m2}.*Tr));
    end
end
theta1=-90:1:90;

%%% DOA Estimation
P_tf_music_adtfd = tf_music(D_s, n_sources, N_sensors, 2,1, theta1);

close all;

plot(theta1,P_tf_music_ckd)
hold on;
plot(theta1,P_tf_music_adtfd,'r:')


