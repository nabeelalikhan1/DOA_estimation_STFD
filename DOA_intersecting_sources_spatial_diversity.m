close all;
clear all;
n=0:127;

%addpath('D:\D\win64_bin\win64_bin');
addpath('D:\tfsa_5-5\windows\win64_bin');
%addpath('C:\Users\DR. Nabeel\Desktop\IF estimation using connectivity\DOA-estimation-of-intersecting-components-master\DOA-estimation-of-intersecting-components-master');


%+exp(2*pi*1i*(0.45*n+0*0.1*n.^2/(2*128)-0*0.5*n.^3/(128*128*3)));
s2=exp(2*pi*1i*(0.11*n+1*0.3*n.^2/(128*3)));

s1=exp(2*pi*1i*(0.05*n+0.45*n.^3/(128*128*3)));

s2=exp(2*pi*1i*(1*0.125*n+1.15*0.2*n.^3/(128*128*3)));

s1=exp(2*pi*1i*(0.05*n+0.45*n.^3/(128*128*3)));

s2=exp(2*pi*1i*(1*0.125*n+1.15*0.2*n.^3/(128*128*3)));

s1=exp(2*pi*1i*(0.05*n+0.45*n.^3/(128*128*3)));

s = [(s1.') (s2.') ];

perc=0.4;
IF_O(1,1:128)=0.05+3*0.45*n.^2/(128*128*3);
IF_O(2,1:128)=0.125+1.15*3*0.2*n.^2/(128*128*3);
n_sources=2;
N_sensors=2;
N_C=2;
s_orig=s;

% set mixing matrix A
%theta = [15,30,50]*pi/180;   % sensor separation angles in radians
theta = [-30,30]*pi/180;   % sensor separation angles in radians

A = exp(1j*pi*[0:N_sensors-1].'*sin(theta));  % mixing matrix A
X = A*s.';                             % mixed source
theta9=round(theta *180/pi);
% generate noise
SNR=5;

sigma = 10^(-SNR/20);
w = sigma*(randn(N_sensors,length(n)) + 1j*(randn(N_sensors,length(n))))/sqrt(2); % noise

X=X+w;


%[ss,IFF] = multi_sensor_source_separation(X, N_C, 3,N_sensors);

%figure; plot(IFF');
 % hold on; plot(IF_O','r')
 %xlabel('time')
 %ylabel('frequency');
%for iii=1:n_sources
 %   for jjj=1:N_sensors
  %      a(jjj,:)=ss(jjj,iii,:);
   % end
    %theta1=-90:1:90;
    
    %p=TMMUSIC(cov(a'), 2, N_sensors, 1, 1, theta1');
    %[x,y]=max(p);
    %P(iii,:)=p;
    %end
%figure;
%plot(theta1,P')
%aa=zeros(1,length(theta1));
%aa(theta9+91)=1;
%hold on; stem(theta1,aa)
%xlabel('DOA Estimation (degrees)');
%ylabel('Spatial Spectrum');
%title('Viterbi Algorithm based on direction of ridges');


% IF estimation
[ss,IFF] = multi_sensor_source_separation_spatial_diversity(X, 3,N_sensors,1);
figure; plot(IFF');
hold on; plot(IF_O','r:')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Estimated IF vs Original IF using CKD');
figure;
[aa,bb,~]=size(ss);
for iii=1:bb%n_sources
    for jjj=1:aa%N_sensors
        a(jjj,:)=ss(jjj,iii,:);
    end
    theta1=-90:1:90;
    
    p=TMMUSIC(cov(a'), 2, N_sensors, 1, 1, theta1');
    [x,y]=max(p);
    P(iii,:)=p;
    
end

plot(theta1,P')
aa=zeros(1,length(theta1));
aa(theta9+91)=1;
hold on; stem(theta1,aa)
xlabel('DOA Estimation (degrees)');
ylabel('Spatial Spectrum');
title('CKD based DOA estimation');

%ss= multi_sensor_source_separation_viterbi(X, N_C, 3,N_sensors);

ss= multi_sensor_source_separation(X, N_C, 3,N_sensors);

%%DOA estimation
clear P;
IP=1;
figure;
for iii=1:n_sources
    for jjj=1:N_sensors
        a(jjj,:)=ss(jjj,iii,:);
    end
    theta1=-90:1:90;
    
   p=TMMUSIC(cov(a'), 2, N_sensors, 1, 1, theta1');
   [x,y]=max(p);
   % y1(iii)=y;
    P(iii,:)=p; 
  % p=music((a), 1, N_sensors, 1, 1, theta1');
   % P(iii,:)=p;
end

   plot(P')
aa=zeros(1,length(theta1));

aa(theta9+90)=1;

hold on; stem(aa)
xlabel('DOA Estimation (degrees)');
ylabel('Spatial Spectrum');

%
% [ss,IF_out] = multi_sensor_source_separation_spatial_TF_direction(X, N_C, 3,N_sensors);
%  figure; plot(IF_out');
%  hold on; plot(IF_O','r:')
%  xlabel('Time (s)')
%  ylabel('Frequency (Hz)')
%  title('Estimated IF vs Original IF using Spatial Adaptive TFD');
% figure;
% for iii=1:n_sources
%     for jjj=1:N_sensors
%         a(jjj,:)=ss(jjj,iii,:);
%     end
%     theta1=-90:1:90;
%     
%     p=TMMUSIC(cov(a'), 2, N_sensors, 1, 1, theta1');
%     [x,y]=max(p);
%     P(iii,:)=p;
% end
% 
% plot(theta1,P')
% aa=zeros(1,length(theta1));
% aa(theta9+91)=1;
% hold on; stem(theta1,aa)
% xlabel('DOA Estimation (degrees)');
% ylabel('Spatial Spectrum');
% title('ADTFD based DOA estimation using Spatial TFD');


% [ss,IFF] = multi_sensor_source_separation_spatial_diversity_post_proc(X, N_C, 3,N_sensors);
% figure; plot(IFF');
% hold on; plot(IF_O','r:')
% xlabel('Time (s)')
% ylabel('Frequency (Hz)')
% title('Estimated IF vs Original IF using Spatial Adaptive TFD');
% 
% figure;
% for iii=1:n_sources
%     for jjj=1:N_sensors
%         a(jjj,:)=ss(jjj,iii,:);
%     end
%     theta1=-90:1:90;
%     
%     p=TMMUSIC(cov(a'), 2, N_sensors, 1, 1, theta1');
%     [x,y]=max(p);
%     P(iii,:)=p;
% end
% 
% plot(theta1,P')
% aa=zeros(1,length(theta1));
% aa(theta9+91)=1;
% hold on; stem(theta1,aa)
% xlabel('DOA Estimation (degrees)');
% ylabel('Spatial Spectrum');
% title('SADTFD based DOA estimation');

%%DOA estimation
D   = mtfd(X, 'ckd',1, 0.05, 0.05, length(X));
%%% Averaged Auto-TFD
D_avg = zeros(length(X), length(X));
for mm = 1:N_sensors, D_avg = D{mm,mm} + D_avg; end
D_avg = D_avg./N_sensors;
%%% Selection of high-energy (t,f) points
thr = 0.4*max(max(D_avg));
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



%P=TMMUSIC(cov(X'), 2, N_sensors, n_sources, 1, theta1');
P = tf_music(D_s, n_sources, N_sensors, 2,1, theta1);

figure;
plot(theta1,P','r:')

hold on; stem(theta1,aa)
xlabel('DOA Estimation (degrees)');
ylabel('Spatial Spectrum');
title('Time domain method');

