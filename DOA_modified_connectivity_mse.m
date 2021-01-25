%% computing MSE curve vs. SNR in under-determined scenario
%In this code we assume that the user add in this path the TFSAP toolboox.
%Please cite the following reference:

close all;
clear all;
addpath('D:\tfsa_5-5\windows\win64_bin');
%addpath('C:\Users\DR. Nabeel\Desktop\IF estimation using connectivity\DOA-estimation-of-intersecting-components-master\DOA-estimation-of-intersecting-components-master');

% NUmber of simulation runs
LL=500;

% Signal model%
%addpath 'D:\tfsa_5-5\windows\win64_bin\'
index=0;
n=0:127;
% Number of SOurces
% NUmber of components
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

for SNR=0:2:10
    
    for ii=1:LL
        M=N_sensors;% Number of sensors
        
        % Signal generation
        s = [(s1.') (s2.') ];
        s_orig=s;
        % set mixing matrix A
        A = exp(1j*pi*[0:M-1].'*sin(theta));  % mixing matrix A
        X = A*s.';                             % mixed source
        theta9=round(theta *180/pi);
        % generate noise
        sigma = 10^(-SNR/20);
        w = sigma*(randn(M,length(n)) + 1j*(randn(M,length(n))))/sqrt(2); % noise
        
        X=X+w;
        ss= multi_sensor_source_separation(X, N_C, 2,N_sensors);
        %%%%%%%%  BSS code  ends
        %%DOA estimation
        IP=1;
        %figure;
        clear y1;
        for iii=1:n_sources
            for jjj=1:N_sensors
                a(jjj,:)=ss(jjj,iii,:);
            end
            theta1=-90:1:90;
            
            p=TMMUSIC(cov(a'), 2, N_sensors, 1, 1, theta1');
            [x,y]=max(p);
            y1(iii)=y-91;
            % p=music((a), 1, N_sensors, 1, 1, theta1');
            % P(iii,:)=p;
        end
        
        if length(y1)>length(theta)
            y1=y1(1:2);
        elseif length(y1)<length(theta)
            y1(length(y1):length(theta))=0;
        end
        theta=3*theta/pi;
        y1=3*y1/180;
        
        mse_r(ii)=mean((sort(y1)-sort(theta)).^2);
        
        [ss,IFF] = multi_sensor_source_separation_spatial_diversity(X, 3,N_sensors,0);
        [aa,bb,~]=size(ss);
                clear y1;
                clear a;

        for iii=1:bb%n_sources
            for jjj=1:aa%N_sensors
                a(jjj,:)=ss(jjj,iii,:);
            end
            
            p=TMMUSIC(cov(a'), 2, N_sensors, 1, 1, theta1');
            [x,y]=max(p);
                y1(iii)=y-91;
            
        end
        
        
        if length(y1)>length(theta)
            y1=y1(1:n_sources);
        elseif length(y1)<length(theta)
            y1(length(y1):length(theta))=0;
        end
        theta=3*theta/pi;
        y1=3*y1/180;
        
        mse_c(ii)=mean((sort(y1)-sort(theta)).^2);
        
        
    end
    
    index=index+1;
    mse_connectivity(index)=mean(mse_c)
    mse_ride_direction(index)=mean(mse_r)
end
index=0;
SNR=0:2:10;

plot(SNR,10*(log10(mse_ride_direction)),'--md','linewidth',2);
hold on;
plot(SNR,10*(log10(mse_connectivity)),'r','linewidth',2);

xlabel('Signal to Noise Ratio');
ylabel('Mean Square Error (dB)');
