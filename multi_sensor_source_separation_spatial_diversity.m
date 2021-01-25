function [ss_out,IF_out] = multi_sensor_source_separation_spatial_diversity(X, delta,n_sensors,dis)
% Extract ridges for multi-component signals.
% In each iteration,the signal component associated with the extrated ridge is
% reconstructed by the ICCD and then removed from the original signal so
% that the ridge curves of other signal components with smaller energies
% can be extracted in the subsequent iterations.
%%%%%%%%%%%%%%%%%%%%%%%    input      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sig£ºmeasured signal,a row vector
% SampFreq: sampling frequency
% num: the number of the signal components
% delta£ºmaximum allowable frequency variation between two consecutive points
% orderIF: the order of the Fourier model used for smoothing the extracted ridge curves
% bw£ºthe bandwidth of the ICCD (unit£ºHz); herein the ICCD can be regarded as a time-frequency filtering technique
% Nfrebin,window are two parameters for implementing the STFT
% alpha£ºTikhonov regularization parameter for ICCD.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  output   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fidexmult: obtained ridge index for multi-component signals
  
%D = mtfd(X, 'CKD',0.2,0.2);
D = mtfd(X, 'CKD',0.15,0.15);
D_avg = zeros(length(X), length(X));
for mm = 1:n_sensors, D_avg = D{mm,mm} + D_avg; end
D_avg=real(D_avg);

%[D_avg,D,~]=SADTFD_new(X,2,20,length(X)/2);


%    [D_avg,~]=post_processing_directional(D_avg,2,20,84);
  D_avg(D_avg<0)=0;
%figure;imagesc(D_avg);
if dis==1
    
figure;imagesc(D_avg)
set(gca,'YDir','normal');
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
%title('(a)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);
end
[IF1,p, o,orient1] = component_linking_spatial(D_avg,D,0.15,5,5);
% [D_avg,D,~]= SADTFD(X,2,20,84);
%     D_avg(D_avg<0)=0;
%     figure;imagesc(D_avg)
% [IF1,p, o,orient1] = component_linking_spatial(D_avg,D,0.15,8,4);
if dis==1
figure;imagesc(p)
set(gca,'YDir','normal');
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'FontSize',20);

figure;imagesc(orient1)
set(gca,'YDir','normal');
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'FontSize',20);

figure;imagesc(o)
set(gca,'YDir','normal');
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'FontSize',20);
end

[IF1]= merge_IFs(IF1,orient1,30,60,length(X)/2);%[IF,out, peaks] = component_linking_neww_direc(tfd,orient,0.2,length(Sig)/4,30,5);



IF_est=fill_zeros(IF1,length(X));

IF_image=zeros(size(o));
IF_est=round(IF_est);
[num,L]=size(IF_est);

%[IF_est,~] = RPRG(IF_est,20);
for ii=1:num
    for jj=1:L
        if IF_est(ii,jj)~=0
        IF_image(IF_est(ii,jj),jj)=1;
        end
    end
end
if dis==1
figure;imagesc(IF_image)
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
set(gca,'FontSize',20);
end
IF_est=IF_est/(2*L);


for i = 1:num
IF=IF_est(i,:);
   % c = findridges_new_viterbi_mtfd(D_avg,orient);
    Phase=2*pi*filter(1,[1 -1],IF);
    s_dechirp=exp(-1i*Phase);
    
    %im_label2=bwmorph(im_label2,'dilate',3);
    
    % For each sensor do the following steps
    
    L=delta;
    %TF filtering for each sensor
    for iii=1:n_sensors
        
        s=(X(iii,:));
        %TF filtering for each sensor
        s1 = s.*(s_dechirp);
        s2=fftshift(fft(s1));
        s3=zeros(1,length(s));
        s3(length(s)/2-L:length(s)/2+L)=s2(length(s)/2-L:length(s)/2+L);
        s2(length(s)/2-L:length(s)/2+L)=0;
        s1=ifft(ifftshift((s3))).*conj(s_dechirp);
        s4=ifft(ifftshift((s2))).*conj(s_dechirp);
        X(iii,:)=s4;
        ss(iii,i,:)=s1;
        % nly for visualization    TFD of each extracted components
        %        l=l+tfrwv(conj(s1'));
    end
        ss_out=ss;
        IF_out=IF_est;
    
end

end