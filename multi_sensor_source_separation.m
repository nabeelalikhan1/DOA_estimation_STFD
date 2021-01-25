function [ss] = multi_sensor_source_separation(X, num, delta,n_sensors)
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
for i = 1:num
    %Spec=tfr_stft_high(Sig);
    I2=quadtfd(X(1,:),length(X)/2-1,1,'wvd',length(X));
    for jj=2:n_sensors
        I2=I2+quadtfd(X(jj,:),length(X)/2-1,1,'wvd',length(X));
    end
    [ADTFD1,orienttfd1]=post_processing_directional(I2,2,15,84);
    [ADTFD2,orienttfd2]=post_processing_directional(I2,2,30,84);
    
    ADTFD=min(ADTFD2,ADTFD1);
    for ii=1:length(ADTFD2)
        for jj=1:length(ADTFD2)
            value=min(ADTFD2(ii,jj),ADTFD1(ii,jj));
            ADTFD(ii,jj)=value;
            if ADTFD2(ii,jj)==value
                orienttfd(ii,jj)=orienttfd1(ii,jj);
            else
                orienttfd(ii,jj)=orienttfd2(ii,jj);
                
            end
            
            
        end
    end
    c = findridges_new1(ADTFD,orienttfd,delta);
%    figure;
 %   imagesc(ADTFD);
    %c = findridges_neww(Spec,orienttfd,delta);
    %figure;
    %tfsapl(X(1,:),ADTFD,'TFfontSize',30,'Grayscale', 'on');
    IF=(c)/(2*length(X));
    
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
    end
    
end

end