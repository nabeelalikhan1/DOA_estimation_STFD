function [fidexmult] = extridge_mult_new_modified1(Sig, num, delta)
%%%%%%%%%%%%%%%%%%%%%%%    input      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sig£ºmeasured signal,a row vector
% num: the number of the signal components
% delta£ºmaximum allowable frequency variation between two consecutive points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  output   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fidexmult: obtained ridge index for multi-component signals

for i = 1:num
%Spec=tfr_stft_high(Sig);

                [Spec2,orienttfd2]=HTFD_neww(Sig,2,15,84);
        [Spec1,orienttfd1]=HTFD_neww(Sig,2,30,84);

Spec=min(Spec1,Spec2);
   %     [Spec3,orienttfd3]=HTFD_neww(Sig,3,8,64);

%Spec=min(Spec,Spec3);

%figure; imagesc(Spec);
for ii=1:length(Spec1)
    for jj=1:length(Spec2)
            value=min(Spec1(ii,jj),Spec2(ii,jj));
            Spec(ii,jj)=value;
            if Spec1(ii,jj)==value
            orienttfd(ii,jj)=orienttfd1(ii,jj);
            else
            orienttfd(ii,jj)=orienttfd2(ii,jj);

            end
            
         
    end
end
c = findridges_new1(Spec,orienttfd,delta);
%c = findridges_neww(Spec,orienttfd,delta);


 IF=(c)/(2*length(Sig));

                Phase=2*pi*filter(1,[1 -1],IF);
                s_dechirp=exp(-1i*Phase);

                
                %im_label2=bwmorph(im_label2,'dilate',3);
                
                % For each sensor do the following steps
                
                L=2;
                %TF filtering for each sensor
                s1 = Sig.*(s_dechirp);
                s2=fftshift(fft(s1));
                s3=zeros(1,length(Sig));
                s3(128-L:128+L)=s2(128-L:128+L);
                s2(128-L:128+L)=0;
                extr_Sig=ifft(ifftshift(s3)).*conj(s_dechirp);
                s2=ifft(ifftshift(s2)).*conj(s_dechirp);
                
                %Sig(iii)=Sig(iii)-extr_Sig(iii);
                Sig=s2;%-extr_Sig(iii);

           % extr_Sig1=extr_Sig1+extr_Sig;
%hold on; quiver(1:256,1:256,cos(orienttfd*pi/180),sin(orienttfd*pi/180));
fidexmult(i,:) = c;

end

end
