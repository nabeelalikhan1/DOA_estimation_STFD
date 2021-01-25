function [IF_est]= fill_zeros(IF,L)
if isempty(IF)
    IF_est=-1;
return 
end
[aa,~]=size(IF);
IF_est=zeros(aa,L);
 for i=1:length(IF)
     A=IF{i};
     t_axis=A(:,2);
     f_axis=A(:,1);
     for tt=1:length(A)
         IF_est(i,t_axis(tt))=f_axis(tt);
     end
 end

for i=1:length(IF)
    jjj=1;

while jjj<length(IF_est)
    a=find(IF_est(i,jjj:end)==0,1);
    if isempty(a)
        break;
    elseif a==1
     jjj=jjj+1;
     continue;
    else
            a=a+jjj-1;

        b=find(IF_est(i,a:end)~=0,1);
        if isempty(b)
            break;
        else
            b=b+a-1;
            x=1:b-a;
            m=(IF_est(i,b)-IF_est(i,a-1))/(b-a);
            y=m*x+IF_est(i,a-1);
            IF_est(i,a:b-1)=y;
                    jjj=b;

        end
    end
end
end
%figure; plot(IF_est.')
end

