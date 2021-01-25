function [c]=check_merge_spatial(IF,orient,i,j,ang,distance)
A=IF{i};
A=sortrows(A,2);
a=A(1,:);
y_si=a(1);
x_si=a(2);
a=A(end,:);
y_ei=a(1);
x_ei=a(2);

A=IF{j};
A=sortrows(A,2);
a=A(1,:);
y_sj=a(1);
x_sj=a(2);
a=A(end,:);
y_ej=a(1);
x_ej=a(2);

%A=IF{j};
%A=sortrows(A,2);

%[y_sj,x_sj]=A(1,:);
%[y_ej,x_ej]=A(end,:);

if x_si>x_ej
    
    dis=sqrt(abs(x_ej-x_si)^2+abs(y_ej-y_si)^2);
   d1=abs(x_ej-x_si);
   d2=abs(y_ej-y_si);
   
    if dis<distance
        %if abs(orient(x_ej,y_ej)-orient(x_si,y_si))<ang
        ang_dif=abs(orient(y_ej,x_ej)-orient(y_si,x_si));
        if ang_dif>90
        ang_dif=180-ang_dif;
        end
        if ang_dif<ang
                if and(d2/d1>2>2, dis>max(distance/2,5))
                    c=0;
                    return;
                end
                c=1;
            return;
        end
    end
    

elseif x_sj>x_ei
        dis=sqrt(abs(x_ei-x_sj)^2+abs(y_ei-y_sj)^2);
   d1=abs(x_ei-x_sj);
   d2=abs(y_ei-y_sj);
        
    if dis<distance
        ang_dif=abs(orient(y_sj,x_sj)-orient(y_ei,x_ei));
        if ang_dif>90
        ang_dif=180-ang_dif;
        end
        %if ang_dif<ang
            
        if ang_dif<ang
            if and(d2/d1>2, dis>max(distance/2,5))
                    c=0;
                    return;
            end

            
            c=1;
            return;
        end
    end
    
    
    else
        c=0;
        return;
end

c=0;
end
