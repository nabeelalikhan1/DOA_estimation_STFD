function [IF]= merge_IFs(IF,orient,ang,dist,min_length)

run_again=1;
NC=length(IF);
while run_again==1
    if NC==1
        break;
    end
    for i=1:NC-1
        
        for j=i+1:NC
            if check_merge_spatial(IF,orient,i,j,ang,dist)==1
                IF=merge_cells(IF,i,j);
                NC=NC-1;
                run_again=1;
                break;
            else
                run_again=0;
            end
        end
        if   run_again==1
            break;
        end

    end
end
jj=1;
IF1={};
for ii=1:length(IF)
    
    if length(IF{ii})>min_length
        IF1{jj}=IF{ii};
        jj=jj+1;
    end
end
if isempty(IF1)==0
IF=IF1;
end
end