function [ IF_new,NC ] = merge_cells( IF,i,j )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
A=IF{i};
B=IF{j};
A=[A;B];

IF{i}=sortrows(A,1);
IF{j}=[];
NC=0;
for jj=1:length(IF)
    if jj==j
    continue;
    else
        NC=NC+1;
        IF_new{NC}=IF{jj};
        
    end

end
end