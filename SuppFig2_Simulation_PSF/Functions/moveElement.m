% Written by Ye Ma, 
% Department of Optical Engineering, Zhejiang University
% Department of Biomedical Engineering, Johns Hopkins University

function [ B ] = moveElement( A,x,y)

B=rowmove(A,x);
B=colmove(B,y);


end

