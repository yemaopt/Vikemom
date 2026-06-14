% Written by Ye Ma, 
% Department of Optical Engineering, Zhejiang University
% Department of Biomedical Engineering, Johns Hopkins University

function B=colmove(A,x)
[m,n]=size(A);
if x>0
   C=A(1:m-x,:);
   D=zeros(x,n);
   B=[D' C']';
end
if x<0
   C=A(abs(x)+1:m,:);
   D=zeros(abs(x),n);
   B=[C' D']';
end
if x==0
   B=A;   
end
end
