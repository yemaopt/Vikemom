function B=rowmove(A,x)
[m,n]=size(A);
if x>0
   C=A(:,1:n-x);
   D=zeros(m,x);
   B=cat(2,D,C);
end
if x<0
   C=A(:,abs(x)+1:n);
   D=zeros(m,abs(x));
   B=cat(2,C,D);
end
if x==0
   B=A;   
end
end