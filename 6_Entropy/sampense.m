function [se,e]=sampense(y,M,r)
%function [e,A,B]=sampense(y,M,r);
%
%Input
%
%y input data
%M maximum template length
%r matching tolerance
%
%Output
%
%se standard error estimates for m=0,1,...,M-1
%e  sample entropy estimates for m=0,1,...,M-1

[F1,R1,F2,R2]=makerun(y,M,r);
F=F1+F2;
n=length(y);
dd=1;
K0=sum(F.*(F-1));
K(1,:)=K0;
for m=1:M
   for d=1:min(m+1,M)
      i1=(d+1):n;
      i2=i1-d;     
      nm1=F1(i1,m);
      nm2=F2(i1,m);
      nm3=F1(i2,m);
      nm4=F2(i2,m);         
      nm1=nm1-sum(R1(i1,1:(dd-1))>=m,2);
      nm2=nm2-sum(R2(i1,1:(2*d))>=m,2);
      nm3=nm3-sum(R1(i2,1:(2*d-1))>=m,2);
      nm4=nm4-sum(R2(i2,1:(dd-1))>=m,2);         
      K(d+1,m)=2*sum((nm1+nm2).*(nm3+nm4));
   end     
end

n1=zeros(M,1);
n2=zeros(M,1);
n1(1)=n*(n-1)*(n-2);
for m=1:(M-1),
   n1(m+1)=sum(K(1:(m+1),m));
end
for m=1:M,
   n2(m)=sum(K(1:m,m));
end

A=sum(F)'/2;
N=n*(n-1)/2;
B=[N;A(1:(M-1))];
p=A./B;
e=-log(p);

vp=p.*(1-p)./B+max((n2-n1.*p.^2)./B.^2,0);
sp=sqrt(vp);
se=sp./p;

function [F1,R1,F2,R2]=makerun(y,M,r);
%function [F1,R1,F2,R2]=makerun(y,M,r);
%
%Input
%
%y input data
%M maximum template length
%r matching tolerance
%
%Output
%
%F1 matches with future points
%R1 runs with future points
%F2 matches with past points
%R2 runs with past points

n=length(y);
run1=zeros(1,n);
MM=2*M;
R1=zeros(n,MM);
R2=zeros(n,MM);
F=zeros(n,M);
F1=zeros(n,M);
for i=1:(n-1)
   j=(i+1):n;
   match=abs(y(j)-y(i))<r;
   k=find(match);
   nj=length(j);
   run=zeros(1,length(j));
   run(k)=run1(k)+1;
   for m=1:M
      k=find(run>=m);      
      nm=length(k);
      F1(i,m)=nm;
      F(i,m)=F(i,m)+nm;
      F(i+k,m)=F(i+k,m)+1;
   end
   nj=min(nj,MM);
   k=(1:nj);
   R1(i,k)=run(k);
   run1=run;   
end
for i=1:n
   nj=min(MM,i-1);
   for j=1:nj
      R2(i,j)=R1(i-j,j);
   end
end
F2=F-F1;
