function [v,errEnd,gamma,err]=cyclicSubCoordDescent(M,b,d)
%J=1/2v'Mv-b'v+c+d'|x|, M=M'>0
%dJ=Mv-b+gamma=0
n=size(M,1); %number of variables
D=diag(diag(M));
iD=inv(D);
k=1;
maxIter=1000*n;
err=zeros(0,1);
err(1)=inf;
v=M\b;
gamma=zeros(n,1);
lambda=iD*d;
while err(k)>1e-3 && k<=maxIter
    v0=v;
    gradF=M*v-b;
    i=mod(k,n)+1;
    z=max(0,abs(v(i)-1/M(i,i)*gradF(i))-lambda(i))*sign(v(i)-1/M(i,i)*gradF(i));
    p=min(d(i),abs(M(i,i)*v(i)-gradF(i)))*sign(M(i,i)*v(i)-gradF(i));
    v(i)=z;
    gamma(i)=p;
    k=k+1;
    err(k)=norm(M*v-b+gamma);
%     err(k)=norm(v-v0);
end
% if k>=maxIter
%     warning('Maximum number of iteration!')
% end
errEnd=err(end);