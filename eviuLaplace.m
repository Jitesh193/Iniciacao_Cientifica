function [v,gamma,err]=eviuLaplace(v0,Lambda,Delta,H,iR,z)
% Return the EVIU variation estimate
[N,n]=size(H);
M=Lambda+H'*iR*H;
b=H'*iR*(z+H*v0);
if sum(diag(Delta)>=0)==n %check local convexity condition
[v,err,gamma]=cyclicSubCoordDescent(M,b,diag(Delta));
else
 v=M\H'*iR*(z+H*v0);
 gamma=0;
 err=0;
end
end