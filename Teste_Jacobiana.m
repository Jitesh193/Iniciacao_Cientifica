close all
clear all
clc
rng(101)
  
%% Real data
load pop_data_intRic_const

%% Parameters
T=max(data(:,2)); %time horizon
time=1:T;
M=66; %total of individuals
V0=1; %initial turmor size
n=4; %numer of parameters to estimate
m=1; %dimension of the measurement channel

%% Richards Model
%r => growth rate (theta(:,2))
%K => carrying capacity (mm^3) (theta(:,1))
%v => affects near which asymptote maximum growth occurs (theta(:,3))
%V => volume (mm^3)
h=@(theta,DeltaT)(theta(:,1).*V0)./(V0.^(theta(:,3))+(theta(:,1).^(theta(:,3))-V0.^(theta(:,3))).*exp(-theta(:,2).*(theta(:,3)).*(theta(:,4)-DeltaT))).^(1/theta(:,3));
ht=@(theta,t)(theta(:,1).*V0)./(V0.^(theta(:,3))+(theta(:,1).^(theta(:,3))-V0.^(theta(:,3))).*exp(-theta(:,2).*(theta(:,3)).*(t))).^(1./(theta(:,3)));
jacobian=@(x,DeltaT)[exp(x(1)).*(((exp(x(1))./V0).^(exp(x(3))) - 1).*exp(-exp(x(2) + x(3)).*(exp(x(4)) - DeltaT)) + 1).^(-exp(-x(3))) - ((exp(x(1))./V0).^(exp(x(3)) - 1).*exp(2.*x(1) - exp(x(2) + x(3)).*(exp(x(4)) - DeltaT)).*(((exp(x(1))./V0).^(exp(x(3))) - 1).*exp(-exp(x(2) + x(3)).*(exp(x(4)) - DeltaT)) + 1).^(-exp(-x(3)) - 1))./V0;...
    (exp(x(4)) - DeltaT).*((exp(x(1))./V0).^(exp(x(3))) - 1).*exp(-exp(x(2) + x(3)).*(exp(x(4)) - DeltaT) + x(2) + x(1)).*(((exp(x(1))./V0).^(exp(x(3))) - 1).*exp(-exp(x(2) + x(3)).*(exp(x(4)) - DeltaT)) + 1).^(-exp(-x(3)) - 1);...
    h(exp(x),DeltaT).*(exp(-x(3)).*log(((exp(x(1))./V0).^(exp(x(3))) - 1).*exp(-exp(x(2) + x(3)).*(exp(x(4)) - DeltaT)) + 1) - (exp(-x(3)).*((exp(x(1))./V0).^(exp(x(3))).*log(exp(x(1))./V0).*exp(x(3) - exp(x(2) + x(3)).*(exp(x(4)) - DeltaT)) - (exp(x(4)) - DeltaT).*((exp(x(1))./V0).^(exp(x(3))) - 1).*exp(-exp(x(2) + x(3)).*(exp(x(4)) - DeltaT) + x(2) + x(3))))./(((exp(x(1))./V0).^(exp(x(3))) - 1).*exp(-exp(x(2) + x(3)).*(exp(x(4)) - DeltaT)) + 1));...
    ((exp(x(1))./V0).^(exp(x(3))) - 1).*exp(-exp(x(2) + x(3)).*(exp(x(4)) - DeltaT) + x(2) + x(1) + x(4)).*(((exp(x(1))./V0).^(exp(x(3))) - 1).*exp(-exp(x(2) + x(3)).*(exp(x(4)) - DeltaT)) + 1).^(-exp(-x(3)) - 1)]';
%% Population prior
i=1; %prior without the first individue
htheta=thetaPop(:,i);
PSI=Qpop(:,i);
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
ht=@(theta,t)(theta(:,1).*V0)./(V0.^(theta(:,3))+(theta(:,1).^(theta(:,3))-V0.^(theta(:,3))).*exp(-theta(:,2).*(theta(:,3)).*(t))).^(1./(theta(:,3)));
M=1000; %Number of experiments
time=0:T;
y=zeros(M,length(time));
for j=1:M
    theta0=exp(mvnrnd(log(htheta),diag(PSI)));
    y(j,:)=ht(theta0,time);
end
ciY=CIFcn(y,95);
ym=mean(y);
% plot
figure
trainData=data(data(:,1)~=i,:);
plot(trainData(:,2),trainData(:,3),'k.','MarkerSize',10)
hold on
plot(time,ym,'-','color',colors('blue'),'linewidth',2)
xlabel('Time (days)')
ylabel('Volume (mm^3)')
fill([time flip(time)],[ciY(1,:) flip(ciY(2,:))],[0, 0.4470, 0.7410], 'FaceAlpha', 0.2,'linestyle','none');
% axis([0 time(end) 0 500])
datai=data(data(:,1)==i,:);
plot(datai(:,2),datai(:,3),'sr--')
legend('Data','Gompertz Model','95% CI','Individual','location','best')
%saveas(gcf, '../results/fig5.png')
%% Individual Estimation (Leave-one-out Strategy)
T=max(data(:,2)); %time horizon
time=1:T;
M=66; %total of individuals
V0=1; %initial turmor size
n=4; %numer of parameters to estimate
m=1; %dimension of the measurement channel
ki=1;
% for alpha=0:1:10 %
% alpha=1;
%     kj=1;
%     for beta=0:1:10    
% beta=0;
    age=zeros(M,1);
for i=1:M
    %% Prior from Population
    theta0=thetaPop(:,i);
    Q=diag(Qpop(:,i));
    gamma0=log(theta0);
    P0=Q;
    iP0=inv(P0);
    %% individual measurements
    data0=data(data(:,1)==i,:);
    y=data0(:,3); %measurements
    age(i)=data0(end,2);
    DeltaT=data0(:,4);
    N=length(y); 
    sigmaR=17;
    R=sigmaR^2;
    iR=1/R;
    W=eye(N)*iR;
    iterMax=30;
    %% MAP (Gauss-Newton)
    gamma2=zeros(4,0);
    gamma2(:,1)=gamma0;
    err2=zeros(0,1);
    err2(1)=inf;
    iter2=1;
    while err2(iter2)>1e-4 && iter2<iterMax
        z2=zeros(N,1);
        H2=zeros(N,n);
        for k=1:N
            z2(k)=y(k)-h(exp(gamma2(:,iter2))',DeltaT(k));
            H2(k,:)=jacobian(gamma2(:,iter2)',DeltaT(k));
        end
        Jrls(iter2)=0.5*z2'*W*z2+0.5*(gamma2(:,iter2)-gamma0)'*iP0*(gamma2(:,iter2)-gamma0);
        gamma2(:,iter2+1)=gamma0+(iP0+H2'*W*H2)\H2'*W*(z2+H2*(gamma2(:,iter2)-gamma0));
        err2(iter2+1)=norm(gamma2(:,iter2+1)-gamma2(:,iter2));
        iter2=iter2+1;
    end
    # Subgradient of the cost function 
    dJ(:,i) = (iP0 + H2'*W*H2)*(gamma2(:,end)-gamma0) - H2'*W*(z2+H2*(gamma2(:,end-1)-gamma0));
    Prls(:,:,i)=(iP0+H2'*W*H2)\eye(n); % approx. error variance
    %% results
    hthetaRls(:,i)=exp(gamma2(:,iter2)); %RLS
    hat_age2=hthetaRls(4,i);
    %% residual
    ageErrRls(i)=age(i)-hat_age2; %RLS
    
    
end

%% Jacobian Evaluation - using the expression (7) given in the article [1] in READme.md
% Expression's variables:
% Lambda = inv(diag(Qpop))
% H = jacobian of h
% R = variance of the measurer (given by matrix W)
% v = solution (hthetaRls)
% z = Y - H.log(u) -> measurements - H.thetaPop

dJ
