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
%a => growth rate -> theta(:,2)
%K => carrying capacity (mm^3) -> theta(:,1)
%V => turmor size (mm^3)
%t => time -> theta(:,4)
%v => affects near which asymptote maximum growth occurs -> theta(:,3)

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
legend('Data','Richards Model','95% CI','Individual','location','best')
% saveas(gcf, '../results/fig5.png')

%% Individual Estimation (Leave-one-out Strategy)
T=max(data(:,2)); %time horizon
time=1:T;
M=66; %total of individuals
V0=1; %initial turmor size
n=4; %numer of parameters to estimate
m=1; %dimension of the measurement channel
ki=1;
% for alpha=1:1:11
alpha = 3.4;
   kj=1;
%    for beta=0:0.1:1   

beta=0.1;
% This 2 values was used to test the EVIU estimator, in order that EVIU results is equal to RLS results with this 2 values of extra parameters
% alpha=1;
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
    sigmaR=17; % Variancia minima do estimador que pegou as medidas!
    R=sigmaR^2;
    iR=1/R;
    W=eye(N)*iR;
    iterMax=30;
    %% EVIU
    gamma=zeros(4,0);
    gamma(:,1)=gamma0;
    err=zeros(0,1);
    err(1)=inf;
    iter=1;
    while err(iter)>1e-6 && iter<iterMax
        z=zeros(N,1);
        H=zeros(N,n);
        for k=1:N
            z(k)=y(k)-h(exp(gamma(:,iter))',DeltaT(k));
            H(k,:)=jacobian(gamma(:,iter)',DeltaT(k)); 
        end
        v0=gamma(:,iter)-gamma0;
        Alpha=alpha*1./diag(Q);
        Beta=beta*1./diag(Q);
        Lambda=diag(Alpha);
        Delta=diag(Beta);
        bar_v=eviuLaplace(v0,Lambda,Delta,H,iR,z);
        gamma(:,iter+1)=gamma0+bar_v;
        err(iter+1)=norm(gamma(:,iter+1)-gamma(:,iter));
        iter=iter+1;
    end
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
    Prls(:,:,i)=(iP0+H2'*W*H2)\eye(n); % approx. error variance
    %% results
    htheta(:,i)=exp(gamma(:,iter)); %EVIU-JENSEN
    hthetaRls(:,i)=exp(gamma2(:,iter2)); %RLS
    for k=1:N
        H(k,:)=jacobian(gamma(:,end)',DeltaT(k));
    end
    Peviu(:,:,i)=(iP0+H'*W*H+Lambda)\eye(n); % approx. error variance 
    hat_age=htheta(4,i); 
    hat_age2=hthetaRls(4,i);
    %% residual
    ageErr(i)=age(i)-hat_age; %EVIU-JENSEN
    ageErrRls(i)=age(i)-hat_age2; %RLS
 
end


% save ind_data_Richards htheta Peviu 

%%
figure
histogram(ageErrRls,'normalization','pdf')
err=-20:0.1:35;
hold on
histogram(ageErr,'normalization','pdf')
plot(err,normpdf(err,mean(ageErr),std(ageErr)),'Color',colors('red'),...
    'LineWidth',1.5)
plot(err,normpdf(err,mean(ageErrRls),std(ageErrRls)),'Color',colors('blue'),...
    'LineWidth',1.5)
legend('RLS','EVIU')
xlabel('Age Error (days)')
ylabel('pdf')
drawnow
% saveas(gcf, '../results/histogram.png')

%%
figure,plot(ageErr,'s-')
hold on
plot(ageErrRls,'o-')
grid on
xlabel('Individual')
ylabel('Age Error (days)')
legend('EVIU','RLS')
xlim([1,66])
axes('pos',[.5 .5 .3 .3])
boxplot([ageErr;ageErrRls]','Labels',{'EVIU','RLS'})
h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),colors('red'),'FaceAlpha',.5);
patch(get(h(2),'XData'),get(h(2),'YData'),colors('blue'),'FaceAlpha',.5);
ylabel('Age Error (days)')
% saveas(gcf, '../results/fig7.png')

%%
rmseJeviu=sqrt(mse(ageErr));
rmseRls=sqrt(mse(ageErrRls));
gain=(1-rmseJeviu/rmseRls)*100
% rmseJeviu(ki,kj)=sqrt(mse(ageErr));
% rmseRls(ki,kj)=sqrt(mse(ageErrRls));
% kj=kj+1;
%      end
%      ki=ki+1;
%      fprintf("ok %d%d\n",ki,kj)
%      save experiment034 rmseJeviu rmseRls
%      if ki == 112
%          break
%      end
% end
%%
id=[1 13 23 50];
%ageErrRls(23)=ageErrRls(23)-1;
for ij=1:4
subplot(2,2,ij)
i=id(ij);
htheta0=thetaPop(1:4,i);
PSI=Qpop(1:4,i);
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
M=1000; %Number of experiments
if i==23
    time=0:50;
end
yt=zeros(M,length(time));
for j=1:M
    theta0=exp(mvnrnd(log(htheta0),diag(PSI)));
    yt(j,:)=ht(theta0,time);
end
ciY=CIFcn(yt,95);
ym=mean(yt);
% plot
%figure
hold on
p1=plot(time,ym,'-','color',colors('yellow'),'linewidth',2);
fill([time flip(time)],[ciY(1,:) flip(ciY(2,:))],colors('yellow'), 'FaceAlpha', 0.2,'linestyle','none');
datai=data(data(:,1)==i,:);

xlabel('Time (days)')
ylabel('Volume (mm^3)')

yt=zeros(M,length(time));
for j=1:M
    theta0=exp(mvnrnd(log(htheta(1:4,i)),diag(diag(Peviu(1:4,1:4)))));
    yt(j,:)=ht(theta0,time);
end
ciY=CIFcn(yt,95);
ym=mean(yt);
p2=plot(time+ageErr(i),ym,'-','color',colors('blue'),'linewidth',2);
fill([time flip(time)]+ageErr(i),[ciY(1,:) flip(ciY(2,:))],colors('blue'), 'FaceAlpha', 0.2,'linestyle','none');


yt=zeros(M,length(time));
for j=1:M
    theta0=exp(mvnrnd(log(hthetaRls(1:4,i)),diag(diag(Prls(1:4,1:4)))));
    yt(j,:)=ht(theta0,time);
end
ciY=CIFcn(yt,95);
ym=mean(yt);
p3=plot(time+ageErrRls(i),ym,'-','color',colors('red'),'linewidth',2);
fill([time flip(time)]+ageErrRls(i),[ciY(1,:) flip(ciY(2,:))],colors('red'), 'FaceAlpha', 0.2,'linestyle','none');
title(['#' num2str(id(ij))])
xlim([0 30])
if i==23
xlim([0 30])
end
pd=plot(datai(:,2),datai(:,3),'k.','MarkerSize',10);
end
legend([p1,pd,p2,p3],'Prior','Individual Data','EVIU','RLS','Orientation','horizontal')
% saveas(gcf, '../results/fig8.png')

