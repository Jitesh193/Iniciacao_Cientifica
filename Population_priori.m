close all
clear all
clc
rng(101)

%% Real data
%data=[id day volume]
data = load('LM2-4LUC.txt'); %load raw data
data(data(:,3)==0,:)=[]; %remove zeros
data(:,1)=data(:,1)+1; %convert non-zero-index
%% Parameters
T=max(data(:,2)); %time horizon = max time from the population
N=66; %total of individuals
V0=1; %initial turmor size (mm^3)

%% Richards Model 
%a => growth rate -> theta(:,2)
%K => carrying capacity (mm^3) -> theta(:,1)
%V => turmor size (mm^3)
%t => time -> theta(:,4)
%v => affects near which asymptote maximum growth occurs -> theta(:,3)
h=@(theta,DeltaT)(theta(:,1).*V0)./(V0.^(theta(:,3))+(theta(:,1).^(theta(:,3))-V0.^(theta(:,3))).*exp(-theta(:,2).*(theta(:,3)).*(theta(:,4)-DeltaT))).^(1/theta(:,3));

% h=@(theta,DeltaT)(theta(:,1))./(1+(-1+(theta(:,1)./V0).^(theta(:,3))).*exp(-theta(:,2).*theta(:,3).*(theta(:,4)-DeltaT))).^(1./theta(:,3));
%% Time Interval
k=1;
DeltaT=zeros(size(data,1),1); % Vector that stores the quantity. of total days passed from the last sample
for i=1:N
    data0=data(data(:,1)==i,:);
    diffTime=diff(data0(:,2));
    for j=1:size(data0,1)
        DeltaT(k)=sum(diffTime(j:end));
        k=k+1;
    end
end
%% Population Estimation (Leave-one-out Strategy)
thetaPop=zeros(4,N);
Qpop=zeros(4,N);
data=[data DeltaT];
for i=1:N
data0=data(data(:,1)~=i,:);
[thetaPop(:,i),Qpop(:,i)]=estimatePopulationIntRic(data0(:,4),data0(:,3),data0(:,1),h,log([2500 1.2 0.03 30]));
%% Confidence Interval
if i==1 || i==22
    htheta=thetaPop(1:3,i);
    PSI=Qpop(1:3,i);
    CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
    ht=@(theta,t)(theta(:,1).*V0)./(V0.^(theta(:,3))+(theta(:,1).^(theta(:,3))-V0.^(theta(:,3))).*exp(-theta(:,2).*(theta(:,3)).*(t))).^(1./(theta(:,3)));
    M=1000; %Number of experiments
    time=1:max(data(:,2));
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
    legend('Data','Richards Model','95% CI','Individual')
end

% if i==34
%     break
% end
end
%% save
save pop_data_intRic_teste thetaPop Qpop data
