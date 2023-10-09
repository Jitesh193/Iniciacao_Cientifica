close all
clear all
clc
rng(101)

%% Carregar os dados populacionais
load pop_data_intRic_const

V0=1;
%% Carregar os dados obtidos do Gompertz
load ind_data_Gompertz

hthetaGom = htheta;
PeviuGom = Peviu;
EMQeviu_Gom = EMQeviu;
htGom=@(theta,t)theta(:,1).*exp(log(V0./theta(:,1)).*exp(-theta(:,2).*t));

%% Carregar os dados obtidos do Richards
load ind_data_Richards

hthetaRic = htheta;
PeviuRic = Peviu;
EMQeviu_Ric = EMQeviu;
htRic=@(theta,t)(theta(:,1).*V0)./(V0.^(theta(:,3))+(theta(:,1).^(theta(:,3))-V0.^(theta(:,3))).*exp(-theta(:,2).*(theta(:,3)).*(t))).^(1./(theta(:,3)));

%% Calculo dos erros

T=max(data(:,2)); %time horizon
time=1:T;
M=66;
age=zeros(M,1);
for i=1:M
    data0=data(data(:,1)==i,:);
    age(i)=data0(end,2);
    hat_age=hthetaGom(3,i); 
    hat_age2=hthetaRic(4,i);
    % hat_age3=hthetaLog(3,i);
    %% residual
    ageErr(i)=age(i)-hat_age; %Gompertz
    ageErrRic(i)=age(i)-hat_age2; %Richards
end

%% Graficos e ganho 

% Histograma de erro de idade 
figure
histogram(ageErrRic,'normalization','pdf')
err=-20:0.1:35;
hold on
histogram(ageErr,'normalization','pdf')
plot(err,normpdf(err,mean(ageErr),std(ageErr)),'Color',colors('red'),...
    'LineWidth',1.5)
plot(err,normpdf(err,mean(ageErrRic),std(ageErrRic)),'Color',colors('blue'),...
    'LineWidth',1.5)
% plot(err,normpdf(err,mean(ageErrLog),std(ageErrLog)),'Color',colors('green'),...
%     'LineWidth',1.5)
legend('Richards','Gompertz')
% legend('Richards','Gompertz','Logistic')
xlabel('Age Error (days)')
ylabel('pdf')
drawnow
% saveas(gcf, '../histogram_Models.png')

% Histograma de erro entre medidas
figure
histogram(EMQeviu_Gom,10,'normalization','pdf')
hold on
histogram(EMQeviu_Ric,10,'normalization','pdf')
legend('Gompertz','Richards')
xlabel('Measurements error')
ylabel('pdf')
drawnow

% Erro entre medidas para cada individuo
figure,plot(EMQeviu_Gom,'s-')
hold on
plot(EMQeviu_Ric,'o-')
grid on
xlabel('Individual')
ylabel('Erro Medidas (mm3)')
legend('Gompertz','Richards')
xlim([1,66])
axes('pos',[.5 .5 .3 .3])
boxplot([EMQeviu_Gom;EMQeviu_Ric]','Labels',{'Gompertz','Richards'})
h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),colors('red'),'FaceAlpha',.5);
patch(get(h(2),'XData'),get(h(2),'YData'),colors('blue'),'FaceAlpha',.5);
ylabel('Erro Medidas (mm3)')

% Erro de idade
figure,plot(ageErr,'s-')
hold on
plot(ageErrRic,'o-')
% plot(ageErrLog,'o-')
grid on
xlabel('Individual')
ylabel('Age Error (days)')
legend('Gompertz','Richards')
% legend('Gompertz','Richards','Logistic')
xlim([1,66])
axes('pos',[.5 .5 .3 .3])
boxplot([ageErr;ageErrRic]','Labels',{'Gompertz','Richards'})
% boxplot([ageErr;ageErrRic;ageErrLog]','Labels',{'Gompertz','Richards','Logistic'})
h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),colors('red'),'FaceAlpha',.5);
patch(get(h(2),'XData'),get(h(2),'YData'),colors('blue'),'FaceAlpha',.5);
ylabel('Age Error (days)')
% saveas(gcf, '../figAgeErr.png')

% Ganho
rmseJeviuGom=sqrt(mse(ageErr));
rmseJeviuRic=sqrt(mse(ageErrRic));
% rmseJeviuLog=sqrt(mse(ageErrLog));
gain=(1-rmseJeviuGom/rmseJeviuRic)*100
% gain2=(1-rmseJeviuRic/rmseJeviuLog)*100

% Comparativo com a tendencia populacional
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
    theta0=exp(mvnrnd(log(hthetaGom(1:3,i)),diag(diag(PeviuGom(1:3,1:3)))));
    yt(j,:)=htGom(theta0,time);
end
ciY=CIFcn(yt,95);
ym=mean(yt);
hold on
p2=plot(time+ageErr(i),ym,'-','color',colors('blue'),'linewidth',2);
fill([time flip(time)]+ageErr(i),[ciY(1,:) flip(ciY(2,:))],colors('blue'), 'FaceAlpha', 0.2,'linestyle','none');
datai=data(data(:,1)==i,:);

xlabel('Time (days)')
ylabel('Volume (mm^3)')

yt=zeros(M,length(time));
for j=1:M
    theta0=exp(mvnrnd(log(hthetaRic(1:4,i)),diag(diag(PeviuRic(1:4,1:4)))));
    yt(j,:)=htRic(theta0,time);
end
ciY=CIFcn(yt,95);
ym=mean(yt);
p3=plot(time+ageErrRic(i),ym,'-','color',colors('red'),'linewidth',2);
fill([time flip(time)]+ageErrRic(i),[ciY(1,:) flip(ciY(2,:))],colors('red'), 'FaceAlpha', 0.2,'linestyle','none');


title(['#' num2str(id(ij))])
xlim([0 30])
if i==23
xlim([0 40])
end
pd=plot(datai(:,2),datai(:,3),'k.','MarkerSize',10);
end
legend([p2,p3,pd],'Gompertz','Richards','Individual Data','Orientation','horizontal')
% legend([p2,p3,p4,pd],'Gompertz','Richards','Logistic','Individual Data','Orientation','horizontal')
% saveas(gcf, '../figComparativo.png')


