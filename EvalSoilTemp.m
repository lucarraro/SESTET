
atst_tmp=xlsread('utilities\data_aT_sT.xls');

for i=1:8
    ind_t=(i-1)*3+1; ind_a=(i-1)*3+2; ind_s=(i-1)*3+3;
    AirTempMeas_s(:,i)=atst_tmp(:,ind_a); SoilTempMeas(:,i)=atst_tmp(:,ind_s);
end
clear atst_tmp

%% find relationship air-to-soil temperature
RMSE=zeros(1,80); coeff=zeros(80,2);
for lag=1:80
    AirT_MA_long=[]; SoilT_long=[];
    for i=1:8
        subset=find(not(isnan(SoilTempMeas(:,i))));
        AirT_MA = movavg(AirTempMeas_s(subset,i),'simple',lag);
        subset2=find(not(isnan(AirT_MA))); AirT_MA=AirT_MA(subset2);
        SoilT=SoilTempMeas(subset,i); SoilT=SoilT(subset2);
        AirT_MA_long=[AirT_MA_long; AirT_MA];
        SoilT_long=[SoilT_long; SoilT];
    end
    mdl=fitlm(AirT_MA_long,SoilT_long);
    RMSE(lag)=mdl.RMSE;
    coeff(lag,:)=mdl.Coefficients{:,1};
end
BestLag=find(RMSE==min(RMSE));
BestCoeff=coeff(BestLag,:);

%% apply soil temperature model

tmp=movavg(AirTemp','simple',BestLag);
% for i=1:BestLag-1
%     tmp(:,i)=mean(AirTemp(:,1:i),2);
% end
SoilTemp=BestCoeff(1)+tmp'*BestCoeff(2);

clear tmp coeff RMSE

%%
tmp=movavg(AirTempMeas_s,'simple',BestLag);
for i=1:BestLag-1
    tmp(i,:)=mean(AirTempMeas_s(1:i,:),1);
end
SoilTempMod=BestCoeff(1)+tmp*BestCoeff(2);

bin=[]; quants=[]; midpoints=[];
for i=1:8
    ind_subset=find(not(isnan(SoilTempMeas(:,i))));
    bin(i,:)=quantile(SoilTempMeas(ind_subset,i),[0:0.05:1]);
    midpoints(i,:)=[diff(bin(i,:))/2+bin(i,1:end-1)];
    for j=1:20
        low=bin(i,j); up=bin(i,j+1);
        tmp=find(and(SoilTempMeas(:,i)>=low,SoilTempMeas(:,i)<up));
        subset=SoilTempMod(tmp,i);
        tmp=quantile(subset,[0.025 0.25 0.75 0.975]);
        for k=1:4; quants(i,j,k)=tmp(k); end
    end
end

if show_fig
    figure('units','centimeters','position',[0 0 18 10])
    for i=1:8
        subplot(2,4,i);
        RMSE_stat(i)=sqrt(nanmean((SoilTempMeas(53:end,i)-SoilTempMod(53:end,i)).^2));
        hold on;
        patch([midpoints(i,:) flip(midpoints(i,:))],[quants(i,:,1) flip(quants(i,:,4))],[0.75 0.75 0.75],'edgecolor','none')
        patch([midpoints(i,:) flip(midpoints(i,:))],[quants(i,:,2) flip(quants(i,:,3))],[0.5 0.5 0.5],'edgecolor','none')
        hold on; plot([0 23],[0 23],'k'); title(['RMSE = ',num2str(RMSE_stat(i))])
        set(gca,'tickdir','out')
        axis([0 21 0 21]); box off
    end
end
RMSE_stat=[];