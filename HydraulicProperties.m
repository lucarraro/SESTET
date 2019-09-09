Cp=4181.3;  % [J Kg^-1 K^-1]
g=9.81;    % [m s^-2]

N_obs_partial=sum(not(isnan(TempMeas)));
N_obs=sum(N_obs_partial);

%% Find coefficients beta for each subcatchment
beta1=zeros(N_reach);
for i=1:N_reach
    tmp=find(AD(:,i));
    if numel(tmp)>0; for j=1:numel(tmp)
            beta1(i,tmp(j))=area_upstream(tmp(j))/area_upstream(i);
        end; end
end
beta1tot=sum(beta1,2);
beta2=1-beta1tot; beta3=1-0.5*beta2;

%% stage discharge relationships
stage_ZOF=sd_data(:,1)-sd_data(1,1); discharge_ZOF=sd_data(:,2); 
clear ndata text alldata

%% Add discharge
TimeQ=(datenum(2014,01,01):datenum(2018,06,30))';
subset=ismember(TimeQ,TimeAir);
Q_ZOF=Q_ZOF(subset);
d_ZOF=interp1(discharge_ZOF,stage_ZOF,Q_ZOF);
Q_outlet=area_upstream(1)/area_upstream(13)*Q_ZOF;
Ks=30; % [m^(1/3) s^-1]
Q_all=(area_upstream(:)/max(area_upstream))*Q_outlet';
dQdt=[zeros(N_reach,1) diff(Q_all,1,2)];

reach_slope(reach_slope==0)=0.001;
u=zeros(N_reach,length(TimeAir));  
for t=1:length(TimeAir)
     for i=1:N_reach  
        u_manning(i,t) = (Ks^1.5*reach_slope(i)^0.75*area_upstream(i)/area_upstream(1)*Q_outlet(t)/reach_width(i))^0.4;
        reach_depth_manning(i,t) = area_upstream(i)/area_upstream(1)*Q_outlet(t)/reach_width(i)/u_manning(i,t);
        reach_volume_manning(i,t)=reach_depth_manning(i,t)*reach_width(i)*length_reach(i);
     end
end

for i=1:N_reach
    one(i,:)=diff(reach_volume_manning(i,:))/86400; two(i,:)=beta2(i)*Q_all(i,2:end);
QL(i,:)=one(i,:)+two(i,:);
end
mean_u=mean(u_manning,2);