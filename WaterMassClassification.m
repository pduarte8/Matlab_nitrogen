
clear all

%FileName = 'European_Sector.mat';
%Sector = 'EuropeanSector7080N'
%Sector = 'EuropeanSector8090N'

%FileName = 'Siberian_Sector.mat';
%Sector = 'SiberianSector7080N'
%Sector = 'SiberianSector8090N'

%FileName = 'Pacific_Sector.mat';
%Sector = 'PacificSector7080N'
%Sector = 'PacificSector8090N'

FileName = 'Canadian_Sector.mat';
%Sector = 'CanadianSector7080N'
Sector = 'CanadianSector8090N'

load(FileName)

long = S.(Sector).longitude.Value;
lat  = S.(Sector).latitude.Value;
SIGMA0=S.(Sector).SIGMA0.Value;     
SP = S.(Sector).SALNTY.Value;
p  = S.(Sector).PRESSURE.Value;
pt = S.(Sector).THETA.Value;
SA = gsw_SA_from_SP(SP,p,long,lat); %Salinity aboslute
CT = gsw_CT_from_pt(SA,pt);         %Conservative temperature 
SIGMA05 = gsw_rho(SA,CT,500)-1000; % Potential density anomaly referred to surface

%% Water masses flagging
% Atlantic water 27.70<?0<27.97 and ?>2oC
[AW_flag(:,1) AW_flag(:,2)]=find(SIGMA0>27.70 & SIGMA0<27.97 & CT>2);
% Polar Surface Water: ?0<27.70 and ?<0oC
[PSW_flag(:,1) PSW_flag(:,2)]=find(SIGMA0<27.70 & CT<0);
% warm Polar Surface Water: ?0<27.70 and ?>0oC
[PSWw_flag(:,1) PSWw_flag(:,2)]=find(SIGMA0<27.70 & CT>0);
% Modified Atlantic Wwater: 27.70<?0<27.97 and ?<2oC
[MAW_flag(:,1) MAW_flag(:,2)]=find(SIGMA0>27.70 & SIGMA0<27.97 & CT<2);
% Intermediate Water: 27.97<?0, ?0.5 <30.444 and ?<0oC
[IW_flag(:,1) IW_flag(:,2)]=find(SIGMA0>27.97 & SIGMA05<30.444 & CT<0);
% MAW/AW: 27.97<?0, ?0.5 <30.444 and ?>0oC
[MAWAW_flag(:,1) MAWAW_flag(:,2)]=find(SIGMA0>27.97 & SIGMA05<30.444 & CT>0);
% Nordic Deep Water: ?0.5 >30.444
[NDW_flag(:,1) NDW_flag(:,2)]=find(SIGMA05>30.444);

%% Water mass matrix
watermass=NaN(size(SIGMA0));
for i=1:length(PSW_flag); watermass(PSW_flag(i,1),PSW_flag(i,2))=70;
end
for i=1:length(PSWw_flag); watermass(PSWw_flag(i,1),PSWw_flag(i,2))=60;
end
for i=1:length(MAW_flag); watermass(MAW_flag(i,1),MAW_flag(i,2))=50;
end
for i=1:length(AW_flag); watermass(AW_flag(i,1),AW_flag(i,2))=40;
end
for i=1:length(MAWAW_flag); watermass(MAWAW_flag(i,1),MAWAW_flag(i,2))=30;
end
for i=1:length(IW_flag); watermass(IW_flag(i,1),IW_flag(i,2))=20;
end
for i=1:length(NDW_flag);watermass(NDW_flag(i,1),NDW_flag(i,2))=10;
end
S.(Sector).watermass.Value = watermass;
save(FileName,'S','-append')
