function [STEC, VTEC, Meantec, Stdtec, VtecStation]=GnssTimeSeriesTec(ELV, DCB, P4_Smooth,F,LEN)
global Re h A Cutoff c f1G f2G f1C f2C GPS_flag GLO_flag BDS_flag
%% GPS
if GPS_flag ==1
Elevation_Angles_gps=ELV.gps*180/pi;
nSat_gps=32;%LEN(1,1);
STEC.gps           = nan(2880,nSat_gps); % Vertical Total Electron Content(VTEC)
VTEC.gps           = nan(2880,nSat_gps); % Slant Total Electron Content(STEC)
% DCB Stellite & Reciver
bs = DCB.GPS_dcb_s*c*1e-9;
br = DCB.GPS_dcb_r*c*1e-9;

ELEV.gps=zeros(size(STEC.gps,1),size(STEC.gps,2));
for i=1:size(Elevation_Angles_gps,1)
    for j=1:size(Elevation_Angles_gps,2)
        ELEV.gps(i,j)=Elevation_Angles_gps(i,j);
    end
end

ELEV.gps(ELEV.gps < Cutoff(1,1)) = nan;
M.gps = (1-(Re*cosd(ELEV.gps)./(Re+h)).^2).^(-1/2);
BS = repmat(bs,size(P4_Smooth.gps,1),1);
[id_m]=find(P4_Smooth.gps==0);
P4_Smooth.gps(id_m)=nan;
P4_Smooth.gps=P4_Smooth.gps-br- BS;
% STEC & VTEC
STEC.gps = (P4_Smooth.gps/A)*(f1G^2*f2G^2)/(f2G^2-f1G^2);
STEC.gps(STEC.gps<0) = nan;
STEC.gps=STEC.gps*1e-16;
MM.gps=zeros(size(STEC.gps,1),size(STEC.gps,2));
for i=1:size(M.gps,1)
    for j=1:size(M.gps,2)
        MM.gps(i,j)=M.gps(i,j);
    end
end

VTEC.gps = (STEC.gps)./(MM.gps);
VTEC.gps(VTEC.gps<0) = nan;

% Compute VTEC for Station
timerange = 80;
VtecStation.gps = zeros(size(VTEC.gps,1),1);

for i = 1:size(VTEC.gps,1)
    i1 = max(i-timerange,1);
    i2 = min(i+timerange,size(VtecStation.gps,1));
    vtecs = VTEC.gps(i1:i2,:);
    vtecs = vtecs(:);
    sindelv=sind(ELEV.gps(i1:i2,:));
    sindelv_v=sindelv(:);
    sindelv_v(isnan(vtecs))=[];
    vtecs(isnan(vtecs))=[];
    
    Meantec.gps (i,1)= nanmean(vtecs);
    Stdtec.gps (i,1)= nanstd(vtecs);
    
    id1=find(vtecs<Meantec.gps(i,1)-3*Stdtec.gps(i,1));
    id2=find(vtecs>Meantec.gps(i,1)+3*Stdtec.gps(i,1));
    id=[id1;id2];
    sindelv_v(id,:) = [];
    vtecs(id,:) = [];
    
    VtecStation.gps(i,1) = sum(vtecs.*sindelv_v)/sum(sindelv_v);
    clear id1 id2 id;
end
end
clear id1 id2 id vtecs sindelv_v i i1 i2 id_m BS br bs Elevation_Angles_gps nSat_gps;


%% GLONASS
if GLO_flag==1
Elevation_Angles_glo=ELV.glo*180/pi;
nSat_glo=24;%LEN(3,1);
STEC.glo           = nan(2880,nSat_glo); % Vertical Total Electron Content(VTEC)
VTEC.glo           = nan(2880,nSat_glo); % Slant Total Electron Content(STEC)
% DCB Stellite & Reciver
bs =  DCB.GLO_dcb_s*c*1e-9;
br =  DCB.GLO_dcb_r*c*1e-9;

Elevation_Angles_glo(Elevation_Angles_glo < Cutoff(3,1)) = nan;
M.glo = (1-(Re*cosd(Elevation_Angles_glo)./(Re+h)).^2).^(-1/2);
BS = repmat(bs,size(P4_Smooth.glo,1),1);
[id_m]=find(P4_Smooth.glo==0);
P4_Smooth.glo(id_m)=nan;
P4_Smooth.glo=P4_Smooth.glo-br- BS;
% STEC & VTEC
P4_Smooth_glo=P4_Smooth.glo;
for i=1:size(P4_Smooth.glo,2)

    if size(F,1)==24 || 70
        fi=F(1:24,:);
    elseif size(F,1)==56 || 102
        fi=F(33:56,:);
    end

    STEC.glo(:,i) = (P4_Smooth_glo(:,i)/A)*(fi(i,1)^2*fi(i,2)^2)/(fi(i,2)^2-fi(i,1)^2);
    fi=[];
end
STEC.glo(STEC.glo<0) = nan;
STEC.glo=STEC.glo*1e-16;
VTEC.glo = (STEC.glo)./(M.glo);
VTEC.glo(VTEC.glo<0) = nan;

% Compute VTEC for Station
timerange = 180;
VtecStation.glo = zeros(size(VTEC.glo,1),1);

for i = 1:size(VTEC.glo,1)
    i1 = max(i-timerange,1);
    i2 = min(i+timerange,size(VtecStation.glo,1));
    vtecs = VTEC.glo(i1:i2,:);
    vtecs = vtecs(:);
    sindelv=sind(Elevation_Angles_glo(i1:i2,:));
    sindelv_v=sindelv(:);
    sindelv_v(isnan(vtecs))=[];
    vtecs(isnan(vtecs))=[];
    
    Meantec.glo(i,1) = nanmean(vtecs);
    Stdtec.glo(i,1) = nanstd(vtecs);
    
    id1=find(vtecs<Meantec.glo(i,1)-3*Stdtec.glo(i,1));
    id2=find(vtecs>Meantec.glo(i,1)+3*Stdtec.glo(i,1));
    id=[id1;id2];
    sindelv_v(id,:) = [];
    vtecs(id,:) = [];
    
    VtecStation.glo(i,1) = sum(vtecs.*sindelv_v)/sum(sindelv_v);
    clear id1 id2 id;
end
end
clear id1 id2 id vtecs sindelv_v i i1 i2 id_m BS br bs Elevation_Angles_glo;

%% BEIDOU
if BDS_flag==1
Elevation_Angles_bds=ELV.bds*180/pi;
nSat_bds=46;%LEN(2,1);
STEC.bds           = nan(2880,nSat_bds); % Vertical Total Electron Content(VTEC)
VTEC.bds           = nan(2880,nSat_bds); % Slant Total Electron Content(STEC)
% DCB Stellite & Reciver
bs =  DCB.BDS_dcb_s*c*1e-9;
br =  DCB.BDS_dcb_r*c*1e-9;

% ELEV.bds=zeros(size(STEC.bds,1),size(STEC.bds,2));
for i=1:size(Elevation_Angles_bds,1)
    for j=1:size(Elevation_Angles_bds,2)
        ELEV.bds(i,j)=Elevation_Angles_bds(i,j);
    end
end

ELEV.bds(ELEV.bds < Cutoff(2,1)) = nan;
M.bds = (1-(Re*cosd(ELEV.bds)./(Re+h)).^2).^(-1/2);
BS = repmat(bs,size(P4_Smooth.bds,1),1);
[id_m]=find(P4_Smooth.bds==0);
P4_Smooth.bds(id_m)=nan;
P4_Smooth.bds=P4_Smooth.bds-br- BS;
% STEC & VTEC
STEC.bds = (P4_Smooth.bds/A)*(f1C^2*f2C^2)/(f2C^2-f1C^2);
STEC.bds(STEC.bds<0) = nan;
STEC.bds=STEC.bds*1e-16;


% MM.bds=zeros(size(STEC.bds,1),size(STEC.bds,2));
for i=1:size(M.bds,1)
    for j=1:size(M.bds,2)
        MM.bds(i,j)=M.bds(i,j);
    end
end

VTEC.bds = (STEC.bds)./(MM.bds);
VTEC.bds(VTEC.bds<0) = nan;

% Compute VTEC for Station
timerange = 180;
VtecStation.bds = zeros(size(VTEC.bds,1),1);

for i = 1:size(VTEC.bds,1)
    i1 = max(i-timerange,1);
    i2 = min(i+timerange,size(VtecStation.bds,1));
    vtecs = VTEC.bds(i1:i2,:);
    vtecs = vtecs(:);
    sindelv=sind(ELEV.bds(i1:i2,:));
    sindelv_v=sindelv(:);
    sindelv_v(isnan(vtecs))=[];
    vtecs(isnan(vtecs))=[];
    
    Meantec.bds(i,1) = nanmean(vtecs);
    Stdtec.bds (i,1)= nanstd(vtecs);
    
    id1=find(vtecs<Meantec.bds(i,1)-3*Stdtec.bds(i,1));
    id2=find(vtecs>Meantec.bds(i,1)+3*Stdtec.bds(i,1));
    id=[id1;id2];
    sindelv_v(id,:) = [];
    vtecs(id,:) = [];
    
    VtecStation.bds(i,1) = sum(vtecs.*sindelv_v)/sum(sindelv_v);
    clear id1 id2 id;
end
end
clear id1 id2 id vtecs sindelv_v i i1 i2 id_m BS br bs Elevation_Angles_bds;


end

