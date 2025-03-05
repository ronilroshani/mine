function [G_R, G_S, IONC, m0, NN] = Get_nonSH_G(fig,doy ,Sites_Info,sate,SDCB_REF,K,M,PG)
%%  estimate satellite and receiver DCBs, ionospheric parameters
% INPUT:
%     fig: group number of products
%     doy: year and doy of year
%     Sites_Info: name and coordinate information of the stations
%     sate: precise coordinates of the satellites
%     K,M: order and degree of nonintegral-SH model
%     PG,PC,PR: weight of different system observations
% OUTPUT:
%     G_R, G_S, C_R, C_S, R_R, R_S: estimated receiver and satellite DCBs
%     IONC: ionospheric parameters
%     m0: standard deviation
%     NN: covariance matrices

%% --------------------------------------------------------------------------
global stationname c

Coor=Sites_Info.coor;
stations=Sites_Info.name;
doys=Sites_Info.doy;
%% check the gps data
Gx=sate.gpsx;Gy=sate.gpsy;Gz=sate.gpsz;
[i1,j1]=find(isnan(Gx));[i2,j2]=find(isnan(Gy));[i3,j3]=find(isnan(Gz));
Gx(i1,j1)=0;Gy(i2,j2)=0;Gz(i3,j3)=0;

path_GPS=['M_P4/GPS/' doy];
list_obs_GPS=[stationname doy 'P4.mat'];
load([path_GPS '/' list_obs_GPS]);
G_n_r=1;%length(list_gps);%the number of receivers
Gsat=32;%size(GPSP4,2);
%--check the number of each satellite's observations
G_PRN=linspace(0,0,Gsat);
G_S=linspace(0,0,Gsat);
for i=1:G_n_r
    for j=1:Gsat
        for k=1:2880
            if GPSP4(k,j)~=0
                G_PRN(j)=G_PRN(j)+1;
            end
        end
    end
    clear GPSP4;
end
gps_d_sat=find(G_PRN==0);
if isempty(gps_d_sat)
    G_n_s=Gsat;
else
    G_n_s=Gsat-length(gps_d_sat);%the number of satellites
    disp(['doy ', doy ,' GPS PRN ',num2str(gps_d_sat) ,' have no observations.']);
    %     for k=length(gps_d_sat):-1:1
    %         Gx(:,gps_d_sat(k))=[];Gy(:,gps_d_sat(k))=[];Gz(:,gps_d_sat(k))=[];
    %     end
end

if G_n_s==Gsat
    G_Wx=0;
else
    %Satellites DCB values must be exsist in related ionox files
    index= SDCB_REF.doy==str2double(doy);
    G_Wx=-sum(SDCB_REF.value.gps(index,gps_d_sat));
    G_S(gps_d_sat)=SDCB_REF.value.gps(index,gps_d_sat);
end


%% --chose the order of spheric harmonic function
%--LS estimate
n_m=((K+1)^2-(K-M)*(K-M+1))*fig;
num=n_m+G_n_s+G_n_r;
N=zeros(num,num);
U=zeros(num,1);
L=0; sizel=0;
C_GPS=linspace(0,0,num);

C_GPS(G_n_r+1:G_n_r+G_n_s)=ones(1,G_n_s);

for i=1:G_n_r
    load([path_GPS '/' list_obs_GPS]);
    %     if ~isempty(gps_d_sat)
    %         for k=length(gps_d_sat):-1:1
    %             GPSP4(:,k)=[];
    %         end
    %     end
    site=stationname;
    indices=doys==str2double(doy);
    index=find(strcmpi(site,stations(indices)), 1);
    sx=Coor(index,1);
    sy=Coor(index,2);
    sz=Coor(index,3);
    [sN,sl,Az_gps,elv_gps,lat_IPP_gps,lon_IPP_gps,lat_sta,lon_sta,IPPZ_gps]=Get_GPSMatrix(fig,GPSP4,Gx,Gy,Gz,sx,sy,sz,G_n_r,G_n_s,i,K,M);
    N=N+sN'*sN*PG;
    U=U+sN'*sl*PG;
    %--RMS
    sizel=sizel+length(sl);
    L=L+sl'*sl*PG;
    %-----------
    AZ.gps=Az_gps;
    ELV.gps=elv_gps;
    LAT_IPP.lat_IPP_gps=lat_IPP_gps;
    LON_IPP.lon_IPP_gps=lon_IPP_gps;
    LAT_LON_STA.lat_sta=lat_sta;
    LAT_LON_STA.lon_sta=lon_sta;
    IPPZ.IPPZ_gps=IPPZ_gps;
    clear GPSP4;
    disp(['1.----- [ ',num2str(i),' / ',num2str(G_n_r),' ] ',num2str(i/G_n_r*100),'% GPS data has constructed !']);
end



N=N+C_GPS'*C_GPS;
U=U+C_GPS'*G_Wx;
L=L+G_Wx'*G_Wx;
R=pinv(N)*U;
G_R=R(1:G_n_r)*10^9/c;
temp_gps=linspace(1,Gsat,Gsat);
temp_gps(gps_d_sat)=[];
G_S(temp_gps)=R(G_n_r+1:G_n_r+G_n_s)*10^9/c;

IONC=R(G_n_r+G_n_s+1:end);
%--RMS
V=L-R'*U;
f=sizel-num;
m0=sqrt(V/f);
NN=N(G_n_r+G_n_s+1:end,G_n_r+G_n_s+1:end);
%------
end

%% ------------------------------sub_function--------------------------------
function [MC,l,Az_dcb,elv_dcb,lat_IPP,lon_IPP,sb,sl,IPPZ]=Get_GPSMatrix(fig,GPSP4,x,y,z,sx,sy,sz,gps_n_r,gps_n_s,ith,K,M)
global Re h
MC=[];l=[];


Az_dcb=zeros(size(GPSP4,1),size(GPSP4,2));
elv_dcb=zeros(size(GPSP4,1),size(GPSP4,2));
lat_IPP=zeros(size(GPSP4,1),size(GPSP4,2));
lon_IPP=zeros(size(GPSP4,1),size(GPSP4,2));
IPPZ=zeros(size(GPSP4,1),size(GPSP4,2));

num=(K+1)^2-(K-M)*(K-M+1);
[sb,sl]=XYZtoBLH(sx,sy,sz);
figt=2880/fig;
for i=1:fig
    for j=1:gps_n_s                %----j is satellite number
        for k=figt*i-(figt-1):figt*i %----k is epoch number
            if GPSP4(k,j)==0
                continue;
            end
            M_col=linspace(0,0,num*fig+gps_n_s+gps_n_r);
            [E,A]=Get_EA(sx,sy,sz,x(k,j),y(k,j),z(k,j)); %----sx,sy,sz station coordinate ; x,y,z satellite coordinate
            IPPz=asin(Re*sin(pi/2-E)/(Re+h));%-----SLM
            t_r=30*(k-1)*pi/43200;
            [b,s,s2]=Get_IPP(E,A,sb,sl,IPPz,t_r);%
            M_col(ith)=(-9.52437)*cos(IPPz);   %-----station dcb coefficient
            M_col(gps_n_r+j)=(-9.52437)*cos(IPPz); %----satallite dcb coefficient
            st=num*(i-1)+gps_n_r+gps_n_s+1;
            ed=num*i+gps_n_r+gps_n_s;
            M_col(st:ed)=Get_nonSH(b,s,K,M);
            M_scol=sparse(M_col);
            MC=[MC;M_scol];
            l=[l;GPSP4(k,j)*(-9.52437)*cos(IPPz)];
            elv_dcb(k,j)= E;
            Az_dcb(k,j) = A;
            lat_IPP(k,j)= b;
            lon_IPP(k,j)= s2;
            IPPZ(k,j)=IPPz;
        end
    end
end
end
%% ------------------------------sub_function--------------------------------
function cof_capP=Get_nonSH(b,s,K,M)

cof_capP=linspace(0,0,(K+1)^2-(K-M)*(K-M+1));
ms=linspace(s,M*s,M);
i=1;
MM=M;
x=cos(b);
for k=0:K
    P=legendre(k,x);
    if k<M
        M=k;
    end
    for m=0:M
        if m==0
            cof_capP(i)=P(m+1)*norm(k,m);                    %------------an0
        else
            cof_capP(i)=P(m+1)*norm(k,m)*cos(ms(m));         %------------anm
            i=i+1;
            cof_capP(i)=P(m+1)*norm(k,m)*sin(ms(m));         %------------bnm
        end
        i=i+1;
    end
    M=MM;
end
end

%% ------------------------------sub_function--------------------------------
function N=norm(n,m)
if m==0
    N=sqrt(factorial(n-m)*(2*n+1)/factorial(n+m));
else
    N=sqrt(factorial(n-m)*(4*n+2)/factorial(n+m));
end
end
