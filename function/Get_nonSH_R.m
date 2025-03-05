function [R_R, R_S, IONC, m0, NN] = Get_nonSH_R(fig,doy ,Sites_Info,sate,SDCB_REF,K,M,PR)
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

%% check the glonass data
Rx=sate.glox;Ry=sate.gloy;Rz=sate.gloz;
[i1,j1]=find(isnan(Rx));[i2,j2]=find(isnan(Ry));[i3,j3]=find(isnan(Rz));
Rx(i1,j1)=0;Ry(i2,j2)=0;Rz(i3,j3)=0;

path_GLO=['M_P4/GLO/' (doy)];
list_obs_GLO=[stationname doy 'P4.mat'];
load([path_GLO '/' list_obs_GLO]);
R_n_r=1;%the number of receivers
Rsat=24;%size(GLOP4,2);

%--check the number of each satellite's observations
R_PRN=linspace(0,0,Rsat);
R_S=linspace(0,0,Rsat);
for i=1:R_n_r
    for j=1:Rsat
        for k=1:2880
            if GLOP4(k,j)~=0
                R_PRN(j)=R_PRN(j)+1;
            end
        end
    end
    clear GLOP4;
end
glo_d_sat=find(R_PRN==0);
if isempty(glo_d_sat)
    R_n_s=Rsat;
else
    R_n_s=Rsat-length(glo_d_sat);%the number of satellites
    disp(['doy ', doy ,' GLO PRN ',num2str(glo_d_sat) ,' have no observations.']);
    %     for k=length(glo_d_sat):-1:1
    %         Rx(:,glo_d_sat(k))=[];Ry(:,glo_d_sat(k))=[];Rz(:,glo_d_sat(k))=[];
    %     end
end
if R_n_s==Rsat
    R_Wx=0;
else
    %Satellites DCB values must be exsist in related ionox files
    index= SDCB_REF.doy==str2double(doy);
    R_Wx=-sum(SDCB_REF.value.glo(index,glo_d_sat));
    R_S(glo_d_sat)=SDCB_REF.value.glo(index,glo_d_sat);
end

%% --chose the order of spheric harmonic function
%--LS estimate
n_m=((K+1)^2-(K-M)*(K-M+1))*fig;
num=n_m+R_n_s+R_n_r;
N=zeros(num,num);
U=zeros(num,1);
L=0; sizel=0;
C_GLO=linspace(0,0,num);
C_GLO(R_n_r+1:R_n_r+R_n_s)=ones(1,R_n_s);

for i=1:R_n_r
    load([path_GLO '/' list_obs_GLO]);
    %     if ~isempty(glo_d_sat)
    %         for k=length(glo_d_sat):-1:1
    %             GLOP4(:,k)=[];
    %         end
    %     end
    site=stationname;
    indices=doys==str2double(doy);
    index=find(strcmpi(site,stations(indices)), 1);
    sx=Coor(index,1);
    sy=Coor(index,2);
    sz=Coor(index,3);
    [sN,sl,Az_glo,elv_glo,lat_IPP_glo,lon_IPP_glo,lat_sta,lon_sta,IPPZ_glo]=Get_GLOMatrix(fig,GLOP4,Rx,Ry,Rz,sx,sy,sz,R_n_r,R_n_s,i,K,M);
    N=N+sN'*sN*PR;
    U=U+sN'*sl*PR;
    %--RMS
    sizel=sizel+length(sl);
    L=L+sl'*sl*PR;
    %-----------
    AZ.glo=Az_glo;
    ELV.glo=elv_glo;
    LAT_IPP.lat_IPP_glo=lat_IPP_glo;
    LON_IPP.lon_IPP_glo=lon_IPP_glo;
    LAT_LON_STA.lat_sta=lat_sta;
    LAT_LON_STA.lon_sta=lon_sta;
    IPPZ.IPPZ_glo=IPPZ_glo;
    clear GLOP4;
    disp(['3.----- [ ',num2str(i),' / ',num2str(R_n_r),' ] ',num2str(i/R_n_r*100),'% GLONASS data has constructed !']);
end

N=N+C_GLO'*C_GLO;
U=U+C_GLO'*R_Wx;
L=L+R_Wx'*R_Wx;
R=pinv(N)*U;


R_R=R(1:R_n_r)*10^9/c;
temp_glo=linspace(1,Rsat,Rsat);
temp_glo(glo_d_sat)=[];
R_S(temp_glo)=R(R_n_r+1:R_n_r+R_n_s)*10^9/c;

IONC=R(R_n_r+R_n_s+1:end);
%--RMS
V=L-R'*U;
f=sizel-num;
m0=sqrt(V/f);
NN=N(R_n_r+R_n_s+1:end,R_n_r+R_n_s+1:end);
%------
end

%% ------------------------------sub_function--------------------------------
function [MC,l,Az_dcb,elv_dcb,lat_IPP,lon_IPP,sb,sl,IPPZ]=Get_GLOMatrix(fig,GLOP4,x,y,z,sx,sy,sz,glo_n_r,glo_n_s,ith,K,M)
global Re h

MC=[];l=[];

Az_dcb=zeros(size(GLOP4,1),size(GLOP4,2));
elv_dcb=zeros(size(GLOP4,1),size(GLOP4,2));
lat_IPP=zeros(size(GLOP4,1),size(GLOP4,2));
lon_IPP=zeros(size(GLOP4,1),size(GLOP4,2));
IPPZ=zeros(size(GLOP4,1),size(GLOP4,2));

R=[-9.76307424,-9.72883589,-9.79050823,-9.76307424,-9.79050823,...
    -9.79737274,-9.74252401,-9.75622176,-9.74937168,-9.74252401,...
    -9.70832174,-9.75622176,-9.74937168,-9.78364612,-9.73567875,...
    -9.77678642,-9.76992913,-9.78364612,-9.73567878,-9.77678642,...
    -9.76992913,-9.76992913,-9.76992913,-9.76992913,-9.76992913];
num=(K+1)^2-(K-M)*(K-M+1);
[sb,sl]=XYZtoBLH(sx,sy,sz);
figt=2880/fig;
for i=1:fig
    for j=1:glo_n_s                %----j is satellite number
        for k=figt*i-(figt-1):figt*i %----k is epoch number
            if GLOP4(k,j)==0
                continue;
            end
            M_col=linspace(0,0,num*fig+glo_n_s+glo_n_r);
            [E,A]=Get_EA(sx,sy,sz,x(k,j),y(k,j),z(k,j)); %----sx,sy,sz station coordinate ; x,y,z satellite coordinate
            IPPz=asin(Re*sin(pi/2-E)/(Re+h));%-----SLM
            t_r=30*(k-1)*pi/43200;
            [b,s,s2]=Get_IPP(E,A,sb,sl,IPPz,t_r);%
            M_col(ith)=R(j)*cos(IPPz);   %----station dcb coefficient
            M_col(glo_n_r+j)=R(j)*cos(IPPz); %----satallite dcb coefficient
            st=num*(i-1)+glo_n_r+glo_n_s+1;
            ed=num*i+glo_n_r+glo_n_s;
            M_col(st:ed)=Get_nonSH(b,s,K,M);
            M_scol=sparse(M_col);
            MC=[MC;M_scol];
            l=[l;GLOP4(k,j)*R(j)*cos(IPPz)];
            elv_dcb(k,j)= E;
            Az_dcb(k,j) = A;
            lat_IPP(k,j)= b;
            lon_IPP(k,j)= s2;
            IPPZ(k,j)=IPPz;
        end
    end
end
end
%% -----------------------------------------------------------------------------
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
