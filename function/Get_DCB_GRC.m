function [ DCB,AZ,ELV,LAT_IPP,LON_IPP,LAT_LON_STA,IPPZ] = Get_DCB_GRC(doy ,Sites_Info,sate,SDCB_REF,order)

global stationname c

    Coor=Sites_Info.coor;
    stations=Sites_Info.name;
    doys=Sites_Info.doy;
    n_r=1;                     %length(list_obs);%the number of receivers
        
    

     %% GPS       
    Gx=sate.gpsx;Gy=sate.gpsy;Gz=sate.gpsz;
    [i1,j1]=find(isnan(Gx));[i2,j2]=find(isnan(Gy));[i3,j3]=find(isnan(Gz));
    Gx(i1,j1)=0;Gy(i2,j2)=0;Gz(i3,j3)=0;
    
    path_GPS=['M_P4/GPS/' (doy)];
    list_obs_GPS=[stationname doy 'P4.mat'];
    load([path_GPS '/' list_obs_GPS]);
    Gsat=32;%size(GPSP4,2);
    
    %--check the number of each satellite's observations
    PRN=linspace(0,0,Gsat);
    DCB.GPS_dcb_s=linspace(0,0,Gsat);
    for i=1:n_r
        for j=1:Gsat
            for k=1:size(GPSP4,1)%2880
                if GPSP4(k,j)~=0
                    PRN(j)=PRN(j)+1;
                end
            end
        end
        clear GPSP4
    end
    
    d_sat=find(PRN==0);
    
    if isempty(d_sat)
        n_s=Gsat;
    else
        n_s=Gsat-length(d_sat);%the number of satellites
        display([doy ' PRN ',num2str(d_sat) ,' have no observations.']);
    end
    
    %--chose the order of spheric harmonic function
    %order=str2double(input('Please input the order of spheric harmonic function (4 order is recommended):','s'));
    %%--LS estimate
    
    
    B=[];l=[];
    C=linspace(0,0,(order+1)^2*12+n_s+n_r);
    C(n_r+1:n_r+n_s)=ones(1,n_s);
    if n_s==Gsat
        Wx=0;
    else
        %Satellites DCB values must be exsist in related ionox files
%         idd=find(d_sat>size(SDCB_REF.value.gps,2));
%         d_sat(idd)=[];        
        index= SDCB_REF.doy==(str2double(doy));
        Wx=-sum(SDCB_REF.value.gps(index,d_sat));
        DCB.GPS_dcb_s(d_sat)=SDCB_REF.value.gps(index,d_sat);
    end
    for i=1:n_r
        load([path_GPS '/' list_obs_GPS]);
        site=stationname;
        indices=doys==(str2double(doy));
        index=find(strcmpi(site,stations(indices)), 1);
        sx=Coor(index,1);
        sy=Coor(index,2);
        sz=Coor(index,3);
        ith=i;
        [sN,sl,Az_gps,elv_gps,lat_IPP_gps,lon_IPP_gps,lat_sta,lon_sta,IPPZ_gps]=Get_Matrix(GPSP4,Gx,Gy,Gz,sx,sy,sz,n_r,ith,order,Gsat,'G');
        B=[B;sN];
        l=[l;sl];
        AZ.gps=Az_gps;
        ELV.gps=elv_gps;
        LAT_IPP.lat_IPP_gps=lat_IPP_gps;
        LON_IPP.lon_IPP_gps=lon_IPP_gps;
        LAT_LON_STA.lat_sta=lat_sta;
        LAT_LON_STA.lon_sta=lon_sta;
        IPPZ.IPPZ_gps=IPPZ_gps;
        
        clear GPSP4;
    end
    
    if ~isempty(d_sat)
        B(:,d_sat+n_r)=[];
    end
    BB=[B;C];
    L=[l;Wx];
    R=BB\L;
    DCB.GPS_dcb_r=R(1:n_r)*10^9/c;
    temp=linspace(1,Gsat,Gsat);
    temp(d_sat)=[];
    DCB.GPS_dcb_s(temp)=R(n_r+1:n_r+n_s)*10^9/c;
%     IONC=R(n_r+n_s+1:end);
    
   clear temp R BB L B l ith index temp sx sy sz DCB_S DCB_R Wx C n_s Gsat PRN d_sat i j k i1 i2 i3 j1 j2 j3 Gx Gy Gz
   clear sN sl Az_gps elv_gps lat_IPP_gps  lon_IPP_gps lat_sta lon_sta IPPZ_gps idd
    %% GLONASS
    
    Rx=sate.glox;Ry=sate.gloy;Rz=sate.gloz;
    [i1,j1]=find(isnan(Rx));[i2,j2]=find(isnan(Ry));[i3,j3]=find(isnan(Rz));
    Rx(i1,j1)=0;Ry(i2,j2)=0;Rz(i3,j3)=0;
    
    path_GLO=['M_P4/GLO/' (doy)];
    list_obs_GLO=[stationname doy 'P4.mat'];
    load([path_GLO '/' list_obs_GLO]);
    Rsat=24;%size(GLOP4,2);
    %--check the number of each satellite's observations
    PRN=linspace(0,0,Rsat);
    DCB.GLO_dcb_s=linspace(0,0,Rsat);
    for i=1:n_r
        for j=1:Rsat
            for k=1:size(GLOP4,1)%2880
                if GLOP4(k,j)~=0
                    PRN(j)=PRN(j)+1;
                end
            end
        end
        clear GLOP4
    end
    d_sat=find(PRN==0);
    if isempty(d_sat)
        n_s=Rsat;
    else
        n_s=Rsat-length(d_sat);%the number of satellites
        display([doy ' PRN ',num2str(d_sat) ,' have no observations.']);
    end
    
    %--chose the order of spheric harmonic function
    %order=str2double(input('Please input the order of spheric harmonic function (4 order is recommended):','s'));
    %--LS estimate
    B=[];l=[];
    C=linspace(0,0,(order+1)^2*12+n_s+n_r);
    C(n_r+1:n_r+n_s)=ones(1,n_s);
    if n_s==Rsat
        Wx=0;
    else
        %Satellites DCB values must be exsist in related ionox files
%         idd=find(d_sat>size(SDCB_REF.value.glo,2));
%         d_sat(idd)=[];
        index= SDCB_REF.doy==(str2double(doy));
        Wx=sum(SDCB_REF.value.glo(index,d_sat));
        DCB.GLO_dcb_s(d_sat)=SDCB_REF.value.glo(index,d_sat);
    end
    for i=1:n_r
        load([path_GLO '/' list_obs_GLO]);
        site=stationname;
        indices=doys==(str2double(doy));
        index=find(strcmpi(site,stations(indices)), 1);
        sx=Coor(index,1);
        sy=Coor(index,2);
        sz=Coor(index,3);
        ith=i;
        [sN,sl,Az_glo,elv_glo,lat_IPP_glo,lon_IPP_glo,lat_sta,lon_sta,IPPZ_glo]=Get_GLOMatrix(GLOP4,Rx,Ry,Rz,sx,sy,sz,n_r,ith,order,Rsat);
        B=[B;sN];
        l=[l;sl];
        AZ.glo=Az_glo;
        ELV.glo=elv_glo;
        LAT_IPP.lat_IPP_glo=lat_IPP_glo;
        LON_IPP.lon_IPP_glo=lon_IPP_glo;
        LAT_LON_STA.lat_sta=lat_sta;
        LAT_LON_STA.lon_sta=lon_sta;
        IPPZ.IPPZ_glo=IPPZ_glo;
        clear GLOP4;
    end
    
    if ~isempty(d_sat)
        B(:,d_sat+n_r)=[];
    end
    BB=[B;C];
    L=[l;Wx];
    R=BB\L;
    DCB.GLO_dcb_r=R(1:n_r)*10^9/c;
    temp=linspace(1,Rsat,Rsat);
    temp(d_sat)=[];
    DCB.GLO_dcb_s(temp)=R(n_r+1:n_r+n_s)*10^9/c;
%     IONC=R(n_r+n_s+1:end);
    
       clear temp R BB L B l ith index temp sx sy sz DCB_S DCB_R Wx C n_s Rsat PRN d_sat i j k i1 i2 i3 j1 j2 j3 Rx Ry Rz
       clear sN sl Az_glo elv_glo lat_IPP_glo lon_IPP_glo lat_sta lon_sta IPPZ_glo idd


    %% BEIDOU

    
    Cx=sate.bdsx;Cy=sate.bdsy;Cz=sate.bdsz;
    [i1,j1]=find(isnan(Cx));[i2,j2]=find(isnan(Cy));[i3,j3]=find(isnan(Cz));
    Cx(i1,j1)=0;Cy(i2,j2)=0;Cz(i3,j3)=0;
    
    path_BDS=['M_P4/BDS/' (doy)];
    list_obs_BDS=[stationname doy 'P4.mat'];
    load([path_BDS '/' list_obs_BDS]);
    Csat=46;%size(BDSP4,2);
    %--check the number of each satellite's observations
    PRN=linspace(0,0,Csat);
    DCB.BDS_dcb_s=linspace(0,0,Csat);
    for i=1:n_r
        for j=1:Csat
            for k=1:size(BDSP4,1)%2880
                if BDSP4(k,j)~=0
                    PRN(j)=PRN(j)+1;
                end
            end
        end
        clear BDSP4
    end
    d_sat=find(PRN==0);
    if isempty(d_sat)
        n_s=Csat;
    else
        n_s=Csat-length(d_sat);%the number of satellites
        display([doy ' PRN ',num2str(d_sat) ,' have no observations.']);
    end
    
    %--chose the order of spheric harmonic function
    %order=str2double(input('Please input the order of spheric harmonic function (4 order is recommended):','s'));
    %--LS estimate
    B=[];l=[];
    C=linspace(0,0,(order+1)^2*12+n_s+n_r);
    C(n_r+1:n_r+n_s)=ones(1,n_s);
    if n_s==Csat
        Wx=0;
    else
        %Satellites DCB values must be exsist in related ionox files
%         idd=find(d_sat>size(SDCB_REF.value.bds,2));
%         d_sat(idd)=[];
        index= SDCB_REF.doy==(str2double(doy));
        Wx=-sum(SDCB_REF.value.bds(index,d_sat));
        DCB.BDS_dcb_s(d_sat)=SDCB_REF.value.bds(index,d_sat);
    end
    for i=1:n_r
        load([path_BDS '/' list_obs_BDS]);
        site=stationname;
        indices=doys==(str2double(doy));
        index=find(strcmpi(site,stations(indices)), 1);
        sx=Coor(index,1);
        sy=Coor(index,2);
        sz=Coor(index,3);
        ith=i;
        [sN,sl,Az_bds,elv_bds,lat_IPP_bds,lon_IPP_bds,lat_sta,lon_sta,IPPZ_bds]=Get_Matrix(BDSP4,Cx,Cy,Cz,sx,sy,sz,n_r,ith,order,Csat,'C');
        B=[B;sN];
        l=[l;sl];
        AZ.bds=Az_bds;
        ELV.bds=elv_bds;
        LAT_IPP.lat_IPP_bds=lat_IPP_bds;
        LON_IPP.lon_IPP_bds=lon_IPP_bds;
        LAT_LON_STA.lat_sta=lat_sta;
        LAT_LON_STA.lon_sta=lon_sta;
        IPPZ.IPPZ_bds=IPPZ_bds;
        clear BDSP4;
    end
    
    if ~isempty(d_sat)
        B(:,d_sat+n_r)=[];
    end
    BB=[B;C];
    L=[l;Wx];
    R=BB\L;
    DCB.BDS_dcb_r=R(1:n_r)*10^9/c;
    temp=linspace(1,Csat,Csat);
    temp(d_sat)=[];
    DCB.BDS_dcb_s(temp)=R(n_r+1:n_r+n_s)*10^9/c;
%     IONC=R(n_r+n_s+1:end);
    
       clear temp R BB L B l ith index temp sx sy sz DCB_S DCB_R Wx C n_s Csat PRN d_sat i j k i1 i2 i3 j1 j2 j3 Cx Cy Cz
       clear sN sl Az_bds elv_bds lat_IPP_bds lon_IPP_bds lat_sta lon_sta IPPZ_bds idd
    
end

%------------------------------sub_function--------------------------------
function [M,l,Az_dcb,elv_dcb,lat_IPP,lon_IPP,sb,sl,IPPZ]=Get_Matrix(P4,x,y,z,sx,sy,sz,n_r,ith,order,nSatTot,cons)

% nSatTot=size(x,2);
M=[];
%P=[];
l=[];
[sb,sl]=XYZtoBLH(sx,sy,sz);
Az_dcb=zeros(size(P4,1),size(P4,2));
elv_dcb=zeros(size(P4,1),size(P4,2));
lat_IPP=zeros(size(P4,1),size(P4,2));
lon_IPP=zeros(size(P4,1),size(P4,2));
IPPZ=zeros(size(P4,1),size(P4,2));

if     cons=='G'
    kfix=(-9.52437);

elseif cons=='C'
    kfix= (-8.99768938);
end


for i=1:12
    for j=1:nSatTot                %-----------------------j is satellite number
        for k=240*i-239:240*i %-------------------------k is epoch number
            if P4(k,j)==0
                continue;
            end
            M_col=linspace(0,0,(order+1)^2*12+nSatTot+n_r);
            [E,A]=Get_EA(sx,sy,sz,x(k,j),y(k,j),z(k,j));
%             IPPz=asin(6378137*sin(0.9782*(pi/2-E))/(6378137+506700));%-----MSLM
            IPPz=asin(6371000*sin(pi/2-E)/(6371000+450000));%-----SLM 
            t_r=30*(k-1)*pi/43200;
            [b,s,s1]=Get_IPP(E,A,sb,sl,IPPz,t_r);
            M_col(ith)=kfix*cos(IPPz);   %-----station dcb coefficient
            M_col(n_r+j)=kfix*cos(IPPz); %---satallite dcb coefficient
            st=(order+1)^2*(i-1)+n_r+(nSatTot+1);
            ed=(order+1)^2*i+n_r+nSatTot;
            M_col(st:ed)=Get_coef(b,s,order);
            M_scol=sparse(M_col);
            M=[M;M_scol];
            l=[l;P4(k,j)*kfix*cos(IPPz)];
            elv_dcb(k,j)= E;
            Az_dcb(k,j) = A;
            lat_IPP(k,j)= b;
            lon_IPP(k,j)= s1;
            IPPZ(k,j)=IPPz;
        end
    end
end
end
function [M,l,Az_dcb,elv_dcb,lat_IPP,lon_IPP,sb,sl,IPPZ]=Get_GLOMatrix(P4,x,y,z,sx,sy,sz,n_r,ith,order,nSatTot)

% nSatTot=size(x,2);
M=[];
%P=[];
l=[];
[sb,sl]=XYZtoBLH(sx,sy,sz);
Az_dcb=zeros(size(P4,1),size(P4,2));
elv_dcb=zeros(size(P4,1),size(P4,2));
lat_IPP=zeros(size(P4,1),size(P4,2));
lon_IPP=zeros(size(P4,1),size(P4,2));
IPPZ=zeros(size(P4,1),size(P4,2));
R=[-9.76307424,-9.72883589,-9.79050823,-9.76307424,-9.79050823,...
      -9.79737274,-9.74252401,-9.75622176,-9.74937168,-9.74252401,...
      -9.70832174,-9.75622176,-9.74937168,-9.78364612,-9.73567875,...
      -9.77678642,-9.76992913,-9.78364612,-9.73567878,-9.77678642,...
      -9.76992913,-9.76992913,-9.76992913,-9.76992913];

for i=1:12
    for j=1:nSatTot                %-----------------------j is satellite number
        for k=240*i-239:240*i %-------------------------k is epoch number
            if P4(k,j)==0
                continue;
            end
            M_col=linspace(0,0,(order+1)^2*12+nSatTot+n_r);
            [E,A]=Get_EA(sx,sy,sz,x(k,j),y(k,j),z(k,j));
%             IPPz=asin(6378137*sin(0.9782*(pi/2-E))/(6378137+506700));%-----MSLM
            IPPz=asin(6371000*sin(pi/2-E)/(6371000+450000));%-----SLM 
            t_r=30*(k-1)*pi/43200;
            [b,s,s1]=Get_IPP(E,A,sb,sl,IPPz,t_r);
            M_col(ith)=R(j)*cos(IPPz);   %-----station dcb coefficient
            M_col(n_r+j)=R(j)*cos(IPPz); %---satallite dcb coefficient
            st=(order+1)^2*(i-1)+n_r+(nSatTot+1);
            ed=(order+1)^2*i+n_r+nSatTot;
            M_col(st:ed)=Get_coef(b,s,order);
            M_scol=sparse(M_col);
            M=[M;M_scol];
            l=[l;P4(k,j)*R(j)*cos(IPPz)];
            elv_dcb(k,j)= E;
            Az_dcb(k,j) = A;
            lat_IPP(k,j)= b;
            lon_IPP(k,j)= s1;
            IPPZ(k,j)=IPPz;
        end
    end
end
end
%------------------------------sub_function--------------------------------
function cof_P=Get_coef(b,s,order)
cof_P=linspace(0,0,(order+1)^2);
% ms=[s,2*s,3*s,4*s,5*s,6*s,7*s,8*s,9*s,10*s,11*s,12*s,13*s,14*s,15*s];
ms=linspace(s,order*s,order);
i=1;
x=sin(b);
for n=0:order
    P=legendre(n,x);
    for m=0:n
        if m==0
            cof_P(i)=P(m+1)*norm(n,m);                    %------------an0
        else
            cof_P(i)=P(m+1)*norm(n,m)*cos(ms(m));         %------------anm
            i=i+1;
            cof_P(i)=P(m+1)*norm(n,m)*sin(ms(m));         %------------bnm
        end
        i=i+1;
    end
end
end
%------------------------------sub_function--------------------------------
function N=norm(n,m)
if m==0
    N=sqrt(factorial(n-m)*(2*n+1)/factorial(n+m));
else
    N=sqrt(factorial(n-m)*(4*n+2)/factorial(n+m));
end
end



