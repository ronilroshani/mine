function[]=Modelling_VTEC_Bspline_shapefile(S_path,tend)

addpath(S_path);
year=(S_path(:,end-8:end-5));
doy=(S_path(1,end-3:end-1));
[dateV] = doy2date(str2double(doy),str2double(year));
dateVV=datevec(dateV);
data_dir= dir(S_path);
num_file=length(data_dir)-2;

if tend==60
    Hend=24;
elseif tend==45
    Hend=64;
elseif tend==30
    Hend=96;
elseif tend==15
    Hend=192;
end

for H =1:Hend
    t_int =tend/60;
    t1=(H-1)*(t_int*3600/30)+1;
    
    if H < Hend
        t2=(H+1)*(t_int*3600/30);
    else
        t2=(H)*(t_int*3600/30);
    end
    
    lat_pole=80.65;
    lon_pole=-72.68;
    
    disp(['Hour: ' num2str((t1*30/3600)) '-' num2str((t2*30/3600)) ' UT'])
    %% Read results
    VTEC_full=[];lat_lon_full=[];ELV=[];lat_lon_geo=[];
    for i=1:num_file
        
        if data_dir(i+2).name(1:3) == 'loc'
            sta_name=data_dir(i+2).name(10:13);
            load(data_dir(i+2).name);
            Lat_IPP_sta=location.lat_IPP(t1:t2,:)*180/pi;
            Lon_IPP_sta=location.lon_IPP(t1:t2,:)*180/pi;
            
            if size(Lat_IPP_sta,2)<102  | size(Lon_IPP_sta,2)<102
                s = size(Lat_IPP_sta,2);
                P1 = cat(2,Lat_IPP_sta,zeros(size(Lat_IPP_sta,1),102-s));
                P2 = cat(2,Lon_IPP_sta,zeros(size(Lon_IPP_sta,1),102-s));    
                
                Lat_IPP_sta=P1;
                Lon_IPP_sta=P2;
                
            end
            
            
            
            Elv_sta=location.ELV(t1:t2,:);
            
            Lat_IPP_sta_sunfix=asin(sind(Lat_IPP_sta)*sind(lat_pole)+...
                cosd(Lat_IPP_sta)*cosd(lat_pole).*cosd(Lon_IPP_sta-lon_pole));
            Lat_IPP_sta_sunfix=Lat_IPP_sta_sunfix*180/pi;
            
            LL_sta(:,1)=location.lat_sta*180/pi;
            LL_sta(:,2)=location.lon_sta*180/pi;
            
            for j=1:num_file
                if length(data_dir(j+2).name)==24 & data_dir(j+2).name(1:4) == 'VTEC'
                    if data_dir(j+2).name(6:9) == sta_name
                        load(data_dir(j+2).name);
                        VTEC_sta=VTEC(t1:t2,:);
                    end
                end
            end
            
            %             if size(location.lat_IPP,2)<93 && size(VTEC,2)<93
            %                 continue;
            %             end
            
            VTEC_sta=movmean(VTEC_sta,80);
            
            [LAT_TKG,LON_TKG]=read_shp_iran_tehran_karaj;
            LAT_TK=[];LON_TK=[];
            kk=1;
            for kt=t1:t2
                t_r=30*(kt-1)*pi/43200;
                LAT_TK1=asin(sind(LAT_TKG)*sind(lat_pole)+...
                    cosd(LAT_TKG)*cosd(lat_pole).*cosd(LON_TKG-lon_pole));
                LAT_TK1=LAT_TK1*180/pi;
                
                LON_TK1=LON_TKG+t_r*180/pi -180;
                Lon_IPP_sta_sunfix(kk,:)=Lon_IPP_sta(kk,:)+t_r*180/pi -180;
                
                VTEC_sta1=reshape(VTEC_sta(kk,:),[],1);
                id_nan=find(isnan(VTEC_sta1)==1);
                VTEC_sta1(id_nan)=[];
                Elv_sta1=reshape(Elv_sta(kk,:),[],1);
                Elv_sta1(id_nan)=[];
                Lat_IPP_sta_sunfix1=reshape(Lat_IPP_sta_sunfix(kk,:),[],1);
                Lon_IPP_sta_sunfix1=reshape(Lon_IPP_sta_sunfix(kk,:),[],1);
                Lat_IPP_sta_sunfix1(id_nan)=[];
                Lon_IPP_sta_sunfix1(id_nan)=[];
                Lat_IPP_sta1=reshape(Lat_IPP_sta(kk,:),[],1);
                Lon_IPP_sta1=reshape(Lon_IPP_sta(kk,:),[],1);
                Lat_IPP_sta1(id_nan)=[];
                Lon_IPP_sta1(id_nan)=[];
                
                
                id_nan=[];
                id_nan=find((Elv_sta1)<15*pi/180);
                VTEC_sta1(id_nan)=[];
                Elv_sta1(id_nan)=[];
                Lat_IPP_sta_sunfix1(id_nan)=[];
                Lon_IPP_sta_sunfix1(id_nan)=[];
                Lat_IPP_sta1(id_nan)=[];
                Lon_IPP_sta1(id_nan)=[];
                
                
                [IN,ON] = inpolygon(Lat_IPP_sta_sunfix1,Lon_IPP_sta_sunfix1,LAT_TK1,LON_TK1) ;
                Lat_IPP_sta_sunfix1=Lat_IPP_sta_sunfix1(IN);
                Lon_IPP_sta_sunfix1=Lon_IPP_sta_sunfix1(IN);
                Lat_IPP_sta1=Lat_IPP_sta1(IN);
                Lon_IPP_sta1=Lon_IPP_sta1(IN);
                VTEC_sta1=VTEC_sta1(IN);
                Elv_sta1=Elv_sta1(IN);
                VTEC_full=[VTEC_full;VTEC_sta1];
                lat_lon_full=[lat_lon_full;Lat_IPP_sta_sunfix1 Lon_IPP_sta_sunfix1];
                lat_lon_geo =[lat_lon_geo;Lat_IPP_sta1 Lon_IPP_sta1];
                ELV=[ELV;1./sin(Elv_sta1).^2];
                
                kk=kk+1;
            end
            
            clear VTEC Elv_sta1 location sta_name VTEC_sta1 Lat_IPP_sta Lon_IPP_sta TF id_out Lat_IPP_sta_sunfix1 Lon_IPP_sta_sunfix1;
        end
        
        
    end
    if size(VTEC_full,1)<=1
        continue
    end
    %% show IPP
    figure;
    worldmap([34.5 36.5],[50 53.5]);
    %     worldmap ([min(LAT_TKG)-2 max(LAT_TKG)+2],[min(LON_TKG)-2 max(LON_TKG)+2]);
    plotm(LAT_TKG,LON_TKG,'.','Color',[51 211 240]/256);hold on;
    scatterm(lat_lon_geo(:,1),lat_lon_geo(:,2),2,'filled');
    title(['IPP [' num2str((t1*30/3600)) '-' num2str((t2*30/3600)) ' UT] ' '   date : ' num2str(dateVV(1,1)) '/' num2str(dateVV(1,2)) '/' num2str(dateVV(1,3))],'FontSize', 10,'FontAngle','italic','FontName','Times New Roman')
    name_IPP=['IPP[' num2str((t1*30/3600)) '-' num2str((t2*30/3600)) 'UT]' '-' num2str(dateVV(1,1)) '-' num2str(dateVV(1,2)) '-' num2str(dateVV(1,3))];
    cd .\plot\IPP;
    mkdir([year])
    cd(year);
    mkdir(doy);
    cd (doy);
    export_fig(name_IPP,'-jpg','-r600');
    cd .\..\..\..\..;
    close all;
    %% Bspline
    max_lat=max(lat_lon_full(:,1));
    min_lat=min(lat_lon_full(:,1));
    max_lon=max(lat_lon_full(:,2));
    min_lon=min(lat_lon_full(:,2));
    
    j =1;
    normal_LAT=(lat_lon_full(:,1)-min_lat)/(max_lat-min_lat);
    normal_LON=(lat_lon_full(:,2)-min_lon)/(max_lon-min_lon);
    
    [N_LAT]=BS(normal_LAT,j);
    [N_LON]=BS(normal_LON,j);
    %     shift=ones(size(N_LAT,1),1);
    for i=1:size(N_LAT,1)
        A(i,:)=kron(N_LAT(i,:),N_LON(i,:));
    end
    
    X_had_tikh=regular(A,VTEC_full,'GCV','tikh');
    
    %% ---------
    % [U,S,V]=svd(A);
    % s=diag(S);
    % [reg_corner,rho,eta,reg_param] = l_curve(U,s,VTEC_full,'tikh');
    % l=reg_corner;
    % X_HAT=inv((A'*A)+(l.*(eye(100,100))))*A'*VTEC_full;
    
    %% make map
    
    t_r=H*15;
    
    LAT_TKM=LAT_TKG;
    LON_TKM=LON_TKG;
    
    lat_map=[min(LAT_TKM):0.01:max(LAT_TKM)]';
    lon_map=[min(LON_TKM):0.01:max(LON_TKM)];
    lat_map=repmat(lat_map,1,length(lon_map));
    lon_map=repmat(lon_map,size(lat_map,1),1);
    
    lat_T=reshape(lat_map,[],1);
    lon_T=reshape(lon_map,[],1);
    
    max_lat_T=max(lat_T);
    min_lat_T=min(lat_T);
    max_lon_T=max(lon_T);
    min_lon_T=min(lon_T);
    
    normal_LAT_T=(lat_T-min_lat_T)/(max_lat_T-min_lat_T);
    normal_LON_T=(lon_T-min_lon_T)/(max_lon_T-min_lon_T);
    
    [N_LAT_T]=BS(normal_LAT_T,j);
    [N_LON_T]=BS(normal_LON_T,j);
    
    for i=1:size(N_LAT_T,1)
        A_T(i,:)=kron(N_LAT_T(i,:),N_LON_T(i,:));
    end
    
    VTEC_CAP=A_T*X_had_tikh;
    VTEC_map=reshape(VTEC_CAP,size(lat_map,1),[]);
    
    figure_geoshow(lon_map, lat_map, VTEC_map,['MAP_ TEC : [' num2str((t1*30/3600)) '-' num2str((t2*30/3600))  ' UT]' '   date : ' num2str(dateVV(1,1)) '/' num2str(dateVV(1,2)) '/' num2str(dateVV(1,3))],'[TECU]');
    name_MAP=['VTECMAP[' num2str((t1*30/3600)) '-' num2str((t2*30/3600)) 'UT]-' num2str(dateVV(1,1)) '-' num2str(dateVV(1,2)) '-' num2str(dateVV(1,3))];
    cd .\plot\VTEC_MAP;
    mkdir([year])
    cd(year);
    mkdir(doy);
    cd (doy);
    export_fig(name_MAP,'-jpg','-r600');
    cd .\..\..\..\..;
    close all;
    
    
    Time1(H,:)=H;
    RES_TECMAP(:,:,H)=VTEC_map;
    RES_LAT(:,:,H)=lat_map;
    RES_LON(:,:,H)=lon_map;

   
    
    VTEC_map=[];VTEC_CAP=[];N_LAT_T=[];N_LON_T=[];X_had_tikh=[];[N_LAT]=[];
    [N_LON]=[];A=[];IN=[];VTEC_full=[];lat_lon_full=[];ELV=[];
    
end
% save RES_TECMAP RES_TECMAP;
Time_Map=[];
Tec_Map=[];
Lat_Map=[];
Lon_Map=[];

for i=1:size(RES_TECMAP,3)
    Time2=repmat(Time1(i,:),size(RES_TECMAP,1),1);
    map_tec=RES_TECMAP(:,:,i);
    map_lat=RES_LAT(:,:,i);
    map_lon=RES_LON(:,:,i);
    
    Time_Map=[Time_Map;Time2] ;
    Tec_Map=[Tec_Map;map_tec];
    Lat_Map=[Lat_Map;map_lat];
    Lon_Map=[Lon_Map;map_lon];
end
T = table(Time_Map,Tec_Map,Lat_Map,Lon_Map,'VariableNames',{'UT','TEC','Lat','Lon'});

cd .\Results\Map_text;
mkdir([year])
cd(year);
mkdir(doy);
cd (doy);
writetable(T,[num2str(dateVV(1,1)) '-' num2str(dateVV(1,2)) '-' num2str(dateVV(1,3)) '.txt'],'Delimiter','\t');
cd .\..\..\..\..;

end
