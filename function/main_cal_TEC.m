function main_cal_TEC(POSItion,ObsGNSS,Sites_Info,SDCB_REF,F,date_RM,Doy)
global GPS_flag GLO_flag  BDS_flag Cutoff order PG PC PR fig stationname p_path
%% GRC
year  = num2str(date_RM(1,1));
month = num2str(date_RM(1,2),'%.2d');
day  = num2str(date_RM(1,3),'%.2d');

if [GPS_flag,GLO_flag ,BDS_flag]==[1,1,1]

    sys='GPS_GLO_BDS';
    
    [GPSP4 ,GLOP4, BDSP4]= Get_P4GRC(POSItion,ObsGNSS,Sites_Info,Cutoff);
    [ DCB,AZ,ELV,LAT_IPP,LON_IPP,LAT_LON_STA,IPPZ] = Get_DCB_GRC(num2str(SDCB_REF.doy) ,Sites_Info,POSItion,SDCB_REF,order);
    [G_R, ~, C_R, ~, R_R, ~, ~, ~, ~] = Get_nonSH_GCR(fig,num2str(SDCB_REF.doy) ,Sites_Info,POSItion,SDCB_REF,order,order,PG,PC,PR);
    DCB.GPS_dcb_r  =G_R;
    DCB.BDS_dcb_r  =C_R;
    DCB.GLO_dcb_r  =R_R;
    P4_Smooth.gps=GPSP4;
    P4_Smooth.glo=GLOP4;
    P4_Smooth.bds=BDSP4;
    
    % latitude & longtitude IPP
    location.lat_IPP=[LAT_IPP.lat_IPP_gps,LAT_IPP.lat_IPP_glo,LAT_IPP.lat_IPP_bds];
    location.lon_IPP=[LON_IPP.lon_IPP_gps,LON_IPP.lon_IPP_glo,LON_IPP.lon_IPP_bds];
    location.AZ  = [AZ.gps,AZ.glo,AZ.bds];
    location.ELV = [ELV.gps,ELV.glo,ELV.bds];
    location.lat_sta=LAT_LON_STA.lat_sta;
    location.lon_sta=LAT_LON_STA.lon_sta;
    location.IPPZ=[IPPZ.IPPZ_gps,IPPZ.IPPZ_glo,IPPZ.IPPZ_bds];
    
    %% ------------------------------- Compute STEC & VTEC --------------------------------
    [STEC, VTEC, ~, ~, VtecStation]=GnssTimeSeriesTec(ELV, DCB, P4_Smooth,F);
    
    num=3;
    VS=[VtecStation.gps,VtecStation.glo,VtecStation.bds];
    VSS=[STEC.gps,STEC.glo,STEC.bds];
    VSV=[VTEC.gps,VTEC.glo,VTEC.bds];
    
    
    VS(isnan(VS))=0;
    
    for i=1:length(VS)
        VV(i,:)=sum(VS(i,:));
    end
    VVV=VV/num;
    Vtec_GNSS =VVV;
    % movavg(VVV,'linear',100);
    %% --------------------------------- Result TXT file ----------------------------------
    % dateT=datetime(date,'Format','yyyy/MM/d');
    dateT=datenum(date_RM);
    TEC__=Vtec_GNSS(1:length(dateT),1);
    Time = [dateT(:,1)-dateT(1,1)]*24;
    TEC_=[num2str(TEC__,3)];
    T = table(num2str(Time,6),TEC_,'VariableNames',{'Time','TEC_'});
    
    YYDOY=num2str(Sites_Info.doy);
    
    
    cd .\Results\TimeSeries_text;
    mkdir([year '_' YYDOY(3:5)])
    result_f=[year '_' YYDOY(3:5)];
    cd (result_f);
    writetable(T,[stationname year YYDOY(3:5) '_' sys '.txt'],'Delimiter','\t');
    cd ..\..\..;
    
    %% -------------------------------- SAVE FILE -----------------------------------
    S_path = [p_path 'Results\Mat_result\' num2str(year) '\' YYDOY(3:5) '\'];
    Save_Result(date_RM, S_path, stationname,Doy,Sites_Info.coor, VSS, VSV, Vtec_GNSS,location);
    
    %% ----------------------------------- PLOT -------------------------------------
    Plot(year,month,day,VtecStation,YYDOY,sys);
    Plot_TimeSeries_TEC(year,month,day,Vtec_GNSS,YYDOY,sys);
    
    %% GR
elseif [GPS_flag,GLO_flag ,BDS_flag]==[1,1,0]

    sys='GPS_GLO';

    
    [GPSP4 ,GLOP4]= Get_P4GR(POSItion,ObsGNSS,Sites_Info,Cutoff);
    [ DCB,AZ,ELV,LAT_IPP,LON_IPP,LAT_LON_STA,IPPZ] = Get_DCB_GR(num2str(SDCB_REF.doy) ,Sites_Info,POSItion,SDCB_REF,order);
    [G_R, ~, R_R, ~, ~, ~, ~] = Get_nonSH_GR(fig,num2str(SDCB_REF.doy) ,Sites_Info,POSItion,SDCB_REF,order,order,PG,PR);
    DCB.GPS_dcb_r  =G_R;
    DCB.GLO_dcb_r  =R_R;
    P4_Smooth.gps=GPSP4;
    P4_Smooth.glo=GLOP4;
    
    % latitude & longtitude IPP
    location.lat_IPP=[LAT_IPP.lat_IPP_gps,LAT_IPP.lat_IPP_glo];
    location.lon_IPP=[LON_IPP.lon_IPP_gps,LON_IPP.lon_IPP_glo];
    location.AZ  = [AZ.gps,AZ.glo];
    location.ELV = [ELV.gps,ELV.glo];
    location.lat_sta=LAT_LON_STA.lat_sta;
    location.lon_sta=LAT_LON_STA.lon_sta;
    location.IPPZ=[IPPZ.IPPZ_gps,IPPZ.IPPZ_glo];
    %% ------------------------------- Compute STEC & VTEC --------------------------------
    [STEC, VTEC, ~, ~, VtecStation]=GnssTimeSeriesTec(ELV, DCB, P4_Smooth,F);
    
    num=2;
    VS=[VtecStation.gps,VtecStation.glo];
    VSS=[STEC.gps,STEC.glo];
    [VSV]=[VTEC.gps,VTEC.glo];
    
    
    VS(isnan(VS))=0;
    
    for i=1:length(VS)
        VV(i,:)=sum(VS(i,:));
    end
    VVV=VV/num;
    Vtec_GNSS =VVV;
    % movavg(VVV,'linear',100);
    %% --------------------------------- Result TXT file ----------------------------------
    % dateT=datetime(date,'Format','yyyy/MM/d');
    dateT=datenum(date_RM);
    TEC__=Vtec_GNSS(1:length(dateT),1);
    Time = [dateT(:,1)-dateT(1,1)]*24;
    TEC_=[num2str(TEC__,3)];
    T = table(num2str(Time,6),TEC_,'VariableNames',{'Time','TEC_'});
    
    YYDOY=num2str(Sites_Info.doy);
    
    
    cd .\Results\TimeSeries_text;
    mkdir([year '_' YYDOY(3:5)])
    result_f=[year '_' YYDOY(3:5)];
    cd (result_f);
    writetable(T,[stationname year YYDOY(3:5) '_' sys '.txt'],'Delimiter','\t');
    cd ..\..\..;
    
    %% -------------------------------- SAVE FILE -----------------------------------
    S_path = [p_path 'Results\Mat_result\' num2str(year) '\' YYDOY(3:5) '\'];
    Save_Result(date_RM, S_path, stationname,Doy,Sites_Info.coor, VSS, VSV, Vtec_GNSS,location);
    
    %% ----------------------------------- PLOT -------------------------------------
    Plot(year,month,day,VtecStation,YYDOY,sys);
    Plot_TimeSeries_TEC(year,month,day,Vtec_GNSS,YYDOY,sys);
    
    %% GC
elseif [GPS_flag,GLO_flag ,BDS_flag]==[1,0,1]
    
    sys='GPS_BDS';
    
    [GPSP4 , BDSP4]= Get_P4GC(POSItion,ObsGNSS,Sites_Info,Cutoff);
    [ DCB,AZ,ELV,LAT_IPP,LON_IPP,LAT_LON_STA,IPPZ] = Get_DCB_GC(num2str(SDCB_REF.doy) ,Sites_Info,POSItion,SDCB_REF,order);
    [G_R, ~, C_R, ~, ~, ~, ~] = Get_nonSH_GC(fig,num2str(SDCB_REF.doy) ,Sites_Info,POSItion,SDCB_REF,order,order,PG,PC);
    DCB.GPS_dcb_r  =G_R;
    DCB.BDS_dcb_r  =C_R;
    P4_Smooth.gps=GPSP4;
    P4_Smooth.bds=BDSP4;
    
    % latitude & longtitude IPP
    location.lat_IPP=[LAT_IPP.lat_IPP_gps,LAT_IPP.lat_IPP_bds];
    location.lon_IPP=[LON_IPP.lon_IPP_gps,LON_IPP.lon_IPP_bds];
    location.AZ  = [AZ.gps,AZ.bds];
    location.ELV = [ELV.gps,ELV.bds];
    location.lat_sta=LAT_LON_STA.lat_sta;
    location.lon_sta=LAT_LON_STA.lon_sta;
    location.IPPZ=[IPPZ.IPPZ_gps,IPPZ.IPPZ_bds];
    
    
    [STEC, VTEC, ~, ~, VtecStation]=GnssTimeSeriesTec(ELV, DCB, P4_Smooth,F);
    
    num=2;
    VS=[VtecStation.gps,VtecStation.bds];
    VSS=[STEC.gps,STEC.bds];
    VSV=[VTEC.gps,VTEC.bds];
    
    
    
    VS(isnan(VS))=0;
    
    for i=1:length(VS)
        VV(i,:)=sum(VS(i,:));
    end
    VVV=VV/num;
    Vtec_GNSS =VVV;
    % movavg(VVV,'linear',100);
    %% --------------------------------- Result TXT file ----------------------------------
    % dateT=datetime(date,'Format','yyyy/MM/d');
    dateT=datenum(date_RM);
    TEC__=Vtec_GNSS(1:length(dateT),1);
    Time = [dateT(:,1)-dateT(1,1)]*24;
    TEC_=[num2str(TEC__,3)];
    T = table(num2str(Time,6),TEC_,'VariableNames',{'Time','TEC_'});
    
    YYDOY=num2str(Sites_Info.doy);
    
    
    cd .\Results\TimeSeries_text;
    mkdir([year '_' YYDOY(3:5)])
    result_f=[year '_' YYDOY(3:5)];
    cd (result_f);
    writetable(T,[stationname year YYDOY(3:5) '_' sys '.txt'],'Delimiter','\t');
    cd ..\..\..;
    
    %% -------------------------------- SAVE FILE -----------------------------------
    S_path = [p_path 'Results\Mat_result\' num2str(year) '\' YYDOY(3:5) '\'];
    Save_Result(date_RM, S_path, stationname,Doy,Sites_Info.coor, VSS, VSV, Vtec_GNSS,location);
    
    %% ----------------------------------- PLOT -------------------------------------
    Plot(year,month,day,VtecStation,YYDOY,sys);
    Plot_TimeSeries_TEC(year,month,day,Vtec_GNSS,YYDOY,sys);
    %% RC
elseif [GPS_flag,GLO_flag ,BDS_flag]==[0,1,1]
    
    sys ='GLO_BDS';
    
    [GLOP4, BDSP4]= Get_P4RC(POSItion,ObsGNSS,Sites_Info,Cutoff);
    [ DCB,AZ,ELV,LAT_IPP,LON_IPP,LAT_LON_STA,IPPZ] = Get_DCB_RC(num2str(SDCB_REF.doy) ,Sites_Info,POSItion,SDCB_REF,order);
    [C_R, ~, R_R, ~, ~, ~, ~] = Get_nonSH_CR(fig,num2str(SDCB_REF.doy) ,Sites_Info,POSItion,SDCB_REF,order,order,PC,PR);
    DCB.BDS_dcb_r  =C_R;
    DCB.GLO_dcb_r  =R_R;
    P4_Smooth.glo=GLOP4;
    P4_Smooth.bds=BDSP4;
    
    % latitude & longtitude IPP
    location.lat_IPP=[LAT_IPP.lat_IPP_glo,LAT_IPP.lat_IPP_bds];
    location.lon_IPP=[LON_IPP.lon_IPP_glo,LON_IPP.lon_IPP_bds];
    location.AZ  = [AZ.glo,AZ.bds];
    location.ELV = [ELV.glo,ELV.bds];
    location.lat_sta=LAT_LON_STA.lat_sta;
    location.lon_sta=LAT_LON_STA.lon_sta;
    location.IPPZ=[IPPZ.IPPZ_glo,IPPZ.IPPZ_bds];
    
    %% ------------------------------- Compute STEC & VTEC --------------------------------
    [STEC, VTEC, ~, ~, VtecStation]=GnssTimeSeriesTec(ELV, DCB, P4_Smooth,F);
    
    num=3;
    VS=[VtecStation.glo,VtecStation.bds];
    VSS=[STEC.glo,STEC.bds];
    VSV=[VTEC.glo,VTEC.bds];
    
    
    VS(isnan(VS))=0;
    
    for i=1:length(VS)
        VV(i,:)=sum(VS(i,:));
    end
    VVV=VV/num;
    Vtec_GNSS =VVV;
    % movavg(VVV,'linear',100);
    %% --------------------------------- Result TXT file ----------------------------------
    % dateT=datetime(date,'Format','yyyy/MM/d');
    dateT=datenum(date_RM);
    TEC__=Vtec_GNSS(1:length(dateT),1);
    Time = [dateT(:,1)-dateT(1,1)]*24;
    TEC_=[num2str(TEC__,3)];
    T = table(num2str(Time,6),TEC_,'VariableNames',{'Time','TEC_'});
    
    YYDOY=num2str(Sites_Info.doy);
    
    
    cd .\Results\TimeSeries_text;
    mkdir([year '_' YYDOY(3:5)])
    result_f=[year '_' YYDOY(3:5)];
    cd (result_f);
    writetable(T,[stationname year YYDOY(3:5) '_' sys '.txt'],'Delimiter','\t');
    cd ..\..\..;
    
    %% -------------------------------- SAVE FILE -----------------------------------
    S_path = [p_path 'Results\Mat_result\' num2str(year) '\' YYDOY(3:5) '\'];
    Save_Result(date_RM, S_path, stationname,Doy,Sites_Info.coor, VSS, VSV, Vtec_GNSS,location);
    
    %% ----------------------------------- PLOT -------------------------------------
    Plot(year,month,day,VtecStation,YYDOY,sys);
    Plot_TimeSeries_TEC(year,month,day,Vtec_GNSS,YYDOY,sys);
    
    %% G
elseif [GPS_flag,GLO_flag ,BDS_flag]==[1,0,0]
    
    sys='GPS';
    
    [GPSP4]= Get_P4G(POSItion,ObsGNSS,Sites_Info,Cutoff);
    [ DCB,AZ,ELV,LAT_IPP,LON_IPP,LAT_LON_STA,IPPZ] = Get_DCB_G(num2str(SDCB_REF.doy) ,Sites_Info,POSItion,SDCB_REF,order);
    [G_R, ~, ~, ~, ~] = Get_nonSH_G(fig,num2str(SDCB_REF.doy) ,Sites_Info,POSItion,SDCB_REF,order,order,PG);
    DCB.GPS_dcb_r  =G_R;
    P4_Smooth.gps=GPSP4;
    
    % latitude & longtitude IPP
    location.lat_IPP=[LAT_IPP.lat_IPP_gps];
    location.lon_IPP=[LON_IPP.lon_IPP_gps];
    location.AZ  = [AZ.gps];
    location.ELV = [ELV.gps];
    location.lat_sta=LAT_LON_STA.lat_sta;
    location.lon_sta=LAT_LON_STA.lon_sta;
    location.IPPZ=[IPPZ.IPPZ_gps];
    
    %% ------------------------------- Compute STEC & VTEC --------------------------------
    [STEC, VTEC, ~, ~, VtecStation]=GnssTimeSeriesTec(ELV, DCB, P4_Smooth,F);
    num=1;
    VS=[VtecStation.gps];
    VSS=[STEC.gps];
    VSV=[VTEC.gps];
    
    
    VS(isnan(VS))=0;
    
    for i=1:length(VS)
        VV(i,:)=sum(VS(i,:));
    end
    VVV=VV/num;
    Vtec_GNSS =VVV;
    % movavg(VVV,'linear',100);
    %% --------------------------------- Result TXT file ----------------------------------
    % dateT=datetime(date,'Format','yyyy/MM/d');
    dateT=datenum(date_RM);
    TEC__=Vtec_GNSS(1:length(dateT),1);
    Time = [dateT(:,1)-dateT(1,1)]*24;
    TEC_=[num2str(TEC__,3)];
    T = table(num2str(Time,6),TEC_,'VariableNames',{'Time','TEC_'});
    
    YYDOY=num2str(Sites_Info.doy);
    
    
    cd .\Results\TimeSeries_text;
    mkdir([year '_' YYDOY(3:5)]);
    result_f=[year '_' YYDOY(3:5)];
    cd (result_f);
    writetable(T,[stationname year YYDOY(3:5) '_' sys '.txt'],'Delimiter','\t');
    cd ..\..\..;
    
    %% -------------------------------- SAVE FILE -----------------------------------
    S_path = [p_path 'Results\Mat_result\' num2str(year) '\' YYDOY(3:5) '\'];
    Save_Result(date_RM, S_path, stationname,Doy,Sites_Info.coor, VSS, VSV, Vtec_GNSS,location);
    
    %% ----------------------------------- PLOT -------------------------------------
    Plot(year,month,day,VtecStation,YYDOY,sys);
    Plot_TimeSeries_TEC(year,month,day,Vtec_GNSS,YYDOY,sys);
    
    
        %% R
elseif [GPS_flag,GLO_flag ,BDS_flag]==[0,1,0]
    
    sys='GLO';
    
    [GLOP4]= Get_P4R(POSItion,ObsGNSS,Sites_Info,Cutoff);
    [ DCB,AZ,ELV,LAT_IPP,LON_IPP,LAT_LON_STA,IPPZ] = Get_DCB_R(num2str(SDCB_REF.doy) ,Sites_Info,POSItion,SDCB_REF,order);
    [R_R, ~, ~, ~, ~] = Get_nonSH_R(fig,num2str(SDCB_REF.doy) ,Sites_Info,POSItion,SDCB_REF,order,order,PR);
    
    DCB.GLO_dcb_r  =R_R;
    P4_Smooth.glo=GLOP4;
    
    % latitude & longtitude IPP
    location.lat_IPP=[LAT_IPP.lat_IPP_glo];
    location.lon_IPP=[LON_IPP.lon_IPP_glo];
    location.AZ  = [AZ.glo];
    location.ELV = [ELV.glo];
    location.lat_sta=LAT_LON_STA.lat_sta;
    location.lon_sta=LAT_LON_STA.lon_sta;
    location.IPPZ=[IPPZ.IPPZ_glo];
    %% ------------------------------- Compute STEC & VTEC --------------------------------
    [STEC, VTEC, ~, ~, VtecStation]=GnssTimeSeriesTec(ELV, DCB, P4_Smooth,F);
    num=1;
    VS=[VtecStation.glo];
    VSS=[STEC.glo];
    VSV=[VTEC.glo];
    
    
    VS(isnan(VS))=0;
    
    for i=1:length(VS)
        VV(i,:)=sum(VS(i,:));
    end
    VVV=VV/num;
    Vtec_GNSS =VVV;
    % movavg(VVV,'linear',100);
    %% --------------------------------- Result TXT file ----------------------------------
    % dateT=datetime(date,'Format','yyyy/MM/d');
    dateT=datenum(date_RM);
    TEC__=Vtec_GNSS(1:length(dateT),1);
    Time = [dateT(:,1)-dateT(1,1)]*24;
    TEC_=[num2str(TEC__,3)];
    T = table(num2str(Time,6),TEC_,'VariableNames',{'Time','TEC_'});
    
    YYDOY=num2str(Sites_Info.doy);
    
    
    cd .\Results\TimeSeries_text;
    mkdir([year '_' YYDOY(3:5)]);
    result_f=[year '_' YYDOY(3:5)];
    cd (result_f);
    writetable(T,[stationname year YYDOY(3:5) '_' sys '.txt'],'Delimiter','\t');
    cd ..\..\..;
    
    %% -------------------------------- SAVE FILE -----------------------------------
    S_path = [p_path 'Results\Mat_result\' num2str(year) '\' YYDOY(3:5) '\'];
    Save_Result(date_RM, S_path, stationname,Doy,Sites_Info.coor, VSS, VSV, Vtec_GNSS,location);
    
    %% ----------------------------------- PLOT -------------------------------------
    Plot(year,month,day,VtecStation,YYDOY,sys);
    Plot_TimeSeries_TEC(year,month,day,Vtec_GNSS,YYDOY,sys);
        %% C
elseif [GPS_flag,GLO_flag ,BDS_flag]==[0,0,1]
    
    sys='BDS';
    
    [BDSP4]= Get_P4C(POSItion,ObsGNSS,Sites_Info,Cutoff);
    [ DCB,AZ,ELV,LAT_IPP,LON_IPP,LAT_LON_STA,IPPZ] = Get_DCB_C(num2str(SDCB_REF.doy) ,Sites_Info,POSItion,SDCB_REF,order);
    [C_R, ~, ~, ~, ~] = Get_nonSH_C(fig,num2str(SDCB_REF.doy) ,Sites_Info,POSItion,SDCB_REF,order,order,PC);
    DCB.GPS_dcb_r  =C_R;
    P4_Smooth.bds=BDSP4;
    
    % latitude & longtitude IPP
    location.lat_IPP=[LAT_IPP.lat_IPP_bds];
    location.lon_IPP=[LON_IPP.lon_IPP_bds];
    location.AZ  = [AZ.bds];
    location.ELV = [ELV.bds];
    location.lat_sta=LAT_LON_STA.lat_sta;
    location.lon_sta=LAT_LON_STA.lon_sta;
    location.IPPZ=[IPPZ.IPPZ_bds];
    
    %% ------------------------------- Compute STEC & VTEC --------------------------------
    [STEC, VTEC, ~, ~, VtecStation]=GnssTimeSeriesTec(ELV, DCB, P4_Smooth,F);
    num=1;
    VS=[VtecStation.bds];
    VSS=[STEC.bds];
    VSV=[VTEC.bds];
    
    
    VS(isnan(VS))=0;
    
    for i=1:length(VS)
        VV(i,:)=sum(VS(i,:));
    end
    VVV=VV/num;
    Vtec_GNSS =VVV;
    % movavg(VVV,'linear',100);
    %% --------------------------------- Result TXT file ----------------------------------
    % dateT=datetime(date,'Format','yyyy/MM/d');
    dateT=datenum(date_RM);
    TEC__=Vtec_GNSS(1:length(dateT),1);
    Time = [dateT(:,1)-dateT(1,1)]*24;
    TEC_=[num2str(TEC__,3)];
    T = table(num2str(Time,6),TEC_,'VariableNames',{'Time','TEC_'});
    
    YYDOY=num2str(Sites_Info.doy);
    
    
    cd .\Results\TimeSeries_text;
    mkdir([year '_' YYDOY(3:5)]);
    result_f=[year '_' YYDOY(3:5)];
    cd (result_f);
    writetable(T,[stationname year YYDOY(3:5) '_' sys '.txt'],'Delimiter','\t');
    cd ..\..\..;
    
    %% -------------------------------- SAVE FILE -----------------------------------
    S_path = [p_path 'Results\Mat_result\' num2str(year) '\' YYDOY(3:5) '\'];
    Save_Result(date_RM, S_path, stationname,Doy,Sites_Info.coor, VSS, VSV, Vtec_GNSS,location);
    
    %% ----------------------------------- PLOT -------------------------------------
    Plot(year,month,day,VtecStation,YYDOY,sys);
    Plot_TimeSeries_TEC(year,month,day,Vtec_GNSS,YYDOY,sys);
    
    
end
end

