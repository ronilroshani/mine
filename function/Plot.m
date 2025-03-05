function Plot(year,month,day,VtecStation,YYDOY,sys)
global stationname GPS_flag GLO_flag BDS_flag
Time_TEC  = (0:30:24*3600-30)/3600;      %   Time rate 30 second

figure1 = figure('Color',[1 1 1]);
axes1 = axes('Parent',figure1);
hold(axes1,'on');


if [GPS_flag,GLO_flag ,BDS_flag]==[1,1,1]
plot(Time_TEC,VtecStation.gps,'-','LineWidth',1,'color','b');
plot(Time_TEC,VtecStation.glo,'-','LineWidth',1,'color','r');
plot(Time_TEC,VtecStation.bds,'-','LineWidth',1,'color','g');
legend('Gps','Glonass','Beidou');
elseif [GPS_flag,GLO_flag ,BDS_flag]==[1,1,0]
plot(Time_TEC,VtecStation.gps,'-','LineWidth',1,'color','b');
plot(Time_TEC,VtecStation.glo,'-','LineWidth',1,'color','r');
legend('Gps','Glonass');
elseif [GPS_flag,GLO_flag ,BDS_flag]==[1,0,1]
plot(Time_TEC,VtecStation.gps,'-','LineWidth',1,'color','b');
plot(Time_TEC,VtecStation.bds,'-','LineWidth',1,'color','g');
legend('Gps','Beidou');
elseif [GPS_flag,GLO_flag ,BDS_flag]==[1,0,0]
plot(Time_TEC,VtecStation.gps,'-','LineWidth',1,'color','b');
legend('Gps');
elseif [GPS_flag,GLO_flag ,BDS_flag]==[0,0,1]
plot(Time_TEC,VtecStation.bds,'-','LineWidth',1,'color','g');
legend('Beidou');
elseif [GPS_flag,GLO_flag ,BDS_flag]==[0,1,0]
plot(Time_TEC,VtecStation.glo,'-','LineWidth',1,'color','r');
legend('Glonass');
elseif [GPS_flag,GLO_flag ,BDS_flag]==[0,1,1]
plot(Time_TEC,VtecStation.glo,'-','LineWidth',1,'color','r');
plot(Time_TEC,VtecStation.bds,'-','LineWidth',1,'color','g');
legend('Glonass','Beidou');
end



set(axes1,'PlotBoxAspectRatio',[12 5 1],'XTick',...
    [0 2 4 6 8 10 12 14 16 18 20 22 24],'XTickLabel',...
    {'00:00','2:00','4:00','6:00','8:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'},...
    'XTickLabelRotation',45);
box(axes1,'on');
grid(axes1,'on');
title([' VTEC at ' stationname ' station date: ' year '/' month '/' day],'FontAngle','italic','FontName','Times Ten LT Std Roman');

names=['VTEC (TECU)'];
ylabel(names,'FontSize',10,'FontAngle','italic','FontName','Times Ten LT Std Roman','FontWeight','bold')
hold off
pbaspect([12,5,1]);

name2=[ stationname year YYDOY(3:5) '_' sys];


% cd .\plot\VTEC_Station\;

if exist(['plot\VTEC_Station\' num2str(year) '\' YYDOY(3:5) ],'dir')==0
    mkdir(['plot\VTEC_Station\' num2str(year) '\' YYDOY(3:5) ]);
end
cd(['plot\VTEC_Station\' num2str(year) '\' YYDOY(3:5) ])
export_fig(name2,'-jpg','-r600');
cd .\..\..\..\..;
close all;

end