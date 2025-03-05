function Plot_TimeSeries_TEC(year,month,day,VTEC_station,YYDOY,sys)
global stationname S_path
addpath(['Results\Mat_result\' num2str(year) '\' num2str(YYDOY(3:5)) '\'])
Time_TEC  = (0:30:24*3600-30)/3600;                 %   Time rate 1 second
Time_TEC_noon  = (12*3600:30:14*3600-30)/3600;      %   Time rate 12h:14h

%% Figure#1 STEC & VTEC
figure1 = figure('Color',[1 1 1]);
axes1 = axes('Parent',figure1);
hold(axes1,'on');

sgtitle(['STEC and VTEC at ' stationname ' station date: ' year '/' month '/' day],'FontAngle','italic','FontName','Times Ten LT Std Roman');

% load the results from save path
filename  = [S_path 'STEC_' stationname '_' year '_' month '_' day];
filename1 = [S_path 'VTEC_' stationname '_' year '_' month '_' day];
filename2 = [S_path 'VTEC_Station_' stationname '_' year '_' month '_' day];

load(filename)
load(filename1)
load(filename2)


name1 = ['STEC_' year '_' month '_' day];
name2 = ['VTEC_' year '_' month '_' day];
name3 = ['VTEC_Station_' year '_' month '_' day];

subplot1 = subplot(2,1,1);
hold(subplot1,'on');

plot(Time_TEC,STEC,'-');
set(subplot1,'XTick',[0 2 4 6 8 10 12 14 16 18 20 22 24],'XTickLabel',...
    {'0:00','2:00','4:00','6:00','8:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'},...
    'XTickLabelRotation',45);
box(subplot1,'on');
grid(subplot1,'on');

ylabel('STEC (TECU)','FontSize',10,'FontAngle','italic','FontName','Times Ten LT Std Roman','FontWeight','bold')
title('Slant Total Electron Content (STEC)','FontSize', 10,'FontAngle','italic','FontName','Times Ten LT Std Roman')

subplot2 = subplot(2,1,2);
hold(subplot2,'on');

plot(Time_TEC,VTEC,'-','LineWidth',1,'Color',[0.50,0.50,0.50]);
hold on;
plot(Time_TEC,VTEC_station,'-','LineWidth',1,'color','b');

set(subplot2,'XTick',[0 2 4 6 8 10 12 14 16 18 20 22 24],'XTickLabel',...
    {'0:00','2:00','4:00','6:00','8:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'},...
    'XTickLabelRotation',45);
box(subplot2,'on');
grid(subplot2,'on');

ylabel('VTEC (TECU)','FontSize',10,'FontAngle','italic','FontName','Times Ten LT Std Roman','FontWeight','bold')
title('Vertical Total Electron Content (VTEC)','FontSize', 10,'FontAngle','italic','FontName','Times Ten LT Std Roman')

hold off

name1=['VTEC & STEC ' stationname year YYDOY(3:5) '_' sys];
% cd .\plot\VTEC_STEC_timeseries\;

if exist(['plot\VTEC_STEC_timeseries\' num2str(year) '\' YYDOY(3:5) ],'dir')==0
    mkdir(['plot\VTEC_STEC_timeseries\' num2str(year) '\' YYDOY(3:5) ]);
end
cd(['plot\VTEC_STEC_timeseries\' num2str(year) '\' YYDOY(3:5) ])
export_fig(name1,'-jpg','-r600');
cd .\..\..\..\..;
close all;

end
