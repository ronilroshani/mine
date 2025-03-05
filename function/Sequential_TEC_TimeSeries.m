function []=Sequential_TEC_TimeSeries (fileobs,seqpath)

% This function calculates the consecutive time series of TEC.
% seqpath : path of the input data to calculate the consecutive time series
% of TEC.
%% ------------------------------------------------------------------------

addpath(seqpath);

data_dir= dir(seqpath);
num_file=length(data_dir)-2;


for i =1:num_file
    if i==1
        Sequential_TEC_TimeSeries=[];DOY=[];
    end
    
    if length(data_dir(i+2).name)==32 & data_dir(i+2).name(1:12) == 'VTEC_Station'
        filemat=load(data_dir(i+2).name);
        
        station=data_dir(i+2).name(14:17);
        year= data_dir(i+2).name (19:22);
        month= data_dir(i+2).name (24:25);
        day= data_dir(i+2).name (27:28);
        
        
        if i==1
            datestr=[year month day];
        elseif i==num_file
            dateend=[year month day];
        end
        
        t1=[str2double(year),str2double(month),str2double(day),0,0,0];
        t2=repmat(t1,size(filemat.(['VTEC_Station_' station '_' year '_' month '_' day]),1),1);
        [~, ~, doy, ~] = greg2gps(t2);
        
        Sequential_TEC_TimeSeries=[Sequential_TEC_TimeSeries;filemat.(['VTEC_Station_' station '_' year '_' month '_' day])];
        DOY=[DOY;(doy)];
        
    end
end

cd ./Results/Sequential_TimeSeries/Mat_file
name=(['Sequential_TEC_TimeSeries_' station '_' datestr '_' dateend]);
eval([name '= Sequential_TEC_TimeSeries']);
save (name,name);
cd .\..\..\..;

cd ./Results/Sequential_TimeSeries/text_file;
T = table(num2str(DOY),num2str(Sequential_TEC_TimeSeries),'VariableNames',{'DOY','Sequential_TEC_TimeSeries'});
writetable(T,[name '.txt'],'Delimiter','\t');
cd .\..\..\..;

%% plot

cd ./plot/Sequential_TimeSeries

figure1 = figure('Color',[1 1 1]);
axes1 = axes('Parent',figure1);
hold(axes1,'on');

plot(Sequential_TEC_TimeSeries,'-','LineWidth',1,'color','b');

tickLocations = DOY; 
tickLabels    = (unique(DOY)); 

for i1=1:length(tickLabels)
    Xtic(i1,1)=i1*2880;
end
set(gca,'XTick',Xtic,'XTickLabel',tickLabels,'XTickLabelRotation',45)


box(axes1,'on');
grid(axes1,'on');
title([' Sequential_ TEC_ TimeSeries at ' station ],'FontAngle','italic','FontName','Times Ten LT Std Roman');

ylabel('TEC [TECU] ','FontSize',10,'FontAngle','italic','FontName','Times Ten LT Std Roman','FontWeight','bold')
hold off
pbaspect([12,5,1]);

export_fig(name,'-jpg','-r600');

cd .\..\..;
warning off;
close;
end
