function tec=figure_geoshow_manual(lon1, lat1,V,title1,titlee,LAT_TK,LON_TK)
% figure('units','normalized','Color',[1 1 1]);
load coast;
coast(:,2)=long;coast(:,1)=lat;
Parent1=figure('units','normalized','Color',[1 1 1]);
axes1 = axes('Parent',Parent1);
hold(axes1,'on');
%  worldmap iran
% worldmap world
% [LAT_TK,LON_TK]=read_shp_iran_tehran_karaj_map;
% [LAT_TK,LON_TK]=read_shp_iran_tehran_karaj;

% LON_TK=[];LAT_TK=[];
% LAT_TK=[37;37;34;34];
% LON_TK=[50;54;54;50];

% worldmap([34.5 36.5],[50 53.5]);
worldmap([min(LAT_TK) max(LAT_TK)],[min(LON_TK) max(LON_TK)]);

% worldmap ([min(LAT_TK)-2 max(LAT_TK)+2],[min(LON_TK)-2 max(LON_TK)+2]);

lat_v=reshape(lat1,[],1);
lon_v=reshape(lon1,[],1);
[IN,ON] = inpolygon(lat_v,lon_v,LAT_TK,LON_TK) ;
tec=nan(size(V,1),size(V,2));

 tec(IN)=V(IN);           
            
%worldmap([min(lat1(:,1))-0.2 max(lat1(:,1))+0.2],[min(lon1(1,:))-0.2 max(lon1(1,:))+0.2]);
%% 
setm(axes1,'FontSize',5,'FontWeight','bold')
geoshow(lat1,lon1,tec,...
    'DisplayType','surface');       %Texturemap                    % show the matrix V
 hold on;
 title(title1,'FontSize', 10,'FontAngle','italic','FontName','Times New Roman');

%  contourm(repmat(lat1,1,length(lon1)),...
%  repmat(lon1,length(lat1),1),V','LevelStep',C,'LineWidth',0.5,'LineColor',[127 127 127]/256);
hold on
% load worldlo.mat;
% a=POline(1);
% latg=a.lat;
% longg=a.long;
% h=plotm(latg,longg,'k-');
%  plotm(LAT_TKz,LON_TKz,'.','Color',[.15 .15 .15]);
% setm(gca,'mapprojection','mollweid','ParallelLabel','on','MeridianLabel','on')
% hold on
% land=shaperead('landareas','UseGeocoords',true);
% geoshow([land.Lat],[land.Lon],'color','black')
%colormap('jet')
%colormap(brewermap(11,'RdBu'))
cmap = [flip(brewermap(256,'spectral'))];
colormap(jet(4096));

hcb = colorbar('FontWeight','bold','FontSize',6); 
set(get(hcb,'Xlabel'),'String',cmap,'FontSize', 6)       % color-bar
set(get(hcb,'Xlabel'),'String',titlee,'FontSize', 6)       % color-bar

vv=reshape(V,[],1);
 %caxis([mean(vv)-2*mean(vv)/3 mean(vv)+2*mean(vv)/3]);
caxis([0 100]);   
% range of the color-bar
% range of the color-bar
hold on; 
plotm(coast(:,1),coast(:,2),'.k');
%maximum=max(abs(V));
% caxis([-max(maximum) max(maximum)]);                      % range of the color-bar
box on;                                                     % adding box to surrounding the box