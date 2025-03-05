 function[lat_tk,lon_tk]=read_shp_iran_tehran_karaj_map
addpath('.\SHAPEFILE');
name='IRN_adm1.shp';warning off; 
S = shaperead(name);

lon_krj=S(1).X;lon_krj=lon_krj';
lat_krj=S(1).Y;lat_krj=lat_krj';

lon_tehn=S(28).X;lon_tehn=lon_tehn';
lat_tehn=S(28).Y;lat_tehn=lat_tehn';

lat_tk=[lat_tehn(1:end-1,:);lat_krj(1:end-1,:)];
lon_tk=[lon_tehn(1:end-1,:);lon_krj(1:end-1,:)];
