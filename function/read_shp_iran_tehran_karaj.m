function [LAT_T,LON_T]=read_shp_iran_tehran_karaj
addpath('.\SHAPEFILE');
name='IRN_adm1.shp';warning off; 
S = shaperead(name);
% LAT_T=[];LON_T=[];
% for i=1:size(S,1)
% LON_T1=S(i).X;
% LAT_T1=S(i).Y;
% LON_T1=LON_T1';
% LAT_T1=LAT_T1';
% LAT_T=[LAT_T;LAT_T1];
% LON_T=[LON_T;LON_T1];
% LON_T1=[];LAT_T1=[];
% end
LON_T1=S(1).X;LON_T1=LON_T1';
LAT_T1=S(1).Y;LAT_T1=LAT_T1';


% LON_T2=S(19).X;LON_T2=LON_T2';
% LAT_T2=S(19).Y;LAT_T2=LAT_T2';
% 
% 
% LON_T3=S(20).X;LON_T3=LON_T3';
% LAT_T3=S(20).Y;LAT_T3=LAT_T3';
% 
% 
% LON_T4=S(22).X;LON_T4=LON_T4';
% LAT_T4=S(22).Y;LAT_T4=LAT_T4';
% 
% 
% LON_T5=S(23).X;LON_T5=LON_T5';
% LAT_T5=S(23).Y;LAT_T5=LAT_T5';
% 
% 
% LON_T6=S(25).X;LON_T6=LON_T6';
% LAT_T6=S(25).Y;LAT_T6=LAT_T6';


LON_T7=S(28).X;LON_T7=LON_T7';
LAT_T7=S(28).Y;LAT_T7=LAT_T7';

LAT_T=[LAT_T1(1:end-1,:);LAT_T7(1:end-1,:)];
LON_T=[LON_T1(1:end-1,:);LON_T7(1:end-1,:)];


% LAT_T=[LAT_T3(1:end,:);LAT_T4(1:end,:);LAT_T1(1:end,:);...
%     LAT_T7(1:end,:);LAT_T6(1:end,:);LAT_T2(1:end,:);LAT_T5(1:end,:)];
% 
% LON_T=[LON_T3(1:end,:);LON_T4(1:end,:);LON_T1(1:end,:);...
%     LON_T7(1:end,:);LON_T6(1:end,:);LON_T2(1:end,:);LON_T5(1:end,:)];







% 
% LON_TK=[];LAT_TK=[];
% LAT_TK=[37;37;34;34];
% LON_TK=[50;54;54;50];


% for i=1:size(S,1)
%     lon_st=S(i).X;
%     lat_st=S(i).Y;
%     
%     LON_TK=[LON_TK;lon_st'];
%     LAT_TK=[LAT_TK;lat_st'];
% end


