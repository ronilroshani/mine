function [E,A]= Get_EA(sx,sy,sz,x,y,z)
%GET_EL Summary of this function goes here
%   Detailed explanation goes here
%% ------------------------------------------------------------------------

[sb,sl]=XYZtoBLH(sx,sy,sz);
llaR = ecef2lla([sx,sy,sz],'WGS84');
llaR(1) = deg2rad(llaR(1)); llaR(2) = deg2rad(llaR(2));
T=[-sin(sb)*cos(sl) -sin(sb)*sin(sl) cos(sb);
    -sin(sl)               cos(sl)         0;
    cos(sb)*cos(sl) cos(sb)*sin(sl)  sin(sb)];%transition matrix(XYZ to NEU)
deta_xyz=[x,y,z]-[sx,sy,sz];
NEU=T*(deta_xyz)';
neu1 = xyz2neu(deta_xyz,llaR);

E=atan(NEU(3)/sqrt(NEU(1)*NEU(1)+NEU(2)*NEU(2)));
A=atan(abs(NEU(2)/NEU(1)));

[h,az] = enu2altaz(neu1(1),neu1(2),neu1(3));

if NEU(1)>0
    if NEU(2)>0
    else
        A=2*pi-A;
    end
else
    if NEU(2)>0
        A=pi-A;
    else
        A=pi+A;
    end 
end
end

