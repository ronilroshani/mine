function neu = xyz2neu(xyz,plh)
%XYZ2NEU dX, dY, dZ (Cartesian) to dN, dE, dU (North, East, Up). 
%        Converts delta quantities dX, dY and dZ (Cartesian) into
%        dN, dE, dU (North, East, Up):
%
%           neu = xyz2neu(xyz,plh)
%
%        plh is a vector wtih Lattitude and Longitude of the point.

%        H. van der Marel, LGR, 07-05-95
%        (c) Geodetic Computing Centre, TU Delft

%[m,n]=size(xyz);
%if n~=3 & m==3, xyz=xyz';plh=plh';, end
%[m1,n1]=size(xyz);
%[m2,n2]=size(plh);
%if m1~=m2 & m2 ~=1, disp(['error XYZ2NEU']);, end

neu = [ -sin(plh(:,1)).*cos(plh(:,2)).*xyz(:,1) - sin(plh(:,1)).*sin(plh(:,2)).*xyz(:,2) + cos(plh(:,1)).*xyz(:,3) ...
                       -sin(plh(:,2)).*xyz(:,1) +                cos(plh(:,2)).*xyz(:,2)                           ...
         cos(plh(:,1)).*cos(plh(:,2)).*xyz(:,1) + cos(plh(:,1)).*sin(plh(:,2)).*xyz(:,2) + sin(plh(:,1)).*xyz(:,3) ];
%if n~=3 & m==3, neu=neu';plh=plh';, end
    


