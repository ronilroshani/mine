function [N2]=BS(u,j);
% clear all;
% clc;

% j=input('level(j)=');

mj=2^j+2;
d=1/(mj-2);

U=[0,0,0:d:1,1,1];        % U is a knote vector

% u=0:0.005:1;
%----------------------ZERO degree------------------------------------------
N0=zeros(length(u),length(U)-1);
for i=1:length(U)-1
    for t=1:length(u)
        if u(t)>=U(i) & u(t)<U(i+1)
            N0(t,i)=1;
        elseif u(t)==1 & i==length(U)-3
            N0(t,i)=1;
        else
            N0(t,i)=0;
        end
    end
end

%-------------------------1-st degree---------------------------------------
N1=zeros(length(u),length(U)-2);
for i=1:length(U)-2
    for t=1:length(u)

        if U(i)==U(i+1)
            a=0;
        else
            a=(u(t)-U(i))/(U(i+1)-U(i));
        end

        if U(i+1)==U(i+2)
            b=0;
        else
            b=(U(i+2)-u(t))/(U(i+2)-U(i+1));
        end

        N1(t,i)=a*N0(t,i)+b*N0(t,i+1);
    end
end

%------------------------2-nd degree---------------------------------------
N2=zeros(length(u),length(U)-3);
for i=1:length(U)-3
    for t=1:length(u)

        if U(i)==U(i+2)
            c=0;
        else
            c=(u(t)-U(i))/(U(i+2)-U(i));
        end

        if U(i+1)==U(i+3)
            d=0;
        else
            d=(U(i+3)-u(t))/(U(i+3)-U(i+1));
        end

        N2(t,i)=c*N1(t,i)+d*N1(t,i+1);
    end
end

% color=['r','g','b','c','m','y','k'];
% color=[color,color,color];
% 
% figure(1);
% for j=1:min(size(N2))
%     plot(u,N2(:,j),color(j),'linewidth',3);
%     hold on;
% end
% hold off;
% 
% 
% % -------------------2D-Bspline---------------------------------------------
% [x,y]=meshgrid(u);
% BS_2D=N2(:,6)*N2(:,5)';
% figure(2);
% surf(x,y,BS_2D);
% shading interp;
% colormap jet;
% colorbar;