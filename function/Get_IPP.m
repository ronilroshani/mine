function [b,s,s1] = Get_IPP(E,A,B,L,z,t_r)
%GET_IPP Summary of this function goes here
%   Detailed explanation goes here
%% ------------------------------------------------------------------------

t=pi/2-E-z;
b=asin(sin(B)*cos(t)+cos(B)*sin(t)*cos(A));
s1=L+asin(sin(t)*sin(A)/cos(b));
s=s1+t_r-pi;
end