function [T] = Georgd2yearday(year,month,day)
% year = 2007;
% month = 7;
% day= 17;
Y = (1940:4:2040)';
temp = find (year ==Y);
    
if isempty(temp)==1
if month ==1
    T = 0+day;
elseif month ==2
    T = 31+day;
elseif month ==3
    T = 31+28+day;
elseif month ==4
    T = 31+28+31+day;
elseif month ==5
    T = 31+28+31+30+day;
elseif month ==6
    T = 31+28+31+30+31+day;
elseif month ==7
    T = 31+28+31+30+31+30+day;
elseif month ==8
    T = 31+28+31+30+31+30+31+day;
elseif month ==9
    T = 31+28+31+30+31+30+31+31+day;
elseif month ==10
    T = 31+28+31+30+31+30+31+31+30+day;
elseif month ==11
    T = 31+28+31+30+31+30+31+31+30+31+day;
elseif month ==12
    T = 31+28+31+30+31+30+31+31+30+31+30+day;
end
else
    if month ==1
    T = 0+day;
elseif month ==2
    T = 31+day;
elseif month ==3
    T = 31+29+day;
elseif month ==4
    T = 31+29+31+day;
elseif month ==5
    T = 31+29+31+30+day;
elseif month ==6
    T = 31+29+31+30+31+day;
elseif month ==7
    T = 31+29+31+30+31+30+day;
elseif month ==8
    T = 31+29+31+30+31+30+31+day;
elseif month ==9
    T = 31+29+31+30+31+30+31+31+day;
elseif month ==10
    T = 31+29+31+30+31+30+31+31+30+day;
elseif month ==11
    T = 31+29+31+30+31+30+31+31+30+31+day;
elseif month ==12
    T = 31+29+31+30+31+30+31+31+30+31+30+day;
    end
end