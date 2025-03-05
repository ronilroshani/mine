function Save_Result(date, S_path,stationname,Doy,pos, STEC, VTEC, VTEC_Station,location)

%% Save file
year  = num2str(date(1,1));
month = num2str(date(1,2),'%.2d');
day  = num2str(date(1,3),'%.2d');
% mkdir(S_path);
   
if (Doy(1)) < 10
    if exist(['Results\Mat_result\' num2str(year) '\00' num2str(Doy(1))],'dir')==0
        mkdir(['Results\Mat_result\' num2str(year) '\00' num2str(Doy(1))]);
    end
elseif (Doy(1)) > 9 & (Doy(1)) < 100
    if exist(['Results\Mat_result\' num2str(year) '\0' num2str(Doy(1))],'dir')==0
        mkdir(['Results\Mat_result\' num2str(year) '\0' num2str(Doy(1))]);
    end
elseif (Doy(1)) > 99
    if exist(['Results\Mat_result\' num2str(year) '\' num2str(Doy(1))],'dir')==0
        mkdir(['Results\Mat_result\' num2str(year) '\' num2str(Doy(1))]);
    end
end

name1 = ['STEC'];
name2 = ['VTEC'];
name3 = ['VTEC_Station_' stationname '_' year '_' month '_' day];
name4 = ['location'];

eval([name1 '= STEC;']);
eval([name2 '= VTEC;']);
eval([name3 '= VTEC_Station']);
eval([name4 '= location;']);

filename  = [S_path 'STEC_' stationname '_' year '_' month '_' day];
filename1 = [S_path 'VTEC_' stationname '_' year '_' month '_' day];
filename2 = [S_path 'VTEC_Station_' stationname '_' year '_' month '_' day];
filename3 = [S_path 'location_' stationname '_' year '_' month '_' day];

save(filename,name1,'pos');
save(filename1,name2,'pos');
save(filename2,name3);%,'pos')
save(filename3,name4);

disp(['Complete to Calculate STEC & VTEC & VTEC_station at ' stationname ' station'])
end