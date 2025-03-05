addpath(genpath(pwd));
global stationname S_path

prompt = {' Year:',' Month:','Start Day:','End Day:'};
title = 'Date';
dims = [1 50];
def = [{'2022','6','1','1'}];

answer = inputdlg(prompt,title,dims,def);
YEAR = answer{1,1};YEAR=str2double(YEAR);
MONTH = answer{2,1};MONTH=str2double(MONTH);
DAYStart = answer{3,1};DAYStart=str2double(DAYStart);
DAYEnd = answer{4,1};DAYEnd=str2double(DAYEnd);
DAY = linspace(DAYStart,DAYEnd,DAYEnd-DAYStart+1);


for i=1:1:length(DAY)
    doy(i,1) = Georgd2yearday((YEAR),(MONTH),(DAY(1,i)));
end

DOY=num2str(doy);


for i=1:length(doy)
    if length(DOY(i,1:end)) < 2
        DOY(:,:) = ['00' DOY(i,1:end)];
    elseif length(DOY(i,1:end)) < 3
        DOY(:,:) = ['0' DOY(i,1:end)];
    elseif length(DOY(i,1:end)) == 3
        DOY(:,:) = [DOY(i,1:end)];
    end
    
    
    %% ----------------------------------------------------------------------------
    folder=['RINEX\' num2str(YEAR) '\' DOY(i,:)];
    
    addpath(folder);
    directory = dir(folder);
    num_file=length(directory)-2;
    
    for j=1%:num_file
        RINEX_filename=directory(j+2).name;
        TEC(RINEX_filename,YEAR,(DOY(i,:)));
        
        disp([num2str(j) ':' 'Station ' stationname  ' was processed.'])
        
    end
    %% ----------------------------------- VTEC MAP ---------------------------------
    
    Modelling_VTEC_Bspline(S_path)
end
rmpath(genpath(pwd));
toc;


% fileID = fopen([ 'Run_VTEC_Estimate_Stations.txt'],'w');
% fList = matlab.codetools.requiredFilesAndProducts('Run_VTEC_Estimate_Stations.m');
% for i=1:length(fList)
%     fprintf(fileID,['%s\n'],fList{1,i}(1,1:end));
% end
% fclose(fileID);
