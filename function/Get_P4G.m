function  [ GPSP4 ] = Get_P4G(POSItion,OBSER,Sites_Info,Cutoff)
%% get smoothed P4 observations
% INPUT:
%     path_obs: storage path of rinex files
%     Sites_Info: name and coordinate information of the stations
%     lim: cut-off angle

%% ---------------------------------------------------------------------
global stationname
Coor=Sites_Info.coor;
stations=Sites_Info.name;
doys=Sites_Info.doy;
% current_doy=num2str(unique(doys));

if exist('M_P4','dir')==0
    mkdir('M_P4');
end


len=1;%length(list_obs);
for i=1:len
    site=stationname;
    doy=doys;
    indices=doys==(doy);
    index=find(strcmpi(site,stations(indices)), 1);
    sx=Coor(index,1);
    sy=Coor(index,2);
    sz=Coor(index,3);
    if sx==0 && sy==0 && sz==0
        continue;
    end
    %get rid of sites without receivers coordinates.
    obs=cutobs(POSItion,sx,sy,sz,OBSER,Cutoff);
    fields=fieldnames(obs);
    
    if ~isnan(find(strcmp(fields, 'GPSP1' )))
        if ~(all(all(obs.GPSP1==0)) || all(all(obs.GPSP2==0)))
            [GPSP4]=GPS_prepro(obs);
            if all(all(GPSP4==0))
                continue;
            end
        else
            GPSP4=nan(size(obs.GPSP1,1),size(obs.GPSP1,2));
        end
        if exist(['M_P4/GPS/' num2str(doy)],'dir')==0
            mkdir(['M_P4/GPS/' num2str(doy)]);
        end
        filenameP4=['M_P4/GPS/' num2str(doy) '/' site num2str(doy) 'P4.mat'];
        save(filenameP4,'GPSP4','-mat');        
    end
    
    %     clear GPSP4
    clear POSItion;
    
end

end
%% ----------------subfunction-----------------
function obs=cutobs(POSItion,sx,sy,sz,obs,Cutoff)
%% get the observations under specified cut-off angle
% INPUT:
%     sate: precise coordinates of the satellites
%     sx: X coordinate of the station
%     sy: Y coordinate of the station
%     sz: Z coordinate of the station
%     obs: original observation structs
%     lim: cut-off angle
% OUTPUT:
%      obs: updated observation structs
fields=fieldnames(obs);
%% cut gps obs
if ~isnan(find(strcmp(fields, 'GPSP1' )))
    gpsline=size(obs.GPSP2,1);
    gpssate=size(POSItion.gpsx,2);
    if gpsline<2880
        obs.GPSL1(gpsline+1:2880,:)=0;obs.GPSL2(gpsline+1:2880,:)=0;
        obs.GPSP1(gpsline+1:2880,:)=0;obs.GPSP2(gpsline+1:2880,:)=0;
    end
    gpsl=size(obs.GPSP2,2);
    if gpsl<gpssate
        obs.GPSL1(:,gpsl+1:gpssate)=0;obs.GPSL2(:,gpsl+1:gpssate)=0;
        obs.GPSP1(:,gpsl+1:gpssate)=0;obs.GPSP2(:,gpsl+1:gpssate)=0;
    end
%     if gpsl>gpssate
%         POSItion.gpsx(:,gpssate+1:gpsl)=0;POSItion.gpsy(:,gpssate+1:gpsl)=0;
%         POSItion.gpsz(:,gpssate+1:gpsl)=0;
%     end
    %     gpsdelete=find(sate_mark.gps==0);
    gpsnum=size(POSItion.gpsx,2);
    if gpsnum<size(obs.GPSP1,2)
        %         if ~isnan(gpsdelete)
        %             k=length(gpsdelete);
        %             for i=k:-1:1
        %                 obs.GPSP1(:,gpsdelete(i))=[];
        %                 obs.GPSP2(:,gpsdelete(i))=[];
        %                 obs.GPSL1(:,gpsdelete(i))=[];
        %                 obs.GPSL2(:,gpsdelete(i))=[];
        %             end
        %         end
        obs.GPSP1=obs.GPSP1(:,1:gpsnum);
        obs.GPSP2=obs.GPSP2(:,1:gpsnum);
        obs.GPSL1=obs.GPSL1(:,1:gpsnum);
        obs.GPSL2=obs.GPSL2(:,1:gpsnum);
    end
    obs.GPSP1(isnan(obs.GPSP1))=0;
    obs.GPSP2(isnan(obs.GPSP2))=0;
    obs.GPSL1(isnan(obs.GPSL1))=0;
    obs.GPSL2(isnan(obs.GPSL2))=0;
    
    gpsx=POSItion.gpsx;gpsy=POSItion.gpsy;gpsz=POSItion.gpsz;
    
    for i=1:gpsnum
        for j=1:2880
            if obs.GPSL1(j,i)==0 || obs.GPSL2(j,i)==0 || obs.GPSP1(j,i)==0 || obs.GPSP2(j,i)==0
                obs.GPSL1(j,i)=0;obs.GPSL2(j,i)=0;obs.GPSP1(j,i)=0;obs.GPSP2(j,i)=0;
                continue;
            end
            [el,aaa]=Get_EA(sx,sy,sz,gpsx(j,i),gpsy(j,i),gpsz(j,i));
            if el<Cutoff(1,1)
                obs.GPSP1(j,i)=0;obs.GPSP2(j,i)=0;obs.GPSL1(j,i)=0;obs.GPSL2(j,i)=0;
                continue;
            end
        end
    end
end
end
%% ---------------------------subfunction----------------------------------
function [GPSPP4]=GPS_prepro(obs)
%% get smoothed GPS P4 observations
% INPUT:
%      obs: struct of rinex files
% OUTPUT:
%     GPSPP4: smoothed GPS P4 observation
%% --------------------------------------------------------------------------
%____delete incomplete epoch_____
global c
size2=size(obs.GPSP1,2);
LenSatG=32;%LEN(1,1);
if size2<LenSatG
    obs.GPSP1(:,size2+1:LenSatG)=0;obs.GPSP2(:,size2+1:LenSatG)=0;
    obs.GPSL1(:,size2+1:LenSatG)=0;obs.GPSL2(:,size2+1:LenSatG)=0;
end
GPSP4=zeros(2880,LenSatG);
GPSL4=zeros(2880,LenSatG);
GPSPP4=zeros(2880,LenSatG);
GPS_f1=1575.42*10^6;                 %_________________________________unit:Hz
GPS_f2=1227.6*10^6;                  %_________________________________unit:Hz
lamda_w=c/(GPS_f1-GPS_f2);       %____________________wide lane wavelength
L6=lamda_w*(obs.GPSL1-obs.GPSL2)-(GPS_f1*obs.GPSP1+GPS_f2*obs.GPSP2)/(GPS_f1+GPS_f2); %__MW observable
Li=obs.GPSL1-GPS_f1*obs.GPSL2/GPS_f2;
Nw=-L6;                   %_____________________wide lane ambiguity
for i=1:LenSatG  %i is PRN number
    %------divide arc---------------------------
    arc=Get_arc(L6(:,i));
    [arc_n,aaa]=size(arc);
    %----delete arc less than 10 epoches-------
    arc_d=[];
    for j=1:arc_n
        n_epoch=arc(j,2)-arc(j,1);
        if n_epoch<20
            for k=arc(j,1):arc(j,2)
                obs.GPSP1(k,i)=0;obs.GPSP2(k,i)=0;obs.GPSL1(k,i)=0;obs.GPSL2(k,i)=0;
                L6(k,i)=0;Nw(k,i)=0;Li(k,i)=0;
            end
            arc_d=[arc_d,j];
        end
    end
    arc(arc_d,:)=[];
    %----mw detect cycle slip------------------
    [arc_n,aaa]=size(arc);
    j=1;
    while j<arc_n+1 %j is arc number
        %----first epoch check----------
        e=arc(j,1);
        while 1
            if e+1==arc(j,2) || e==arc(j,2)
                break;
            end
            fir=Nw(e,i);sec=Nw(e+1,i);thi=Nw(e+2,i);
            firl=Li(e,i);secl=Li(e+1,i);thil=Li(e+2,i);
            sub=abs(fir-sec);sub2=abs(sec-thi);
            subl=abs(firl-secl);subl2=abs(secl-thil);
            if sub>1||sub2>1||subl>1||subl2>1
                L6(e,i)=0;obs.GPSL1(e,i)=0;obs.GPSL2(e,i)=0;obs.GPSP1(e,i)=0;obs.GPSP2(e,i)=0;Nw(e,i)=0;Li(e,i)=0;
                e=e+1;
                arc(j,1)=e;
            else
                arc(j,1)=e;
                break;
            end
        end
        %----detect------------------
        if arc(j,2)-arc(j,1)<20
            for k=arc(j,1):arc(j,2)
                obs.GPSP1(k,i)=0;obs.GPSP2(k,i)=0;obs.GPSL1(k,i)=0;obs.GPSL2(k,i)=0;
                L6(k,i)=0;Nw(e,i)=0;Li(k,i)=0;
            end
            arc(j,:)=[];
            arc_n=arc_n-1;
            continue;
        end
        ave_N=linspace(0,0,arc(j,2)-arc(j,1)+1);
        sigma2=linspace(0,0,arc(j,2)-arc(j,1)+1);
        sigma=linspace(0,0,arc(j,2)-arc(j,1)+1);
        ave_N(1)=Nw(arc(j,1),i);
        sigma2(1)=0;
        sigma(1)=0;
        count=2;
        for k=arc(j,1)+1:arc(j,2)-1 %----------------------check epoch k+1
            ave_N(count)=ave_N(count-1)+(Nw(k,i)-ave_N(count-1))/count;
            sigma2(count)=sigma2(count-1)+((Nw(k,i)-ave_N(count-1))^2-sigma2(count-1))/count;
            sigma(count)=sqrt(sigma2(count));
            T=abs(Nw(k+1,i)-ave_N(count));
            I1=abs(Li(k+1,i)-Li(k,i));
            if T<4*sigma(count)&&I1<0.28   %-------------------------no cycle slip
                count=count+1;
                continue;
            else
                if k+1==arc(j,2)            %---------------------arc end
                    if k+1-arc(j,1)>20
                        L6(k+1,i)=0;obs.GPSP1(k+1,i)=0;obs.GPSP2(k+1,i)=0;
                        obs.GPSL1(k+1,i)=0;obs.GPSL2(k+1,i)=0;Nw(k+1,i)=0;
                        Li(k,i)=0;arc(j,2)=k;
                    else                     %------delete scatter epoches
                        for l=arc(j,1):k+1
                            L6(l,i)=0;obs.GPSP1(l,i)=0;obs.GPSP2(l,i)=0;
                            obs.GPSL1(l,i)=0;obs.GPSL2(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,:)=[];
                        j=j-1;
                        arc_n=arc_n-1;
                    end
                    break;
                end
                I2=abs(Li(k+2,i)-Li(k+1,i));
                if abs(Nw(k+2,i)-Nw(k+1,i))<1&&I2<1%-----------cycle slip
                    if k+1-arc(j,1)>20
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+1,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else                    %------delete scatter epoches
                        for l=arc(j,1):k
                            L6(l,i)=0;obs.GPSP1(l,i)=0;obs.GPSP2(l,i)=0;
                            obs.GPSL1(l,i)=0;obs.GPSL2(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,1)=k+1;
                        j=j-1;
                    end
                else                        %-----------------gross error
                    if k+1-arc(j,1)>20
                        L6(k+1,i)=0;obs.GPSP1(k+1,i)=0;obs.GPSP2(k+1,i)=0;
                        obs.GPSL1(k+1,i)=0;obs.GPSL2(k+1,i)=0;Nw(k+1,i)=0;Li(k,i)=0;
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+2,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else
                        for l=arc(j,1):k+1  %------delete scatter epoches
                            L6(l,i)=0;obs.GPSP1(l,i)=0;obs.GPSP2(l,i)=0;
                            obs.GPSL1(l,i)=0;obs.GPSL2(l,i)=0;Nw(l,i)=0;Li(k,i)=0;
                        end
                        arc(j,1)=k+2;
                        j=j-1;
                    end
                end
                break;
            end
        end
        j=j+1;
    end
    GPSP4(:,i)=obs.GPSP1(:,i)-obs.GPSP2(:,i);
    GPSL4(:,i)=(c/GPS_f1)*obs.GPSL1(:,i)-(c/GPS_f2)*obs.GPSL2(:,i);
    
    %--------------------------smoothing-------------------------
    GPSPP4(:,i)=GPSP4(:,i);
    for j=1:arc_n
        t=2;
        for k=arc(j,1)+1:arc(j,2)
            %P4(k,i)=mean(P4(arc(j,1):k,i))-L4(k,i)+mean(L4(arc(j,1):k,i));
            GPSPP4(k,i)=GPSPP4(k,i)/t+(GPSPP4(k-1,i)+GPSL4(k-1,i)-GPSL4(k,i))*(t-1)/t;
            t=t+1;
        end
        GPSPP4(arc(j,1):arc(j,1)+4,i)=0;
    end
    %--------remove bad P4---------------------
    arc=Get_arc(GPSPP4(:,i));
    [arc_n,aaa]=size(arc);
    for j=1:arc_n
        ave=mean(GPSPP4(arc(j,1):arc(j,2),i));
        if abs(ave)>10
            for kk=arc(j,1):arc(j,2)
                GPSPP4(kk,i)=0;
            end
        end
    end
end
end

function arc = Get_arc(array)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
len=length(array);
arc=[];
for i=1:len
    if i==len
        if array(i)~=0
            arc=[arc,i];
        end
        continue;
    end
    if i==1&&array(i)~=0
        arc=[arc,i];
    end
    if array(i)==0&&array(i+1)~=0
        arc=[arc,i+1];
        continue;
    end
    if array(i)~=0&&array(i+1)==0
        arc=[arc,i];
        continue;
    end
end
if len==0,return,end
arc=reshape(arc,2,[]);
arc=arc';
end