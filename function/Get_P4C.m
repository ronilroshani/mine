function  [BDSP4] = Get_P4C(POSItion,OBSER,Sites_Info,Cutoff)
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
    
    if ~isnan(find(strcmp(fields, 'BDSP1' )))
        if ~(all(all(obs.BDSP1==0)) || all(all(obs.BDSP2==0)))
            [BDSP4]=BDS_prepro(obs);
            if all(all(BDSP4==0))
                continue;
            end
        else
            BDSP4=nan(size(obs.BDSP1,1),size(obs.BDSP1,2));
        end
        if exist(['M_P4/BDS/' num2str(doy)],'dir')==0
            mkdir(['M_P4/BDS/' num2str(doy)]);
        end
        filenameP4=['M_P4/BDS/' num2str(doy) '/' site num2str(doy) 'P4.mat'];
        save(filenameP4,'BDSP4','-mat');     
    end
    
    %     clearBDSP4;
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
%% cut bds obs
if ~isnan(find(strcmp(fields, 'BDSP1' )))
    bdsline=size(obs.BDSP1,1);
    bdssate=size(POSItion.bdsx,2);
    if bdsline<2880
        obs.BDSL1(bdsline+1:2880,:)=0;obs.BDSL2(bdsline+1:2880,:)=0;
        obs.BDSP1(bdsline+1:2880,:)=0;obs.BDSP2(bdsline+1:2880,:)=0;
    end
    bdsl=size(obs.BDSP1,2);
    if bdsl<bdssate
        obs.BDSL1(:,bdsl+1:bdssate)=0;obs.BDSL2(:,bdsl+1:bdssate)=0;
        obs.BDSP1(:,bdsl+1:bdssate)=0;obs.BDSP2(:,bdsl+1:bdssate)=0;
    end
%     if bdsl>bdssate
%         POSItion.glox(:,bdssate+1:bdsl)=0;POSItion.gloy(:,bdssate+1:bdsl)=0;
%         POSItion.gloz(:,bdssate+1:bdsl)=0;
%     end
    
    %     bdsdelete=find(sate_mark.bds==0);
    bdsnum=size(POSItion.bdsx,2);
    if bdsnum<size(obs.BDSP1,2)
        %         if ~isnan(bdsdelete)
        %             k=length(bdsdelete);
        %             for i=k:-1:1
        %                 obs.BDSP1(:,bdsdelete(i))=[];
        %                 obs.BDSP2(:,bdsdelete(i))=[];
        %                 obs.BDSL1(:,bdsdelete(i))=[];
        %                 obs.BDSL2(:,bdsdelete(i))=[];
        %             end
        %         end
        obs.BDSP1=obs.BDSP1(:,1:bdsnum);
        obs.BDSP2=obs.BDSP2(:,1:bdsnum);
        obs.BDSL1=obs.BDSL1(:,1:bdsnum);
        obs.BDSL2=obs.BDSL2(:,1:bdsnum);
    end
    obs.BDSP1(isnan(obs.BDSP1))=0;
    obs.BDSP2(isnan(obs.BDSP2))=0;
    obs.BDSL1(isnan(obs.BDSL1))=0;
    obs.BDSL2(isnan(obs.BDSL2))=0;
    
    bdsx=POSItion.bdsx;bdsy=POSItion.bdsy;bdsz=POSItion.bdsz;
    
    for i=1:bdsnum
        for j=1:2880
            if obs.BDSL1(j,i)==0 || obs.BDSL2(j,i)==0 || obs.BDSP1(j,i)==0 || obs.BDSP2(j,i)==0
                obs.BDSL1(j,i)=0;obs.BDSL2(j,i)=0;obs.BDSP1(j,i)=0;obs.BDSP2(j,i)=0;
                continue;
            end
            [el,aaa]=Get_EA(sx,sy,sz,bdsx(j,i),bdsy(j,i),bdsz(j,i));
            if el<Cutoff(3,1)
                obs.BDSP1(j,i)=0;obs.BDSP2(j,i)=0;obs.BDSL1(j,i)=0;obs.BDSL2(j,i)=0;
                continue;
            end
        end
    end
end

end
%% ---------------------------subfunction----------------------------------

%% ---------------------------subfunction-------------------------
function [BDSPP4]=BDS_prepro(obs)
%% get smoothed BDS P4 observations
% INPUT:
%      obs: struct of rinex files
% OUTPUT:
%     BDSPP4: smoothed BDS P4 observation
%% -----------------------------------------------------------------
%____delete incomplete epoch_____
global c
size2=size(obs.BDSP1,2);
LenSatC=46;%LEN(2,1);
if size2<LenSatC
    obs.BDSP1(:,size2+1:LenSatC)=0;obs.BDSP2(:,size2+1:LenSatC)=0;
    obs.BDSL1(:,size2+1:LenSatC)=0;obs.BDSL2(:,size2+1:LenSatC)=0;
end
BDSP4=zeros(2880,LenSatC);
BDSL4=zeros(2880,LenSatC);
BDSPP4=zeros(2880,LenSatC);
BDS_f2=1561.098*10^6;                 %_________________________________unit:Hz
BDS_f7=1207.140*10^6;                  %_________________________________unit:Hz
lamda_w=c/(BDS_f2-BDS_f7);       %____________________wide lane wavelength
L6=lamda_w*(obs.BDSL1-obs.BDSL2)-(BDS_f2*obs.BDSP1+BDS_f7*obs.BDSP2)/(BDS_f2+BDS_f7); %__MW observable
Li=obs.BDSL1-BDS_f2*obs.BDSL2/BDS_f7;
Nw=-L6;                   %_____________________wide lane ambiguity
for i=1:LenSatC  %i is PRN number
    %------divide arc---------------------------
    arc=Get_arc(L6(:,i));
    [arc_n,aaa]=size(arc);
    %----delete arc less than 10 epoches-------
    arc_d=[];
    for j=1:arc_n
        n_epoch=arc(j,2)-arc(j,1);
        if n_epoch<20
            for k=arc(j,1):arc(j,2)
                obs.BDSP1(k,i)=0;obs.BDSP2(k,i)=0;obs.BDSL1(k,i)=0;obs.BDSL2(k,i)=0;
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
                L6(e,i)=0;obs.BDSL1(e,i)=0;obs.BDSL2(e,i)=0;obs.BDSP1(e,i)=0;obs.BDSP2(e,i)=0;Nw(e,i)=0;Li(e,i)=0;
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
                obs.BDSP1(k,i)=0;obs.BDSP2(k,i)=0;obs.BDSL1(k,i)=0;obs.BDSL2(k,i)=0;
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
                        L6(k+1,i)=0;obs.BDSP1(k+1,i)=0;obs.BDSP2(k+1,i)=0;
                        obs.BDSL1(k+1,i)=0;obs.BDSL2(k+1,i)=0;Nw(k+1,i)=0;
                        Li(k,i)=0;arc(j,2)=k;
                    else                     %------delete scatter epoches
                        for l=arc(j,1):k+1
                            L6(l,i)=0;obs.BDSP1(l,i)=0;obs.BDSP2(l,i)=0;
                            obs.BDSL1(l,i)=0;obs.BDSL2(l,i)=0;Nw(l,i)=0;
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
                            L6(l,i)=0;obs.BDSP1(l,i)=0;obs.BDSP2(l,i)=0;
                            obs.BDSL1(l,i)=0;obs.BDSL2(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,1)=k+1;
                        j=j-1;
                    end
                else                        %-----------------gross error
                    if k+1-arc(j,1)>20
                        L6(k+1,i)=0;obs.BDSP1(k+1,i)=0;obs.BDSP2(k+1,i)=0;
                        obs.BDSL1(k+1,i)=0;obs.BDSL2(k+1,i)=0;Nw(k+1,i)=0;Li(k,i)=0;
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+2,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else
                        for l=arc(j,1):k+1  %------delete scatter epoches
                            L6(l,i)=0;obs.BDSP1(l,i)=0;obs.BDSP2(l,i)=0;
                            obs.BDSL1(l,i)=0;obs.BDSL2(l,i)=0;Nw(l,i)=0;Li(k,i)=0;
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
    BDSP4(:,i)=obs.BDSP1(:,i)-obs.BDSP2(:,i);
    BDSL4(:,i)=(c/BDS_f2)*obs.BDSL1(:,i)-(c/BDS_f7)*obs.BDSL2(:,i);
    
    %--------------------------smoothing-------------------------
    BDSPP4(:,i)=BDSP4(:,i);
    for j=1:arc_n
        t=2;
        for k=arc(j,1)+1:arc(j,2)
            %P4(k,i)=mean(P4(arc(j,1):k,i))-L4(k,i)+mean(L4(arc(j,1):k,i));
            BDSPP4(k,i)=BDSPP4(k,i)/t+(BDSPP4(k-1,i)+BDSL4(k-1,i)-BDSL4(k,i))*(t-1)/t;
            t=t+1;
        end
        BDSPP4(arc(j,1):arc(j,1)+4,i)=0;
    end
    %--------remove bad P4---------------------
    arc=Get_arc(BDSPP4(:,i));
    [arc_n,aaa]=size(arc);
    for j=1:arc_n
        ave=mean(BDSPP4(arc(j,1):arc(j,2),i));
        if abs(ave)>10
            for kk=arc(j,1):arc(j,2)
                BDSPP4(kk,i)=0;
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