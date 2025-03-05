function  [ GLOP4 ] = Get_P4R(POSItion,OBSER,Sites_Info,Cutoff)
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
    
    if ~isnan(find(strcmp(fields, 'GLOP1' )))
        if ~(all(all(obs.GLOP1==0)) || all(all(obs.GLOP2==0)));
            [GLOP4]=GLO_prepro(obs);
            if all(all(GLOP4==0))
                continue;
            end
        else
            GLOP4=nan(size(obs.GLOP1,1),size(obs.GLOP1,2));
        end
        if exist(['M_P4/GLO/' num2str(doy)],'dir')==0
            mkdir(['M_P4/GLO/' num2str(doy)]);
        end
        filenameP4=['M_P4/GLO/' num2str(doy) '/' site num2str(doy) 'P4.mat'];
        save(filenameP4,'GLOP4','-mat');      
    end
    %     clear GPSP4 GLOP4
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

%% cut glonass obs
if ~isnan(find(strcmp(fields, 'GLOP1' )))
    gloline=size(obs.GLOP1,1);
    glosate=size(POSItion.glox,2);
    if gloline<2880
        obs.GLOL1(gloline+1:2880,:)=0;obs.GLOL2(gloline+1:2880,:)=0;
        obs.GLOP1(gloline+1:2880,:)=0;obs.GLOP2(gloline+1:2880,:)=0;
    end
    glol=size(obs.GLOP1,2);
    if glol<glosate
        obs.GLOL1(:,glol+1:glosate)=0;obs.GLOL2(:,glol+1:glosate)=0;
        obs.GLOP1(:,glol+1:glosate)=0;obs.GLOP2(:,glol+1:glosate)=0;
    end
%     if glol>glosate
%         POSItion.glox(:,glosate+1:glol)=0;POSItion.gloy(:,glosate+1:glol)=0;
%         POSItion.gloz(:,glosate+1:glol)=0;
%     end
    %     glodelete=find(sate_mark.glo==0);
    glonum=size(POSItion.glox,2);
    if glonum<size(obs.GLOP1,2)
        %         if ~isnan(glodelete)
        %             k=length(glodelete);
        %             for i=k:-1:1
        %                 obs.GLOP1(:,glodelete(i))=[];
        %                 obs.GLOP2(:,glodelete(i))=[];
        %                 obs.GLOL1(:,glodelete(i))=[];
        %                 obs.GLOL2(:,glodelete(i))=[];
        %             end
        %         end
        obs.GLOP1=obs.GLOP1(:,1:glonum);
        obs.GLOP2=obs.GLOP2(:,1:glonum);
        obs.GLOL1=obs.GLOL1(:,1:glonum);
        obs.GLOL2=obs.GLOL2(:,1:glonum);
    end
    obs.GLOP1(isnan(obs.GLOP1))=0;
    obs.GLOP2(isnan(obs.GLOP2))=0;
    obs.GLOL1(isnan(obs.GLOL1))=0;
    obs.GLOL2(isnan(obs.GLOL2))=0;
    
    glox=POSItion.glox;gloy=POSItion.gloy;gloz=POSItion.gloz;
    
    for i=1:glonum
        for j=1:2880
            if obs.GLOL1(j,i)==0 || obs.GLOL2(j,i)==0 || obs.GLOP1(j,i)==0 || obs.GLOP2(j,i)==0
                obs.GLOL1(j,i)=0;obs.GLOL2(j,i)=0;obs.GLOP1(j,i)=0;obs.GLOP2(j,i)=0;
                continue;
            end
            [el,aaa]=Get_EA(sx,sy,sz,glox(j,i),gloy(j,i),gloz(j,i));
            if el<Cutoff(2,1)
                obs.GLOP1(j,i)=0;obs.GLOP2(j,i)=0;obs.GLOL1(j,i)=0;obs.GLOL2(j,i)=0;
                continue;
            end
        end
    end
end
end
%% ------------------------subfunction----------------------------
function [GLOPP4]=GLO_prepro(obs)
%% get smoothed GLONASS P4 observations
% INPUT:
%      obs: struct of rinex files
% OUTPUT:
%     GLOPP4: smoothed GLONASS P4 observation
%% ------------------------------------------------------------------
%____delete incomplete epoch_____
global c
size2=size(obs.GLOP1,2);
LenSatR=24;%LEN(3,1);
if size2<LenSatR
    obs.GLOP1(:,size2+1:LenSatR)=0;obs.GLOP2(:,size2+1:LenSatR)=0;
    obs.GLOL1(:,size2+1:LenSatR)=0;obs.GLOL2(:,size2+1:LenSatR)=0;
end
GLOP4=zeros(2880,LenSatR);
GLOL4=zeros(2880,LenSatR);
GLOPP4=zeros(2880,LenSatR);
for i=1:LenSatR
    %     Fre=[1,-4, 5, 1, 5, 6, -2, 0, -1, -2, -7, 0, -1, 4, -3, 3, 2, 4, -3, 3, 2];
    Fre=[1,-4, 5, 6, 1, -4, 5, 6, -2, -7, 0, -1, -2, -7, 0, -1, 4, -3, 3, 2, 4, -3, 3, 2];
%         Fre=[1 -4 5 6 1 -4 5 6 -2 -7 0 -1 -2 -7 0 -1 4 -3 3 2 4 -3 3 2 0 0];
    GLO_f1(i)=(1602+Fre(i)*0.5625)*10^6;                 %_________________________________unit:Hz
    GLO_f2(i)=(1246+Fre(i)*0.4375)*10^6;                  %_________________________________unit:Hz
    lamda_w(i)=c/(GLO_f1(i)-GLO_f2(i));       %____________________wide lane wavelength
    L6(:,i)=lamda_w(i)*(obs.GLOL1(:,i)-obs.GLOL2(:,i))-(GLO_f1(i)*obs.GLOP1(:,i)+GLO_f2(i)*obs.GLOP2(:,i))/(GLO_f1(i)+GLO_f2(i)); %__MW observable
    Li(:,i)=obs.GLOL1(:,i)-GLO_f1(i)*obs.GLOL2(:,i)/GLO_f2(i);
    Nw(:,i)=-L6(:,i);               %_____________________wide lane ambiguity
end
for i=1:LenSatR  %i is PRN number
    %------divide arc---------------------------
    arc=Get_arc(L6(:,i));
    [arc_n,aaa]=size(arc);
    %----delete arc less than 10 epoches-------
    arc_d=[];
    for j=1:arc_n
        n_epoch=arc(j,2)-arc(j,1);
        if n_epoch<20
            for k=arc(j,1):arc(j,2)
                obs.GLOP1(k,i)=0;obs.GLOP2(k,i)=0;obs.GLOL1(k,i)=0;obs.GLOL2(k,i)=0;
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
                L6(e,i)=0;obs.GLOL1(e,i)=0;obs.GLOL2(e,i)=0;obs.GLOP1(e,i)=0;obs.GLOP2(e,i)=0;Nw(e,i)=0;Li(e,i)=0;
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
                obs.GLOP1(k,i)=0;obs.GLOP2(k,i)=0;obs.GLOL1(k,i)=0;obs.GLOL2(k,i)=0;
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
                        L6(k+1,i)=0;obs.GLOP1(k+1,i)=0;obs.GLOP2(k+1,i)=0;
                        obs.GLOL1(k+1,i)=0;obs.GLOL2(k+1,i)=0;Nw(k+1,i)=0;
                        Li(k,i)=0;arc(j,2)=k;
                    else                     %------delete scatter epoches
                        for l=arc(j,1):k+1
                            L6(l,i)=0;obs.GLOP1(l,i)=0;obs.GLOP2(l,i)=0;
                            obs.GLOL1(l,i)=0;obs.GLOL2(l,i)=0;Nw(l,i)=0;
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
                            L6(l,i)=0;obs.GLOP1(l,i)=0;obs.GLOP2(l,i)=0;
                            obs.GLOL1(l,i)=0;obs.GLOL2(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,1)=k+1;
                        j=j-1;
                    end
                else                        %-----------------gross error
                    if k+1-arc(j,1)>20
                        L6(k+1,i)=0;obs.GLOP1(k+1,i)=0;obs.GLOP2(k+1,i)=0;
                        obs.GLOL1(k+1,i)=0;obs.GLOL2(k+1,i)=0;Nw(k+1,i)=0;Li(k,i)=0;
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+2,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else
                        for l=arc(j,1):k+1  %------delete scatter epoches
                            L6(l,i)=0;obs.GLOP1(l,i)=0;obs.GLOP2(l,i)=0;
                            obs.GLOL1(l,i)=0;obs.GLOL2(l,i)=0;Nw(l,i)=0;Li(k,i)=0;
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
    GLOP4(:,i)=obs.GLOP1(:,i)-obs.GLOP2(:,i);
    GLOL4(:,i)=(c/GLO_f1(i))*obs.GLOL1(:,i)-(c/GLO_f2(i))*obs.GLOL2(:,i);
    
    %--------------------------smoothing-------------------------
    GLOPP4(:,i)=GLOP4(:,i);
    for j=1:arc_n
        t=2;
        for k=arc(j,1)+1:arc(j,2)
            %P4(k,i)=mean(P4(arc(j,1):k,i))-L4(k,i)+mean(L4(arc(j,1):k,i));
            GLOPP4(k,i)=GLOPP4(k,i)/t+(GLOPP4(k-1,i)+GLOL4(k-1,i)-GLOL4(k,i))*(t-1)/t;
            t=t+1;
        end
        GLOPP4(arc(j,1):arc(j,1)+4,i)=0;
    end
    %--------remove bad P4---------------------
    arc=Get_arc(GLOPP4(:,i));
    [arc_n,aaa]=size(arc);
    for j=1:arc_n
        ave=mean(GLOPP4(arc(j,1):arc(j,2),i));
        if abs(ave)>10
            for kk=arc(j,1):arc(j,2)
                GLOPP4(kk,i)=0;
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