function [Alpha,Beta,F,LAMBEDA1,LAMBEDA2]=info_FreWaveAlphaBeta(selct_system,LengthSAT)
global c

%% GRC
if strcmp(selct_system,'GRC')==1
    
    for i=1:LengthSAT(1)
        freq_gps (i,1)   = 10.23*10^6*154;     % GPS L1 frequency [Hz]
        freq_gps (i,2)   = 10.23*10^6*120;     % GPS L2 frequency [Hz]
        wavl_gps (i,1)   = c/(10.23*10^6*154); % Corresponding wavelenghts L1 [m]
        wavl_gps (i,2)   = c/(10.23*10^6*120); % Corresponding wavelenghts L2 [m]
        alpha_gps(i,1)   =  ((freq_gps(i,1)^2)/(freq_gps(i,1)^2 - freq_gps(i,2)^2)); % Ion-free combinations coefficient alpha
        beta_gps (i,1)   = -((freq_gps(i,2)^2)/(freq_gps(i,1)^2 - freq_gps(i,2)^2)); % Ion-free combinations coefficient beta
    end
    alpha_gps   = alpha_gps(1,1);
    beta_gps    = beta_gps(1,1);
    f1_gps      = freq_gps(1,1);
    f2_gps      = freq_gps(1,2);
    Lambda1_gps = wavl_gps(1,1);
    Lambda2_gps = wavl_gps(1,2);
    
    glok = [1 -4 5 6 1 -4 5 6 -2 -7 0 -1 -2 -7 0 -1 4 -3 3 2 4 3 2 -6];
    for i=1:LengthSAT(2)
        freq_glonass(i,1) = (1602 + 0.5625*glok(i))*10^6;     % GLONASS G1 frequency [Hz]
        freq_glonass(i,2) = (1246 + 0.4375*glok(i))*10^6;     % GLONASS G2 frequency [Hz]
        wavl_glonass(i,1) = c/((1602 + 0.5625*glok(i))*10^6); % Corresponding wavelenghts G1 [m]
        wavl_glonass(i,2) = c/((1246 + 0.4375*glok(i))*10^6); % Corresponding wavelenghts G2 [m]
        alpha_glonass(i,1)   =  ((freq_glonass(i,1)^2)/(freq_glonass(i,1)^2 - freq_glonass(i,2)^2)); % Ion-free combinations coefficient alpha
        beta_glonass (i,1)   = -((freq_glonass(i,2)^2)/(freq_glonass(i,1)^2 - freq_glonass(i,2)^2)); % Ion-free combinations coefficient beta
    end
    alpha_glonass   = alpha_glonass(1,1);
    beta_glonass    = beta_glonass(1,1);
    f1_glonass      = freq_glonass(:,1);
    f2_glonass      = freq_glonass(:,2);
    Lambda1_glonass = wavl_glonass(:,1);
    Lambda2_glonass = wavl_glonass(:,2);
    
    for i=1:LengthSAT(3)
        freq_beidou (i,1)   = 10.23*10^6*152.6;     % BeiDou C2 frequency [Hz]
        freq_beidou (i,2)   = 10.23*10^6*118;       % BeiDou C7 frequency [Hz]
        wavl_beidou (i,1)   = c/(10.23*10^6*152.6); % Corresponding wavelenghts E1  [m]
        wavl_beidou (i,2)   = c/(10.23*10^6*118);   % Corresponding wavelenghts E5a [m]
        alpha_beidou(i,1)   =  ((freq_beidou(i,1)^2)/(freq_beidou(i,1)^2 - freq_beidou(i,2)^2)); % Ion-free combinations coefficient alpha
        beta_beidou (i,1)   = -((freq_beidou(i,2)^2)/(freq_beidou(i,1)^2 - freq_beidou(i,2)^2)); % Ion-free combinations coefficient beta
    end
    alpha_beidou   = alpha_beidou(1,1);
    beta_beidou    = beta_beidou(1,1);
    f1_beidou      = freq_beidou(1,1);
    f2_beidou      = freq_beidou(1,2);
    Lambda1_beidou = wavl_beidou(1,1);
    Lambda2_beidou = wavl_beidou(1,2);
    
    Alpha.gps=alpha_gps;
    Alpha.glonass=alpha_glonass;
    Alpha.beidou=alpha_beidou;
    
    Beta.gps=beta_gps;
    Beta.glonass=beta_glonass;
    Beta.beidou=beta_beidou;
    
    f1_gps=repmat(f1_gps,LengthSAT(1),1);f2_gps=repmat(f2_gps,LengthSAT(1),1);
    f1_beidou=repmat(f1_beidou,LengthSAT(3),1);f2_beidou=repmat(f2_beidou,LengthSAT(3),1);
    
    F1=[f1_gps;f1_glonass;f1_beidou];
    F2=[f2_gps;f2_glonass;f2_beidou];
    F=[F1,F2];
    
    Lambda1_gps=repmat(Lambda1_gps,LengthSAT(1),1);Lambda2_gps=repmat(Lambda2_gps,LengthSAT(1),1);
    Lambda1_beidou=repmat(Lambda1_beidou,LengthSAT(3),1);Lambda2_beidou=repmat(Lambda2_beidou,LengthSAT(3),1);
    
    LAMBEDA1.gps=Lambda1_gps;
    LAMBEDA1.glonass=Lambda1_glonass;
    LAMBEDA1.beidou=Lambda1_beidou;
    LAMBEDA2.gps=Lambda2_gps;
    LAMBEDA2.glonass=Lambda2_glonass;
    LAMBEDA2.beidou=Lambda2_beidou;
    
    %% GR
elseif strcmp(selct_system,'GR')==1
    
    for i=1:LengthSAT(1)
        freq_gps (i,1)   = 10.23*10^6*154;     % GPS L1 frequency [Hz]
        freq_gps (i,2)   = 10.23*10^6*120;     % GPS L2 frequency [Hz]
        wavl_gps (i,1)   = c/(10.23*10^6*154); % Corresponding wavelenghts L1 [m]
        wavl_gps (i,2)   = c/(10.23*10^6*120); % Corresponding wavelenghts L2 [m]
        alpha_gps(i,1)   =  ((freq_gps(i,1)^2)/(freq_gps(i,1)^2 - freq_gps(i,2)^2)); % Ion-free combinations coefficient alpha
        beta_gps (i,1)   = -((freq_gps(i,2)^2)/(freq_gps(i,1)^2 - freq_gps(i,2)^2)); % Ion-free combinations coefficient beta
    end
    alpha_gps   = alpha_gps(1,1);
    beta_gps    = beta_gps(1,1);
    f1_gps      = freq_gps(1,1);
    f2_gps      = freq_gps(1,2);
    Lambda1_gps = wavl_gps(1,1);
    Lambda2_gps = wavl_gps(1,2);
    
    glok = [1 -4 5 6 1 -4 5 6 -2 -7 0 -1 -2 -7 0 -1 4 -3 3 2 4 -3 3 2 0 0];
    for i=1:LengthSAT(2)
        freq_glonass(i,1) = (1602 + 0.5625*glok(i))*10^6;     % GLONASS G1 frequency [Hz]
        freq_glonass(i,2) = (1246 + 0.4375*glok(i))*10^6;     % GLONASS G2 frequency [Hz]
        wavl_glonass(i,1) = c/((1602 + 0.5625*glok(i))*10^6); % Corresponding wavelenghts G1 [m]
        wavl_glonass(i,2) = c/((1246 + 0.4375*glok(i))*10^6); % Corresponding wavelenghts G2 [m]
        alpha_glonass(i,1)   =  ((freq_glonass(i,1)^2)/(freq_glonass(i,1)^2 - freq_glonass(i,2)^2)); % Ion-free combinations coefficient alpha
        beta_glonass (i,1)   = -((freq_glonass(i,2)^2)/(freq_glonass(i,1)^2 - freq_glonass(i,2)^2)); % Ion-free combinations coefficient beta
    end
    alpha_glonass   = alpha_glonass(1,1);
    beta_glonass    = beta_glonass(1,1);
    f1_glonass      = freq_glonass(:,1);
    f2_glonass      = freq_glonass(:,2);
    Lambda1_glonass = wavl_glonass(:,1);
    Lambda2_glonass = wavl_glonass(:,2);
    

    Alpha.gps=alpha_gps;
    Alpha.glonass=alpha_glonass;
    
    Beta.gps=beta_gps;
    Beta.glonass=beta_glonass;
    
    f1_gps=repmat(f1_gps,LengthSAT(1),1);f2_gps=repmat(f2_gps,LengthSAT(1),1);
    
    F1=[f1_gps;f1_glonass];
    F2=[f2_gps;f2_glonass];
    F=[F1,F2];
    
    Lambda1_gps=repmat(Lambda1_gps,LengthSAT(1),1);Lambda2_gps=repmat(Lambda2_gps,LengthSAT(1),1);
    
    LAMBEDA1.gps=Lambda1_gps;
    LAMBEDA1.glonass=Lambda1_glonass;
    LAMBEDA2.gps=Lambda2_gps;
    LAMBEDA2.glonass=Lambda2_glonass;

    %% GC
elseif strcmp(selct_system,'GC')==1
    
    for i=1:LengthSAT(1)
        freq_gps (i,1)   = 10.23*10^6*154;     % GPS L1 frequency [Hz]
        freq_gps (i,2)   = 10.23*10^6*120;     % GPS L2 frequency [Hz]
        wavl_gps (i,1)   = c/(10.23*10^6*154); % Corresponding wavelenghts L1 [m]
        wavl_gps (i,2)   = c/(10.23*10^6*120); % Corresponding wavelenghts L2 [m]
        alpha_gps(i,1)   =  ((freq_gps(i,1)^2)/(freq_gps(i,1)^2 - freq_gps(i,2)^2)); % Ion-free combinations coefficient alpha
        beta_gps (i,1)   = -((freq_gps(i,2)^2)/(freq_gps(i,1)^2 - freq_gps(i,2)^2)); % Ion-free combinations coefficient beta
    end
    alpha_gps   = alpha_gps(1,1);
    beta_gps    = beta_gps(1,1);
    f1_gps      = freq_gps(1,1);
    f2_gps      = freq_gps(1,2);
    Lambda1_gps = wavl_gps(1,1);
    Lambda2_gps = wavl_gps(1,2);
        
    for i=1:LengthSAT(2)
        freq_beidou (i,1)   = 10.23*10^6*152.6;     % BeiDou C2 frequency [Hz]
        freq_beidou (i,2)   = 10.23*10^6*118;       % BeiDou C7 frequency [Hz]
        wavl_beidou (i,1)   = c/(10.23*10^6*152.6); % Corresponding wavelenghts E1  [m]
        wavl_beidou (i,2)   = c/(10.23*10^6*118);   % Corresponding wavelenghts E5a [m]
        alpha_beidou(i,1)   =  ((freq_beidou(i,1)^2)/(freq_beidou(i,1)^2 - freq_beidou(i,2)^2)); % Ion-free combinations coefficient alpha
        beta_beidou (i,1)   = -((freq_beidou(i,2)^2)/(freq_beidou(i,1)^2 - freq_beidou(i,2)^2)); % Ion-free combinations coefficient beta
    end
    alpha_beidou   = alpha_beidou(1,1);
    beta_beidou    = beta_beidou(1,1);
    f1_beidou      = freq_beidou(1,1);
    f2_beidou      = freq_beidou(1,2);
    Lambda1_beidou = wavl_beidou(1,1);
    Lambda2_beidou = wavl_beidou(1,2);
    
    Alpha.gps=alpha_gps;
    Alpha.beidou=alpha_beidou;
    
    Beta.gps=beta_gps;
    Beta.beidou=beta_beidou;
    
    f1_gps=repmat(f1_gps,LengthSAT(1),1);f2_gps=repmat(f2_gps,LengthSAT(1),1);
    f1_beidou=repmat(f1_beidou,LengthSAT(2),1);f2_beidou=repmat(f2_beidou,LengthSAT(2),1);
    
    F1=[f1_gps;f1_beidou];
    F2=[f2_gps;f2_beidou];
    F=[F1,F2];
    
    Lambda1_gps=repmat(Lambda1_gps,LengthSAT(1),1);Lambda2_gps=repmat(Lambda2_gps,LengthSAT(1),1);
    Lambda1_beidou=repmat(Lambda1_beidou,LengthSAT(2),1);Lambda2_beidou=repmat(Lambda2_beidou,LengthSAT(2),1);
    
    LAMBEDA1.gps=Lambda1_gps;
    LAMBEDA1.beidou=Lambda1_beidou;
    LAMBEDA2.gps=Lambda2_gps;
    LAMBEDA2.beidou=Lambda2_beidou;
   %% RC
elseif strcmp(selct_system,'RC')==1
    
    glok = [1 -4 5 6 1 -4 5 6 -2 -7 0 -1 -2 -7 0 -1 4 -3 3 2 4 -3 3 2 0 0];
    for i=1:LengthSAT(1)
        freq_glonass(i,1) = (1602 + 0.5625*glok(i))*10^6;     % GLONASS G1 frequency [Hz]
        freq_glonass(i,2) = (1246 + 0.4375*glok(i))*10^6;     % GLONASS G2 frequency [Hz]
        wavl_glonass(i,1) = c/((1602 + 0.5625*glok(i))*10^6); % Corresponding wavelenghts G1 [m]
        wavl_glonass(i,2) = c/((1246 + 0.4375*glok(i))*10^6); % Corresponding wavelenghts G2 [m]
        alpha_glonass(i,1)   =  ((freq_glonass(i,1)^2)/(freq_glonass(i,1)^2 - freq_glonass(i,2)^2)); % Ion-free combinations coefficient alpha
        beta_glonass (i,1)   = -((freq_glonass(i,2)^2)/(freq_glonass(i,1)^2 - freq_glonass(i,2)^2)); % Ion-free combinations coefficient beta
    end
    alpha_glonass   = alpha_glonass(1,1);
    beta_glonass    = beta_glonass(1,1);
    f1_glonass      = freq_glonass(:,1);
    f2_glonass      = freq_glonass(:,2);
    Lambda1_glonass = wavl_glonass(:,1);
    Lambda2_glonass = wavl_glonass(:,2);
    
    for i=1:LengthSAT(2)
        freq_beidou (i,1)   = 10.23*10^6*152.6;     % BeiDou C2 frequency [Hz]
        freq_beidou (i,2)   = 10.23*10^6*118;       % BeiDou C7 frequency [Hz]
        wavl_beidou (i,1)   = c/(10.23*10^6*152.6); % Corresponding wavelenghts E1  [m]
        wavl_beidou (i,2)   = c/(10.23*10^6*118);   % Corresponding wavelenghts E5a [m]
        alpha_beidou(i,1)   =  ((freq_beidou(i,1)^2)/(freq_beidou(i,1)^2 - freq_beidou(i,2)^2)); % Ion-free combinations coefficient alpha
        beta_beidou (i,1)   = -((freq_beidou(i,2)^2)/(freq_beidou(i,1)^2 - freq_beidou(i,2)^2)); % Ion-free combinations coefficient beta
    end
    alpha_beidou   = alpha_beidou(1,1);
    beta_beidou    = beta_beidou(1,1);
    f1_beidou      = freq_beidou(1,1);
    f2_beidou      = freq_beidou(1,2);
    Lambda1_beidou = wavl_beidou(1,1);
    Lambda2_beidou = wavl_beidou(1,2);
    
    Alpha.glonass=alpha_glonass;
    Alpha.beidou=alpha_beidou;
    
    Beta.glonass=beta_glonass;
    Beta.beidou=beta_beidou;
    
    f1_beidou=repmat(f1_beidou,LengthSAT(2),1);f2_beidou=repmat(f2_beidou,LengthSAT(2),1);
    
    F1=[f1_glonass;f1_beidou];
    F2=[f2_glonass;f2_beidou];
    F=[F1,F2];
    
    Lambda1_beidou=repmat(Lambda1_beidou,LengthSAT(1),1);Lambda2_beidou=repmat(Lambda2_beidou,LengthSAT(1),1);
    
    LAMBEDA1.glonass=Lambda1_glonass;
    LAMBEDA1.beidou=Lambda1_beidou;
    LAMBEDA2.glonass=Lambda2_glonass;
    LAMBEDA2.beidou=Lambda2_beidou;
     
    %% G
elseif strcmp(selct_system,'G')==1
    
    for i=1:LengthSAT(1)
        freq_gps (i,1)   = 10.23*10^6*154;     % GPS L1 frequency [Hz]
        freq_gps (i,2)   = 10.23*10^6*120;     % GPS L2 frequency [Hz]
        wavl_gps (i,1)   = c/(10.23*10^6*154); % Corresponding wavelenghts L1 [m]
        wavl_gps (i,2)   = c/(10.23*10^6*120); % Corresponding wavelenghts L2 [m]
        alpha_gps(i,1)   =  ((freq_gps(i,1)^2)/(freq_gps(i,1)^2 - freq_gps(i,2)^2)); % Ion-free combinations coefficient alpha
        beta_gps (i,1)   = -((freq_gps(i,2)^2)/(freq_gps(i,1)^2 - freq_gps(i,2)^2)); % Ion-free combinations coefficient beta
    end
    alpha_gps   = alpha_gps(1,1);
    beta_gps    = beta_gps(1,1);
    f1_gps      = freq_gps(1,1);
    f2_gps      = freq_gps(1,2);
    Lambda1_gps = wavl_gps(1,1);
    Lambda2_gps = wavl_gps(1,2);
    
    Alpha.gps=alpha_gps;
    
    Beta.gps=beta_gps;
    
    f1_gps=repmat(f1_gps,LengthSAT(1),1);f2_gps=repmat(f2_gps,LengthSAT(1),1);
    
    F1=[f1_gps];
    F2=[f2_gps];
    F=[F1,F2];
    
    Lambda1_gps=repmat(Lambda1_gps,LengthSAT(1),1);Lambda2_gps=repmat(Lambda2_gps,LengthSAT(1),1);
    
    LAMBEDA1.gps=Lambda1_gps;
    LAMBEDA2.gps=Lambda2_gps;

  %% R
elseif strcmp(selct_system,'R')==1
    
    glok = [1 -4 5 6 1 -4 5 6 -2 -7 0 -1 -2 -7 0 -1 4 -3 3 2 4 -3 3 2 0 0];
    for i=1:LengthSAT(1)
        freq_glonass(i,1) = (1602 + 0.5625*glok(i))*10^6;     % GLONASS G1 frequency [Hz]
        freq_glonass(i,2) = (1246 + 0.4375*glok(i))*10^6;     % GLONASS G2 frequency [Hz]
        wavl_glonass(i,1) = c/((1602 + 0.5625*glok(i))*10^6); % Corresponding wavelenghts G1 [m]
        wavl_glonass(i,2) = c/((1246 + 0.4375*glok(i))*10^6); % Corresponding wavelenghts G2 [m]
        alpha_glonass(i,1)   =  ((freq_glonass(i,1)^2)/(freq_glonass(i,1)^2 - freq_glonass(i,2)^2)); % Ion-free combinations coefficient alpha
        beta_glonass (i,1)   = -((freq_glonass(i,2)^2)/(freq_glonass(i,1)^2 - freq_glonass(i,2)^2)); % Ion-free combinations coefficient beta
    end
    alpha_glonass   = alpha_glonass(1,1);
    beta_glonass    = beta_glonass(1,1);
    f1_glonass      = freq_glonass(:,1);
    f2_glonass      = freq_glonass(:,2);
    Lambda1_glonass = wavl_glonass(:,1);
    Lambda2_glonass = wavl_glonass(:,2);

    Alpha.glonass=alpha_glonass;
    Beta.glonass=beta_glonass;

    F1=[f1_glonass];
    F2=[f2_glonass];
    F=[F1,F2];
    
    LAMBEDA1.glonass=Lambda1_glonass;
    LAMBEDA2.glonass=Lambda2_glonass;
 
    %% C
elseif strcmp(selct_system,'C')==1
    
    for i=1:LengthSAT(1)
        freq_beidou (i,1)   = 10.23*10^6*152.6;     % BeiDou C2 frequency [Hz]
        freq_beidou (i,2)   = 10.23*10^6*118;       % BeiDou C7 frequency [Hz]
        wavl_beidou (i,1)   = c/(10.23*10^6*152.6); % Corresponding wavelenghts E1  [m]
        wavl_beidou (i,2)   = c/(10.23*10^6*118);   % Corresponding wavelenghts E5a [m]
        alpha_beidou(i,1)   =  ((freq_beidou(i,1)^2)/(freq_beidou(i,1)^2 - freq_beidou(i,2)^2)); % Ion-free combinations coefficient alpha
        beta_beidou (i,1)   = -((freq_beidou(i,2)^2)/(freq_beidou(i,1)^2 - freq_beidou(i,2)^2)); % Ion-free combinations coefficient beta
    end
    alpha_beidou   = alpha_beidou(1,1);
    beta_beidou    = beta_beidou(1,1);
    f1_beidou      = freq_beidou(1,1);
    f2_beidou      = freq_beidou(1,2);
    Lambda1_beidou = wavl_beidou(1,1);
    Lambda2_beidou = wavl_beidou(1,2);

    Alpha.beidou=alpha_beidou;
    Beta.beidou=beta_beidou;
    
    f1_beidou=repmat(f1_beidou,LengthSAT(1),1);f2_beidou=repmat(f2_beidou,LengthSAT(1),1);
    
    F1=[f1_beidou];
    F2=[f2_beidou];
    F=[F1,F2];
    
    Lambda1_beidou=repmat(Lambda1_beidou,LengthSAT(1),1);Lambda2_beidou=repmat(Lambda2_beidou,LengthSAT(1),1);
    LAMBEDA1.beidou=Lambda1_beidou;
    LAMBEDA2.beidou=Lambda2_beidou;
    
end
end