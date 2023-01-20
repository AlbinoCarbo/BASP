% Astrometric_filter: orbital best fit script for single satellite astrometric observations 
% using find_orb by Bill Gay. Delete observations that are more than DELTA arcsec away 
% from the best fit value..
%
% ALGORITMO:
% 1-Cancella le osservazioni con tempi uguali. Questa è una funzione per il
% filtro delle osservazioni di satelliti singoli: se i tempi sono ripetuti 
% vuol dire che ASTRiDE ha trovato più di una traccia nella stessa immagine,
% probabilmente perché la traccia del satellite è debole e l'ha divisa in
% più troconi oppure ha trovato una galassia o un raggio cosmico. I tempi
% ripetuti vanno cancellati, altrimenti introducono del rumore nella
% determinazione dell'orbita geocentrica. Questo problema si può attenuare
% aumentando il cut sul valore dell'area della traccia nel file SST_astride_TDM.py.
%
% 2-Cancella le osservazioni astrometriche che cadono dentro la fascia
% della Via lattea, dove le misure sono maggiormente alterate dal fitto fondo stellare.
% La parte di filtro per la Via Lattea invernale si può lasciare
% disattivato visto la minore densità di stelle rispetto alla Via Lattea
% estiva.
%
% 3-Conteggio del numero di satelliti diversi osservati.
%
% 4-Separazione dei dati astrometrici per ciascun satellite.
%
% 5-Filtro delle migliori osservazioni per i satelliti con traccia singola 
% con il best fit orbitale di find_orb.
%
% 6-Riporta le osservazioni astrometriche filtrate allo script principale.
%
% 7-Se le osservazioni astrometriche sono pari o inferiori a tre salta il
% best fit orbitale per l'eliminazione degli outliers.
%
% INPUT
% data_pathx = path delle immagini SST
% DELTA = massimo residuo di best fit orbitale delle osservazioni astrometriche (arcsec).
% DMW = numero per la cancellazione delle osservazioni astrometriche che
% cadono nella fascia della Via Lattea (se vale 1 le cancella, altrimenti no).
% YYYYx = vettore degli anni
% MMx = vettore dei mesi (in numero)
% DDx = vettore dei giorni
% hhx = vettore delle ore
% mmx = vettore dei minuti
% ssx = vettore dei secondi e decimali
% AR0x = vettore delle AR osservate al J2000 (gradi)
% DEC0x = vettore delle declinazioni osservate al J2000 (gradi)
% NORADx = vettore dei numeri Norad dei satelliti
%
% OUTPUT
% YYYY_ord = vettore degli anni
% MM_ord = vettore dei mesi (in numero)
% DD_ord = vettore dei giorni
% hi_ord = vettore delle ore
% mi_ord = vettore dei minuti
% si_ord = vettore dei secondi e decimali
% AR_ord = vettore delle best AR osservate al J2000 (gradi)
% DEC_ord = vettore delle best declinazioni osservate al J2000 (gradi)
% NORAD_ord = vettore dei numeri Norad dei satelliti
%
% Funzioni utilizzate:
% equat2galactic.m, funzione per la selezione delle sole posizioni astrometriche
% dei satelliti al di fuori della fascia della Via Lattea. 
%
% Fit_orb: funzione per il best fit dell'orbita di un satellite artificiale 
% eliminando le osservazioni astrometriche con un residuo maggiore di DELTA.
%
% Albino Carbognani, INAF-OAS
% Versione del 27 ottobre 2021

function [YYYY_ord, MM_ord, DD_ord, hi_ord, mi_ord, si_ord, AR_ord, DEC_ord, NORAD_ord]=Astrometric_filter(data_pathx, DELTA, DMW, YYYYx, MMx, DDx, hhx, mmx, ssx, AR0x, DEC0x, NORADx)

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%         Filter the satellites worst observations - single trace     %')
disp('%                             Oct 2021 version                        %') 
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('  ')

% Numero di osservazioni in input
Nobs0=length(AR0x);
disp(strcat('Astrometric filter: number of starting astrometric observations:', {' '}, num2str(Nobs0)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DELETE THE OBSERVATIONS WITH THE REPEATED TIMES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrice di tutte le osservazioni astrometriche di input
MO1=[YYYYx MMx DDx hhx mmx ssx AR0x DEC0x NORADx];
% format long   % ===> Per il controllo dell'output
% disp(MO1)     % ===> Per il controllo dell'output

% Calcolo vettore MJD (tempi crescenti)
MJDx=Mjday(YYYYx, MMx, DDx, hhx, mmx, ssx);

% Estrazione vettore MJD con i tempi unici
k=1;
for i=1:length(MJDx)  
   
    MJD_aux=abs(MJDx-MJDx(i));
    numero_zeri=sum(MJD_aux == 0); % Conteggio numero di zeri, se mumero_zeri=1 l'elemento è unico
    if numero_zeri == 1
        MJD_single(k)=MJDx(i); % Vettore MJD con i tempi unici
        k=k+1;
    end
    
end

% Salva le righe della matrice delle osservazioni aventi tempi osservativi
% che compaiono una sola volta
for i=1:length(MJD_single)
    for j=1:length(MJDx)
        if MJD_single(i)==MJDx(j)
         MO2(i,:) = MO1(j,:);
        end
    end
end

% disp(MO2) % ===> Per il controllo dell'output

Nobs1=length(MO2(:, 9));
deleted_astrometry=Nobs0-Nobs1;
disp(strcat('Astrometric filter: deleted', {' '}, num2str(deleted_astrometry), {' '},'astrometric observations with repeated time'))

% Istruzione if-else-end per cancellare eventualmente le osservazioni
% astrometriche che ricadono nella fascia della Via Lattea. Se DMW=1 le
% cancella, altrimenti le tiene tutte.

if DMW==1
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DELETE THE ASTROMETRIC OBSERVATIONS THAT FALL WITHIN THE MILKY WAY BELT %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('     Astrometric filter: deleting astrometry inside Milky Way')
disp('   ')

[YYYY, MM, DD, hh, mm, ss, AR0, DEC0, NORAD] = equat2galactic(MO2(:, 1), MO2(:, 2), MO2(:, 3), MO2(:, 4), MO2(:, 5), MO2(:, 6), MO2(:, 7), MO2(:, 8), MO2(:, 9));

% Numero delle osservazioni astrometriche con tempi diversi e fuori dalla Via Lattea
Nobs=length(NORAD);
disp(strcat('Astrometric filter: number of single astrometric observations:', {' '}, num2str(Nobs)))

else
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOLDS ALL ASTROMETRIC OBSERVATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
YYYY = MO2(:, 1);
MM = MO2(:, 2);
DD = MO2(:, 3);
hh = MO2(:, 4);
mm = MO2(:, 5);
ss = MO2(:, 6);
AR0 = MO2(:, 7);
DEC0 = MO2(:, 8);
NORAD = MO2(:, 9);

% Numero delle osservazioni astrometriche con tempi diversi
Nobs=length(NORAD);
disp('     Astrometric filter: Milky Way filter non-active')
disp('   ')
disp(strcat('Astrometric filter: number of single astrometric observations:', {' '}, num2str(Nobs)))

end

% Calcolo vettore MJD
MJD=Mjday(YYYY, MM, DD, hh, mm, ss);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COUNT THE NUMBER OF DIFFERENT SATELLITES OBSERVED %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=1;
SAT(n)=NORAD(1);
for i=1:length(NORAD)
    if abs(SAT(n)-NORAD(i))>0
        SAT(n+1)=NORAD(i);    % Vettore che contiene il numero NORAD dei satelliti non ripetuti
        n=n+1;
    end
end

Nsat=length(SAT); % Numero di satelliti diversi
disp(strcat('Astrometric filter: number of different satellites to analize:', {' '}, num2str(n)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ciclo per la separazione dei dati astrometrici per ciascun satellite
% usando il vettore che contiene il numero NORAD dei satelliti non ripetuti
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:Nsat
    k=1;
    for j=1:Nobs
        if abs(SAT(i)-NORAD(j))==0
        ARsat(k)=AR0(j);
        DECsat(k)=DEC0(j);
        MJDsat(k)=MJD(j);
        YYYYsat(k)=YYYY(j);
        MMsat(k)=MM(j);
        DDsat(k)=DD(j);
        hhsat(k)=hh(j);
        mmsat(k)=mm(j);
        sssat(k)=ss(j);
        NORADsat(k)=NORAD(j);
        k=k+1;
        end
    end
    
    AR_sat{i}=ARsat;     % Vettore AR osservata i-esimo satellite
    DEC_sat{i}=DECsat;   % Vettore DEC osservata i-esimo satellite
    MJD_sat{i}=MJDsat;   % Vettore MJD osservato i-esimo satellite
    YYYY_sat{i}=YYYYsat;
    MM_sat{i}=MMsat;
    DD_sat{i}=DDsat;
    hh_sat{i}=hhsat;
    mm_sat{i}=mmsat;
    ss_sat{i}=sssat;
    NORAD_sat{i}=NORADsat;
    
    % Clear dei vettori temporanei
    clear ARsat DECsat MJDsat YYYYsat MMsat DDsat hhsat mmsat sssat NORADsat
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtro delle migliori osservazioni per i satelliti con traccia singola %
% con il best fit orbitale con find_orb                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:Nsat
    obs=length(MJD_sat{i});            % Numero osservazioni i-esimo satellite
    name_sat=num2str(NORAD_sat{i}(1)); % Numero Norad i-esimo satellite
    
    if obs >= 3 % Numero minimo di osservazioni necessarie per il best fit orbitale
        
        % Selezione coordinate di best fit orbitale usando find_orb
        disp(strcat('Astrometric filter: select BEST FIT ORBIT astrometry for NORAD', {' '}, name_sat))
        disp('   ')
        [YYYY3, MM3, DD3, hh3, mm3, ss3, AR3, DEC3, S]=Fit_orb(data_pathx, DELTA, name_sat, YYYY_sat{i}, MM_sat{i}, DD_sat{i}, hh_sat{i}, mm_sat{i}, ss_sat{i}, AR_sat{i}, DEC_sat{i});
        disp('   ')
        disp(strcat('Astrometric filter: angle traveled by the satellite', {' '}, name_sat, {' '},'(degree):', {' '}, num2str(S)))
        disp('   ')
        
        % Calcolo vettore MJD delle best osservazioni del satellite
        MJD3=Mjday(YYYY3, MM3, DD3, hh3, mm3, ss3);
        
        % Costruzione vettore con il numero Norad del satellite 
        % (ogni osservazione astrometrica deve avere il numero Norad 
        % del satellite cui si riferisce)
        NORAD3=zeros(length(1),length(AR3)) + str2double(name_sat);
        
        AR_ord1{i}=AR3;       % Vettore best AR osservata i-esimo satellite
        DEC_ord1{i}=DEC3;     % Vettore best DEC osservata i-esimo satellite
        MJD_ord1{i}=MJD3;     % Vettore best MJD osservato i-esimo satellite
        YYYY_ord1{i}=YYYY3;
        MM_ord1{i}=MM3;
        DD_ord1{i}=DD3;
        hi_ord1{i}=hh3;
        mi_ord1{i}=mm3;
        si_ord1{i}=ss3;
        NORAD_ord1{i}=NORAD3;
        
        % Clear dei vettori temporanei per l'uso con il satellite
        % successivo
        clear AR3 DEC3 MJD3 YYYY3 MM3 DD3 hh3 mm3 ss3 NORAD3
    else
        % Se le osservazioni non sono sufficienti per il best fit orbitale
        % si lasciano come sono
        disp(strcat('Astrometric filter: satellite', {' '}, name_sat, {' '},'NOT FILTERED with find_orb, TOO FEW observations'))
        disp('   ')
        
        AR_ord1{i}=AR_sat{i};       % Vettore AR osservata i-esimo satellite
        DEC_ord1{i}=DEC_sat{i};     % Vettore DEC osservata i-esimo satellite
        MJD_ord1{i}=MJD_sat{i};     % Vettore MJD osservato i-esimo satellite
        YYYY_ord1{i}=YYYY_sat{i};
        MM_ord1{i}=MM_sat{i};
        DD_ord1{i}=DD_sat{i};
        hi_ord1{i}=hh_sat{i};
        mi_ord1{i}=mm_sat{i};
        si_ord1{i}=ss_sat{i};
        NORAD_ord1{i}=NORAD_sat{i}; 
    end
    
    
end

% Concatenazione cell array per l'output
AR_ord=[AR_ord1{:}];
DEC_ord=[DEC_ord1{:}];
YYYY_ord=[YYYY_ord1{:}];
MM_ord=[MM_ord1{:}];
DD_ord=[DD_ord1{:}];
hi_ord=[hi_ord1{:}];
mi_ord=[mi_ord1{:}];
si_ord=[si_ord1{:}];
NORAD_ord=[NORAD_ord1{:}];

        disp(strcat('Astrometric filter: filtering single', {' '}, 'traces satellites completed'))
        disp('   ')
end
