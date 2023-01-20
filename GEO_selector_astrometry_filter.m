% Script for choosing a GEO satellite that is among a group of other GEOs in the 
% same image using the SGP4 model with the current TLE from celestrak.com. 
% After identification of the desired GEO satellite, observations can be filtered 
% if they fall within the Milky Way and analyzed with Fit_orb.m for best orbit fit 
% and to delete the worst astrometric measurements.
%
% ALGORITMO:
% Una volta calcolate le posizioni topocentriche (ossia per 598) teoriche del 
% satellite con un dato numero NORAD di riferimento, lo script considera solo 
% le posizioni astrometriche che sono più vicine a quelle calcolate per il satellite 
% ed elimina tutte le altre osservazioni astrometriche. Per eliminare le immagini dove
% non tutto il gruppo di satelliti è stato rilevato (e fra questi ci potrebbe essere 
% il satellite che ci interessa) fa una statistica del numero di satelliti per ciascuna 
% immagine e prende come numero di satelliti quello che compare con una frequenza maggiore.
%
% Se non esiste il TLE online del satellite di cui si vuole l'effemeride, 
% la funzione impone che i valori di RA e DEC della effemeride siano quelli 
% osservati e passa al satellite successivo della lista. In questo caso 
% tutte le osservazioni passano senza essere filtrate.
%
% Le osservazioni che si riferiscono a un dato satellite possono essere
% filtrate per le posizioni che cadono nella Via Lattea, infine vengono
% mandate allo script Fit_orb.m, per il best fit orbitale e l'eliminazione
% delle posizioni astrometriche con i residui maggiori di DELTA arcsec.
%
% INPUT
% data_path = path dove si trova la funzione, ci salva il file del TLE corrente
% Lat_O = latitudine osservatore (gradi)
% Long_O = longitudine osservatore (gradi)
% h_O = quota osservatore (metri) 
% DELTA = massimo residuo di best fit orbitale delle osservazioni astrometriche (arcsec).
% DMW = numero per la cancellazione delle osservazioni astrometriche che
% cadono nella fascia della Via Lattea (se vale 1 le cancella, altrimenti no).
% YYYY = vettore degli anni
% MM = vettore dei mesi (in numero)
% DD = vettore dei giorni
% hh = vettore delle ore
% mm = vettore dei minuti
% ss = vettore dei secondi e decimali
% AR0 = vettore delle AR osservate al J2000 (gradi)
% DEC0 = vettore delle declinazioni osservate al J2000 (gradi)
% NORAD = vettore dei numeri Norad dei satelliti
%
% OUTPUT
% YYYY_ord = vettore degli anni
% MM_ord = vettore dei mesi (in numero)
% DD_ord = vettore dei giorni
% hi_ord = vettore delle ore
% mi_ord = vettore dei minuti
% si_ord = vettore dei secondi e decimali
% AR_ord = vettore delle AR osservate al J2000 entro 10" dalle effemeridi (gradi)
% DEC_ord = vettore delle declinazioni osservate al J2000 entro 10" dalle effemeridi (gradi)
% NORAD_ord = vettore dei numeri Norad dei satelliti
%
% Funzioni utilizzate:
% SGP4_Ephemeris.m, calcola AR e DEC topocentrici alla data del satellite (gradi)
%
% EQ2EQ2000.m, trasforma le coordinate equatoriali alla data al J2000 (radianti).
% 
% equat2galactic.m, funzione per la selezione delle sole posizioni astrometriche
% dei satelliti al di fuori della fascia della Via Lattea.
%
% Fit_orb.m, funzione per il best fit dell'orbita di un satellite artificiale 
% eliminando le osservazioni astrometriche con un residuo maggiore di DELTA.
%
% Albino Carbognani, OAS
% Version 21 ottobre 2021

function [YYYY_ord, MM_ord, DD_ord, hi_ord, mi_ord, si_ord, AR_ord, DEC_ord, NORAD_ord]=GEO_selector_astrometry_filter(data_path, Long_O, Lat_O, h_O, DELTA, DMW, YYYY, MM, DD, hh, mm, ss, AR0, DEC0, NORAD)

C=180/pi; % Trasformazione gradi <---> radianti

% File temporaneo per il salvataggio del TLE corrente
temp_tle=strcat(data_path, 'tle.txt');     

% Allocazione vettori delle effemeridi dei satelliti osservati
AR_eph=zeros(1, length(NORAD));
DEC_eph=zeros(1, length(NORAD));

% Calcolo vettore MJD
MJD=Mjday(YYYY, MM, DD, hh, mm, ss); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TLE download cycle and daily ephemeris calculation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('                   GEO selector and astrometry filter                  ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('  ')

disp(strcat('==> GEO selector: total astrometric observations', {' '}, num2str(length(AR0))))

%%% Calcolo effemeride per la prima posizione del primo satellite della lista

% Numero NORAD del satellite di cui si vogliono i TLE
Stringa_NORAD=num2str(NORAD(1));
        
% Acquisizione TLE da celestrak.com
fullURL = strcat('https://celestrak.com/satcat/tle.php?CATNR=', Stringa_NORAD);
disp('           ')
disp(strcat('==> GEO selector: TLE acquisition for satellite', {' '}, Stringa_NORAD))
        
% Stringa html con i TLE del satellite con numero NORAD(i)
str = webread(fullURL);

% Se non esiste il TLE online cerca un file locale con il TLE in
% data_path avente il nome "TLE_NORAD.txt" e lo copia nel file temporaneo "tle.txt"
% Se il file locale con il TLE non c'è mantiene i valori di AR0 e DEC0 e passa 
% al satellite successivo della lista.
        
tfx = strcmp(str,'No TLE found'); % Confronto fra stringhe

if tfx == 1
     disp(strcat('GEO selector, warning: no TLE online for', " ", Stringa_NORAD, " ", 'satellite'))
     disp('   ')
            ExTLE=strcat(data_path, '/TLE_', Stringa_NORAD, '.txt');
            file_mancanti=isfile(ExTLE); % Variabile logica, vale 1 se il file esiste, 0 altrimenti
            if file_mancanti==1
                answer=strcat('File', " ",'TLE_', Stringa_NORAD, '.txt', " ", 'exist');
                disp(answer) 
                disp('   ')
                copyfile(ExTLE, temp_tle) % Copia il TLE nel file temporaneo "tle.txt"
            else
                answer=strcat('File', " ",'TLE_', Stringa_NORAD, '.txt', " ", 'does not exist');
                disp(answer)
                disp('   ')
                AR_eph(1)=AR0(1); 
                DEC_eph(1)=DEC0(1);
                disp('GEO selector: skip to next satellite')
                disp('   ')
            end   
     
else
    
    % Salvataggio del TLE corrente nel file "tle.txt"  
    fid2 = fopen(temp_tle,'w');
    fprintf(fid2, '%s \n', str);
    fclose(fid2);
      
    disp(strcat('GEO selector: compute ephemeris (first position) for satellite', {' '}, Stringa_NORAD))
        
    % Funzione per il calcolo dell'effemeride del satellite nel
    % modello SGP4 in gradi per la prima posizione (equinozio alla data)
            
    [AR_eph(1), DEC_eph(1)]=SGP4_Ephemeris(temp_tle, Lat_O, Long_O, h_O, MJD(1));
            
end    
   
% Calcolo effemeridi per i satelliti rimanenti (i > 1). Se il satellite resta sempre lo stesso
% viene usato il TLE corrente senza scaricarne uno nuovo
            
for i=2:length(NORAD) % Ciclo per il calcolo delle effemeridi dei rimanenti satelliti osservati
        
    if abs(NORAD(i-1)-NORAD(i)) > 0
        % Numero NORAD del satellite di cui si vogliono i TLE
        Stringa_NORAD=num2str(NORAD(i));
        
        % Acquisizione TLE da celestrak.com
        fullURL = strcat('https://celestrak.com/satcat/tle.php?CATNR=', Stringa_NORAD);
        disp('           ')
        disp(strcat('==> GEO selector: TLE acquisition for satellite', {' '}, Stringa_NORAD))
        
        % Stringa html con i TLE del satellite con numero NORAD(i)
        str = webread(fullURL);
        
        % Se non esiste il TLE online cerca un file locale con il TLE in
        % data_path avente il nome "TLE_NORAD.txt" e lo copia nel file temporaneo "tle.txt"
        % Se il file locale con il TLE non c'è mantiene i valori di AR0 e DEC0 e passa 
        % al satellite successivo della lista.       
        
        tfx = strcmp(str,'No TLE found'); % Confronto fra stringhe
        if tfx == 1
            disp(strcat('GEO selector, warning: no TLE online for', " ", Stringa_NORAD, " ", 'satellite'))
            disp('   ')
            ExTLE=strcat(data_path, '/TLE_', Stringa_NORAD, '.txt');
            file_mancanti=isfile(ExTLE); % Variabile logica, vale 1 se il file esiste, 0 altrimenti
            if file_mancanti==1
                answer=strcat('File', " ",'TLE_', Stringa_NORAD, '.txt', " ", 'exist');
                disp(answer) 
                disp('   ')
                copyfile(ExTLE, temp_tle) % Copia il TLE nel file temporaneo "tle.txt"
            else
                answer=strcat('File', " ",'TLE_', Stringa_NORAD, '.txt', " ", 'does not exist');
                disp(answer)
                disp('   ')
                AR_eph(i)=AR0(i); 
                DEC_eph(i)=DEC0(i);
                disp('GEO selector: skip to next satellite')
                disp('   ')
                continue
            end   
            
        end    
        
        % Salvataggio del TLE corrente nel file "tle.txt"  
        fid2 = fopen(temp_tle,'w');
        fprintf(fid2, '%s \n', str);
        fclose(fid2);
    end
    
       
        disp(strcat('GEO selector: compute ephemeris for satellite', {' '}, Stringa_NORAD))
        
        % Funzione per il calcolo della effemeride del satellite nel
        % modello SGP4 in gradi (equinozio alla data)
            
        [AR_eph(i), DEC_eph(i)]=SGP4_Ephemeris(temp_tle, Lat_O, Long_O, h_O, MJD(i));
                    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transformation of calculated equatorial coordinates for satellites from %
% equinox to date at J2000 (degrees)                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JDin=2400000.5+mean(MJD); % giurno giuliano medio dell'epoca alla data
[AR_eph2000, DEC_eph2000] = EQ2EQ2000_B(AR_eph, DEC_eph, JDin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Offset observed-computed %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('     ')
disp('     ')
disp(strcat('GEO selector: COMPUTE OFFSET', {' '}, 'observed-ephemeris for every satellites'))

% Inizializzazione vettori
delta_AR=zeros(1, length(NORAD));
delta_DEC=zeros(1, length(NORAD));
offset=zeros(1, length(NORAD));

for i=1:length(NORAD)

        % Calcolo differenze AR e DEC osservato-calcolato in arcsec
        delta_AR(i)=3600*C*angdiff(AR_eph2000(i)/C, AR0(i)/C)*cos(DEC0(i)/C);
        delta_DEC(i)=3600*C*angdiff(DEC_eph2000(i)/C, DEC0(i)/C);
        
        % Distanza fra punto osservato e calcolato (arcsec)
        offset(i)=sqrt(delta_AR(i)^2+delta_DEC(i)^2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extraction of the observation closest to the satellite ephemeris %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=1;                       % Indice progressivo delle osservazioni originali
k=1;                       % Indice del vettore di osservazioni selezionate per un solo satellite
Number_block=length(MJD);  % Numero massimo dei blocchi di osservazioni dei GEO

% Ciclo preliminare per il conteggio del numero di satelliti nelle immagini
% Serve per selezionare le immagini che mostrano sempre lo stesso numero di
% satelliti
while Number_block>=1  
   
    MJD_aux0=abs(MJD-MJD(i));
    numero(k)=sum(MJD_aux0 == 0);        % Conteggio numero di tempi uguali, se numero=1 il satellite è unico
    i=i+numero(k);                       % Aggiorna indice per esaminare il blocco di osservazioni astrometriche successive
    Number_block=Number_block-numero(k); % Aggiorna numero delle osservazioni rimanenti da esaminare
    k=k+1;    
end

% Numero maggiormente ricorrente di satelliti per ogni immagine
Sat_number=mode(numero); 

i=1;                       % Indice progressivo delle osservazioni originali
k=1;                       % Indice del vettore di osservazioni selezionate per un solo satellite
start=1;                   % Indice dell'elemento inziale del blocco di osservazioni dei satelliti GEO
Number_block=length(MJD);  % Numero massimo dei blocchi di osservazioni dei GEO

% Ciclo per la selezione del satellite che si trova più vicino all'effemeride calcolata
while Number_block>=1  
   
    MJD_aux=abs(MJD-MJD(i));
    numero_zeri=sum(MJD_aux == 0); % Conteggio numero di tempi uguali, se numero_zeri=1 il satellite è unico
    
        [Min, I]=min(offset(start:start+numero_zeri-1)); % Estrazione minimo dal sottoinsieme del vettore offset con i tempi uguali
        I=I+start-1; % Trasforma l'indice del minimo del sottoinsieme dell'offset nel corrispondente indice delle osservazioni astrometriche
        
        disp(strcat('GEO selector: min offset (block i=', num2str(start),')', {' '},'for satellite', {' '}, num2str(NORAD(start)),{' '}, '(arcsec):',{' '}, num2str(offset(I))))
        
        if numero_zeri >= Sat_number % Estrae solo le osservazioni nella cui immagine c'è un numero maggiore o uguale a quello più frequente dei satelliti
        
             disp(strcat('GEO selector: extract observed astrometry (block i=', num2str(start),')', {' '}, 'nearest the ephemeris for satellite', {' '}, num2str(NORAD(start))))
        
             % Estrae le osservazioni astrometriche associabili al numero Norad
             % del satellite nell'header delle immagini
             anno(k)=YYYY(I); mese(k)=MM(I); giorno(k)=DD(I); ora(k)=hh(I); minuto(k)=mm(I); secondo(k)=ss(I); 
             AR1(k)=AR0(I); DEC1(k)=DEC0(I); NORAD1(k)=NORAD(I); MJD1(k)=MJD(i);
             
             k=k+1;                            % Aggiorna indice del sottoinsieme estratto delle misure astrometriche
        end
        
        start=start+numero_zeri;               % Aggiorna indice dell'elemento iniziale del blocco da esaminare
        i=i+numero_zeri;                       % Aggiorna indice per esaminare il blocco di osservazioni astrometriche successive
        Number_block=Number_block-numero_zeri; % Aggiorna numero delle osservazioni rimanenti da esaminare
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If-else-end statement to delete astrometric observations that fall    %
% within the belt of the Milky Way. If DMW=1 it deletes them, otherwise %
% it keeps them all.                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if DMW==1
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DELETE THE ASTROMETRIC OBSERVATIONS THAT FALL WITHIN THE MILKY WAY BELT %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp('   ')
disp('     GEO selector: deleting astrometry inside Milky Way')
disp('   ')

[YYYY2, MM2, DD2, hh2, mm2, ss2, AR2, DEC2, NORAD2] = equat2galactic(anno, mese, giorno, ora, minuto, secondo, AR1, DEC1, NORAD1);

% Numero delle osservazioni astrometriche con tempi diversi e fuori dalla Via Lattea
Nobs=length(NORAD2);
disp(strcat('GEO selector: number of single astrometric observations:', {' '}, num2str(Nobs)))

else
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOLDS ALL ASTROMETRIC OBSERVATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
YYYY2 = anno;
MM2 = mese;
DD2 = giorno;
hh2 = ora;
mm2 = minuto;
ss2 = secondo;
AR2 = AR1;
DEC2 = DEC1;
NORAD2 = NORAD1;

% Numero delle osservazioni astrometriche con tempi diversi
Nobs=length(NORAD2);
disp('   ')
disp('   ')
disp('     GEO selector: Milky Way filter non-active')
disp('   ')
disp(strcat('==> GEO selector - number of single astrometric observations:', {' '}, num2str(Nobs)))

end

% Calcolo vettore MJD2 delle sole osservazioni che interessano
MJD2=Mjday(YYYY2, MM2, DD2, hh2, mm2, ss2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COUNTING THE NUMBER OF DIFFERENT SATELLITES OBSERVED %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=1;
SAT(n)=NORAD2(1);
for i=1:length(NORAD2)
    if abs(SAT(n)-NORAD2(i))>0
        SAT(n+1)=NORAD2(i);    % Vettore che contiene il numero NORAD dei satelliti non ripetuti
        n=n+1;
    end
end

Nsat=length(SAT); % Numero di satelliti diversi
disp(strcat('==> GEO selector - number of different satellites to analize:', {' '}, num2str(n)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cycle for separating astrometric data for each satellite using the vector %
% containing the NORAD number of non-repeating satellites                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:Nsat
    k=1;
    for j=1:Nobs
        if abs(SAT(i)-NORAD2(j))==0
        ARsat(k)=AR2(j);
        DECsat(k)=DEC2(j);
        MJDsat(k)=MJD2(j);
        YYYYsat(k)=YYYY2(j);
        MMsat(k)=MM2(j);
        DDsat(k)=DD2(j);
        hhsat(k)=hh2(j);
        mmsat(k)=mm2(j);
        sssat(k)=ss2(j);
        NORADsat(k)=NORAD2(j);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter of the best observations for single track satellites using orbital best fit with find_orb %                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:Nsat
    obs=length(MJD_sat{i});            % Numero osservazioni i-esimo satellite
    name_sat=num2str(NORAD_sat{i}(1)); % Numero Norad i-esimo satellite
    
    if obs >= 3 % Numero minimo di osservazioni necessarie per il best fit orbitale
        
        % Selezione coordinate di best fit orbitale usando find_orb
        disp(strcat('GEO selector: select best fit orbit astrometry for Norad', {' '}, name_sat))
        disp('   ')
        [YYYY3, MM3, DD3, hh3, mm3, ss3, AR3, DEC3, S]=Fit_orb(data_path, DELTA, name_sat, YYYY_sat{i}, MM_sat{i}, DD_sat{i}, hh_sat{i}, mm_sat{i}, ss_sat{i}, AR_sat{i}, DEC_sat{i});
        disp('   ')
        disp(strcat('GEO selector: angle traveled by the satellite', {' '}, name_sat, {' '},'(degree):', {' '}, num2str(S)))
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
        disp(strcat('GEO selector: satellite', {' '}, name_sat, {' '},'not filtered with find_orb, too few observations'))
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

        disp(strcat('GEO selector: filtering multiple', {' '}, 'GEO satellite traces completed'))
        disp('   ')

end
