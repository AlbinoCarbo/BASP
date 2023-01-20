% BASP main script for plate solve, tracks extraction, astrometric measurements and TDM saving from SST images in fits format. 
% The coordinates of the center of the satellite tracks reported in the TDM are topocentric and referred to the J2000 equinox.
% 
% ALGORITHM AND HISTORY:
%
% Legge i parametri iniziali dal file "Settings_astrometry.txt" (Matlab)
%
% Verifica che esistano i file fit raw SST da analizzare. Se i file non
% esistono torna in home ed esce dal programma. Se ne manca solo qualcuno
% prosegue l'esecuzione (Matlab).
%
% Verifica l'esistenza dei file raw dei bias da cui ricavare il master bias. Se
% mancano tutti questi file esce dall'esecuzione, se ne manca solo qualcuno prosegue. 
% I file dei bias vengono rinominati con l'aggiunta del suffisso "BIAS" per distinguerli dalle
% immagini SST evitando così che vengano analizzati da Astrometry.net come immagini normali.
% Se le immagini con suffisso "BIAS" sono già presenti, mentre mancano quelle originali,
% le conta normalmente e prosegue l'esecuzione (Matlab).
%
% Creazione di un master bias frame mediano dai singoli bias frame rinominati. 
% I bias possono non avere una numerazione continua, se ne manca qualcuno lo script
% li salta e passa al successivo (Python3)
%
% Calibrazione delle immagini SST con sottrazione del master bias frame (Python3)
%
% Lettura header keys delle immagini calibrate con il master bias e salvataggio nel file "Data_keys.txt" 
% per avere RA e DEC del centro immagine. Questo dato facilita la calibrazione con astrometry.net. 
% Dalla versione del 26 febbraio 2021 se l'heade del file non ha le key richieste dalla pipeline cancella 
% il file e passa al successivo. Se il file da cui leggere le key non esiste passa a quello successivo (Python 3)
%
% Calibrazione astrometrica delle immagini SST (master bias subtracted), con Astrometry.net locale. 
% Questo processo crea le immagini con le costanti WCS nell'header (Linux)
%
% Verifica che esistano i file di tipo WCS. Se non sono stati creati file WCS esce dall'esecuzione (Matlab)
%
% Lettura header keys delle immagini WCS e salvataggio in "Data_keys.txt" (Python 3)
%
% Estrazione sia delle header keys importanti sia delle tracce dei satelliti dalle immagini aventi suffisso di calibrazione WCS con ASTRiDE.
% Dalla versione del 4 febbraio 2021 in poi, se una immagine CONTIENE DIVERSE TRACCE GEO vengono estratte tutte, ma avranno lo stesso numero Norad del
% satellite osservato. Starà poi alla funzione GEO_selector_astrometry_filter.m scegliere le osservazioni corrette associabili al satellite
% con il numero Norad dell'header (Python3)
% 
% Dalla versione del 6 ottobre 2021 in poi è attiva la funzione "Astrometric_filter", applicabile ai casi 
% in cui c'è un solo satellite per ogni immagine. La funzione sceglie le osservazioni astrometriche migliori
% in base al fit orbitale fatto con Find_orb attivato dallo script Fit_orb.m. 
% Astrometric filter taglia anche le osservazioni multiple con i tempi uguali, indice di un cattivo 
% riconoscimento della traccia del singolo satellite da parte di ASTRiDE. In Astrometric filter si può attivare 
% anche la funzione per cancellare le misure astrometriche che ricadono nella fascia della Via Lattea prima 
% di fare il best fit orbitale.
%
% Dalla versione del 13 ottobre 2021 in poi è attiva la funzione "GEO_selector_astrometry_filter.m", applicabile ai casi 
% in cui ci sono diversi satelliti GEO per ogni immagine. La funzione sceglie le osservazioni astrometriche associabili 
% al satellite Norad indicato nell'header delle immagini rigettando l'astrometria delle altre tracce satellitari. 
% Dopo di che, come in "Astrometric_filter" che vale per i singoli satelliti, alimina le osservazioni con i residui più 
% elevati in base al fit orbitale fatto con find_orb attivato dallo script Fit_orb.m. In Astrometric filter si può attivare 
% anche la funzione per cancellare le misure astrometriche che ricadono nella fascia della Via Lattea prima 
% di fare il best fit orbitale. 
%
% Salvataggio astrometria satelliti nel file "MPC_YYYYMMDD.txt" formato MPC. Attenzione che il numero Norad del satellite 
% non può avere più di 7 cifre, altrimenti non rientra nel formato del Minor Planet Center (Matlab)
%
% Salvataggio del TDM, per ogni satellite osservato nella sessione, in file del tipo "IT_CASSINIYYYYMMDDThhmmss_NORAD.txt".
% WARNING: salva anche i TDM con una sola osservazione anche se non se ne può controllare la qualità astrometrica. 
% Sono necessarie minimo tre osservazioni per avere una stima affidabile del residuo medio del TDM (Matlab)
%
% Dalla versione del 30 marzo 2022 in poi, nel caso di osservazione di satelliti Galileo per la calibrazione, i residui vengono valutati dal
% software python del PoliMI che si trova dentro la cartella "Evaluate_TDM".
%
% Dal 20 gennaio 2023 in poi sono stati introdotti gli ulteriori file di settings: Settings_Astrometry.txt che contiene il limite inferiore e
% superiore in arcsec/pixel della scala dell'immagine e Settings_MPC.txt che contiene i dati per scrivere l'intestazione del file delle
% osservazioni astrometriche nel formato del Minor Planet Center.
%
% N.B. In tutta la pipeline se un file fits non esiste si passa al successivo senza bloccare l'eleborazione.
%
% FUNZIONI ESTERNE DI PYTHON 3 UTILIZZATE:
%
% fits_master_bias.py, script per la creazione di un master bias a partire
% dai singoli fits bias
%
% fits_calibration.py, script per la calibrazione dei file fits usando il
% master_bias creato da "fits_master_bias.py"
%
% fits_keys_reader.py, script per la lettura di alcune keys nell'header dei file fits
%
% SST_Astride_TDM.py, script per l'estrazione sia delle header keys importanti sia delle coordinate dei punti medi delle tracce 
% dei satelliti dai file calibrati fits WCS
%
% evaluate_TDM_measurements.py, script del PoliMI per verificare i residui
% dei TDM dei Galileo usati per le campagne di calibrazione. Si trova
% dentro la cartella "Evaluate_TDM".
%
% FUNZIONI DI MATLAB UTILIZZATE:
% GEO_selector_astrometry_filter.m, funzione per la scelta di un satellite GEO osservato in mezzo a una "flottiglia" 
% di GEO presenti nella stessa immagine con il modello SGP4 usando il TLE corrente da celestrak.com. 
% Dopo l'identificazione delle osservazioni che si riferiscono al satellite GEO dell'header le osservazioni 
% vengono analizzate con Fit_orb.m per l'eliminazione delle misure astrometriche peggiori.
%
% Astrometric_filter.m, funzione analoga a GEO_selector_astrometry_filter.m, ma per il best 
% fit orbitale delle osservazioni astrometriche dei satelliti singoli usando find_orb di Bill Gay. 
% Cancella le osservazioni peggiori, che distano più di DELTA arcsec dal valore di best fit.
%
% Fit_orb.m, funzione richiamata da Astrometric_filter.m per il best fit 
% dell'orbita di un satellite artificiale eliminando le osservazioni astrometriche 
% con un residuo maggiore di DELTA arcsec.
%
% Astrometry_MPC.m, procedura usata da Fit_orb.m per il salvataggio delle osservazioni
% astrometriche in un file di testo nel formato del MPC. In questo modo le osservazioni 
% astrometriche possono essere analizzate anche manualmente da Find_orb di
% Bill Gray.
%
% equat2galactic.m, funzione usata da Astrometric_filter.m per la selezione 
% delle posizioni astrometriche dei satelliti che si trovano al di fuori della fascia 
% della Via Lattea. Utile per evitare di includere l'astrometria dei satelliti fatta 
% nei campi stellari più affollati, dove le misure di posizione accurate sono difficili.
%
% format_seconds.m, funzione per la formattazione dei secondi di tempo
% nel formato XX.XXX. Serve per rispettare lo standard EUSST dei TDM.
%
% A. Carbognani, INAF-OAS
%
% Versione del 20 gennaio 2023

clear all

disp('                                                                       ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('% Bologna Astrometry Satellites Pipeline (BASP)-SST PROJECT INAF-OAS  %')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('                                                                       ')
disp('                    by Albino Carbognani (INAF-OAS)                    ')
disp('                            Linux Version                              ')
disp('                              Mar 2022                                 ')
disp('                                                                       ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('   ')

% Lettura del file dei Settings
wholefile_set = fileread('./Settings_astrometry.txt');
disp('READ SETTINGS FILE')
disp('   ') 

% Split dei Settings e inizializzazione delle variabili stringa
set = regexp(wholefile_set,'\$+','split');

image_prefix=strtrim(set{4});                % Nome comune delle immagini SST nel formato "YYYYMMDD_NNN.fit"
data_path=strtrim(set{7});                   % Path delle immagini SST
home_path=strtrim(set{10});                  % Path della cartella dove si trovano gli script Python
Nmin=str2double(strtrim(set{13}));           % Minimum SST images number
Nmax=str2double(strtrim(set{16}));           % Maximum SST image number
NBmin=str2double(strtrim(set{19}));          % Minimum bias number
NBmax=str2double(strtrim(set{22}));          % Maximum bias number
Long_O=str2double(strtrim(set{25}));         % Longitudine E osservatore (gradi)
Lat_O=str2double(strtrim(set{28}));          % Latitudine osservatore (gradi)
h_O=str2double(strtrim(set{31}));            % Quota osservatore slm (metri)
Ephem_comp=str2double(strtrim(set{34}));     % Filtro per le misure astrometriche basato sul best fit orbitale (1=multiple GEO satellites, 2=single satellite, 3=No, 4=Galileo per calibrazione)
                                             % Passano solo le osservazioni che distano al massimo DELTA arcsec dalla posizione teorica del
                                             % satellite.
DELTA=str2double(strtrim(set{37}));          % Massimo residuo consentito per filtrare l'astrometria peggiore (arcsec)
DMW=str2double(strtrim(set{40}));            % Se DMW=1 il filtro che cancella le osservazioni astrometriche fatte nella fascia della Via Lattea è attivato.
soglia_astride=str2double(strtrim(set{43})); % Se soglia_astride=1 la rilevazione del bordo del satellite è più sensibile, adatto per tracce deboli. Valori 2 o 3 vanno bene per satelliti brillanti

% NOTA: la funzione strtrim cancella gli spazi vuoti ad inizio e fine stringa

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK EXISTENCE OF FITS RAW IMAGES AND BIAS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('CHECK RAW SST IMAGES AND BIAS EXISTENCE')
disp('   ')

file_raw=0;  % Conteggio delle immagini SST
file_rawB=0; % Conteggio delle immagini BIAS

% Matlab va nella cartella dove ci sono le immagini SST e BIAS
cd(data_path)

% Conteggio immagini SST
for i=Nmin:1:Nmax
   
   image_name0=strcat(image_prefix, '_', num2str(i), '.fit');
   file_mancanti=isfile(image_name0); % Variabile logica, vale 1 se il file esiste, 0 altrimenti
   
   % Questa istruzione if fa saltare il conteggio dei file dei bias nel caso i loro numeri d'ordine 
   % siano compresi fra Nmin e Nmax, ossia fra quelli delle immagini da analizzare.
   if (i>= NBmin) && (i<= NBmax)
       continue
   end
   
   if file_mancanti==1
     answer=strcat('File', " ",image_name0, " ", 'exist');
     disp(answer)
     file_raw=file_raw+1;
     
   else
     answer=strcat('File', " ",image_name0, " ", 'does not exist');
     disp(answer)
     continue % Se il file non esiste si passa all'immagine successiva
   end
   
end

disp('   ')

% Conteggio immagini BIAS
for i=NBmin:1:NBmax
   
   image_name0B=strcat(image_prefix, '_', num2str(i), '.fit');
   file_mancanti=isfile(image_name0B); % Variabile logica, vale 1 se il file esiste, 0 altrimenti
   
   if file_mancanti==1
     answer=strcat('BIAS', " ",image_name0B, " ", 'exist');
     disp(answer)
     file_rawB=file_rawB+1;
     
     % Rename dei file BIAS per toglierli dal ciclo di analisi delle immagini
     image_name_BIAS=strcat(image_prefix, '_BIAS_', num2str(i), '.fit');
     rename_command=strcat('mv', " ", image_name0B, " ",image_name_BIAS); 
     system(rename_command);
     
   else
     answer=strcat('BIAS', " ",image_name0B, " ", 'does not exist');
     disp(answer)
     continue % Se il file non esiste si passa all'immagine successiva
   end
   
end

disp('   ')

% Conteggio file dei bias con suffisso BIAS rinominate nel ciclo precedente
for i=NBmin:1:NBmax
   
   image_name_BIAS=strcat(image_prefix, '_BIAS_', num2str(i), '.fit');
   file_mancanti=isfile(image_name_BIAS); % Variabile logica, vale 1 se il file esiste, 0 altrimenti
   
   if file_mancanti==1
     answer=strcat('BIAS', " ",image_name_BIAS, " ", 'exist');
     disp(answer)
     file_rawB=file_rawB+1;
          
   else
     answer=strcat('BIAS', " ",image_name_BIAS, " ", 'does not exist');
     disp(answer)
     continue % Se il file non esiste si passa all'immagine successiva
   end
   
end


disp('   ')
disp('CHECK IMAGES SST/BIAS FILES DONE')
disp('   ')

% Se non ci sono file raw SST torna in home e termina l'esecuzione
if file_raw == 0
    cd(home_path)
    answer1=strcat('SST files does not exist');
     disp(answer1)
     answer2=('Exit from BASP');
     disp(answer2)
     return
end

% Se non ci sono file BIAS torna in home e termina l'esecuzione
if file_rawB == 0
    cd(home_path)
    answer1B=strcat('BIAS files does not exist');
     disp(answer1B)
     answer2B=('Exit from BASP');
     disp(answer2B)
     return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MASTER BIAS FRAME GENERATION AND FITS IMAGE CALIBRATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('COMPUTE MASTER BIAS FRAME')
disp('   ')

% Matlab va nella cartella dove ci sono gli script Python della pipeline
cd(home_path)

% Define command line for fits_master_bias.py
command_line_bias=strcat('python3', " ", 'fits_master_bias.py', " ", data_path, " ", image_prefix, '_BIAS_', " ", '.fit', " ", num2str(NBmin), " ", num2str(NBmax-NBmin+1));
 
% Generate a file "master_bias.fit" in data_path
system(command_line_bias);

disp('RAW IMAGES CALIBRATION (MASTER BIAS SUBTRACTION)')
disp('  ')

% Define command line for fits_calibration.py
command_line_calibration=strcat('python3', " ", 'fits_calibration.py', " ", data_path, " ", image_prefix, '_', " ", '.fit', " ", num2str(Nmin), " ", num2str(Nmax-Nmin+1));

% Images calibration, final format: "YYYYMMDD_NNN_cal.fit"
% N.B. Sono le immagini "*_cal.fit" che vanno calibrate con astrometry.net
system(command_line_calibration);

disp('READ HEADERS KEYS FROM CALIBRATED IMAGES')
disp('  ')

% Define command line for fits_calibration.py
command_line_header=strcat('python3', " ", 'fits_keys_reader.py', " ", data_path, " ", image_prefix, '_', " ", '_cal.fit', " ", num2str(Nmin), " ", num2str(Nmax-Nmin+1));

% Read heder keys: data e ora, nome oggetto, tempo di esposizione, AR e DEC.
% Save in: "Data_keys.txt".
system(command_line_header);

disp('EXIT MASTER BIAS GENERATION, FITS CALIBRATION AND IMPORTANT HEADER KEYS READING')
disp('  ')

% Matlab va nella cartella dove ci sono le immagini SST calibrate
% e legge il file "Data_keys.txt" (generato da fits_keys_reader.py), contenete le 
% coordinate del centro immagine per velocizzare la calibrazione astrometrica con Astrometry.net
cd(data_path)

T0 = readtable('Data_keys.txt'); 
disp('Read RA and DEC center images from Data_keys.txt')
disp('  ')

RA_ast=T0.Var4;    % RA centro immagine
DEC_ast=T0.Var5;   % DEC centro immagine

disp('IMAGES ASTROMETRY CALIBRATION WITH ASTROMETRY.NET')
disp('   ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WCS ASTROMETRIC CALIBRATION WITH LOCAL ASTROMETRY.NET %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matlab va nella cartella dove ci sono gli script
cd(home_path)

% Lettura del file dei Settings
wholefile_set2 = fileread('./Settings_Astrometry.txt');
disp('READ ASTROMETRY SETTINGS FILE')
disp('   ')

% Split dei Settings e inizializzazione delle variabili stringa
set2 = regexp(wholefile_set2,'\$+','split');

% Limite inferiore e superiore della scala dell'immagine in arcsec/pixel
InfScale=strtrim(set2{5});      
SupScale=strtrim(set2{8});

% Matlab va nella cartella dove ci sono le immagini SST e BIAS
cd(data_path)

kk=1; % Indice di RA_ast e DEC_ast (coordinate del centro immagine)
for i=Nmin:1:Nmax
       
   image_name=strcat(image_prefix, '_', num2str(i), '_cal.fit');
   
   % Verifica dell'esistenza del file
   file_mancanti=isfile(image_name); % Variabile logica, vale 1 se il file esiste, 0 altrimenti
   
   if file_mancanti==1
         image_name_WCS=strcat(image_prefix, '_WCS_', num2str(i), '.fit');
         s1 = string(RA_ast(kk),'hh:mm:ss'); % Trasformazione in stringa dell'AR
         s2 = char(DEC_ast(kk));             % Trasformazione in stringa della DEC
   
         % Define command line for Astrometry.net calibration
         %command_line_astrometry=strcat('solve-field ', " ", '-L 0.5', " ", '-H 0.7', " ", '-u arcsecperpix', " ", '--ra', " ", s1, " ", '--dec', " ", s2, " ",'--radius 0.5', " ", image_name, ' -N',  " ", image_name_WCS);
         command_line_astrometry=strcat('solve-field ', " ", InfScale, " ", SupScale, " ", '-u arcsecperpix', " ", '--ra', " ", s1, " ", '--dec', " ", s2, " ",'--radius 0.5', " ", image_name, ' -N',  " ", image_name_WCS);
   
         % Display name image to be processed with astrometry.net
         name_processing=strcat('Astrometry.net processing', " ", image_name);
         disp(name_processing)
         disp('    ')
        
         system(command_line_astrometry); % Calibration image with astrometry.net
   
         kk=kk+1;
   else
       continue % Se il file da calibrare con Astrometry.net non esiste si passa all'immagine successiva
   end
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK EXISTENCE OF WCS IMAGES FITS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('CHECK WCS IMAGES EXISTENCE')
disp('   ')

file_WCS=0; % Conteggio del numero di file WCS

for i=Nmin:1:Nmax
   
   image_name0=strcat(image_prefix, '_WCS_', num2str(i), '.fit');
   file_mancanti=isfile(image_name0); % Variabile logica, vale 1 se il file esiste, 0 altrimenti
   
   if file_mancanti==1
     answer=strcat('File', " ",image_name0, " ", 'exist');
     disp(answer)
     file_WCS=file_WCS+1;
     
   else
     answer=strcat('File', " ",image_name0, " ", 'does not exist');
     disp(answer)
     continue
   end
   
end

% Matlab ritorna nella cartella degli script Python e lancia l'estrazione delle header's keys e delle
% tracce dei satelliti dalle immagini WCS. Se una traccia non viene trovata, nelle coordinate
% AR e DEC mette NaN

cd(home_path)

% Se non ci sono file WCS termina l'esecuzione
if file_WCS == 0
    answer1=strcat('Files WCS does not exist');
     disp(answer1)
     answer2=('Exit from BASP');
     disp(answer2)
     return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRACTING KEYWORDS FROM WCS FITS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('  ')
disp('READ HEADERS KEYS FROM CALIBRATED WCS IMAGES')
disp('  ')

% Define command line for fits_calibration.py
command_line_header=strcat('python3', " ", 'fits_keys_reader.py', " ", data_path, " ", image_prefix, '_WCS_', " ", '.fit', " ", num2str(Nmin), " ", num2str(Nmax-Nmin+1));

% Read heder keys: data e ora, nome oggetto, tempo di esposizione,
% AR e DEC from WCS images. Save in: "Data_keys.txt".

system(command_line_header);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRACTION OF SATELLITE TRACES FROM WCS IMAGES, READING OF HEADER KEYS AND MEASUREMENT OF AR AND DEC COORDINATES OF THE TRACE CENTER %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('START SATELLITES TRACE EXTRACTION WITH ASTRiDE')
disp('   ')

% Define command line for SST_ASTRiDE_TDM.py
command_line_astrid=strcat('python3', " ", 'SST_Astride_TDM.py', " ", data_path, " ", image_prefix, '_WCS_', " ", '.fit', " ", num2str(Nmin), " ", num2str(Nmax-Nmin+1), " ", num2str(soglia_astride));
 
% Satellite headers and trace extraction, save header keys, AR and DEC coordinates in "Data_headers_streaks.txt"
% Warning: a comma is placed after each key to facilitate the recognition of columns with readtable.

system(command_line_astrid);

disp('EXIT ASTRiDE SCRIPT RETURN TO MATLAB')
disp('  ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATION OF A TDM FILES FOR EACH OBSERVED SATELLITE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Torna nella cartella dove ci sono i file di output generati da ASTRiDE
% e li legge per inizializzare la matrice delle osservazioni

cd(data_path)

T1 = readtable('Data_headers_streaks.txt'); 
disp('Read data from Data_headers_streaks.txt (WCS images)')
disp('  ')

data_ora=T1.Var1;   % Data e ora di inizio esposizione
NORAD=T1.Var2;      % Numero NORAD del satellite 
exposure=T1.Var3;   % Tempo di esposizione (s)
AR0=T1.Var4;        % AR del centro della traccia del satellite (gradi)
DEC0=T1.Var5;       % DEC del centro della traccia del satellite (gradi)

% Conversione di data e ora da stringa a formato data
C=datetime(data_ora, 'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS');

% Estrazione anno, mese e giorno dalla data
[YYYY, MM, DD]=ymd(C);

% Estrazione ore, minuti e secondi dalla data
[hi,mi,si] = hms(C);

%%%%%%%%%%%%%%%%%%%%%%%%%
% Mean time computation %
%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Compute mean exposure times')
disp('  ')

si=si+exposure/2; 

% Test per non avere i secondi >= di 60
for i=1:length(si)
   if si(i) >= 60
       si(i)=si(i)/60; 
       integ=floor(si(i));
       si(i)=60*(si(i)-integ);
       mi(i)=mi(i)+1;
   end
end

% Test per non avere i minuti >= di 60
for i=1:length(mi)
   if mi(i) >= 60
       mi(i)=mi(i)/60; 
       integ=floor(mi(i));
       mi(i)=60*(mi(i)-integ);
       hi(i)=hi(i)+1;
   end
end

% Test per non avere le ore >= 24
for i=1:length(hi)
   if hi(i) >= 24
       hi(i)=hi(i)/24; 
       integ=floor(hi(i));
       hi(i)=24*(hi(i)-integ);
       DD(i)=DD(i)+1;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrix initialization of astrometric observations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matlab ritorna nella cartella con gli script Python + Matlab
cd(home_path)

% Controlla se a ogni data corrisponde una sola coordinata AR e DEC
if length(YYYY)~= length(AR0)
disp('   ')
    disp(' EXIT BASP: RA number is different from date number')
    disp('   ')
    return
end

if length(YYYY)~= length(DEC0)
disp('   ')
    disp(' EXIT BASP: DEC number is different from date number')
    disp('   ')
    return
end

% Matrice delle osservazioni astrometriche
MO=[YYYY MM DD hi mi si AR0 DEC0 NORAD];

% Cancella le righe della matrice MO che contengono NaN in AR0 o DEC0
MO(any(isnan(MO), 2),:)=[];

% Ordinamento della matrice MO per ordine crescente del numero NORAD
% N.B. Questo ordinamento permette di raggruppare le osservazioni di un satellite
% su righe contigue e diventa facile contare quanti satelliti diversi ci
% sono per generare i rispettivi TDM.

MO=sortrows(MO,9);

% Se non sono state rilevate tracce satellitari esce da SST2TDM
if length(MO(:,1))<1
disp('   ')
    disp(' EXIT BASP: no satellites traces found!')
    disp('   ')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Filter of the observed astrometric positions using find_orb  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter for images with multiple GEO tracks %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Ephem_comp == 1

% Vettori riga delle osservazioni ordinate secondo il numero crescente del
% satellite (più di uno per ogni immagine), filtrate in base al best fit orbitale
% con find_orb
[YYYY_ord0, MM_ord0, DD_ord0, hi_ord0, mi_ord0, si_ord0, AR_ord0, DEC_ord0, NORAD_ord0]=GEO_selector_astrometry_filter(data_path, Long_O, Lat_O, h_O, DELTA, DMW, MO(:, 1), MO(:, 2), MO(:, 3), MO(:, 4), MO(:, 5), MO(:, 6), MO(:, 7), MO(:, 8), MO(:, 9));    
   
% Trasformazione da vettore riga a vettore colonna

YYYY_ord = YYYY_ord0';
MM_ord = MM_ord0';
DD_ord = DD_ord0';
hi_ord = hi_ord0';
mi_ord = mi_ord0;
si_ord = si_ord0;
AR_ord = AR_ord0;
DEC_ord = DEC_ord0;
NORAD_ord = NORAD_ord0;

    if isempty(YYYY_ord)
        disp(' EXIT BASP: all observations filtered, change DELTA coefficient in settings!')
        disp('   ')
        return
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter for images with traces of single satellites %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Ephem_comp == 2

% Vettori riga delle osservazioni ordinate secondo il numero crescente del
% satellite (uno per ogni immagine), filtrate in base al best fit orbitale
% con find_orb
[YYYY_ord0, MM_ord0, DD_ord0, hi_ord0, mi_ord0, si_ord0, AR_ord0, DEC_ord0, NORAD_ord0]=Astrometric_filter(data_path, DELTA, DMW, MO(:, 1), MO(:, 2), MO(:, 3), MO(:, 4), MO(:, 5), MO(:, 6), MO(:, 7), MO(:, 8), MO(:, 9));    

% Trasformazione da vettore riga a vettore colonna

YYYY_ord = YYYY_ord0';
MM_ord = MM_ord0';
DD_ord = DD_ord0';
hi_ord = hi_ord0';
mi_ord = mi_ord0;
si_ord = si_ord0;
AR_ord = AR_ord0;
DEC_ord = DEC_ord0;
NORAD_ord = NORAD_ord0;

    if isempty(YYYY_ord)
        disp(' EXIT BASP: all observations filtered, change coefficient in settings!')
        disp('   ')
        return
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unfiltered astrometric observations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Ephem_comp == 3
    
% Vettori colonna delle osservazioni astrometriche ordinate secondo il numero crescente del
% satellite senza nessun filtro basato su find_orb

YYYY_ord = MO(:, 1);
MM_ord = MO(:, 2);
DD_ord = MO(:, 3);
hi_ord = MO(:, 4);
mi_ord = MO(:, 5);
si_ord = MO(:, 6);
AR_ord = MO(:, 7);
DEC_ord = MO(:, 8);
NORAD_ord = MO(:, 9);

end


if Ephem_comp == 4
    
% Vettori colonna delle osservazioni astrometriche ordinate secondo il numero crescente del
% satellite senza nessun filtro basato su find_orb

YYYY_ord = MO(:, 1);
MM_ord = MO(:, 2);
DD_ord = MO(:, 3);
hi_ord = MO(:, 4);
mi_ord = MO(:, 5);
si_ord = MO(:, 6);
AR_ord = MO(:, 7);
DEC_ord = MO(:, 8);
NORAD_ord = MO(:, 9);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saving astrometry satellites in the extended MPC format %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lettura del file dei Settings
wholefile_set = fileread('./Settings_MPC.txt');
disp('READ MPC SETTINGS FILE')
disp('   ')

% Split dei Settings e inizializzazione delle variabili stringa
set = regexp(wholefile_set,'\$+','split');

MPC=strtrim(set{5});      
Code=strtrim(set{8});
Address=strtrim(set{11});
Observers=strtrim(set{14});
Measurers=strtrim(set{17});
Telescope=strtrim(set{20});
StarCatalog=strtrim(set{23});
Comments=strtrim(set{26});

disp('Save satellites astrometry in MPC format')
disp('  ')

% Nome file astrometria in formato MPC salvato nella cartella indicata da data_path
output_MPC_file=strcat(data_path, 'MPC_',num2str(YYYY_ord(1)), num2str(MM_ord(1),'%02.f'), num2str(DD_ord(1),'%02.f'), '.txt');
fid0 = fopen(output_MPC_file, 'w');

fprintf(fid0, strcat(Code, ' \n'));
fprintf(fid0, strcat(Address, ' \n'));
fprintf(fid0, strcat(Observers, ' \n'));
fprintf(fid0, strcat(Measurers, ' \n'));
fprintf(fid0, strcat(Telescope, ' \n'));
fprintf(fid0, 'ACK MPCReport file updated %s \n', datetime(now,'ConvertFrom','datenum'));
fprintf(fid0, strcat(StarCatalog, ' \n'));
fprintf(fid0, strcat(Comments, ' \n'));

for i=1:length(NORAD_ord)

% Calcolo giorno e decimali (UT)
Data=DD_ord(i)+(si_ord(i)/86400)+mi_ord(i)/(1440)+hi_ord(i)/24;
    
% Formattazione Data
Data_int=fix(Data);                              % Parte intera Data
Data_dec=fix(1000000*(Data-Data_int));           % Parte con 6 decimali della data
str_Data_int=num2str(Data_int,'%02.f');          % Parte intera data con zero davanti
str_Data_dec=num2str(Data_dec,'%06.f');          % Prima cifra della parte decimale 
Data_s=strcat(str_Data_int, '.',str_Data_dec);   % Data stringa nel formato XX.XXXXXX

%%% Formattazione AR

% Calcolo AR in h m s
AR=degrees2dms(AR_ord(i)/15);    % Ascensione retta J2000

% Formattazione dei secondi di AR
AR3_int=fix(AR(3));                         % Parte intera secondi di tempo
AR3_dec=fix(1000*(AR(3)-AR3_int));          % Parte con 3 cifre decimale secondi di tempo
str_AR3_int=num2str(AR3_int,'%02.f');       % Parte intera secondi di tempo con zero davanti
str_AR3_dec=num2str(AR3_dec,'%03.f');       % Prime 3 cifre della parte decimale

AR3=strcat(str_AR3_int, '.',str_AR3_dec);   % Numero stringa nel formato XX.XXX

%%% Formattazione Dec

if DEC_ord(i) >= 0
    segno='+';
else
    segno='-';
end

DECx=abs(DEC_ord(i)); % Eliminazione del segno della declinazione

% Calcolo Dec in ° ' "
DEC=degrees2dms(DECx); % Declinazione al J2000

% Formattazione dei " di Dec
DEC3_int=fix(DEC(3));                              % Parte intera secondi d'arco
DEC3_dec=fix(100*(DEC(3)-DEC3_int));               % Parte con 2 decimali dei secondi d'arco
str_DEC3_int=num2str(DEC3_int,'%02.f');            % Parte intera secondi d'arco con zero davanti
str_DEC3_dec=num2str(DEC3_dec,'%02.f');            % Prima cifra della parte decimale dei secondi d'arco
DEC3=strcat(str_DEC3_int, '.',str_DEC3_dec);       % Numero stringa nel formato XX.XX

ANNO1=strcat('KC', num2str(YYYY_ord(i)));

L_sat_name = strlength(num2str(NORAD_ord(i)));      % Cifre del numero Norad del satellite

if L_sat_name==1
NameSat=strcat('000000', num2str(NORAD_ord(i)));    % Numero NORAD del satellite nel formato MPC (deve avere 7 cifre)  
end

if L_sat_name==2
NameSat=strcat('00000', num2str(NORAD_ord(i)));     % Numero NORAD del satellite nel formato MPC (deve avere 7 cifre)  
end

if L_sat_name==3
   NameSat=strcat('0000', num2str(NORAD_ord(i)));   % Numero NORAD del satellite nel formato MPC (deve avere 7 cifre)  
end

if L_sat_name==4
NameSat=strcat('000', num2str(NORAD_ord(i)));       % Numero NORAD del satellite nel formato MPC (deve avere 7 cifre)  
end

if L_sat_name==5
NameSat=strcat('00', num2str(NORAD_ord(i)));        % Numero NORAD del satellite nel formato MPC (deve avere 7 cifre)  
end

if L_sat_name==6
NameSat=strcat('0', num2str(NORAD_ord(i)));         % Numero NORAD del satellite nel formato MPC (deve avere 7 cifre)  
end

if L_sat_name==7
NameSat=num2str(NORAD_ord(i));                      % Numero NORAD del satellite nel formato MPC (deve avere 7 cifre)  
end

fprintf(fid0, '     %s %s %02d %s%02d %02d %s%s%02d %02d %s         %s      %s\n', NameSat, ANNO1, MM_ord(i), Data_s, AR(1), AR(2), AR3, segno, DEC(1), DEC(2), DEC3, '15.5 G', MPC);
 
end

fclose(fid0);

%%%%%%%%%%%%%%%%%%
% TDM GENERATION %
%%%%%%%%%%%%%%%%%%

% Estrazione dei numeri NORAD dei satelliti diversi per cui vanno generati i TDM

n=1;
SAT(n)=NORAD_ord(1);
ANNO(n)=YYYY_ord(1);   
MESE(n)=MM_ord(1);
GIORNO(n)=DD_ord(1);
ORA(n)=hi_ord(1);
MINUTI(n)=mi_ord(1);
SECONDI(n)=si_ord(1);
Righe(n)=1;

for i=1:length(NORAD_ord)
    if abs(SAT(n)-NORAD_ord(i))>0
        SAT(n+1)=NORAD_ord(i);    % Vettore dei satelliti (non ripetuti) di cui fare i TDM
        ANNO(n+1)=YYYY_ord(i);    % Anno iniziale del TDM
        MESE(n+1)=MM_ord(i);      % Mese iniziale del TDM
        GIORNO(n+1)=DD_ord(i);    % Giorno iniziale del TDM
        ORA(n+1)=hi_ord(i);       % Ora iniziale del TDM
        MINUTI(n+1)=mi_ord(i);    % Minuti iniziali del TDM
        SECONDI(n+1)=si_ord(i);   % Secondi iniziali del TDM
        Righe(n+1)=i;
        n=n+1;
    else
        Righe(n)=i;      % Numero della riga finale del satellite
                         % Serve per scrivere lo STOP_TIME nel TDM
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TDM saving cycle: one for each different satellite observed %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:length(SAT)

disp(strcat('Save TDM file for satellite number', " ", num2str(SAT(j))))
disp('  ')
    
% Nome file astrometria in formato TDM salvato nella cartella indicata da data_path
output_TDM_file=strcat(data_path, 'SST_',num2str(ANNO(j)), num2str(MESE(j),'%02.f'), num2str(GIORNO(j),'%02.f'), 'T', num2str(ORA(j),'%02.f'), num2str(MINUTI(j),'%02.f'), num2str(fix(SECONDI(j)),'%02.f'), '_', num2str(SAT(j)), '.tdm');

% Cancella il file del TDM se è già presente
if exist(output_TDM_file, 'file')==2
   delete(output_TDM_file);
end

% Formato della data di creazione che viene scritta nel file contenente il TDM
% Se si mette hh il formato è 0-12 h
% Se si mette HH il formato è 0-24 h
datetime.setDefaultFormats('default','yyyy-MM-dd''T''HH:mm:ss.sss')

fid1 = fopen(output_TDM_file, 'at');

fprintf(fid1, 'CCSDS_TDM_VERS = 1.0 \n');
fprintf(fid1, '   \n');
fprintf(fid1, 'CREATION_DATE = %s \n', datetime(now,'ConvertFrom','datenum'));
fprintf(fid1, 'ORIGINATOR = INAF-OAS \n');
fprintf(fid1, '   \n');
fprintf(fid1, 'META_START');
fprintf(fid1, '   \n');
fprintf(fid1, 'TIME_SYSTEM = UTC \n');
fprintf(fid1, 'START_TIME = %02.0f-%02.0f-%02.0fT%02.0f:%02.0f:%s \n', ANNO(j), MESE(j), GIORNO(j), ORA(j), MINUTI(j), format_seconds(SECONDI(j)));
fprintf(fid1, 'STOP_TIME = %02.0f-%02.0f-%02.0fT%02.0f:%02.0f:%s \n', YYYY_ord(Righe(j)), MM_ord(Righe(j)), DD_ord(Righe(j)), hi_ord(Righe(j)), mi_ord(Righe(j)), format_seconds(si_ord(Righe(j))));
fprintf(fid1, 'PARTICIPANT_1 = IT_CASSINI \n');
fprintf(fid1, 'PARTICIPANT_2 = %d \n', SAT(j));
fprintf(fid1, 'MODE = SEQUENTIAL \n');
fprintf(fid1, 'PATH = 2,1 \n');
fprintf(fid1, 'ANGLE_TYPE = RADEC \n');
fprintf(fid1, 'REFERENCE_FRAME = EME2000 \n');
fprintf(fid1, '   \n');
fprintf(fid1, 'META_STOP \n');
fprintf(fid1, '   \n');
fprintf(fid1, 'DATA_START \n');
fprintf(fid1, '   \n');

for i=1:length(NORAD_ord)
    if SAT(j)==NORAD_ord(i)
       fprintf(fid1, 'ANGLE_1 = %02.0f-%02.0f-%02.0fT%02.0f:%02.0f:%s %03.4f \n', YYYY_ord(i), MM_ord(i), DD_ord(i), hi_ord(i), mi_ord(i), format_seconds(si_ord(i)), AR_ord(i)); 
       fprintf(fid1, 'ANGLE_2 = %02.0f-%02.0f-%02.0fT%02.0f:%02.0f:%s %02.4f \n', YYYY_ord(i), MM_ord(i), DD_ord(i), hi_ord(i), mi_ord(i), format_seconds(si_ord(i)), DEC_ord(i));
       fprintf(fid1, '   \n');
    end
end                 

fprintf(fid1, '   \n');
fprintf(fid1, 'DATA_STOP \n');

fclose(fid1);
    
end

% Matlab ritorna nella cartella con gli script Python + Matlab
cd(home_path)

disp('BASP DONE')
disp('  ')
