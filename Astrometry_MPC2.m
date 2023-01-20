% Astrometry_MPC2, alternative procedure to "Astrometry_MPC.m" for saving astrometry observations 
% in a text file in the extended Find Orb format (not compatible with the MPC format).
% In this way the time precision reaches up to 1/1000 second and the astrometric observations 
% can also be analyzed manually by Bill Gray's Find_orb (www.projectpluto.com).
%
% INPUT
% data_pathy: path delle immagini SST
% NORAD: numero NORAD del satellite
% YYYY_ord: vettore anno osservazione
% MM_ord: vettore mese osservazione
% DD_ord: vettore giorno osservazione
% hi_ord: vettore ore osservazione
% mi_ord: vettore minuti osservazione
% si_ord: vettore secondi osservazione
% AR_ord: vettore AR osservata in gradi al J2000
% DEC_ord: vettore DEC osservata in gradi al J2000
% name_astrometry_file: nome del file astrometrico di output (stringa)
%
% OUTPUT
% File di testo astrometrico salvato nella cartella dove si trovano le immagini SST da elaborare. 
%
% Albino Carbognani, INAF-OAS
% Versione del 20 gennaio 2023

function []=Astrometry_MPC2(data_pathy, NORAD, YYYY_ord, MM_ord, DD_ord, hi_ord, mi_ord, si_ord, AR_ord, DEC_ord, name_astrometry_file)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saving astrometry satellites in the extended Find_Orb format using the MJD %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(strcat('Astrometry MPC: save', {' '}, num2str(NORAD), {' '}, 'astrometry in MPC format'))
disp('  ')

% Save file in folder with images SST
fid0 = fopen(strcat(data_pathy, name_astrometry_file), 'w');

% Intestazione del file di output nel formato MPC
fprintf(fid0, strcat(Code, ' \n'));
fprintf(fid0, strcat(Address, ' \n'));
fprintf(fid0, strcat(Observers, ' \n'));
fprintf(fid0, strcat(Measurers, ' \n'));
fprintf(fid0, strcat(Telescope, ' \n'));
fprintf(fid0, 'ACK MPCReport file updated %s \n', datetime(now,'ConvertFrom','datenum'));
fprintf(fid0, strcat(StarCatalog, ' \n'));
fprintf(fid0, strcat(Comments, ' \n'));

for i=1:length(AR_ord)
    
% Formattazione Data e ora
secolo=fix(YYYY_ord(i)/100);
cifre_anno=num2str(floor((YYYY_ord(i)/100-floor(YYYY_ord(i)/100))*100), '%02.f'); % Ultime due cifre dell'anno (valido fino al 2099)
if secolo == 19
    Lettera='J';
else
    Lettera='K';
end

Mese=num2str(MM_ord(i), '%02.f'); 
Giorno=num2str(DD_ord(i), '%02.f'); 
Ora=num2str(hi_ord(i), '%02.f'); 
Minuto=num2str(mi_ord(i), '%02.f');
Secondo_int=fix(si_ord(i));                             % Parte intera secondi
Secondo_dec=fix(1000*(si_ord(i)-Secondo_int));          % Parte con 3 decimali dei secondi
str_Secondo_int=num2str(Secondo_int,'%02.f');           % Parte intera secondi con zero davanti
str_Secondo_dec=num2str(Secondo_dec,'%03.f');           % Prima cifra della parte decimale 
Secondo=strcat(str_Secondo_int,str_Secondo_dec);        % Stringa dei secondi nel formato XX XXX (precisione 1/1000 di secondo) 

Data_s=strcat('C',Lettera,cifre_anno,Mese,Giorno,':',Ora,Minuto,Secondo);

%%% Formattazione AR

% AR in gradi
AR=AR_ord(i);    % Ascensione retta J2000

% Formattazione dei decimali di AR
AR_int=fix(AR);                             % Parte intera AR
AR_dec=fix(10000*(AR-AR_int));              % Parte con 4 cifre decimali
str_AR_int=num2str(AR_int,'%03.f');         % Parte intera gradi con zero davanti
str_AR_dec=num2str(AR_dec,'%04.f');         % Prime 4 cifre della parte decimale

AR=strcat(str_AR_int, '.',str_AR_dec);   % Numero stringa nel formato XXX.XXXXXX

%%% Formattazione Dec

if DEC_ord(i) >= 0
    segno='+';
else
    segno='-';
end

DECx=abs(DEC_ord(i)); % Eliminazione del segno della declinazione

% Formattazione Dec
DEC_int=fix(DECx);                                % Parte intera gradi
DEC_dec=fix(100000*(DECx-DEC_int));               % Parte con 5 decimali 
str_DEC_int=num2str(DEC_int,'%02.f');             % Parte intera secondi d'arco con zero davanti
str_DEC_dec=num2str(DEC_dec,'%05.f');             % Prima cifra della parte decimale dei secondi d'arco
DEC=strcat(str_DEC_int, '.',str_DEC_dec);         % Numero stringa nel formato XX.XXXXX

ANNO1=strcat('KC', num2str(YYYY_ord(i)));

L_sat_name = strlength(num2str(NORAD));      % Cifre del numero Norad del satellite

if L_sat_name==1
NameSat=strcat('000000', num2str(NORAD));    % Numero NORAD del satellite nel formato MPC (deve avere 7 cifre)  
end

if L_sat_name==2
NameSat=strcat('00000', num2str(NORAD));     % Numero NORAD del satellite nel formato MPC (deve avere 7 cifre)  
end

if L_sat_name==3
   NameSat=strcat('0000', num2str(NORAD));   % Numero NORAD del satellite nel formato MPC (deve avere 7 cifre)  
end

if L_sat_name==4
NameSat=strcat('000', num2str(NORAD));       % Numero NORAD del satellite nel formato MPC (deve avere 7 cifre)  
end

if L_sat_name==5
NameSat=strcat('00', num2str(NORAD));        % Numero NORAD del satellite nel formato MPC (deve avere 7 cifre)  
end

if L_sat_name==6
NameSat=strcat('0', num2str(NORAD));         % Numero NORAD del satellite nel formato MPC (deve avere 7 cifre)  
end

if L_sat_name==7
NameSat=num2str(NORAD);                      % Numero NORAD del satellite nel formato MPC (deve avere 7 cifre)  
end

fprintf(fid0, '     %s  %s%s    %s%s            %s      %s\n', NameSat, Data_s, AR, segno, DEC, '15.5 G', MPC);
 
end

fclose(fid0);

end
