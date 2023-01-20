% Astrometry_MPC, procedure for saving astrometric observations in a text 
% file in the extended MPC format which has a time precision of 86.4/1000 seconds. 
% In this way astrometric observations can also be analyzed manually by Bill Gray's 
% Find_orb (www.projectpluto.com).
%
% For a format with a time precision of 1/1000 s you can also use
% the "Astrometry_MPC2.m" script.
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

function []=Astrometry_MPC(data_pathy, NORAD, YYYY_ord, MM_ord, DD_ord, hi_ord, mi_ord, si_ord, AR_ord, DEC_ord, name_astrometry_file)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saving astrometry satellites in the extended MPC format %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(strcat('Astrometry MPC: save', {' '}, num2str(NORAD), {' '}, 'astrometry in MPC format'))
disp('  ')

% MPC='598';                % Codice MPC della stazione

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

fprintf(fid0, '     %s %s %02d %s%02d %02d %s%s%02d %02d %s         %s      %s\n', NameSat, ANNO1, MM_ord(i), Data_s, AR(1), AR(2), AR3, segno, DEC(1), DEC(2), DEC3, '15.5 G', MPC);
 
end

fclose(fid0);

end
