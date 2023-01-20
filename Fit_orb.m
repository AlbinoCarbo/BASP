% Fit_orb: script for the orbital best fit of an artificial satellite using Bill Gray's find_orb 
% and to delete the astrometric observations with a residual greater than DELTA arcsec.

% Using this function it is possible to identify observations that deviate too much from the 
% orbital best fit and eliminate them.
%
% ALGORITMO:
% Ci sono due livelli di best fit orbitale per l'eliminazione delle misure 
% astrometriche bad. 
% Nel primo livello si eliminano le osservazioni astrometriche la cui distanza 
% dal best fit è superiore ai 10". Dopo l'eliminazione di queste osservazioni 
% si fa un nuovo best fit dell'orbita e si eliminano le osservazioni astrometriche 
% che distano più di DELTA arcsec dal best fit. Le osservazioni astrometriche che restano sono
% l'output della funzione, insieme al tempo di ripresa, e vanno a finire
% nel TDM di output della pipeline in cui lavora Fit_orb.
%
% INPUT:
% data_path = path delle immagini SST
% DELTA = residuo massimo dell'osservazione astrometrica (arcsec)
% name_sat = numero Norad satellite 
% YYYY = vettore colonna anno dell'osservazione astrometrica
% MM = vettore colonna mese dell'osservazione astrometrica
% DD = vettore colonna giorno dell'osservazione astrometrica
% hi = vettore colonna ore dell'osservazione astrometrica
% mi = vettore colonna minuti dell'osservazione astrometrica
% AR = vettore colonna colonna ascensione retta osservata (gradi)
% DEC = vettore colonna declinazione osservata (gradi)
%
% OUTPUT:
% YYYY2 = vettore riga anno dell'osservazione astrometrica di best fit
% MM2 = vettore riga mese dell'osservazione astrometrica di best fit
% DD2 = vettore riga giorno dell'osservazione astrometrica di best fit
% hi2 = vettore riga ora dell'osservazione astrometrica di best fit
% mi2 = vettore riga minuti dell'osservazione astrometrica di best fit
% AR2 = vettore riga ascensione retta di best fit (gradi)
% DEC2 = vettore riga declinazione di best fit (gradi)
% S = angolo totale percorso in cielo dal satellite (gradi)
%
% Nella cartella delle immagini SST elaborate vengono salvati anche tre
% file di testo Aux che servono come controllo:
%
% A) Il primo contiene l'astrometria originale di input del satellite, gli elementi
% orbitali e i residui che si ottengono con find_orb. 
%
% B) Il secondo contiene l'astrometria del satellite tolte le osservazioni che
% distano più di 10" dal best fit orbitale, i relativi elementi e i residui.
%
% C) Il terzo contiene l'astrometria del satellite tolte le osservazioni sia che
% distano più di 10" dal primo best fit orbitale, sia più di DELTA arcsec dal 
% secondo best fit orbitale, i relativi elementi e i residui. L'astrometria di 
% quest'ultimo file è quella che viene salvata nel TDM finale.
%
% Nello script esistono degli if-end di controllo. Se con il fit orbitale delle 
% osservazioni astrometriche di input l'orbita ottenuta ha e >= 1 vuol dire che 
% le osservazioni necessitano di find_orb manuale, quindi lo script esce dalla funzione 
% riportando le posizioni astrometriche senza filtrarle.
%
% Lo stesso accade se, dopo avere tolto le osservazioni con residui > di
% 10", restano meno di 3 osservazioni: lo script esce dalla funzione 
% riportando le posizioni astrometriche senza filtrarle automaticamente. Se
% le osservazioni astrometriche superstiti sono >= 3, ma l'orbita ha e >= 1
% lo script esce riportando le osservazioni di input originarie senza
% filtrarle.
%
% Albino Carbognani, INAF-OAS
% Versione del 10 nov 2021

function [YYYY2, MM2, DD2, hi2, mi2, si2, AR2, DEC2, S]=Fit_orb(data_path, DELTA, name_sat, YYYY, MM, DD, hi, mi, si, AR, DEC)

% Def nomi file astrometrici ausiliari

% Nome del file astrometrico originale nel formato MPC da analizzare con find_orb
Name1=strcat('Aux_astrometry_MPC1_', num2str(name_sat), '.txt'); 

% Nome del secondo file astrometrico, tolte le very bad observations, nel formato MPC 
% da analizzare con find_orb
Name2=strcat('Aux_astrometry_MPC2_', num2str(name_sat), '.txt'); 

% Nome del file astrometrico finale nel formato MPC 
Name3=strcat('Aux_astrometry_MPC3_', num2str(name_sat), '.txt'); 

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%                      Start Fit_orb script                       %')
disp('%                            Nov 2021                             %')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(strcat('                    Satellite', {' '}, num2str(name_sat)))
disp('   ')

% Numero di osservazioni fatte
N=length(DEC);

% Calcolo dell'angolo originale percorso in cielo dal satellite (gradi)
S=acosd(sind(DEC(1))*sind(DEC(N))+cosd(DEC(1))*cosd(DEC(N))*cosd(AR(1)-AR(N)));

disp(strcat('Fit_orb: observations number for satellite', {' '}, num2str(name_sat), ':', {' '}, num2str(N)))

% Salvataggio file astrometrico nel formato MPC esteso per l'analisi con find_orb
Astrometry_MPC2(data_path, name_sat, YYYY, MM, DD, hi, mi, si, AR, DEC, Name1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1-Orbital best fit with find_orb to eliminate measurements that are more than 10 arcsec from the orbit %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define command line for find_orb
command_line_fo=strcat('fo', " ", strcat(data_path, Name1), " ",'-h');

% Execute find_orb
disp('     Fit_orb: 1-EXECUTE FIND_ORB NON-INTERACTIVE VERSION')
disp('   ')
system(command_line_fo);
disp('   ')

% Read file total.json genereted by fo (non-interactive find_orb version)
filetext = fileread('~/.find_orb/total.json');

% String conversion
Str=string(filetext);

% Extract orbital elements and append in file Name1
orbit1=extractBetween(Str, '"elements":', '"MOIDs":');
fid1 = fopen(strcat(data_path, Name1), 'a+');
fprintf(fid1, '%s\n', orbit1);
fclose(fid1);

% Extract central body
body1=extractBetween(Str, '"central body":', ',');
tf1 = strcmp(body1, ' "Sun"');
if tf1 == 1
disp(strcat('     Fit_orb WARNING about central body:'," ", body1))
disp('   ')
end

% Extract RA residual values (arcsec)
RA_residual = str2double(extractBetween(Str, '"dRA" :', ','));
 
% Extract DEC residual values (arcsec)
DEC_residual = str2double(extractBetween(Str, '"dDec": ', ','));

% Compute total residual values (arcsec)
Total_residual=sqrt(RA_residual.^2+DEC_residual.^2);

% Save residuals and append in file Name1
fid1 = fopen(strcat(data_path, Name1), 'a+');
fprintf(fid1, 'dRA (")  dDEC (")  Tot res (")\n');

% Esce se il numero dei residui è inferiore alle osservazioni dopo il primo best fit orbitale
if  length(Total_residual) < N 
    disp('     Fit_orb WARNING: Residuals number lower than observations number after the first orbital fit: exit without filtering')
    YYYY2=YYYY; MM2=MM; DD2=DD; hi2=hi; mi2=mi; si2=si; AR2=AR; DEC2=DEC;
    return
end

for i=1:N
fprintf(fid1, '%03.2f       %03.2f       %03.2f\n', abs(RA_residual(i)), abs(DEC_residual(i)), Total_residual(i));
end

fclose(fid1);

mean1=num2str(mean(Total_residual));

disp(strcat('Fit_orb: Mean residual astrometry original set for', {' '}, num2str(name_sat), {' '}, '(arcsec):', {' '}, mean1));

% Esce se l'orbita non è ellittica dopo il primo best fit orbitale
eccentricity=str2double(extractBetween(Str, '"e":  ', ','));
if  eccentricity >= 1 
    disp('     Fit_orb WARNING - Non elliptical orbit after the first orbital fit: exit without filtering')
    YYYY2=YYYY; MM2=MM; DD2=DD; hi2=hi; mi2=mi; si2=si; AR2=AR; DEC2=DEC;
    return
end

% Purge very bad observations
k=1;
for i=1:N

    if Total_residual(i) < 10
        anno(k)=YYYY(i);
        mese(k)=MM(i);
        giorno(k)=DD(i);
        ora(k)=hi(i);
        min(k)=mi(i);
        sec(k)=si(i);
        ar(k)=AR(i);
        dec(k)=DEC(i);
        k=k+1;
    end

end
   
% Esce se ci sono meno di tre osservazioni dopo il primo best fit orbitale
if k < 3
    disp('     Fit_orb WARNING - Too few observations after the first orbital filter: exit without filtering')
    YYYY2=YYYY; MM2=MM; DD2=DD; hi2=hi; mi2=mi; si2=si; AR2=AR; DEC2=DEC;
    return
end

% Calcolo percentuale di immagini eliminate
Percent=round(100*(N-length(dec))/N, 1);
disp(strcat('Fit_orb: deleted observations after the first orbital filter:', {' '}, num2str(Percent),'%'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2-Orbital best fit with find_orb to eliminate measurements that are more than DELTA arcsec away from the orbit %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Numero di osservazioni rimaste dopo il primo filtro
N2=length(dec);

% Salvataggio file astrometrico nel formato MPC esteso per l'analisi con find_orb
Astrometry_MPC2(data_path, name_sat, anno, mese, giorno, ora, min, sec, ar, dec, Name2);

% Define command line for find_orb
command_line_fo2=strcat('fo', " ", strcat(data_path, Name2), " ",'-h');

% Execute find_orb
disp('     Fit_orb: 2-EXECUTE FIND_ORB NON-INTERACTIVE VERSION')
disp('   ')
system(command_line_fo2);
disp('   ')
% Read file total.json genereted by fo (non-interactive find_orb version)
filetext2 = fileread('~/.find_orb/total.json');

% String conversion
Str2=string(filetext2);

% Extract orbital elements and append in file Name2
orbit2=extractBetween(Str2, '"elements":', '"MOIDs":');
fid2 = fopen(strcat(data_path, Name2), 'a+');
fprintf(fid2, '%s\n', orbit2);
fclose(fid2);

% Extract central body
body2=extractBetween(Str2, '"central body":', ',');
tf2 = strcmp(body2, ' "Sun"');
if tf2 == 1
disp(strcat('     Fit_orb WARNING about central body:'," ", body2))
disp('   ')
end

% Extract RA residual values (arcsec)
RA_residual2 = str2double(extractBetween(Str2, '"dRA" :', ','));
 
% Extract DEC residual values (arcsec)
DEC_residual2 = str2double(extractBetween(Str2, '"dDec": ', ','));

% Compute total residual values (arcsec)
Total_residual2=sqrt(RA_residual2.^2+DEC_residual2.^2);

% Save residuals and append in file Name2
fid2 = fopen(strcat(data_path, Name2), 'a+');
fprintf(fid2, 'dRA (")  dDEC (")  Tot res (")\n');

% Esce se il numero dei residui è inferiore alle osservazioni astrometriche dopo il secondo best fit orbitale
if  length(Total_residual2) <  N2
    disp('     Fit_orb WARNING: Residuals number lower than observations number after the second orbital fit: exit without filtering')
    YYYY2=YYYY; MM2=MM; DD2=DD; hi2=hi; mi2=mi; si2=si; AR2=AR; DEC2=DEC;
    return
end

for i=1:N2
fprintf(fid2, '%03.2f       %03.2f       %03.2f\n', abs(RA_residual2(i)), abs(DEC_residual2(i)), Total_residual2(i));
end

fclose(fid2);

mean2=num2str(mean(Total_residual2));

disp(strcat('Fit_orb: mean residual astrometry second set for', {' '}, num2str(name_sat), {' '}, '(arcsec):', {' '}, mean2));

% Esce se l'orbita non è ellittica dopo il secondo best fit orbitale
eccentricity2=str2double(extractBetween(Str2, '"e":  ', ','));
if  eccentricity2 >= 1 
    disp('     Fit_orb WARNING: non elliptical orbit after the second orbital fit: exit without filtering')
    YYYY2=YYYY; MM2=MM; DD2=DD; hi2=hi; mi2=mi; si2=si; AR2=AR; DEC2=DEC;
    return
end

% Purge bad observations
k=1;
for i=1:N2

    if Total_residual2(i) < DELTA
        YYYY2(k)=anno(i);
        MM2(k)=mese(i);
        DD2(k)=giorno(i);
        hi2(k)=ora(i);
        mi2(k)=min(i);
        si2(k)=sec(i);
        AR2(k)=ar(i);
        DEC2(k)=dec(i);
        k=k+1;
    end

end

% Esce se ci sono meno di tre osservazioni dopo il secondo best fit orbitale
if k < 3
    disp('Fit_orb WARNING - Too few observations after the second orbital fit, increase DELTA: exit without filtering')
    YYYY2=YYYY; MM2=MM; DD2=DD; hi2=hi; mi2=mi; si2=si; AR2=AR; DEC2=DEC;
    return
end

% Calcolo percentuale di immagini eliminate
Percent2=round(100*(N-length(DEC2))/N, 1);
disp(strcat('Fit_orb: deleted observations after the second orbital filter:', {' '}, num2str(Percent2), '%'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3-Orbital best fit with find_orb of best astrometric measurements %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Numero di osservazioni rimaste
N3=length(DEC2);

% Salvataggio file astrometrico finale nel formato MPC esteso
Astrometry_MPC2(data_path, name_sat, YYYY2, MM2, DD2, hi2, mi2, si2, AR2, DEC2, Name3);

% Define command line for find_orb
command_line_fo3=strcat('fo', " ", strcat(data_path, Name3), " ",'-h');

% Execute find_orb
disp('     Fit_orb: 3-EXECUTE FIND_ORB NON-INTERACTIVE VERSION')
disp('   ')
system(command_line_fo3);
disp('   ')
% Read file total.json genereted by fo (non-interactive find_orb version)
filetext3 = fileread('~/.find_orb/total.json');

% String conversion
Str3=string(filetext3);

% Extract orbital elements and append in file Name3
orbit3=extractBetween(Str3, '"elements":', '"MOIDs":');
fid3 = fopen(strcat(data_path, Name3), 'a+');
fprintf(fid3, '%s\n', orbit3);
fclose(fid3);

% Extract central body
body3=extractBetween(Str3, '"central body":', ',');
tf3 = strcmp(body3, ' "Sun"');
if tf3 == 1
disp(strcat('     Fit_orb WARNING about central body:'," ", body3))
disp('   ')
end

% Extract RA residual values (arcsec)
RA_residual3 = str2double(extractBetween(Str3, '"dRA" :', ','));
 
% Extract DEC residual values (arcsec)
DEC_residual3 = str2double(extractBetween(Str3, '"dDec": ', ','));

% Compute total residual values (arcsec)
Total_residual3=sqrt(RA_residual3.^2+DEC_residual3.^2);

% Save residuals and append in file Name3
fid3 = fopen(strcat(data_path, Name3), 'a+');
fprintf(fid3, 'dRA (")  dDEC (")  Tot res (")\n');

for i=1:N3
fprintf(fid3, '%03.2f       %03.2f       %03.2f\n', abs(RA_residual3(i)), abs(DEC_residual3(i)), Total_residual3(i));
end

fclose(fid3);

mean3=num2str(mean(Total_residual3));

disp(strcat('Fit_orb: Mean residual astrometry final set for', {' '}, num2str(name_sat), {' '}, '(arcsec)', {' '}, mean3));
disp('   ')
disp(strcat('Fit_orb: total deleted observations:', {' '}, num2str(Percent2), '%'))

% Calcolo dell'angolo percorso in cielo dal satellite (gradi)
S=acosd(sind(DEC2(1))*sind(DEC2(N3))+cosd(DEC2(1))*cosd(DEC2(N3))*cosd(AR2(1)-AR2(N3)));

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%                      End Fit_orb script                         %')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('   ')
end
