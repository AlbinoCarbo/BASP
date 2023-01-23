% Script per la cancellazione automatica delle osservazioni astrometriche "bad"
% all'interno dei TDM. Necessita del numero delle righe di AR e DEC da
% cancellare determinate con find_orb, dopo di che determina le corrispondenti 
% righe di AR e DEC dentro al TDM e le cancella. Salva il nuovi file aggiungendo _purged nel nome.
%
% Albino Carbognani, INAF-OAS
% Versione Linux del 28 settembre 2021

clear all

%%% INIZIO INPUT

% Numero delle righe di osservazioni astrometriche da cancellare determinate con find_orb
%bad_astrometry=[1 8 10 12 14 16 18 20 30 34 38 41 42 47 55 57 59]; 
bad_astrometry=[4 8 11 44 45 50 57];

% Definizione file TDM di input
data_path='/home/albino/Documents/MATLAB/SST_Astrometry_purge/';
data_name='IT_CASSINI_20220328T190103_41859';

%%% FINE INPUT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numero delle righe da cancellare nel file del TDM %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rows_delete1=20+3*bad_astrometry-2; % Righe dell'AR
rows_delete2=20+3*bad_astrometry-1; % Righe della DEC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lettura del TDM iniziale %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = fopen(strcat(data_name, '.tdm'));
data = textscan(f, '%s', 'Delimiter', '\n');
fclose(f);
data = data{1};

% Numero di righe del TDM
N=length(data);

% Definizione del nome del TDM di output 
output_file=strcat(data_path, '/', data_name, '_purged.tdm');

% Apertura del TDM di output 
fid1 = fopen(output_file, 'w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ciclo di eliminazione delle misure astrometriche errate %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:N
    [M1,I1] = min(abs(rows_delete1-i)); % Se c'� un elemento nullo vuol dire che la misura in AR va scartata
    [M2,I2] = min(abs(rows_delete2-i)); % Se c'� un elemento nullo vuol dire che la misura in DEC va scartata
      if (M1 > 0) && (M2 > 0) 
      % Salvataggio della stringa con le misure astrometriche buone nel file di output
      fprintf(fid1, strcat(char(data(i)), '\n'));
      end
end

% Chiusura del nuovo TDM di output ripulito dalla astrometria errata
status = fclose(fid1);

% Elimina le righe bianche multiple lasciate dalla cancellazione delle
% misure errate
new_text = regexprep( fileread(output_file), '(\r?\n)(\r?\n)+', '\r \n');
%write the result to the same file
fid2 = fopen(output_file, 'w');
fwrite(fid2, new_text);
fclose(fid2);