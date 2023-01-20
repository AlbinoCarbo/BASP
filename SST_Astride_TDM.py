# Python script for extracting header and satellite tracks from WCS fits files
# with ASTRiDE software (Automated Streak Detection for Astronomical Images)
# https://github.com/dwkim78/ASTRiDE
#
# Salva i dati dell'header e delle coordinate AR e DEC misurate del centro della traccia del
# satellite nel file "Data_headers_streaks.txt". Dopo ogni dato aggiunge una virgola, in modo che il riconoscimento del
# comando "readtable" di Matlab proceda senza intoppi. Ogni immagine può contenere più tracce satellitari,
# vengono lette tutte e assegnate allo stesso satellite. Sarà poi l'analisi del TDM finale a dover
# eliminare le misure astrometriche che non c'entrano con il numero Norad del satellite indicato nell'header.
# Se non trova un file passa a quello successivo senza interrompere l'eleborazione. 
# 
# A partire dalla versione del 6 luglio 2022 calcola le coordinate RA e DEC del centro 
# della traccia usando una funzione della libreria "Track_best_fit.py". 
#
# Albino Carbognani, INAF-OAS
# Versione del 4 Agosto 2022

# Import astropy.io library.
from astropy.io import fits

# Import astropy.WCS library.
from astropy.wcs import WCS

# Import the ASTRiDE library. 
from astride import Streak 

# Import numeric Python library
import numpy as np

# Importa libreria per input multipli da riga di comando
import sys

# Importa libreria per lavorare con i path dei file
import os.path

# Importa funzione di best fit per le tracce
import Track_best_fit as tbf

# Parametri di input:
#
# Nome script, SST_Astride_TDM.py
# path0, path della cartella con le immagini calibrate WCS (esempio: path0='/home/albino/Test/')
# name = parte comune nome file fit (esempio: SST20201102_WCS_ )
# ext = estensione (esempio: .fit)
# Ni = numero iniziale immagine da analizzare (esempio: 112)
# num_im = numero immagini da analizzare (esempio: 8)
# soglia = soglia di sensibilità di ASTRiDE (1 per i satelliti deboli, 2 o 3 per quelli brillanti)
# Esempio di input da riga di comando: > python3 SST_Astride_TDM.py /home/albino/Test/ SST20201102_WCS_ .fit 112 8 1

print('                                                                      ')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('%               FITS.WCS2HEADERS.ARDEC - SST PROJECT OAS             %')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('%                                                                    %')
print('%                  by Albino Carbognani (INAF-OAS)                   %')
print('%                         Linux Version                              %')
print('%                           Aug 2022                                 %')
print('%                                                                    %')
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('   \n')

# Input dei dati da riga di comando
nome_script, path0, name, ext, Ni, num_im, soglia=sys.argv  

# Estrazione header e tracce dei satelliti dalle immagini WCS
print('HEADERS AND STREAKS SATELLITES EXTRACTION   \n')

with open(path0+'Data_headers_streaks.txt', 'w') as g:
   for j in range(0, int(num_im)):

        num_file=str(j+int(Ni))

        file_to_open=path0+name+num_file+ext

        # Verifica l'esistenza dell'immagine da cui estrarre le tracce
        if os.path.isfile(file_to_open):

              print('Extract satellite streak from ' + file_to_open + '\n')
   
              # Read a fits image and create a Streak instance.
              #streak = Streak(file_to_open, area_cut=50, contour_threshold=3.0, shape_cut=0.07)
              #streak = Streak(file_to_open, area_cut=700, shape_cut=0.40)
              streak = Streak(file_to_open, area_cut=600, contour_threshold=float(soglia))

              # Detect streaks.
              streak.detect()
   
              # Write outputs and save figures.
              streak.write_outputs() 
              streak.plot_figures() 

              # Save in a file the best fit coordinates RA and DEC of all the tracks
              with open(path0 + name + num_file + '/streaks_center.txt', 'w') as ii:
                  ii.write('#    RA (deg)        DEC (deg)   \n')
                  N_tracks=len(streak.streaks) # Numero delle tracce rilevate nell'immagine
                  head = fits.getheader(file_to_open)
                  w = WCS(head)   # Legge costanti WCS nell'header dell'immagine
                  if N_tracks >= 1:                  
                      for jj in range(N_tracks):
                          # Track inclination from ASTRiDE (radiant)
                          theta=(streak.streaks[jj]['slope']) # Serve quando si usa la funzione tbf.track_center(x, y, phi)

                          # Compute and save best fit coordinates of the tracks's center in RA and DEC
                          X, Y = tbf.track_center2(streak.streaks[jj]['x'], streak.streaks[jj]['y'])
                          ra, dec = w.wcs_pix2world(X, Y, 1)          # Trasforma da pixel a RA e DEC (gradi)
                          coordinates=str(ra)+','+" "+str(dec)+'\n'   # Coordinate di best fit del centro tracce rivelate nell'immagine
                          ii.write(coordinates) 
           
              # Estrazione delle coordinate delle tracce di tutti i satelliti trovati da ASTRiDE sull'immagine
              path = path0 + name + num_file + '/streaks_center.txt' # Path del file con i dati del centro delle tracce
              with open(path) as f:
                    dati_streak = f.readline()
                    while dati_streak:
                          dati_streak = f.readline()
                          coordinates2=dati_streak[0:37] # Estrazione AR e DEC (gradi) dalla n-esima riga del file streaks.txt
                          
                          # Se dati_streaks è una stringa vuota perché la traccia non è stata 
                          # riconosciuta da ASTRiDE mette NaN in AR e DEC
                          if not dati_streak:
                                 coordinates2='NaN'+','+'\t'+'NaN' # il carattere \t equivale a un tab

                          hdul = fits.open(file_to_open) # Legge header del file in modo da assegnare le stesse keys alle diverse tracce
                          # Salva data e ora in Data_headers_streaks.txt
                          g.write(hdul[0].header['DATE-OBS']+','+' ')
                          # Salva nome oggetto in Data_headers_streaks.txt
                          g.write(str(hdul[0].header['OBJECT'])+','+' ')
                          # Salva tempo di esposizione (s) in Data_headers_streaks.txt
                          g.write(str(hdul[0].header['EXPTIME'])+','+' ')
                          # Save center coordinates satellite's streak (degree) in Data_headers_streaks.txt
                          g.write(coordinates2 +'\n')
                          

        else:
              print(file_to_open + ' does not exist ' + '\n')
              continue # Se il file non esiste passa a quello successivo
 
# Chiusura del file Data_headers_streaks.txt con keys header e coordinate del centro delle tracce dei satelliti
g.close()



