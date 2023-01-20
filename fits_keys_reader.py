# Python script for reading the header of the fits file.
# Save data: (file name), image date and time, object name, exposure time, image center AR and DEC (J2000)
# in the "Data_Keys.txt" file.
# If the file header doesn't have the keys required by the pipeline, delete the file and move on to the next one
# If the file to read keys from doesn't exist go to the next one
#
# Albino Carbognani, INAF-OAS
# Versione del 26 febbraio 2021

# Import astropy.io library.
from astropy.io import fits

# Importa libreria per input multipli da riga di comando
import sys

# Importa libreria per lavorare con i path dei file
import os.path

# Importa libreria per lavorare sui file
import os

# Parametri di input:
#
# Nome script, fits_keys_reader.py
# path0, path della cartella con le immagini calibrate WCS (esempio: path0='/home/albino/Test/')
# name = parte comune nome file fit (esempio: SST20201102_WCS_ )
# ext = estensione (esempio: .fit)
# Ni = numero iniziale immagine da analizzare (esempio: 112)
# num_im = numero immagini da analizzare (esempio: 8)
# Esempio di input da riga di comando: > python3 fits_keys_reader.py /home/albino/Test/ SST20201102_WCS_ .fit 112 8

# Input dei dati da riga di comando
nome_script, path0, name, ext, Ni, num_im=sys.argv  

# Estrazione key header fits
print('KEYS EXTRACTION FROM HEADER FITS   \n')
with open(path0+'Data_keys.txt', 'w') as f:

  for i in range(0, int(num_im)):

        num_file=str(i+int(Ni))

        file_to_open=path0+name+num_file+ext

        # Verifica l'esistenza del file
        if os.path.isfile(file_to_open):

              print('Extract header keys from ' + file_to_open +'\n')

              hdul = fits.open(file_to_open) # hdul = Header Data Unit List

              hdr = hdul[0].header
              if ('DATE-OBS' in hdr and 'OBJECT' in hdr and 'EXPTIME' in hdr and 'RA' in hdr and 'DEC' in hdr): # Check for existence

                             # Data e ora
                             f.write(hdul[0].header['DATE-OBS']+' ')
        
                             # Nome oggetto
                             f.write(str(hdul[0].header['OBJECT'])+' ')

                             # Tempo di esposizione (s)
                             f.write(str(hdul[0].header['EXPTIME'])+' ')

                             # AR (hh:mm:ss)
                             f.write(str(hdul[0].header['RA'])+' ')

                             # DEC (dd:mm:ss)
                             f.write(str(hdul[0].header['DEC'])+'\n')
              else:
                             print(file_to_open + ' ' + "without correct header!" + '\n')
                             print(file_to_open + ' ' + "will be deleted" + '\n')
                             os.remove(file_to_open) 
                             continue # Se l'header non ha tutte le voci richieste passa a quello successivo

        else:
              continue # Se il file non esiste passa a quello successivo

# Chiusura del file contenente le key dei fits letti
f.close()
