# Python script for calibrating fits file using master_bias created by "fits_master_bias.py".
# Subtract master bias from fits images and save them with the same name by appending 'cal' to the end
# If the file to be calibrated is missing, go to the next one.
#
# Albino Carbognani, INAF-OAS
# Versione del 18 dicembre 2020

# Import astropy.io library
from astropy.io import fits

# Load python library used for working with arrays
import numpy as np

# Importa libreria per input multipli da riga di comando
import sys

# Importa libreria per lavorare con i path dei file
import os.path

# Parametri di input:
#
# Nome script, fits_calibrazione.py
# path0, path della cartella con le immagini dei bias (esempio: path0='/home/albino/Test/')
# name = parte comune nome file fit dei bias (esempio: SST20201102_ )
# ext = estensione (esempio: .fit)
# Ni = numero iniziale immagine da calibrare (esempio: 101)
# num_im = numero immagini da calibrare (esempio: 10)
# Esempio di input da riga di comando: > python3 fits_calibrazione.py /home/albino/Test/ SST20201102_ .fit 101 10


# Input dei dati da riga di comando
nome_script, path0, name, ext, Ni, num_im=sys.argv  

master_bias=fits.getdata(path0+'master_bias.fit', ext=0)

# Ciclo di calibrazione
# Se l'immagine manca passa a quella successiva
for i in range(0, int(num_im)):

     num_file=str(i+int(Ni))

     file_to_open=path0+name+num_file+ext

     # Verifica l'esistenza del file
     if os.path.isfile(file_to_open): 

         data=fits.getdata(file_to_open, ext=0)
         head=fits.getheader(file_to_open, ext=0)

         data_calibrated=data-master_bias

         # Salvataggio immagine calibrata
         fits.writeto(path0+name+num_file+'_cal'+ext, data_calibrated, head, overwrite=True)
     else:
         continue # Se il file non esiste passa a quello successivo

