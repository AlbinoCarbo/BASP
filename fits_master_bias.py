# Python script for creating a master bias with individual fits bias.
# If a bias file doesn't exist go to the next one.
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
# Nome script, fits_master_bias.py
# path0, path della cartella con le immagini dei bias (esempio: path0='/home/albino/Test/')
# name = parte comune nome file fit dei bias (esempio: SST20201102_ )
# ext = estensione (esempio: .fit)
# Ni = numero iniziale immagine da analizzare (esempio: 101)
# num_im = numero immagini da analizzare (esempio: 10)
# Esempio di input da riga di comando: > python3 fits_master_bias.py /home/albino/Test/ SST20201102_ .fit 101 10


# Input dei dati da riga di comando
nome_script, path0, name, ext, Ni, num_im=sys.argv  

# Inizializzazione variabile di tipo lista con gli ADU della prima immagine
data0=fits.getdata(path0+name+Ni+ext, ext=0)
fitslist=[data0]

# Copia header della prima immagine (che sicuramente esiste sempre)
head=fits.getheader(path0+name+Ni+ext, ext=0)

# Incremento variabile di tipo lista con le restanti immagini bias
for i in range(0, int(num_im)-1):

     num_file=str(i+int(Ni)+1)

     file_to_open=path0+name+num_file+ext

     # Verifica l'esistenza del file
     if os.path.isfile(file_to_open):

           data=fits.getdata(path0+name+num_file+ext, ext=0)

           # Aggiunta dell'i-esima immagine nella variabile lista
           fitslist.insert(i+1, data)

     else:
           continue # Se il file non esiste passa a quello successivo

# Creazione master bias
median_image=np.median(fitslist, axis=0)

# Salvataggio master bias
fits.writeto(path0+'master_bias.fit', median_image, head, overwrite=True)

