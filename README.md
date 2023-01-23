# BASP
The main purposes of BASP (Bologna Astrometry Satellites Pipeline), is the automatic astrometric calibration of the images and the extraction of the satellite's traces to obtain a text file with the J2000 astrometric positions of the center of the track in Tracking Data Message. The software was developed under Linux Xubuntu, but can run under any Linux system as long as the following software and their dependencies are installed and running:

1-MATLAB R2019b

2-Python 3 (https://www.python.org/download/releases/3.0/) with Astropy package 

3-Astrometry.net (Plate solve of fits images with 2MASS catalog, http://astrometry.net/)

4-ASTRiDE (Automated Streak Detection for Astronomical Images, https://github.com/dwkim78/ASTRiDE)

5-Find Orb (Geocentric orbit determination from observations, non-interactive version, https://www.projectpluto.com/find\_orb.htm)

BASP requires three settings files: Settings_BASP.txt, is the main configuration file with the most frequently changing settings such as the name and path of the images to be analysed, the type of analysis and the intensity of the satellite track; Settings_Astrometry.txt, is the configuration file for Astrometry.net, here there are the upper and lower limit of the image scale in arcsec/pixel; finally there is Settings_MPC.txt, where there is the header of the output file with the astrometric measurements in the Minor Planet Center format. These last two files contain quantities that are rarely changed so, once you save the settings, you can forget about their existence. The settings files are self-explanatory and with sample settings, so we don't need to go through them in detail.
