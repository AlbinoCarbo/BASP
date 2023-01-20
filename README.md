# BASP
The main purposes of BASP (Bologna Astrometry Satellites Pipeline), is the automatic astrometric calibration of the images and the extraction of the satellite's traces to obtain a text file with the J2000 astrometric positions of the center of the track in Tracking Data Message. The software was developed under Linux Xubuntu, but can run under any Linux system as long as the following software and their dependencies are installed and running:

1-MATLAB R2019b
2-Python 3 (https://www.python.org/download/releases/3.0/) with Astropy package 
3-Astrometry.net (Plate solve of fits images with 2MASS catalog, http://astrometry.net/)
4-ASTRiDE (Automated Streak Detection for Astronomical Images, https://github.com/dwkim78/ASTRiDE)
5-Find Orb (Geocentric orbit determination from observations, https://www.projectpluto.com/find\_orb.htm)

BASP consists of a main script running under MATLAB (MathWorks) that launches both external software as Astrometry.net (Lang et al., 2010) or Find Orb and the Python
3 scripts needed to process the images in fits format taken at the telescope with the CCD camera. The pipeline core is ASTRiDE-based (Dae-Won et al., 2005). ASTRiDE detect the streaks in astronomical images using a contour tracing algorithm for each object and their morphological parameters as area and perimeter. This software is specially designed for the detection of tracks left by fast-moving objects such as meteors, satellites and near-Earth objects (NEOs), so it is very suitable for our case. Finally, a geocentric orbit is automatically computed for each satellite using Find Orb, a popular orbit determination software under Windows or Linux for fitting
orbits of solar system objects and satellites, written by the US American amateur astronomer Bill Gray. For a terrestrial satellite Find Orb takes into account the gravitational perturbations of Sun and Moon, non-sphericity of the Earth, pressure of solar radiation and atmospheric drag. It does not consider relativistic effects except for bodies in heliocentric orbit. As specified later, the non-interactive version of Find Orb - known as fo - was used in the pipeline.
