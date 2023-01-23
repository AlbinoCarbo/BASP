# BASP
The main purposes of **BASP (Bologna Astrometry Satellites Pipeline)**, is the automatic astrometric calibration of the fits images and the extraction of the satellite's traces to obtain a text file with the J2000 astrometric positions of the center of the track in **Tracking Data Message** (TDM) format (but also save in MPC extended format). No correction is made to the extracted equatorial coordinates, such as that for annual aberration. The software was developed under Linux Xubuntu, but can run under any Linux system as long as the following external software and their dependencies are installed and running:

1- [MATLAB R2019b](https://it.mathworks.com/products/new_products/release2019b.html). Note: it is not necessary if you are using the compiled version of BASP.

2- [Python 3](https://www.python.org/download/releases/3.0/) with Astropy package. 

3- [Astrometry.net](http://astrometry.net/) Plate solve of fits images with 2MASS catalog.

4- [ASTRiDE](https://github.com/dwkim78/ASTRiDE) Automated Streak Detection for Astronomical Images. ASTRiDE detect the streaks in astronomical images using a contour tracing  algorithm for each object and their morphological parameters as area and perimeter. This software is specially designed for the detection of tracks left by fast-moving objects such as meteors, satellites and near-Earth objects (NEOs), so it is very suitable for our case.

5- [Find Orb](https://www.projectpluto.com/find\_orb.htm) Geocentric orbit determination from observations, [non-interactive version](https://www.projectpluto.com/find_sou.htm). Find_Orb is a popular orbit determination software under Windows or Linux for fitting orbits of solar system objects and satellites, written by the US American amateur astronomer Bill Gray. For a terrestrial satellite Find_Orb takes into account the gravitational perturbations of Sun and Moon, non-sphericity of the Earth, pressure of solar radiation and atmospheric drag.

BASP consists of a main script, **BASP.m**, that run under Matlab, but there is also the compiled version which does not require Matlab in your system but only the freely available [MATLAB Runtimes](https://www.mathworks.com/products/compiler/mcr/index.html) (see readme.txt). The main script launches both external software as Astrometry.net and the Python 3 scripts needed to process the images in fits format taken at the telescope with the CCD camera.

BASP has three settings files: **Settings_BASP.txt**, is the main configuration file with the most frequently changing settings such as the name and path of the images to be analysed, the type of analysis and the intensity of the satellite track; **Settings_Astrometry.txt**, is the configuration file for Astrometry.net, here there are the upper and lower limit of the image scale in arcsec/pixel; finally there is **Settings_MPC.txt**, where there is the header of the output file with the astrometric measurements in the Minor Planet Center format. These last two files contain quantities that are rarely changed so, once you save the settings, you can forget about their existence. The settings files are self-explanatory and with sample settings, so we don't need to go through them in detail.

### Make a test

In the **Sample_Images** folder there are 19 images that can be used as a test. Image numbers 101 to 110 are calibration biases, 113 to 121 contain satellite images. These images are already selected in Settings_BASP.txt, just change the paths to make them conform to your system and you can run BASP with the console command ./BASP (compiled version), or by running the BASP.m script with Matlab. During the processing, messages appear informing about the execution status.
