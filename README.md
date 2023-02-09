# BASP

Suppose to have a sequence of digital images in **fits format**, taken with a **telescope** and **CCD camera**, of satellites and space debris with a known Norad number. From these images we want the astrometry of the satellite in **Tracking Data Message** (TDM) format, to compute its geocentric orbit. Each image was taken with **sidereal tracking active**, therefore the satellite left a track in the image and we want to measure the equatorial coordinates - right ascension and declination - of all the tracks center. 

To solve this problem it is possible to use **BASP (Bologna Astrometry Satellites Pipeline)**. 
The main purposes of BASP, is the **automatic** astrometric calibration of the images and the extraction of the satellite's traks to obtain a text file with the J2000 astrometric positions of the center of the track in TDM format (but also save in Minor Planet Center extended format). No correction is made to the extracted equatorial coordinates, such as that for annual aberration. The software was developed under **Linux Xubuntu**, but can run under any Linux system as long as the following external software and their dependencies are installed and running:

**1-** [MATLAB R2019b](https://it.mathworks.com/products/new_products/release2019b.html). Note: it is not necessary if you are using the compiled version of BASP.

**2-** [Python 3](https://www.python.org/download/releases/3.0/) with Astropy package. 

**3-** [Astrometry.net](http://astrometry.net/) Plate solve of fits images.The astrometric calibration takes place using the 2MASS infrared star catalog [Two Micron All Sky Survey](https://irsa.ipac.caltech.edu/Missions/2mass.html). This catalog contains J2000 position data over 410 million stars, with uncertainties ranging from 0.15 to 0.30 arcsec.

**4-** [ASTRiDE](https://github.com/dwkim78/ASTRiDE) Automated Streak Detection for Astronomical Images. ASTRiDE detect the streaks in astronomical images using a contour tracing  algorithm for each object and their morphological parameters as area and perimeter. This software is specially designed for the detection of tracks left by fast-moving objects such as meteors, satellites and near-Earth objects (NEOs), so it is very suitable for our case.

**5-** [Find Orb](https://www.projectpluto.com/find\_orb.htm) Geocentric orbit determination from observations, [non-interactive version](https://www.projectpluto.com/find_sou.htm). Find_Orb is a popular orbit determination software under Windows or Linux for fitting orbits of solar system objects and satellites, written by the US American amateur astronomer Bill Gray. For a terrestrial satellite Find_Orb takes into account the gravitational perturbations of Sun and Moon, non-sphericity of the Earth, pressure of solar radiation and atmospheric drag.

BASP consists of a main script, **BASP.m**, that run under Matlab, but there is also the compiled version which does not require Matlab in your system but only the freely available [MATLAB Runtimes](https://www.mathworks.com/products/compiler/mcr/index.html) (see readme.txt). The main script launches both external software as Astrometry.net and the Python 3 scripts needed to process the images in fits format taken at the telescope with the CCD camera. Astrometry.net package saves fully astrometrically calibrated fits images by adding the suffix **WCS** to the name. WCS stands for **World Coordinate Systems** and indicates a set of image constants that describe the geometric transformations necessary to move from the pixel position on the image to the RA and DEC coordinates on the celestial sphere and vice versa. These images, calibrated by Astronomy.net, are used by ASTRiDE for the extraction of the satellite tracks.

BASP has three settings files: **Settings_BASP.txt**, is the main configuration file with the most frequently changing settings such as the name and path of the images to be analysed, the type of analysis and the intensity of the satellite track; **Settings_Astrometry.txt**, is the configuration file for Astrometry.net, here there are the upper and lower limit of the image scale in arcsec/pixel; finally there is **Settings_MPC.txt**, where there is the header of the output file with the astrometric measurements in the Minor Planet Center format. These last two files contain quantities that are rarely changed so, once you save the settings, you can forget about their existence. The settings files are self-explanatory and with sample settings, so we don't need to go through them in detail.

![Satellite_track](/Pictures/Galileo_41550.jpg)

*A typical image for satellite astrometry showing the trace, extract by ASTRiDE, of the Galileo NORAD 41550 taken on the evening of Nov 6, 2020 with the "Cassini" telescope.*

![Perimeter_track](/Pictures/Galileo_41550_track.jpg)

*The perimeter of the NORAD 41550 track, from the previous image, as identified by ASTRiDE.*

### Algorithm

**1-** Read the settings parameters and verify that the fits files to be analyzed exist, both the SST images and the calibration biases. The images can have any name followed by a progressive integer. If some number of the file sequence is missing, the software goes to the next one. If all SST files to be processed are missing or all the bias images are missing, the software shows a warning message and stops running. Note that the presence of a settings files makes it possible to use the compiled version without the need to have MATLAB installed in the system. 

**2-** Bias files are renamed with the addition of the suffix "BIAS", to distinguish them from real SST images, preventing them from being analyzed by Astrometry.net for astrometric calibration. If the images with the suffix "BIAS" are already present (which happens if you run the software several times with the same data), it normally counts them directly as bias and continues the execution.

**3-** Creation of a median master bias starting from the single bias frames.

**4-** Calibration of the SST images with subtraction of the master bias frame.

**5-** Reading of the header fits keys to have approximate RA and DEC of the image center in order to speed up the plate solve operations. The header keys that are important for astrometry (date, time, NORAD satellite number, RA and DEC of image center), are saved in the **Data_keys.txt** file. Date, time and NORAD number contained in this file - updated for WCS images, see next step - will be used later for TDM files.

**6-** Plate solve of the calibrated SST images (i.e., master bias subtracted), with Astrometry.net installed locally. The local installation and the approximate knowledge of RA and DEC of the image center make the astrometric calibration very fast. The astrometric calibration takes place using the 2MASS infrared star catalog (Two Micron All Sky Survey). This catalog contains J2000 position data over 410 million stars, with uncertainties ranging from 0.15 to 0.30 arcsec. Astrometry.net package saves fully astrometrically calibrated fits images by adding the suffix WCS to the name. These images, calibrated by Astronomy.net, are used by ASTRiDE for the extraction of satellite tracks. 

**7-** The main script verifies that WCS-type files exist in the working folder. If no WCS files have been created by Astrometry.net BASP execution stop.

**8-** Reading header keys of WCS images and saving in the **Data_keys.txt** file, overwritting the data previously contained. In this way in the **Data_keys.txt** file there is only the data of the WCS images that ASTRiDE will process. In fact we are not sure that Astrometry.net is able to calibrate all the images, e.g. some may not contain enough comparison stars for plate solve because of cloudy skies.

**9-** Extraction with ASTRiDE of satellite tracks from WCS images. This software saves, in a folder with the same name of the image, some png files containing the extracted traces and a final image with all the identified traces in it. In this way we can visually check if the track recognition process was successful. Besides the png images it also saves the text files: **Streaks.txt**. In this file there are all the numerical data about the extracted tracks, the most interesting values are x_center and y_center (in pixels) and the corresponding RA and DEC values (in degrees). The x_center and y_center values are computed by averaging the coordinates of the detected perimeter of the trace. This is a good approximation of the coordinates of the center of the trace but the presence of background stars incorporated in the contour - which generates an edge asymmetry - could alter them. So we have written a Python library, **Track_best_fit.py**, for computing the best fit center of a satellite track. In this library, the fastest and most effective default method is to delete the outliers points of the upper and lower edges of the track placed in the canonical system (i.e. major axis parallel to the x-axis), then the median of the 4 independent sides of the new track is computed to smooth the perimeter and finally the center of the track is computed as the mean of the coordinates of the median segments. The traces with the contours and the best fit median can be displayed provided that the image plot is enabled in the function **track_center2** inside **Track_best_fit.py**. ASTRiDE in the **Streaks.txt** file also measures the area of the track it has detected and it is possible to put a filter on this value, for example by telling ASTRiDE not to consider tracks below a certain size. This cutoff is very useful to avoid being considered as satellites tracks background galaxies or pairs of stars joined together by the seeing that can mimic a track. The Python script that manage ASTRiDE is **SST_Astride_TDM.py** and saves the RA and DEC coordinates of the track center in the **Data_streaks.txt** file.

![Median_track_filter](/Pictures/Median_track.jpg)
*An example of a median edge filter applied to a track with the perimeter altered by background field stars.*

**10-** After the astrometry of the satellites tracks, the control returns to main script which read the files **Data_keys.txt** and **Data_streaks.txt**, combines them to obtain the observations matrix (OM), having the format **[YYYY MM DD hh mm ss.sss AR DEC NORAD]**. Here YYYY, MM, DD, hh, mm, ss.sss, are the column vectors for year, month, day, hours, minutes, seconds and decimals of UT of the mean images time (starting time + esposure time/2), while AR and DEC are the vectors with the coordinates J2000 of the tracks centers found by ASTRiDE in degree. Finally, NORAD is the vector with the satellite's NORAD numbers in ascending order. If OM contains lines where the coordinates AR and DEC are NaNs (because ASTRiDE has not found any trace), these are automatically deleted. In OM time has the precision of 1/1000 of a second, while RA and DEC reach 1/10,000 of a degree (about 0.36 arcsec).

**11-** Sorting of the OM matrix according to the increasing NORAD number of satellites: in this way all the observations of a given satellite are contiguous and can be easily separated into the various TDMs. For a given NORAD number, the temporal ordering is increasing. At this point the computation flow is different when dealing with images having single track or multiple tracks.

**12-** In the case of a single track the astrometric observations are used to determine a best fit orbit of each satellite using the non-interactive version of Find_Orb. In this way it is possible to automatically compute a geocentric orbit and eliminate the observations that - for various reasons that we will see later - have the largest residuals.

**13-** In the case of images with multiple tracks (typically in the case of GEO satellites), BASP selects the track of the desired satellite by automatically downloading the most updated TLE and computing the ephemeris for the observing session using SGP4 model. The SGP4 model functions in SST_SGP4 folder is from 
[Mahooti 2021](https://www.mathworks.com/matlabcentral/fileexchange/62013-sgp4), with some adjustments to make the output compatible with our pipeline. The trace that is less distant from the ephemeris is chosen. After this selection step, as before the astrometric observations are automatically analyzed with the non-interactive version of Find_Orb for deleting the observations with the largest residuals.

**14-** Saving the astrometry of all satellites in the **MPC_YYYYMMDD.txt** file, with astrometry in the Minor Planet Center extended format. This file is for control purposes only because it is possible to get the geocentric orbital elements.

**15-** Saving astrometry in TDM format for each satellite observed in the session. 

### Setting

There is no installation procedure, just put the software in any folder within the Home and run BASP.m or ./BASP from console. The images with satellite/space debris tracks to be analyzed can be placed in any folder on the hard disk. **Warning: make sure the paths in the configuration files are correct**.

### Make a test

In the **Sample_Images** folder there are 19 images that can be used as a test. Image numbers 101 to 110 are calibration biases, 113 to 121 contain satellite images. These images are already selected in the file **Settings_BASP.txt**, just change the paths to make them conform to your system and you can run BASP with the console command ./BASP (compiled version), or by running the BASP.m script with Matlab. During the processing, messages appear informing about the execution status. Note that in the header fits of the images, in addition to the date of the beginning of the exposure with a precision of the order of a few milliseconds and the exposure time, there is the NORAD number of the object and the approximate equatorial coordinates of the center of the image. These keys are required for images analysis and must **always be present**.
