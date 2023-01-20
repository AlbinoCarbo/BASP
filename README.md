# BASP
The main purposes of BASP (Bologna Astrometry Satellites Pipeline), is the automatic astrometric calibration of the images and the extraction of the satellite's traces to obtain a text file with the J2000 astrometric positions of the center of the track in Tracking Data Message. The software was developed under Linux Xubuntu, but can run under any Linux system as long as the following software and their dependencies are installed and running:

\begin{itemize}
\item MATLAB\textsuperscript{\textregistered} R2019b
\item Python 3\footnote{https://www.python.org/download/releases/3.0/} with Astropy package \citep{astropy2018}
\item Astrometry.net (Plate solve of fits images with 2MASS catalog\footnote{http://astrometry.net/})
\item ASTRiDE (Automated Streak Detection for Astronomical Images\footnote{https://github.com/dwkim78/ASTRiDE})
\item Find\_Orb (Geocentric orbit determination from observations\footnote{https://www.projectpluto.com/find\_orb.htm})
\end{itemize}
