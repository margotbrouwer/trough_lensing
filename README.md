# trough_lensing
Weak lensing of troughs with KiDS and GAMA


During the trough selection, we define circles (with chosen angular size \theta) on our galaxy field.
For each circle, we find the projected density of galaxies within. Circles that are heavily overlapping are flagged.
Underdense circles are called "troughs", while overdense once are defined as "ridges".

The following list shows all the steps that comprise the trough selection. Each item is also mentioned in the code itself, so you can follow the code step by step.

1) Defining the galaxy sample:
  1a) We import the galaxy information from the KiDS/GAMA/mock catalogue.
  1b) We select the sample of galaxies that will be counted to define the troughs.
      This often entails a redshift and/or magnitude cut, because we wish to select troughs that are in front of the "background galaxies" used for lensing.

2) Creating the grid:
  2a) We fill our galaxy field with a cartesian grid of narrowly spaced (1 arcmin) points.
  2b) We remove all points that lie within 1 arcmin of masked areas or the edge of the field.

3) Measuring galaxy density:
  3a) For each grid point, we count the number of sample galaxies within a circle of chosen radius \theta.
      This is done for multiple radii: \theta = [5, 10, 15, 20] arcmin, which allows us to select and compare troughs of different sizes.
  3b) For each grid point, we count the number of other grid points within the same circle of radius \theta.
      Because grid points trace the galaxy field, this allows us to measure the effective area of the circle (taking into account the masks and the field edges).
  3c) For each grid point, the galaxy density is defined as the galaxy count within the circle, divided by the effective area of the circle.

4) Flagging overlapping circles
  4a) We sort the circles by galaxy density, and divide the circles in two samples: underdense and overdense.
      Underdense cirlces have a lower density than the median (50th percentile), overdense circles have a higher density than the median.
  4b) Definition: Two circles are "overlapping" when their centers (= their corresponding grid points) are closer to eachother than 0.5 \theta.
  4c) For the underdense circles: We start from the lowest density circle, and flag all circles that are overlapping with it.
      We then go to the next lowest density circle that is still unflagged, and repeat this process for all underdense cirles.
  4d) For the overdense circles: We start from the highest density circle, and flag all circles that are overlapping with it.
      We then go to the next highest density circle that is still unflagged, and repeat this process for all overdense cirles.

5) Selecting the trough sample:
  5a) For each grid point/circle, calculate the percentage Ptheta of grid points that has a lower galaxy density.
      Following the definition from DES, we usually define the 20% of all circles with the lowest density to be the troughs.
      The 20% of all circles with the highest density are defined as overdensities ("ridges").
      (In the future, it might be possible to optimize this percentage for the S/N of the trough lensing signal.)
  5b) Now we have two trough samples:
    - The overlapping sample, by taking all troughs
    - The non-overlapping sample, by taking only unflagged troughs (unflagged = 1, flagged = 0).
    These two trough samples can be used to create the trough lensing signals.

For each grid point/circle, we save the following information to the catalog:
the location (RA/DEC), number count (Nstheta), effective area in arcmin (grid count Ngtheta), galaxy density (rhotheta), 
percentage of circles below its density (Ptheta), and non-overlapping selection flag (Stheta).
    
Done!
