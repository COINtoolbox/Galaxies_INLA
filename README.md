#  <img  align="right" src="https://github.com/COINtoolbox/DRACULA/blob/master/images/coin.png" width="400"> [Spatial field reconstruction with INLA: Application to IFU galaxy data](https://arxiv.org/abs/1802.06280)


[S. Gonzalez-Gaitan](https://github.com/gongsale), [R. S. de Souza](https://github.com/RafaelSdeSouza), [A. Krone-Martins](https://github.com/algolkm), E. Cameron, P. Coelho, L. Galbany, E.E.O. Ishida


This is one of the products of [COIN Residence Program #4](http://iaacoin.wix.com/crp2017), which took place in August/2017 in Clermont-Ferrand (France). 

We present here: 

- A tutorial to run [R-INLA](http://www.r-inla.org/) on galaxy maps, 
-  A full catalog of INLA reconstruted maps of mass, metallicity, age and extinction obtained with [STARLIGHT](http://www.starlight.ufsc.br/) from the [CALIFA](http://califa.caha.es/) and [PISCO](http://adsabs.harvard.edu/abs/2018arXiv180201589G) surveys. 

### *Code snipppet to reproduce the spatial maps* 

We provide a simple  [example](https://github.com/COINtoolbox/Galaxies_INLA/blob/master/Run_INLA.R) to read the data, perform the INLA analysis and display the results. A simple way to *display* the results in python is shown in a small [notebook](https://github.com/COINtoolbox/Galaxies_INLA/blob/master/read_cube.ipynb).


### CALIFA/PISCO Galactic Property Catalogue 

The current release consists of 667 nearby galaxies from CALIFA and 104 galaxies from PISCO fitted pixel-by-pixel with STARLIGHT obtaining luminosity-weighted age and metallicity, mass-weighted age and metallicity, stellar mass and extinction, each with its respective mean and dispersion from the INLA reconstruction.

Data structure:

A) A single table (divided in 15 [files](https://github.com/COINtoolbox/Galaxies_INLA/blob/master/data/allgalaxies1.dat.gz)) composed of 21 columns:

1. GALAXY NAME
2. x position in pixels
3. y position in pixels
4. STARLIGHT log(mass) in solar masses
5. STARLIGHT luminosity-weighted log(age) in yr
6. STARLIGHT luminosity-weighted log(metallicity)
7. STARLIGHT mass-weighted log(age) in yr
8. STARLIGHT mass-weighted log(metallicity)
9. STARLIGHT visual extinction Av
10. INLA mean of log(mass) in solar masses
11. INLA std deviation of log(mass) in solar masses
12. INLA mean of luminosity-weighted log(age) in yr
13. INLA std deviation of luminosity-weighted log(age) in yr
14. INLA mean of luminosity-weighted log(metallicity)
15. INLA std deviation of luminosity-weighted log(metallicity)
16. INLA mean of mass-weighted log(age) in yr
17. INLA std deviation of mass-weighted log(age) in yr
18. INLA mean of mass-weighted log(metallicity)
19. INLA std deviation of mass-weighted log(metallicity)
20. INLA mean of visual extinction Av
21. INLA std deviation of visual extinction Av

B) One invidual FITS files per galaxy, e.g. [NGC-309](https://github.com/COINtoolbox/Galaxies_INLA/blob/master/data/fits/NGC0309_starlight.fits), for all 6 STARLIGHT properties:
1. STARLIGHT log(mass) in solar masses
2. STARLIGHT luminosity-weighted log(age) in yr
3. STARLIGHT luminosity-weighted log(metallicity)
4. STARLIGHT mass-weighted log(age) in yr
5. STARLIGHT mass-weighted log(metallicity)
6. STARLIGHT visual extinction Av

And another individual FITS file per galaxy, e.g. [NGC-741](https://github.com/COINtoolbox/Galaxies_INLA/blob/master/data/fits/NGC0741_inla.fits), for the respective INLA mean and dispersion:
1. INLA mean of log(mass) in solar masses
2. INLA std deviation of log(mass) in solar masses
3. INLA mean of luminosity-weighted log(age) in yr
4. INLA std deviation of luminosity-weighted log(age) in yr
5. INLA mean of luminosity-weighted log(metallicity)
6. INLA std deviation of luminosity-weighted log(metallicity)
7. INLA mean of mass-weighted log(age) in yr
8. INLA std deviation of mass-weighted log(age) in yr
9. INLA mean of mass-weighted log(metallicity)
10. INLA std deviation of mass-weighted log(metallicity)
11. INLA mean of visual extinction Av
12. INLA std deviation of visual extinction Av

Additionally, 6 png figures are included per galaxy showing input STARLIGHT map and respective INLA mean/standard deviation, e.g. [mass of NGC-309](https://github.com/COINtoolbox/Galaxies_INLA/blob/master/data/plots/NGC0309_mass.png) or [Av of NGC-741]((https://github.com/COINtoolbox/Galaxies_INLA/blob/master/data/plots/NGC0309_Av.png)).

Comments: We note that all galaxies were modeled in the same way and with the same input parameters to R-INLA. For optimum performance however, the user may consider adapting the INLA model to the specific needs of the given problem. We point out for instance that several galaxies have foreground stars for which STARLIGHT erroneously assigns large stellar masses and other properties. These pixels should be removed by the user prior to the INLA reconstruction (or STARLIGHT fit). One such example is KUG 0210-078.
