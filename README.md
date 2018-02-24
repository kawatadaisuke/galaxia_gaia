
# galaxia_gaia

 Adding Gaia errors to Galaxia output.

## References:

Galaxia:

Sharma et al. et al. (2011) http://adsabs.harvard.edu/abs/2011ApJ...730....3S

Gaia errors:

Romero-Gomez et al. (2015) http://adsabs.harvard.edu/abs/2015MNRAS.447..218R

## HOWTO

Generate input files for Gaia Errors estimate code from galaxia output file (default: galaxia/galaxy1.ebf). 
> cd py
> python galaxia2ubgaaiae.py
This generates ubgaiaein.bin.

Estimate Gaia errors.
> cd ubgaiaerrors
> make
> gaia_errors
This reads ../ubgaiaein.bin, and generates ../ubgaiaed-out.bin and ../ubgaiaei-out.bin.

Create a fits file
> cd ..
> python ubgaia2fits.py
This generates galaxia_gaia.fits.

Plot distance errors vs. true distance.
> python plotde_d.py
Plot proper motion (km/s) errors vs. true distance.
> python plotvpe_d.py
Plot radial velocity (km/s) errors vs. true distance.
> python plotvre_d.py
