
# Gaia errors estimate tool 

## Orignal authors:

 M. Romero-Gomez, F. Figueras, T. Antoja (Universitat de Barcelona, ICCUB-IEEC)

## Modified by:

 D. Kawata (MSSL, UCL)

## DR2 calibration

Parallax and proper motion errors are adjusted to https://www.cosmos.esa.int/web/gaia/dr2.  Other errors (including RVS) are not adjusted. 

CAfactor=1
month=22
jflag=1
parallax error increase by sqrt(month/60).
proper motion error increase by (month/60)$^1.5$.
For G<15 mag stars, parallax errors set 0.04 mas. 

Note that for V>15 mag stars, errors are generally than the DR2 website values.  


