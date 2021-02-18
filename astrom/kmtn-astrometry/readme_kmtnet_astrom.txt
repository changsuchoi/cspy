# KMTNET astrometry brief method - 2015/09/17 Changsu Choi
# using the solution made by Dr. Yujin Yang and Dr. Chunguk Lee
# neccessary files : mkfour, kmtnet.sex, kmtnet.param, kmtnet.scamp, kmtnet.swarp, kmtn_makemef.py, kmtn_resetcrval.py, kmtnet_global_ctio.ahead, weight.fits
#


# 1. making list file as 'o.list'
ls kmtc*.fits > o.list

# 2. run perl script for iraf cl sript 'mkfour.cl'
perl mkfour

# 3. open IRAF window and run the script you made
cl < mkfour.cl # iraf

# 4. making MEF file 
python kmtn_makemef.py filename

# 5. give new CRVAL to newfile
python kmtn_resetcrval.py newfile -c ra, dec (deg unit) 

# 6. run sextractor and scamp to solve astrometry : change projection to TPV,

sex 030216.fits -c  kmtnet.sex \
     -CATALOG_NAME    030216.cat \
     -HEADER_SUFFIX   NONE \
     -DETECT_THRESH   50.0 \
     -ANALYSIS_THRESH 50.0 \
     -SATUR_LEVEL     60000.0     \
     -WEIGHT_TYPE     MAP_WEIGHT  \
     -WEIGHT_IMAGE    weight.fits

# scamp
scamp 030216.cat -c  kmtnet.scamp   \
        -ASTREF_CATALOG    2MASS    \
        -POSITION_MAXERR   20.0       \
        -CROSSID_RADIUS    5.0        \
        -DISTORT_DEGREES   3          \
        -PROJECTION_TYPE   TPV        \
        -AHEADER_GLOBAL    kmtnet_global_ctio.ahead \
        -CHECKPLOT_TYPE    ASTR_REFSYSMAP,FGROUPS,DISTORTION,ASTR_REFERROR2D,ASTR_REFERROR1D \
        -CHECKPLOT_NAME    astr_refsysmap,fgroups,distort,astr_referror2d,astr_referror1d \
        -STABILITY_TYPE    INSTRUMENT \

# 7. run swarp to make one image of 4 chips in it
     but in my try, it was hard to do with single swarp run
     so I change the method, 

# 7.1 give each chip WCS info from filename.head (kk,mm,nn,tt order from top)
# 7.2 swarp to coadd 4 images to make one whole image


