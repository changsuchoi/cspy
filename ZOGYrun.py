Suggested steps to get started:

(1) download one of the install scripts, e.g. for mac OS:

wget https://raw.githubusercontent.com/pmvreeswijk/ZOGY/master/InstallScripts/install_zogy_macos.sh

(2) read through the install script to make sure that it won't interfere with your present set-up, e.g. the python version

(3) execute the script in a folder where you'd like the ZOGY folder with python modules and settings/configuration subfolders to be created

(4) edit the settings file located in [some path]/ZOGY/Settings/set_zogy.py to adapt it to your images and their headers. For telescope-dependent parameters, you could add your telescope name to the keys with the corresponding value for that parameter, which will then be used if the telescope input parameter (see below) is set (default: "ML1" for MeerLICHT). If a parameter is not a dictionary but a single value, that value will be used regardless of the telescope input parameter.

(5) check out the main input parameters:

python [some path]/ZOGY/zogy.py -h

(6) some examples how to run it:

- run it on a "new.fits" and "ref.fits" using the default MeerLICHT settings:

python [some path]/ZOGY/zogy.py --new_fits new.fits --ref_fits ref.fits


- instead of MeerLICHT, use the dictionary keys corresponding to "my_tel" defined in the settings file:

python [some path]/ZOGY/zogy.py --new_fits new.fits --ref_fits ref.fits --telescope my_tel


- instead of the default settings file (set_zogy.py), use a copy of it that was adapted to your images
  (depending on whether copy contains dictionaries for some parameters or not, the telescope input
   parameter should be provided or not):

python [some path]/ZOGY/zogy.py --new_fits new.fits --ref_fits ref.fits --set_file mycopy [--telescope my_tel]


- if you have mask images available (settings parameters [transient_mask_max] and [mask_value] should be
  updated if your masks contain different values for the mask pixel type; the keys cannot be changed):

python [some path]/ZOGY/zogy.py --new_fits new.fits --new_fits_mask new_mask.fits --ref_fits ref.fits --ref_fits_mask ref_mask.fits --set_file mycopy [--telescope my_tel]
