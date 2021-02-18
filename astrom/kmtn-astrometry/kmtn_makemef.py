##! /usr/bin/env python

import os, sys, string

try: 
    import pyfits
except ImportError:
    print "Usage: kmtn_makemef.py requires pyfits module!"
    exit(-1)

#----------------------------------------------------------------------
def makemef_batch(files, suffix) :
   import glob

   # determine image root name
   file_lst = read_list(files)

   for file in file_lst:
      root = file.replace('.fits','')
      makemef(root, suffix)

#----------------------------------------------------------------------
def makemef(root, suffix) :
   # Default KMTNet suffix
   if suffix == []:
       suffix = ['.kk','.mm','.nn','.tt']
   # check how many extensions are:
   import glob
   files = [root+sfx+'.fits' for sfx in suffix]

   # Check if all extensions exist:
   all_exist = True
   for file_ext in files:
      if os.path.isfile(file_ext) == False:
         print "File %s not exists!" % file_ext
         all_exist = False
   if all_exist == False:
      print "Some extensions are missing. Exiting..."
      return

   extvec = range(1,len(files)+1)

   # create empty PrimaryHDU
   primary_hdu = pyfits.PrimaryHDU()

   # create new HDUList to fill
   hdulist = pyfits.HDUList()
   
   # fill in with extensions
   hdulist.append(primary_hdu)
   # hdulist.info()

   file_mef = root + '.fits'

   if os.path.isfile(file_mef) == True :
      print "File %s already exists!" % file_mef
      print "Skipping...."
      return

   #---------------------------------------------------------------------- 
   #  FIRST METHOD: directly append/write image each time 
   #                [slow, but less memory intensive]
   #  hdulist.writeto(file_mef)
   #  for (iext, file) in zip(extvec, files):
   #     print iext, file
   #
   #     # open single fits file
   #     hdulist_ext = pyfits.open(file)
   #     hdu = hdulist_ext[0]
   #
   #     # append hdu to file and write to disk
   #     pyfits.append('temp.fits', data=hdu.data)
   #  #  pyfits.append('temp.fits', data=hdu.data, header=hdu.header)
   #---------------------------------------------------------------------- 

   # open single fits file
   hdulist_ext = pyfits.open(files[0])
   hdu = hdulist_ext[0]

   # SECOND METHOD: create mef and write at the end
   print 'Combining extensions: ' 
   for (iext, file) in zip(extvec, files):
      print iext, file

      # open single fits file
      hdulist_ext = pyfits.open(file)
      hdu = hdulist_ext[0]

      hdu.header.set('IMAGEID',  iext)

      # create new ImageHDU
      s = pyfits.ImageHDU(data=hdu.data, header=hdu.header, name=("im%d" % iext))
      s.header.set('EXTNAME', ("im%d" % iext))
   
      # append HDU to HDUList 
      hdulist.append(s)

#  hdulist.info()
   print 'Writing ' + file_mef
   hdulist.writeto(file_mef)

   hdulist_ext.close()
   hdulist.close()

#----------------------------------------------------------------------
# Function Call
#----------------------------------------------------------------------

def usage() :     
   print "Usage: kmtn_makemef.py image_root -s "
   print "       kmtn_makemef.py image_root --suffix=suffix"
   print "       kmtn_makemef.py image --suffix=_1,_2,_3,_4"
   print "         creating image.fits from image_1.fits,image_2.fits,image_3.fits,image_4.fits"
   print "Used package: pyfits"
   print "Input  : images to combine"
   print "Output : MEF file"
   sys.exit(2)

# Function to read the filelist
def read_list(files) :
   file_list = []           # output file list
   for file in files:
      # if argument is @filelist, read the list and extend
      if file[0] == '@' :
         file = file[1:]
         try:
            f = open(file, 'r')
         except IOError:
            print "File %s does not exists!" % file
            print "Exiting...."
            sys.exit(-1)

         list = []
         for line in f.readlines() :
            if line[0] != '#' and line[0] != '\n' :
               # append after remove all spaces
               line = line.replace('\n','')
               line = line.replace('\t','')
               line = line.replace(' ', '')
               list.append(line)
         f.close
         file_list.extend(list)
      # if a single file, just append 
      else:
          file_list.append(file)
   return file_list


if __name__ == "__main__":
   import os, getopt

   suffix=[]
   # Use gnu_getopt for the intermixed options and arguments
   try:
      opts, files = getopt.gnu_getopt(sys.argv[1:], 's:', ['suffix='])
   except getopt.GetoptError:
      usage()

   for opt, arg in opts:
      if opt in ('-s', '--suffix'):
         suffix = arg
         suffix = suffix.split(',')

   if len(files) < 1 : usage()
   # DEBUG: print files, suffix

   makemef_batch(files, suffix)


