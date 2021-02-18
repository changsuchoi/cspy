#! /usr/bin/env python

import os, sys, string

try: 
    import pyfits
except ImportError:
    print "Usage: kmtn_resetcrval.py requires pyfits module!"
    exit(-1)

#----------------------------------------------------------------------
def resetcrval_batch(files, radec) :
   import glob

   # determine image file name
   file_lst = read_list(files)

   for file in file_lst:
      resetcrval(file, radec)

#----------------------------------------------------------------------
def resetcrval(file, radec) :
   if os.path.isfile(file) == False :
      print "File %s not exists!" % file
      print "Skipping...."
      return

   # Open fits file
   hdulist = pyfits.open(file)
   hdulist = pyfits.open(file, mode='update')
   # DEBUG hdulist.info()

   # Update header
   print 'Updating header...' 
   for iext in [1,2,3,4]:
      print file, iext, radec[0], radec[1]
      hdulist[iext].header.set('CRVAL1', radec[0])
      hdulist[iext].header.set('CRVAL2', radec[1])

   hdulist.close()

#----------------------------------------------------------------------
# Function Call
#----------------------------------------------------------------------

def usage() :     
   print "Usage: kmtn_resetcrval.py image_name -c ra,dec"
   print "       kmtn_resetcrval.py image_name --coord=ra,dec"
   print "         reset CRVAL1 & CRVAL2 header"
   print "Used package: pyfits"
   print "Input  : images to correct CRAVL1 & CRVAL2"
   print "Output : the same image with updated CRVAL1 & CRVAL2"
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

   radec=[]
   # Use gnu_getopt for the intermixed options and arguments
   try:
      opts, files = getopt.gnu_getopt(sys.argv[1:], 'c:', ['coord='])
   except getopt.GetoptError:
      usage()

   for opt, arg in opts:
      if opt in ('-c', '--coord'):
         radec = arg
         radec = radec.split(',')
         if len(radec) != 2:
             print '\n\nNeed ra and dec input!\n\n'
             usage()

   if len(files) < 1 : usage()
#  DEBUG: print files, radec

   resetcrval_batch(files, radec)


