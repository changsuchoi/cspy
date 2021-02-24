'''
This is adapted from http://obswww.unige.ch/%7Etewes/alipy/tutorial.html
'''

import alipy
import glob
import os
from multiprocessing import Process,Pool

# images_to_align = sorted(glob.glob("images/*.fits"))
# images_to_align = sorted(glob.glob('Calib-*.fits'))
ref_image = "ref.fits"


def identify_transform(ref_image, images_to_align,rad= 5, nb=500,verbose=False, visual=False):
    identifications = alipy.ident.run(ref_image, images_to_align,r=rad, n=nb, visu=visual)
    # That's it !
    # Put visu=True to get visualizations in form of png files
    #   (nice but much slower)
    # On multi-extension data, you will want to specify the hdu (see API doc).

    if verbose:
        # The output is a list of Identification objects,
        # which contain the transforms :
        for id in identifications:  # list of the same length as images_to_align.
            if id.ok is True:  # i.e., if it worked
                print("%20s : %20s, flux ratio %.2f" % (id.ukn.name,
                                                        id.trans,
                                                        id.medfluxratio))
                # id.trans is a alipy.star.SimpleTransform object.
                #    Instead of printing it out as a string,
                # you can directly access its parameters :
                # print id.trans.v # the raw data,
                #                  # [r*cos(theta)  r*sin(theta)
                #                  #  r*shift_x  r*shift_y]
                # print id.trans.matrixform()
                # print id.trans.inverse() # this returns a
                #                          # new SimpleTransform object

            else:
                print( "%20s : no transformation found !" % (id.ukn.name))

    return identifications


# identifications= identify_transform(ref_image, images_to_align, rad= 5, nb=500, verbose=False, visual=False)


def align_images(ref_image, identifications, iraf=True, outdir='alipy_out'):
    # Minimal example of how to align images :
    outputshape = alipy.align.shape(ref_image)
    # This is simply a tuple (width, height)...
    # you could specify any other shape.
    for id in identifications:
            if id.ok is True:
                if iraf:
                    # Variant 2, using geomap/gregister,
                    # correcting also for distortions :
                    alipy.align.irafalign(id.ukn.filepath, id.uknmatchstars,
                                          id.refmatchstars, shape=outputshape,
                                          makepng=False, outdir=outdir)
                    # id.uknmatchstars and id.refmatchstars are simply lists of
                    # corresponding Star objects.
                else:
                    # Variant 1, using only scipy and the simple affine transorm :
                    alipy.align.affineremap(id.ukn.filepath, id.trans,
                                            shape=outputshape, makepng=True,
                                            outdir=outdir)
                # By default, the aligned images are written into a directory
                # "alipy_out".






# match input images to ref image
# input image will be aligned to ref image, no change in ref image
def alipy_in2ref(ref_image,images_to_align):
    '''
    images_to_align should be a list!
    '''
    identifications= identify_transform(ref_image, images_to_align, rad= 5, nb=500, verbose=False, visual=False)
    align_images(ref_image, identifications, iraf=True, outdir='alipy_out')
    for i in images_to_align :
        newname = 'gr_'+os.path.splitext(i)[0]+os.path.splitext(i)[1]
        os.system('mv alipy_out/'+os.path.splitext(i)[0]+'_gregister'+os.path.splitext(i)[1]+' '+newname)




# match ref image to input image
# ref image will be aligned to input image, no change in input image
# im='Calib-LSGT-NGC3367-20190630-084141-r-180.fits'
# ref_image=['ref.fits']
# output name = gr_Calib-LSGT-NGC3367-20190630-084141-r-180.fits
def alipy_ref2in(im,ref_image=['ref.fits']):
    '''
    ref_image should be a list!
    like ref_image=['ref.fits']
    '''
    identifications= identify_transform(im, ref_image, rad= 5, nb=500, verbose=False, visual=False)
    align_images(im, identifications, iraf=True, outdir='alipy_out')
    newname = 'gr_'+os.path.splitext(im)[0]+os.path.splitext(im)[1]
    os.system('mv alipy_out/'+os.path.splitext(ref_image[0])[0]+'_gregister'+os.path.splitext(ref_image[0])[1]+' '+newname)




# with open('comfiles.txt','r') as file_handle: lines = file_handle.read().splitlines()
def alipy_epoch(lines):
    for iii in lines:
        ii=iii[:-1].split(',')
        print('='*60,'\n')
        print(lines.index(iii), 'of', len(iii))
        if len(ii) == 1 :
            print(ii[0], 'single image, No touch')
            #os.system('mv '+ii[0]+' gr'+ii[0])
            #os.system('mv '+os.path.splitext(ii[0])+'_gregister.fits')
        else :
            ref_image=ii[0]
            images_to_align=ii
            print (len(ii), 'files will be registerd')
            identifications= identify_transform(ref_image, images_to_align, rad= 5, nb=500, verbose=False, visual=False)
            align_images(ref_image, identifications, iraf=True, outdir='alipy_out')
    os.system('mv alipy_out/*gregister.fits .')


#os.system('rename _gregister.fits .fits alipy_out/*gregister*')
#os.system('rename Calib grCalib alipy_out/Calib*.fits')
#os.system('mv alipy_out/grCalib*.fits .')






'''

with open('comfiles.txt','r') as file_handle: lines = file_handle.read().splitlines()

# gregister
for i in lines :
        if len(i.split(',')) <= 2 :
            print(i, 'single image')
        else :
            ref_image=i.split(',')[0]
            images_to_align=i.split(',')[:-1]
            print (i.split(','))
            identifications= identify_transform(ref_image, images_to_align, rad= 5, nb=500, verbose=False, visual=False)
            align_images(ref_image, identifications, iraf=True, outdir='alipy_out')


'''
