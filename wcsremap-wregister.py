# https://github.com/jonnybazookatone/imclass/blob/master/image.py#L602
# imclass/image.py


	def reMap(self, remapFits, outname="input_remapped.fits", verbose=False, wcsregister=False):
	  
		# Create an object class for the output image for fluency of continuing code
		outFits = imFits()
		outFits._Name = outname
		outFits._logger["ReMapping"] = []
	  
		# remapFits = template image
		# this object = source image
		# outFits = output fits image
		
		if wcsregister:
			iraf.wregister(input=self._Name, output=outFits._Name, wcs="world", reference=remapFits._Name, Stdout=1)

		else:
			wcscmd = ["wcsremap","-template", self._Name, "-source", remapFits._Name, "-outIm", outFits._Name]
			wcsremap = subprocess.Popen(wcscmd, stdout=subprocess.PIPE)
			outFits._logger["ReMapping"].append(wcsremap.communicate()[0])
		
		if verbose:
			for log in outFits._logger["ReMapping"]:
				print log
		
		return outFits