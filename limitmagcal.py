import numpy as np

sigma=[8.1, 9.0 , 9.5 ]
seeing=[4.91, 4.35, 4.11]
#seeing=[6.25,6.25,6.25]
zp=[26.11, 26.21, 26.32]
for i in range(len(sigma)) :
	limitmag = zp[i] - 2.5*np.log10( 5 * sigma[i] * np.sqrt( seeing[i] *2 *seeing[i] *2 *np.pi) )
	print 'limit mag = ',limitmag



