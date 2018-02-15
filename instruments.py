path_to_ref = '/home/msiebert/Documents/UCSC/Research/spectral_reduction/QL_ref/'
kast_blue = {'name': 'kast_blue',
			 'read_noise': 3.7,
			 'gain': 1.2,
			 'grism': 'temp',
			 'filter': 'temp',
			 'slit': 'temp',
			 'dispaxis': 1,
			 'biassec': '[1:2048,300:350]',
			 'trimsec': '[1:2048,1:300]',
			 # 'trimsec': '[355:1605,1:300]',
			 'QL_arc': path_to_ref + 'blue_ref.fits',
			 'QL_sens': path_to_ref + 'sens.blue.fits',
			 'QL_arc_ex': path_to_ref + 'tb1003.ms.fits',
			 'section': 'middle line'}


kast_red = { 'name': 'kast_red',
             'read_noise': 3.8,
			 'gain': 1.9,
			 'grism': 'temp',
			 'filter': 'temp',
			 'slit': 'temp',
			 'dispaxis': 2, 
			 'biassec': '[355:405,1:2725]',
			 # 'trimsec': '[45:355,1:2725]',
			 'trimsec': '[45:355,250:2200]',
			 'QL_arc': path_to_ref + 'red_ref.fits',
			 'QL_sens': path_to_ref + 'sens.red.fits', #TODO: Change to yen-chen's file
			 # 'QL_arc_ex': path_to_ref + 'red_ex_ref.ms',
			 'section': 'middle column'}