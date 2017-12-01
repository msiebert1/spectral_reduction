kast_blue = {'name': 'kast_blue',
			 'read_noise': 2.1,#temp
			 'gain': 5.4,#temp
			 'grism': 'temp',
			 'filter': 'temp',
			 'slit': 'temp',
			 'dispaxis': 1,
			 'biassec': '[1:2048,300:350]',
			 'trimsec': '[1:2048,1:300]',
			 'QL_arc': 'b1008_ref.fits',
			 'QL_sens': 'sens.blue.fits',
			 'section': 'middle line'}


kast_red = { 'name': 'kast_red',
             'read_noise': 2.1,#temp
			 'gain': 5.4,#temp
			 'grism': 'temp',
			 'filter': 'temp',
			 'slit': 'temp',
			 'dispaxis': 2, 
			 'biassec': '[355:405,1:2725]',
			 'trimsec': '[45:355,1:2725]',
			 'QL_arc': 'r1007_ref.fits',
			 'QL_sens': 'sens.red.fits', #TODO: Change to yen-chen's file
			 'section': 'middle column'}