import numpy as np 
import matplotlib.pyplot as plt 
from scipy import interpolate
import scipy.optimize as opt

def optimize_scales(blue, red):
    """Scales each unique supernova in SN_Array by minimizing the square residuals
       between the supernova flux and the template flux. This also works for bootstrap
       arrays (can contain repeated spectra) because the objects in SN_Array are not 
       copies. Returns scaled SN_Array and the scales that were used.
    """
    guess = np.average(red[1])/np.average(blue[1])
    scale = opt.minimize(sq_residuals, guess, args = (blue, red), 
                         method = 'Nelder-Mead').x
    return scale
    
    
def sq_residuals(s,blue,red):
    """Calculates the sum of the square residuals between two supernova flux 
       arrays usinig a given scale s. Returns the sum.
    """
#    print s
    s = np.absolute(s)
    new_blue = s*blue[1]
    res = red[1] - new_blue
    sq_res = res*res
    blue_ivar = 1./(blue[2]*blue[2])
    temp_ivar = blue_ivar/(s*s)
    w_res = blue_ivar*sq_res
    return np.sum(w_res)

def combine_blue_red(blue_asci, red_asci, name):
	blue = np.transpose(np.genfromtxt(name + '/' + blue_asci))
	red = np.transpose(np.genfromtxt(name + '/' + red_asci))

	#cut noisy part of red spectrum
	valid_range = np.where(red[0]>5480.)
	red = [red[0][valid_range], red[1][valid_range], red[2][valid_range]]

	#Create red and blue spectra in overlap region
	red_olap = np.where(red[0]<blue[0][-1])
	red_olap_spec = [red[0][red_olap], red[1][red_olap], red[2][red_olap]]
	blue_olap = np.where(blue[0]>red[0][0])
	blue_olap_spec = [blue[0][blue_olap], blue[1][blue_olap], blue[2][blue_olap]]

	#interpolate blue flux and error spectra in overlap region
	blue_spl = interpolate.splrep(blue_olap_spec[0], blue_olap_spec[1])
	blue_err_spl = interpolate.splrep(blue_olap_spec[0], blue_olap_spec[2])
	new_blue_flux = interpolate.splev(red_olap_spec[0], blue_spl)
	new_blue_err = interpolate.splev(red_olap_spec[0], blue_err_spl)
	blue_resampled = [red_olap_spec[0], new_blue_flux, new_blue_err]

	#scale blue in overlap region
	scale = optimize_scales(blue_resampled,red_olap_spec)[0]
	blue_olap_rescaled = [blue_resampled[0], scale*blue_resampled[1], scale*blue_resampled[2]]
	print scale

	#inverse variance weighted average in overlap region
	blue_olap_ivar = 1./blue_olap_rescaled[2]
	red_olap_ivar = 1./red_olap_spec[2]
	olap_flux = np.average([blue_olap_rescaled[1], red_olap_spec[1]], weights = [blue_olap_ivar, red_olap_ivar], axis = 0)
	olap_err = np.average([blue_olap_rescaled[2], red_olap_spec[2]], weights = [blue_olap_ivar, red_olap_ivar], axis = 0)
	olap_arr = [red_olap_spec[0], olap_flux, olap_err]

	#Exclude overlap region from original spectra
	blue_cut = np.where(blue[0]<=red[0][0])
	red_cut = np.where(red[0]>=blue[0][-1])
	blue_short = [blue[0][blue_cut], blue[1][blue_cut], blue[2][blue_cut]]
	red_short = [red[0][red_cut], red[1][red_cut], red[2][red_cut]]

	#scale shortened blue/red spectra to overlap region
	bscale = olap_arr[1][0]/blue_short[1][-1]
	rscale = olap_arr[1][-1]/red_short[1][0]
	blue_short = [blue_short[0], bscale*blue_short[1], bscale*blue_short[2]]
	red_short = [red_short[0], rscale*red_short[1], rscale*red_short[2]]

	#combine into one spectrum
	wavelength = np.concatenate((blue_short[0], olap_arr[0], red_short[0]))
	flux = np.concatenate((blue_short[1], olap_arr[1], red_short[1]))
	err = np.concatenate((blue_short[2], olap_arr[2], red_short[2]))


	plt.plot(blue[0], blue[1])
	plt.plot(red[0], red[1])
	plt.show()
	
	plt.plot(blue_short[0], blue_short[1])
	plt.plot(olap_arr[0], olap_arr[1])
	plt.plot(red_short[0], red_short[1])
	plt.show()


	return [wavelength, flux, err]

# combine_blue_red('feige34_kast_blue_2017-07-03T04:14:09.91_1_f.asci', 'feige34_kast_red_2017-07-03T04:14:10.43_1_f.asci', 'feige34')