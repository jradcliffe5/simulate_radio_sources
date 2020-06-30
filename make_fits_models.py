from astropy.io import fits
import numpy as np
from code_mailer import headless
#from make_fits_models import *
from astropy.modeling.functional_models import Gaussian2D
import pandas as pd
from astropy.convolution import CustomKernel
import os
from astropy.convolution import convolve
from astropy.wcs import WCS
from joblib import Parallel, delayed
import multiprocessing

def convertAIPStoPythonImage(filename,outfilename):
	hdu_list = fits.open(filename)

	head = hdu_list['PRIMARY'].header
	image_data = hdu_list[0].data
	hdu = fits.PrimaryHDU(image_data, header=head)
	hdulist = fits.HDUList([hdu])
	hdulist.writeto(outfilename,overwrite=True)
	return outfilename

def setup_source_pixel_grid(fitsheader, npoint, pc_edge_cut, random):
	RA_pix = fitsheader['NAXIS1']
	DEC_pix = fitsheader['NAXIS2']
	RA_pc_cut = (RA_pix)*pc_edge_cut
	RA = np.array([int(0+RA_pc_cut), int((RA_pix-1)-RA_pc_cut)])
	DEC_pc_cut = (DEC_pix)*pc_edge_cut
	DEC = np.array([int(0+DEC_pc_cut), int((DEC_pix-1)-RA_pc_cut)])
	if random == 'True':
		print('Making random grid')
		RA = np.random.randint(RA[0],RA[1]+1, size=npoint**2)
		DEC = np.random.randint(DEC[0],DEC[1]+1,size=npoint**2)
		pointings = np.array([RA,DEC]).T
		for i in range(len(pointings)):
			while len(np.where(np.equal([True,True], np.isclose(pointings[i],pointings,atol=1e-10,rtol=0)).all(axis=1)==True)[0])>1:
				print('replacing value')
				pointings[i] = np.array([np.random.randint(RA[0],RA[1]+1, size=1)[0], np.random.randint(DEC[0],DEC[1]+1,size=1)[0]])

	else:
		print('Making boringly spaced grid')
		RA = np.linspace(RA[0],RA[1],npoint,endpoint=True).astype(int)
		DEC = np.linspace(DEC[0],DEC[1],npoint,endpoint=True).astype(int)
		pointings = []
		for i in RA:
			for j in DEC:
				pointings = pointings + [[i,j]]
	return pointings

def generate_fits_models_delta_fcn(fitsfile, SN, rms, pixel_grid,type_test,save_model):
	'''
	Generates a delta function grid model sky which can be input into a measure
	ment set using uvsub
	'''
	hdu_list = fits.open(fitsfile)
	image_data_implane = np.copy(hdu_list[0].data)
	image_data0 = np.copy(hdu_list[0].data*0.0)
	head = hdu_list['PRIMARY'].header
	wcs = WCS(head)

	try:
		float(rms)
		print('RMS value found')
		rms_float = True
	except ValueError:
		print('Using supplied rms image to calculate S/N flux')
		rms_float = False
		rms_hdu = fits.open(rms)
		rms_im = rms_hdu[0].data.squeeze()
		rms_hdu.close()


	for i in range(len(SN)):
		image_data = np.copy(image_data0)
		SN_flux=np.empty(len(pixel_grid))
		x = []
		y = []
		for j in range(len(pixel_grid)):
			if rms_float == True:
				SN_flux[j] = float(SN[i]*float(rms))
			else:
				SN_flux[j] = float(SN[i]*rms_im[pixel_grid[j][1],pixel_grid[j][0]])
			image_data[0,0,pixel_grid[j][1],pixel_grid[j][0]] = SN_flux[j]
			x = x + [pixel_grid[j][0]]
			y = y + [pixel_grid[j][1]]
		x_world = wcs.all_pix2world(x,y,1,1,1)[0]
		y_world = wcs.all_pix2world(x,y,1,1,1)[1]
		if type_test == 'implane':
			kernel = make_Gaussian_beam_kernel(head,23)
			image_data_conv = convolve(image_data.squeeze(), kernel, normalize_kernel=False)
			image_data = np.add(image_data_conv,image_data_implane.squeeze())
		hdu = fits.PrimaryHDU(image_data, header=head)
		hdulist = fits.HDUList([hdu])
		if type_test == 'uvplane':
			hdulist.writeto(fitsfile.split('.fits')[0]+'_uv_delt_SN%s.fits' % SN[i] ,overwrite=True)
			os.system('rm %s_uv_delt_SN%s_input_model.csv' % (fitsfile.split('.fits')[0],SN[i]))
			pd.DataFrame({'x':x,'y':y,'x_deg':x_world,'y_deg':y_world,'mode_flux':SN_flux}).to_csv('%s_uv_delt_SN%s_input_model.csv' % (fitsfile.split('.fits')[0],SN[i]))
			hdulist.close()
		elif type_test == 'implane':
			hdulist.writeto(fitsfile.split('.fits')[0]+'_im_delt_SN%s.fits' % SN[i] ,overwrite=True)
			hdulist.close()
			os.system('rm %s_im_delt_SN%s_input_model.csv' % (fitsfile.split('.fits')[0],SN[i]))
			pd.DataFrame({'x':x,'y':y,'x_deg':x_world,'y_deg':y_world,'mode_flux':SN_flux}).to_csv('%s_im_delt_SN%s_input_model.csv' % (fitsfile.split('.fits')[0],SN[i]))
		else:
			print('hello')
		if save_model == True:
			hdu = fits.PrimaryHDU(image_data_conv, header=head)
			hdulist = fits.HDUList([hdu])
			hdulist.writeto(fitsfile.split('.fits')[0]+'_im_delt_SN%s.initmodel.fits' % SN[i] ,overwrite=True)
			hdulist.close()

def makeGaussian(size, fwhm, center):
	""" Make a square gaussian kernel.

	size is the length of a side of the square
	fwhm is full-width-half-maximum, which
	can be thought of as an effective radius.
	"""

	x = np.arange(0, size, 1, float)
	y = x[:,np.newaxis]

	if center is None:
		x0 = y0 = size // 2
	else:
		x0 = center[0]
		y0 = center[1]

	return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)

def makeGaussian_bpa(size,amplitude, std, bpa,center):
	""" Make a non-square gaussian kernel.

	size is the length of a side of the square
	fwhm is full-width-half-maximum, which
	can be thought of as an effective radius.
	"""
	bpa = (np.pi/180.)*bpa ## convert to radians
	x = np.arange(0, size, 1, float)
	y = x[:,np.newaxis]

	if center is None:
		x0 = y0 = size // 2
	else:
		x0 = center[0]
		y0 = center[1]

	a = ((np.cos(bpa)**2)/(2*(std[0]**2))) + ((np.sin(bpa)**2)/(2*(std[1]**2)))
	b = ((np.sin(2*bpa))/(2*(std[0]**2))) - ((np.sin(2*bpa))/(2*(std[1]**2)))
	c = ((np.sin(bpa)**2)/(2*(std[0]**2))) + ((np.cos(bpa)**2)/(2*(std[1]**2)))

	return amplitude*np.exp(-1*a*((x-x0)**2) - b*(x-x0)*(y-y0) - c*((y-y0)**2))

def makeGaussian_astropy(size,amplitude, std, bpa, center):
	""" Make a square gaussian kernel.

	size is the length of a side of the square
	fwhm is full-width-half-maximum, which
	can be thought of as an effective radius.
	"""

	x, y = np.mgrid[0:size:size,0:size:size]
	return Gaussian2D(amplitude=amplitude, x_mean=center[0], y_mean=center[1], x_stddev=std[0], y_stddev=std[1], theta=bpa*-1, cov_matrix=None)(x,y)

def generate_fits_models_gaus(fitsfile, SN, rms, pixel_grid, type_test):
	for i in range(len(SN)):
		print('Making %s fits model for S/N %d' % (type_test,SN[i]))
		hdu_list = fits.open(fitsfile)
		head = hdu_list['PRIMARY'].header
		wcs = WCS(head)
		if type_test == 'uvplane':
			SN_fluxes = np.array(SN)*rms
			image_data0 = hdu_list[0].data*0.0
		elif type_test == 'implane':
			SN_fluxes = np.array(SN)*rms
			image_data0 = hdu_list[0].data*0.0
			rms0 = hdu_list[0].data
			rms1 = rms0
			kernel = make_Gaussian_beam_kernel(head,10)
		else:
			print('type_test can only be \'uvplane\' or \'implane\'')
		hdu_list.close()
		NAXIS1 = head['NAXIS1'] ## Assume size is square
		image_data = image_data0
		### Generate random samples
		bmaj_fwhm = (12. - 1.) * np.random.random_sample(len(pixel_grid)) + 1.
		bmin_fwhm = (12. - 1.) * np.random.random_sample(len(pixel_grid)) + 1.
		for j in range(len(bmaj_fwhm)):
			if bmaj_fwhm[j] < bmin_fwhm[j]:
				bmaj_temp = bmaj_fwhm[j]
				bmin_temp = bmin_fwhm[j]
				bmin_fwhm[j] = bmaj_temp
				bmaj_fwhm[j] = bmin_temp
		for j in range(len(bmaj_fwhm)):
			while bmin_fwhm[j]/bmaj_fwhm[j] < 0.3:
				bmin_fwhm[j] = (12. - 1.) * np.random.random_sample(1) + 1.
		bpa = (360. - 0.) * np.random.random_sample(len(pixel_grid)) + 0.
		x = []
		y = []
		for j in range(len(pixel_grid)):
			gaus = makeGaussian_bpa(size=NAXIS1, amplitude=SN_fluxes[i], std=[bmaj_fwhm[j],bmin_fwhm[j]], center=[pixel_grid[j][0],pixel_grid[j][1]], bpa=bpa[j])
			image_data[:,:] = image_data[:,:] + gaus
			x = x + [pixel_grid[j][0]]
			y = y + [pixel_grid[j][1]]
		x_world = wcs.all_pix2world(x,y,1,1,1)[0]
		y_world = wcs.all_pix2world(x,y,1,1,1)[1]
		if type_test == 'implane':
			image_data_conv = convolve(image_data, kernel)
			print(np.shape(image_data_conv))
			image_data = image_data_conv
			#image_data = np.add(image_data_conv,rms1)
		hdu = fits.PrimaryHDU(image_data, header=head)
		hdulist = fits.HDUList([hdu])
		if type_test == 'uvplane':
			hdulist.writeto(fitsfile.split('.fits')[0]+'_uv_gaus_SN%s.fits' % SN[i] ,overwrite=True)
			hdulist.close()
			os.system('rm %s_uv_gaus_SN%s_input_model.csv' % (fitsfile.split('.fits')[0],SN[i]))
			pd.DataFrame({'x':x,'y':y,'x_deg':x_world,'y_deg':y_world,'x_fwhm':bmaj_fwhm*np.sqrt(8*np.log(2)), 'y_fwhm':bmin_fwhm*np.sqrt(8*np.log(2)), 'theta':bpa}).to_csv('%s_uv_gaus_SN%s_input_model.csv' % (fitsfile.split('.fits')[0],SN[i]))
		elif type_test == 'implane':
			hdulist.writeto(fitsfile.split('.fits')[0]+'_im_gaus_SN%s.fits' % SN[i] ,overwrite=True)
			hdulist.close()
			os.system('rm %s_im_gaus_SN%s_input_model.csv' % (fitsfile.split('.fits')[0],SN[i]))
			pd.DataFrame({'x':x,'y':y,'x_deg':x_world,'y_deg':y_world,'x_fwhm':bmaj_fwhm*np.sqrt(8*np.log(2)), 'y_fwhm':bmin_fwhm*np.sqrt(8*np.log(2)), 'theta':bpa}).to_csv('%s_im_gaus_SN%s_input_model.csv' % (fitsfile.split('.fits')[0],SN[i]))
		else:
			print('hi')

def make_Gaussian_beam_kernel(header,oversampling):
	bmaj = np.abs(header['BMAJ']/header['CDELT1'])/(2*np.sqrt(2*np.log(2)))
	bmin = np.abs(header['BMIN']/header['CDELT2'])/(2*np.sqrt(2*np.log(2)))
	bpa = header['BPA']+90.
	size = int(oversampling*bmaj)
	if size % 2 == 0:## to catch non odd kernels
		size = size +1
	gauss = makeGaussian_bpa(size=size,amplitude=1, std=[bmaj,bmin], bpa=bpa,center=None)
	np.save('Gauss_model.npy',gauss)
	return CustomKernel(gauss)



### Inputs
inputs = headless('inputs.txt')
fitsfile = inputs['fitsfile']
rms = inputs['fits_rms']
if ',' in str(inputs['SN']):
	SN = inputs['SN'].split(',')
	SN = [float(i) for i in SN]
else:
	SN = [float(inputs['SN'])]
npoint = int(inputs['ngrid_points'])
pc_edge_cut = float(inputs['pc_edge_cut'])/100.
random = str(inputs['random'])
save_model = bool(inputs['save_model'])
model = str(inputs['model'])
type_test = str(inputs['type_test'])
parallel = str(inputs['parallel'])

### Convert fits to not have degenerate axes
outfitsname = fitsfile.split('.fits')[0]+'_CASA.fits'
convertAIPStoPythonImage(fitsfile,outfitsname)
header = fits.open(outfitsname)[0].header
while True:
    try:
        pixel_grid = setup_source_pixel_grid(header, npoint=npoint, pc_edge_cut=pc_edge_cut,random=random)
        break
    except:
        pass
if model == 'gaussian':
	generate_fits_models_gaus(fitsfile=outfitsname, SN=SN, rms=rms, pixel_grid=pixel_grid,type_test=type_test)
elif model == 'delta':
	if parallel == 'True':
		inputs = SN
		num_cores = multiprocessing.cpu_count()
		Parallel(n_jobs=num_cores)(delayed(generate_fits_models_delta_fcn)(fitsfile=outfitsname, SN=[i], rms=rms, pixel_grid=pixel_grid,type_test=type_test,save_model=save_model) for i in inputs)
	else:
		generate_fits_models_delta_fcn(fitsfile=outfitsname, SN=SN, rms=rms, pixel_grid=pixel_grid,type_test=type_test,save_model=save_model)
