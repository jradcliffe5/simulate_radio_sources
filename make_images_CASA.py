import numpy as np
#from code_mailer import headless
import pyfits
import os,re, sys

def headless(inputfile):
	''' Parse the list of inputs given in the specified file. (Modified from evn_funcs.py)'''
	INPUTFILE = open(inputfile, "r")
	control = {}
	# a few useful regular expressions
	newline = re.compile(r'\n')
	space = re.compile(r'\s')
	char = re.compile(r'\w')
	comment = re.compile(r'#.*')
	# parse the input file assuming '=' is used to separate names from values
	for line in INPUTFILE:
		if char.match(line):
			line = comment.sub(r'', line)
			line = line.replace("'", '')
			(param, value) = line.split('=')
			param = newline.sub(r'', param)
			param = param.strip()
			param = space.sub(r'', param)
			value = newline.sub(r'', value)
			value = value.replace(' ','').strip()
			valuelist = value.split(',')
			if len(valuelist) == 1:
				if valuelist[0] == '0' or valuelist[0]=='1' or valuelist[0]=='2':
					control[param] = int(valuelist[0])
				else:
					control[param] = str(valuelist[0])
			else:
				control[param] = ','.join(valuelist)
	return control

inputs = headless('inputs.txt')
fitsfile = str(inputs['fitsfile'])
uvfile = inputs['uvfile']
rms = inputs['fits_rms']
model = inputs['model']
wsclean_run = inputs['wsclean_run']
wsclean_beam = inputs['wsclean_beam'].split(',')



if model == 'delta':
	extension = 'delt'
if model == 'gaussian':
	extension = 'gaus'
else:
	print 'Model can only be \'delta\' or \'gaussian\''
	exit()

for i in os.listdir('./'):
	print '%s_CASA_uv_%s_SN' % (fitsfile.split('.fits')[0],extension)
	print (i.startswith('%s_CASA_uv_%s_SN' % (fitsfile.split('.fits')[0],extension)))
	if (i.startswith('%s_CASA_uv_%s_SN' % (fitsfile.split('.fits')[0],extension))) and (i.endswith('.fits')):
		print i
		hdu = pyfits.open(i)
		header = hdu['PRIMARY'].header
		hdu.close()
		CASA_model = i.split('.fits')[0]+'.mod'
		UV_temp = uvfile.split('.ms')[0]+'_temp.ms'

		importfits(fitsimage=i,imagename=CASA_model)
		if os.path.exists(UV_temp) == True:
			clearcal(vis=UV_temp)
		else:
			os.system('cp -r %s %s' % (uvfile, UV_temp))
		#partition(vis=uvfile,outputvis=uvfile.split('.ms')[0]+'_temp.ms', createmms=True)
		ft(vis=UV_temp,model=CASA_model, usescratch=True)
		uvsub(vis=UV_temp,reverse=True)
		if wsclean_run == 'True':
			os.system('python run_wsclean.py %s %d %d %s %d %s %s %s %s' % (i.split('.fits')[0], int(header['NAXIS1']), int(header['NAXIS2']), '0.3arcsec', 30000, wsclean_beam[0], wsclean_beam[1], wsclean_beam[2], UV_temp))
		else:
			tclean(vis=UV_temp,imagename=i.split('.fits')[0],imsize=[header['NAXIS1'],header['NAXIS2']], cell='0.3arcsec', niter=30000,nsigma=2.5, deconvolver='clark',parallel=True, usemask='auto-multithresh')
		#os.system('rm -r %s*' % UV_temp)
		if wsclean_run == 'True':
			os.system('mv %s-image.fits %s.image.fits' % (i.split('.fits')[0],i.split('.fits')[0]))
		else:
			exportfits(imagename=i.split('.fits')[0]+'.image',fitsimage=i.split('.fits')[0]+'.image.fits')
