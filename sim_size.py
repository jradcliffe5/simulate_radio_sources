def convert_fits_to_casa(fitsimage):
    importfits(fitsimage=fitsimage, overwrite=True, imagename=fitsimage.split('.fits')[0]+'.im')

def pull_casa_im_info(casa_image):
    iminfo = {}
    imheader = imhead(imagename=casa_image)
    iminfo['RA'] = imheader['refval'][0]
    iminfo['RA_range'] = np.array([(imheader['refval'][0]-((imheader['refpix'][0]+1)*imheader['incr'][0])),\
    (imheader['shape'][0] - (imheader['refpix'][0]+1))*imheader['incr'][0] + imheader['refval'][0]])
    iminfo['RA_unit'] = imheader['axisunits'][0]
    iminfo['DEC'] = imheader['refval'][1]
    iminfo['DEC_range'] = np.array([(imheader['refval'][1]-((imheader['refpix'][1]+1)*imheader['incr'][1])),\
    (imheader['shape'][1] - (imheader['refpix'][1]+1))*imheader['incr'][1] + imheader['refval'][1]])
    iminfo['DEC_unit'] = imheader['axisunits'][1]
    iminfo['frequency'] = imheader['refval'][2]
    iminfo['frequency_unit'] = imheader['axisunits'][3]
    return iminfo

def setup_source_grid(iminfo, npoint, pc_edge_cut):
    RA = np.sort(iminfo['RA_range'])
    DEC = np.sort(iminfo['DEC_range'])
    RA_pc_cut = (RA[1]-RA[0])*pc_edge_cut
    RA = np.array([RA[0]+RA_pc_cut, RA[1]-RA_pc_cut])
    DEC_pc_cut = (DEC[1]-DEC[0])*pc_edge_cut
    DEC = np.array([DEC[0]+RA_pc_cut, DEC[1]-RA_pc_cut])
    RA = np.linspace(RA[0],RA[1],npoint,endpoint=True)
    DEC = np.linspace(DEC[0],DEC[1],npoint,endpoint=True)
    pointings = []
    for i in RA:
        for j in DEC:
            pointings = pointings + [['J2000','%s%s' % (i,iminfo['RA_unit']), '%s%s' % (j,iminfo['DEC_unit'])]]
    return pointings



### Inputs
fitsimage = 'HDFL0074_PBCOR_NA_IM.fits'
rms = 9.9E-6 ##

imagename = fitsimage.split('.fits')[0]+'.im'
imagename0 = fitsimage.split('.fits')[0]+'0.im'

### 1. First convert fits to casa image
convert_fits_to_casa(fitsimage=fitsimage)
### 2. Pull information from image
iminfo =  pull_casa_im_info(imagename)
pointings = setup_source_grid(iminfo,10,0.1)
### 3. Zero the image to preserve the fits header
os.system('rm -r %s' % imagename0)
immath(imagename=imagename,outfile=imagename0, mode='evalexpr',expr='IM0*0.0')
### 4. Create component list
#cl.open()
for i in pointings:
    print i
    cl.addcomponent(flux=1, fluxunit='mJy', polarization='Stokes',dir=i, shape='gaussian', majoraxis='0.01arcsec',minoraxis='0.01arcsec', positionangle='0deg', freq='1.6GHz',spectrumtype='spectral index', index=-0.8)
os.system('rm -r test.cl')
cl.rename('test.cl')

ia.open(imagename0)
cl2=cltool()
cl2.fromrecord(ia.deconvolvecomponentlist(cl.torecord()))
cl2.rename('deconvolved_sources.cl')

cl.close()
cl.done()
cl2.close()
ia.close()
