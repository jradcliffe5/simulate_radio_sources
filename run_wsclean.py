import os, sys

print(sys.argv)
imagename = str(sys.argv[1])
imsize1 = int(sys.argv[2])
imsize2 = int(sys.argv[3])
cell =  str(sys.argv[4])
niter = int(sys.argv[5])
bmaj = str(sys.argv[6])
bmin = str(sys.argv[7])
bpa = str(sys.argv[8])
vis = str(sys.argv[9])

if (float(bmaj)==0) or (float(bmin)==0):
	rest_beam = ''
else:
	rest_beam = '-beam-shape %s %s %s' % (bmaj,bmin,bpa)


os.system('/Volumes/HD-LXU3/anaconda2/envs/wsclean2.5/bin/wsclean -name %s -size %d %d -scale %s -auto-threshold 0.3 -auto-mask 3 -mgain 0.9 -niter %d -weight natural %s %s' % (imagename,imsize1,imsize2,cell,niter,rest_beam,vis))
