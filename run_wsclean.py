import os, sys

print(sys.argv)
imagename = sys.argv[1]
imsize1 = sys.arg[2]
imsize2 = sys.argv[3]
cell =  sys.argv[4]
niter = sys.argv[5]
vis = sys.argv[6]


os.system('wsclean -name %s -size %d %d -scale %s -auto-threshold 0.3 -auto-mask 3 -mgain 0.9 -niter %d -weight natural %s' % (imagename,imsize1,imsize2,cell,niter,vis))
