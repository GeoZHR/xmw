"""
Demonstrate 3D seismic image processing for faults
Author: Xinming Wu, Colorado School of Mines
Version: 2016.01.22
"""

from utils import *
setupForSubset("cact")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count
# Names and descriptions of image files used below.
gxfile  = "gx" # input image (maybe after bilateral filtering)
fxfile  = "fx" # input image (maybe after bilateral filtering)
gwfile  = "gw" # input image (maybe after bilateral filtering)
hxfile  = "horizon"
gsxfile = "gsx" # image after lsf with fault likelihoods
epfile  = "ep" # eigenvalue-derived planarity
wpfile  = "wp" # eigenvalue-derived planarity
p2file  = "p2" # inline slopes
p3file  = "p3" # crossline slopes
p2kfile = "p2k" # inline slopes (known)
p3kfile = "p3k" # crossline slopes (known)
flfile  = "fl" # fault likelihood
fpfile  = "fp" # fault strike (phi)
ftfile  = "ft" # fault dip (theta)
flvfile  = "flv" # fault likelihood
fpvfile  = "fpv" # fault strike (phi)
ftvfile  = "ftv" # fault dip (theta)
fltfile = "flt" # fault likelihood thinned
fptfile = "fpt" # fault strike thinned
fttfile = "ftt" # fault dip thinned
fs1file = "fs1" # fault slip (1st component)
fs2file = "fs2" # fault slip (2nd component)
fs3file = "fs3" # fault slip (3rd component)
fskbase = "fsk" # fault skin (basename only)
fskr = "fsr" # fault skin (basename only)
fskrs = "fsrs" # fault skin (basename only)
fslbase = "fsl" # fault skin (basename only)
fskgood = "fsg" # fault skin (basename only)
fsktv = "fst" # fault skin (basename only)
fwsfile = "fws" # unfaulted image
sw1file = "sw1" # 1st component of unfaulting shifts
sw2file = "sw2" # 2nd component of unfaulting shifts
sw3file = "sw3" # 3rd component of unfaulting shifts
gufile = "gu" # flattened image
gtfile = "gt" # flattened image
ghfile = "gh" # flattened image
grfile = "gr" # flattened image
x1file = "x1" # horizon volume
u1file = "u1" # first component of normal
u2file = "u2" # second component of normal
u3file = "u3" # third component of normal
smfile = "sm"
cmfile = "cm"
clfile = "cl"

# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
minPhi,maxPhi = 0,360
minTheta,maxTheta = 75,89
sigmaPhi,sigmaTheta = 8,30

# These parameters control the construction of fault skins.
# See the class FaultSkinner for more information.
lowerLikelihood = 0.05
upperLikelihood = 0.4
minSkinSize = 500

# These parameters control the computation of fault dip slips.
# See the class FaultSlipper for more information.
minThrow = 0.0
maxThrow = 20.0

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
pngDir = None
#pngDir = "../../../png/beg/hongliu/"
#pngDir = "../../../png/nwc/"
plotOnly = False

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goDisplay()
  #goSlopes()
  #goScan()
  #goThin()
  #goSkin()
  #goSkinTv()
  #goReskin()
  #goSkinMerge()
  #goSmooth()
  #goSlip()
  #goUnfaultS()
  #goDisplay()
  #goFaultImages()
  #gx = readImage(gxfile)
  #sk = readSkins(fskr)
  #plot3(gx,skins=sk)
  #goTest()
  #goFlatten()
  goRefine()
  '''
  gu1 = readImage(gtfile)
  gu2 = readImage(gufile)
  zm = ZeroMask(0.10,1,1,1,gu2)
  gu1 = gain(gu1)
  gu2 = gain(gu2)
  zero,tiny=0.0,0.01
  zm.setValue(zero,gu2)
  plot3(gu1)
  plot3(gu2)
  '''
  #goResults()
  #goSlices()
  #goChannel()
  #goReskinx()
def goReskinx():
  print n1
  print n2
  print n3
  gx = readImage(gxfile)
  skins = readSkins(fskr)
  skins = [skins[1]]
  fsr = FaultReskin()
  skinr = fsr.applyForSkins(n1,n2,n3,400,skins)
  #writeSkins(fskrs,skinr)
  plot3(gx,skins=skinr)
  plot3(gx,skins=skins)

def goChannel():
  gx = readImage(gufile)
  gc = zerofloat(n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      gc[i3][i2] = gx[i3][i2][160]
  writeImage("gx160",gc)
  '''
  lof = LocalOrientFilterP(2,6);
  u1 = zerofloat(n1,n2,n3)
  u2 = zerofloat(n1,n2,n3)
  u3 = zerofloat(n1,n2,n3)
  lof.applyForNormal(gx,u1,u2,u3)

  gs = zerofloat(n1,n2,n3)
  lsf = LocalSmoothingFilter();
  lof = LocalOrientFilterP(2,6);
  ets = lof.applyForTensors(gx);
  ets.setEigenvalues(0.0001,0.0001,1.0);
  lsf.apply(ets,64,gx,gs);
  '''
  '''
  print min(gx)
  print max(gx)
  cs = ChannelScanner(2,5)
  #cl = cs.scan(4000,gx)
  #writeImage(clfile,cl)
  cl = readImage(clfile)
  print min(cl)
  print max(cl)
  gx = gain(gx)
  plot3(gx)
  #plot3(gs)
  plot3(cl,cmin=0.0,cmax=max(cl)/2)
  cc = cs.setSamples(0.1,cl,u1,u2,u3)
  #print len(cc)
  plot3(gx,cells=cc)
  '''


def goSlices():
  gu = readImage(gufile)
  gu = gain(gu)
  gu = mul(gu,-1)
  plot3(gu,k1=238,cmap=ColorMap.BLUE_WHITE_RED)
  '''
  for k1 in range(440,480,1):
    plot3(gu,k1=k1,png="gu"+str(k1))
  '''
def goResults():
  gx = readImage(gxfile)
  gw = readImage("fws1")
  gu = readImage(gufile)
  sk = readSkins(fskr)
  gx = gain(gx)
  gw = gain(gw)
  gu = gain(gu)
  plot3(gx,png="gx")
  flt = zerofloat(n1,n2,n3)
  fsx = FaultSkinnerX()
  fsx.getFls(sk,flt)
  plot3(gx,skins=sk,clab="Fault likelihood",png="gxf")
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  plot3(gw,png="gw")
  plot3(gu,png="gu")


def goTest():
  gx = readImage(gxfile)
  gw = readImage(fwsfile)
  gw1 = readImage("fws1")
  gw = gain(gw)
  gw1 = gain(gw1)
  zm = ZeroMask(0.10,1,1,1,gx)
  zero,tiny=0.0,0.01
  zm.setValue(zero,gw1)
  plot3(gx)
  plot3(gw)
  plot3(gw1)
def goDisplay():
  gx = readImage(gxfile)
  zm = ZeroMask(0.01,0,0,0,gx)
  zero,tiny=0.0,0.01
  gx = gain(gx)
  zm.setValue(zero,gx)
  writeImage(gxfile,gx)
  plot3(gx)
def goSlopes():
  print "goSlopes ..."
  gx = readImage(gxfile)
  sigma1,sigma2,sigma3,pmax = 16.0,2.0,2.0,5.0
  p2,p3,ep = FaultScanner.slopes(sigma1,sigma2,sigma3,pmax,gx)
  writeImage(p2file,p2)
  writeImage(p3file,p3)
  writeImage(epfile,ep)
  print "p2  min =",min(p2)," max =",max(p2)
  print "p3  min =",min(p3)," max =",max(p3)
  print "ep min =",min(ep)," max =",max(ep)
  plot3(gx,p2, cmin=-1,cmax=1,cmap=bwrNotch(1.0),
        clab="Inline slope (sample/sample)",png="p2")
  plot3(gx,p3, cmin=-1,cmax=1,cmap=bwrNotch(1.0),
        clab="Crossline slope (sample/sample)",png="p3")
  plot3(gx,sub(1,ep),cmin=0,cmax=1,cmap=jetRamp(1.0),
        clab="Planarity")

def goScan():
  print "goScan ..."
  if not plotOnly:
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    gx = readImage(gxfile)
    zm = ZeroMask(0.10,1,1,1,gx)
    zero,tiny=0.0,0.01
    zm.setValue(tiny,gx)
    gx = FaultScanner.taper(10,0,0,gx)
    fs = FaultScanner(sigmaPhi,sigmaTheta)
    fl,fp,ft = fs.scan(minPhi,maxPhi,minTheta,maxTheta,p2,p3,gx)
    print "fl min =",min(fl)," max =",max(fl)
    print "fp min =",min(fp)," max =",max(fp)
    print "ft min =",min(ft)," max =",max(ft)
    writeImage(flfile,fl)
    writeImage(fpfile,fp)
    writeImage(ftfile,ft)
  else:
    gx = readImage(gxfile)
    fl = readImage(flfile)
    fp = readImage(fpfile)
    ft = readImage(ftfile)
  '''
  plot3(gx,fl,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
      clab="Fault likelihood",png="fl")
  plot3(gx,fp,cmin=0,cmax=360,cmap=hueFill(1.0),
      clab="Fault strike (degrees)",cint=45,png="fp")
  plot3(gx,ft,cmin=65,cmax=89,cmap=jetFill(1.0),
      clab="Fault dip (degrees)",png="ft")
  '''

def goThin():
  print "goThin ..."
  gx = readImage(gxfile)
  if not plotOnly:
    fl = readImage(flfile)
    fp = readImage(fpfile)
    ft = readImage(ftfile)
    flt,fpt,ftt = FaultScanner.thin([fl,fp,ft])
    writeImage(fltfile,flt)
    writeImage(fptfile,fpt)
    writeImage(fttfile,ftt)
  else:
    flt = readImage(fltfile)
    fpt = readImage(fptfile)
    ftt = readImage(fttfile)
  plot3(gx,flt,cmin=0.1,cmax=0.7,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  '''
  plot3(gx,ftt,cmin=65,cmax=85,cmap=jetFillExceptMin(1.0),
        clab="Fault dip (degrees)",png="ftt")
  plot3(gx,fpt,cmin=0,cmax=360,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=45,png="fpt")
  '''

def goSkinTv():
  print "go skin..."
  gx = readImage(gxfile)
  if not plotOnly:
    fs = FaultScanner(sigmaPhi,sigmaTheta)
    sp = fs.makePhiSampling(minPhi,maxPhi)
    st = fs.makeThetaSampling(minTheta,maxTheta)
    fsx = FaultSkinnerX()
    fsx.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
    fsx.setMinSkinSize(minSkinSize)
    fsx.setMaxPlanarDistance(0.2)
    fsk = readSkins(fskbase)
    fr = FaultReskin()
    fcs = FaultSkin.getCells(fsk)
    cells = []
    for ic in range(0,len(fcs),4):
      cells.append(fcs[ic])
    fsx.resetCells(cells)
    fsx.setGaussWeights(sp,st)
    skins = fsx.findSkins(n1,n2,n3,cells)
    removeAllSkinFiles(fsktv)
    writeSkins(fsktv,skins)
  else:
    skins = readSkins(fsktv)
  plot3(gx,skins=skins)

def goSkin():
  print "goSkin ..."
  gx = readImage(gxfile)
  if not plotOnly:
    fl = readImage(flfile)
    fp = readImage(fpfile)
    ft = readImage(ftfile)
    fs = FaultSkinner()
    fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
    fs.setMaxDeltaStrike(10)
    fs.setMaxPlanarDistance(0.2)
    fs.setMinSkinSize(minSkinSize)
    cells = fs.findCells([fl,fp,ft])
    skins = fs.findSkins(cells)
    for skin in skins:
      skin.smoothCellNormals(4)
    fd = FaultDisplay()
    #skins = fd.getLowerFaults(350,skins)
    #plot3(gx,skins=skins)
    #sk = fd.reskin(160,350,skins[1])
    #plot3(gx,skins=sk)
    #skins[1] = sk[1]
    print "total number of cells =",len(cells)
    print "total number of skins =",len(skins)
    print "number of cells in skins =",FaultSkin.countCells(skins)
    removeAllSkinFiles(fskbase)
    writeSkins(fskbase,skins)
  else:
    skins = readSkins(fskbase)
  #plot3(gx,cells=cells)
  plot3(gx,skins=skins)
  for iskin,skin in enumerate(skins):
    plot3(gx,skins=[skin],clab="skin"+str(iskin))

  '''
  flt = zerofloat(n1,n2,n3)
  fsx = FaultSkinnerX()
  fsx.getFl(skins,flt)
  plot3(gx,skins=skins)
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  '''
def goSkinMerge():
  gx = readImage(gxfile)
  if not plotOnly:
    skins = readSkins(fskbase)
    fsc = FaultScanner(sigmaPhi,sigmaTheta)
    sp = fsc.makePhiSampling(minPhi,maxPhi)
    st = fsc.makeThetaSampling(minTheta,maxTheta)

    fs = FaultSkinner()
    fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
    fs.setMaxDeltaStrike(10)
    fs.setMaxPlanarDistance(0.2)
    fs.setMinSkinSize(minSkinSize)

    fr = FaultReskin()
    sks1 = [skins[0],skins[2],skins[3],skins[5]] 
    cells = FaultSkin.getCells(sks1)
    fl,fp,ft = fr.faultImagesFromCells(n1,n2,n3,cells)
    div(fl,max(fl),fl)
    cells = fs.findCells([fl,fp,ft])
    skt1 = fs.findSkins(cells)

    cells = FaultSkin.getCells(skt1)
    fl,fp,ft = fr.faultImagesFromCells(n1,n2,n3,cells)
    div(fl,max(fl),fl)
    cells = fs.findCells([fl,fp,ft])
    skt1 = fs.findSkins(cells)


    sks2 = [skins[1]] #[skins[11],skins[6]] #[skins[7],skins[9]]#skins[5]
    cells = FaultSkin.getCells(sks2)
    fl,fp,ft = fr.faultImagesFromCells(n1,n2,n3,cells)
    div(fl,max(fl),fl)
    cells = fs.findCells([fl,fp,ft])
    skt2 = fs.findSkins(cells)

    sks3 = [skins[4]] #[skins[11],skins[6]] #[skins[7],skins[9]]#skins[5]
    cells = FaultSkin.getCells(sks3)
    fl,fp,ft = fr.faultImagesFromCells(n1,n2,n3,cells)
    div(fl,max(fl),fl)
    cells = fs.findCells([fl,fp,ft])
    skt3 = fs.findSkins(cells)


    sks4 = [skins[6]] #[skins[11],skins[6]] #[skins[7],skins[9]]#skins[5]
    cells = FaultSkin.getCells(sks4)
    fl,fp,ft = fr.faultImagesFromCells(n1,n2,n3,cells)
    div(fl,max(fl),fl)
    cells = fs.findCells([fl,fp,ft])
    skt4 = fs.findSkins(cells)

    skrs = [skt1[0],skt2[0],skt3[0],skt4[0]]

    removeAllSkinFiles(fskr)
    writeSkins(fskr,skrs)
  else:
    skrs = readSkins(fskr)
  plot3(gx,skins=skrs)


def goReskin(): 
  gx = readImage(gxfile)
  if not plotOnly:
    skins = readSkins(fskbase)
    fr = FaultReskin()
    sks = fr.reskin(176,0,skins[0])
    skins[0] = sks[0]
    #plot3(gx,skins=skins)
    fs = FaultScanner(sigmaPhi,sigmaTheta)
    sp = fs.makePhiSampling(minPhi,maxPhi)
    st = fs.makeThetaSampling(minTheta,maxTheta)

    fs = FaultSkinner()
    fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
    fs.setMaxDeltaStrike(10)
    fs.setMaxPlanarDistance(0.2)
    fs.setMinSkinSize(minSkinSize)
    cells = FaultSkin.getCells(skins)
    #plot3(gx,cells=cells)
    fr = FaultReskin()
    fl,fp,ft = fr.rescan(n1,n2,n3,sp,st,cells)
    div(fl,max(fl),fl)
    cells = fs.findCells([fl,fp,ft])
    skins = fs.findSkins(cells)
    removeAllSkinFiles(fsktv)
    writeSkins(fsktv,skins)
  else:
    skins = readSkins(fsktv)
  #plot3(gx,skins=skins)
  '''
  fr = FaultReskin()
  sks = [skins[1],skins[3]] #[skins[11],skins[6]] #[skins[7],skins[9]]#skins[5]
  cells = FaultSkin.getCells(sks)
  fs = FaultScanner(sigmaPhi,sigmaTheta)
  sp = fs.makePhiSampling(minPhi,maxPhi)
  st = fs.makeThetaSampling(minTheta,maxTheta)

  fl,fp,ft = fr.faultImagesFromCells(n1,n2,n3,cells)
  div(fl,max(fl),fl)
  fs = FaultSkinner()
  fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  fs.setMaxDeltaStrike(10)
  fs.setMaxPlanarDistance(0.2)
  fs.setMinSkinSize(minSkinSize)
  cells = fs.findCells([fl,fp,ft])
  sks = fs.findSkins(cells)
  removeAllSkinFiles(fskr)
  writeSkins(fskr,sks)
  '''
  plot3(gx,skins=skins)
  k = 0
  for skin in skins:
    plot3(gx,skins=[skins[k]],clab=str(k))
    k = k+1

  '''
  plot3(gx,fl,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
      clab="Fault likelihood",png="fl")
  plot3(gx,fp,cmin=0,cmax=360,cmap=hueFill(1.0),
      clab="Fault strike (degrees)",cint=45,png="fp")
  plot3(gx,convertDips(ft),cmin=15,cmax=55,cmap=jetFill(1.0),
      clab="Fault dip (degrees)",png="ft")
  '''


def goFaultImages():
  gx = readImage(gxfile)
  if not plotOnly:
    fd = FaultDisplay()
    skins = readSkins(fsktv)
    flt = fillfloat(-0.001,n1,n2,n3)
    fpt = fillfloat(-0.001,n1,n2,n3)
    ftt = fillfloat(-0.001,n1,n2,n3)
    fd = FaultDisplay()
    fd.getFlt(skins,gx,flt)
    fd.getFpt(skins,gx,fpt)
    fd.getFtt(skins,gx,ftt)
    writeImage(fltfile,flt)
    writeImage(fptfile,fpt)
    writeImage(fttfile,ftt)
  else:
    #flt = readImage(fltfile)
    fpt = readImage(fptfile)
    #ftt = readImage(fttfile)
  '''
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  plot3(gx,ftt,cmin=65,cmax=85,cmap=jetFillExceptMin(1.0),
        clab="Fault dip (degrees)",png="ftt")
  '''
  plot3(gx,fpt,cmin=0,cmax=180,cmap=jetFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=20,png="fpt")
def goSmooth():
  print "goSmooth ..."
  flstop = 0.1
  fsigma = 8.0
  fl = readImage(flfile)
  gx = readImage(gxfile)
  skins = readSkins(fskr)
  flt = zerofloat(n1,n2,n3)
  fsx = FaultSkinnerX()
  fsx.getFl(skins,flt)
  p2,p3,ep = FaultScanner.slopes(8.0,1.0,1.0,5.0,gx)
  gsx = FaultScanner.smooth(flstop,fsigma,p2,p3,flt,gx)
  writeImage(gsxfile,gsx)
  plot3(gsx,png="gsx")

def goSlip():
  print "goSlip ..."
  gx = readImage(gxfile)
  if not plotOnly:
    gsx = readImage(gsxfile)
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    skins = readSkins(fskr)
    fsl = FaultSlipper(gsx,p2,p3)
    fsl.setOffset(2.0) # the default is 2.0 samples
    fsl.setZeroSlope(False) # True only if we want to show the error
    fsl.computeDipSlips(skins,minThrow,maxThrow)
    print "  dip slips computed, now reskinning ..."
    print "  number of skins before =",len(skins),
    fsk = FaultSkinner() # as in goSkin
    fsk.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
    fsk.setMinSkinSize(minSkinSize)
    fsk.setMinMaxThrow(minThrow,maxThrow)
    #skins = fsk.reskin(skins)
    print ", after =",len(skins)
    removeAllSkinFiles(fslbase)
    writeSkins(fslbase,skins)
    smark = -999.999
    #s1,s2,s3 = fsl.getDipSlips(skins,smark)
    #s1,s2,s3 = fsl.interpolateDipSlips([s1,s2,s3],smark)
    #gw = fsl.unfault([s1,s2,s3],gx)
    #writeImage(gwfile,gw)
    '''
    writeImage(fs1file,s1)
    writeImage(fs2file,s2)
    writeImage(fs3file,s3)
    '''
  else:
    #gw = readImage(gwfile)
    #s1 = readImage(fs1file)
    skins = readSkins(fslbase)
    #skinr = readSkins(fskr)
  plot3(gx,skins=skins,smax=maxThrow)
  #plot3(gx,skins=skinr)
  '''
  plot3(gx,s1,cmin=-10,cmax=10.0,cmap=jetFillExceptMin(1.0),
        clab="Fault throw (samples)",png="gxs1")
  plot3(gx,s1,cmin=0.0,cmax=10.0,cmap=jetFill(0.3),
        clab="Vertical shift (samples)",png="gxs1i")
  plot3(gx,s2,cmin=-2.0,cmax=2.0,cmap=jetFill(0.3),
        clab="Inline shift (samples)",png="gxs2i")
  plot3(gx,s3,cmin=-1.0,cmax=1.0,cmap=jetFill(0.3),
        clab="Crossline shift (samples)",png="gxs3i")
  plot3(gx)
  plot3(gw,png="gw")
  '''

def goUnfaultS():
  print "goUnfault ..."
  gx = readImage(gxfile)
  if not plotOnly:
    fw = zerofloat(n1,n2,n3)
    lof = LocalOrientFilter(8.0,4.0,4.0)
    et = lof.applyForTensors(gx)
    et.setEigenvalues(0.001,1.0,1.0)

    wp = fillfloat(1.0,n1,n2,n3)
    skins = readSkins(fslbase)
    fsc = FaultSlipConstraints(skins)
    sp = fsc.screenPoints(wp)
    mul(sp[3][0],10,sp[3][0])

    uf = UnfaultS(8.0,8.0)
    uf.setIters(100)
    uf.setTensors(et)
    [t1,t2,t3] = uf.findShifts(sp,wp)
    #uf.convertShifts(40,[t1,t2,t3])
    uf.applyShifts([t1,t2,t3],gx,fw)
    writeImage(fwsfile,fw)
    writeImage(sw1file,t1)
    writeImage(sw2file,t2)
    writeImage(sw3file,t3)
  else :
    #gw = readImage(gwfile)
    fw = readImage(fwsfile)
    #fw = readImage("fwt")
  fw = gain(fw)
  plot3(gx,png="gxuf")
  plot3(fw,png="fwuf")
  '''
  plot3(gw,png="fwuf")
  skins = readSkins(fslbase)
  mark = -999.99
  s1 = fillfloat(mark,n1,n2,n3)
  FaultSkin.getThrow(mark,skins,s1)
  plot3(gx,s1,cmin=-10,cmax=10.0,cmap=jetFillExceptMin(1.0),
        clab="Fault throw (samples)",png="gxs1")
  plot3(gx,t1,cmin=-6.0,cmax=6.0,cmap=jetFill(0.3),
        clab="Vertical shift (samples)",png="gxs1i")
  plot3(gx,t2,cmin=-2.0,cmax=2.0,cmap=jetFill(0.3),
        clab="Inline shift (samples)",png="gxs2i")
  plot3(gx,t3,cmin=-1.0,cmax=1.0,cmap=jetFill(0.3),
        clab="Crossline shift (samples)",png="gxs3i")
  '''
def goFlatten():
  fx = readImage(fwsfile)
  p2 = zerofloat(n1,n2,n3)
  p3 = zerofloat(n1,n2,n3)
  ep = zerofloat(n1,n2,n3)
  lsf = LocalSlopeFinder(2.0,1.0)
  lsf.findSlopes(fx,p2,p3,ep);
  ep = pow(ep,4)
  fl = Flattener3()
  fl.setIterations(0.01,300)
  fm = fl.getMappingsFromSlopes(s1,s2,s3,p2,p3,ep)
  gt = fm.flatten(fx)
  writeImage(gtfile,gt)
  #writeImage(gufile,gt)
  gt = readImage(gtfile)
  fx = gain(fx)
  gt = gain(gt)
  plot3(fx)
  plot3(gt)

def goRefine():
  gt = readImage(gtfile)
  gr = zerofloat(n1,n2,n3)
  g0 = gt[0][0]
  for i3 in range(n3):
    for i2 in range(n2):
      gr[i3][i2] = g0
  gt = gain(gt)
  gr = gain(gr)
  plot3(gr)
  dw = DynamicWarping(-10,10)
  dw.setStrainMax(0.25,0.25,0.25)
  dw.setErrorSmoothing(3)
  dw.setShiftSmoothing(1)
  us = dw.findShifts(gr,gt)
  gr = dw.applyShifts(us,gt)
  plot3(gt)
  plot3(gr)


def like(x):
  n3 = len(x)
  n2 = len(x[0])
  n1 = len(x[0][0])
  return zerofloat(n1,n2,n3)

def slog(f):
  return mul(sgn(f),log(add(1.0,abs(f))))

def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(100.0)
  ref.apply1(g,g)
  y = like(x)
  div(x,sqrt(g),y)
  return y
def array(x1,x2,x3=None,x4=None):
  if x3 and x4:
    return jarray.array([x1,x2,x3,x4],Class.forName('[[[F'))
  elif x3:
    return jarray.array([x1,x2,x3],Class.forName('[[[F'))
  else:
    return jarray.array([x1,x2],Class.forName('[[[F'))




#############################################################################
# graphics

def jetFill(alpha):
  return ColorMap.setAlpha(ColorMap.JET,alpha)
def jetFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.JET,a)
def jetRamp(alpha):
  return ColorMap.setAlpha(ColorMap.JET,rampfloat(0.0,alpha/256,256))
def bwrFill(alpha):
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,alpha)
def bwrNotch(alpha):
  a = zerofloat(256)
  for i in range(len(a)):
    if i<128:
      a[i] = alpha*(128.0-i)/128.0
    else:
      a[i] = alpha*(i-127.0)/128.0
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,a)
def hueFill(alpha):
  return ColorMap.getHue(0.0,1.0,alpha)
def hueFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.getHue(0.0,1.0),a)

def addColorBar(frame,clab=None,cint=None):
  cbar = ColorBar(clab)
  if cint:
    cbar.setInterval(cint)
  cbar.setFont(Font("Arial",Font.PLAIN,32)) # size by experimenting
  cbar.setWidthMinimum
  cbar.setBackground(Color.WHITE)
  frame.add(cbar,BorderLayout.EAST)
  return cbar

def convertDips(ft):
  return FaultScanner.convertDips(0.2,ft) # 5:1 vertical exaggeration

def plot3(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          horizon=None,xyz=None,cells=None,skins=None,smax=0.0,slices=None,
          k1=n1/2,links=False,curve=False,trace=False,png=None):
  n3 = len(f)
  n2 = len(f[0])
  n1 = len(f[0][0])
  s1,s2,s3=Sampling(n1),Sampling(n2),Sampling(n3)
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  sf = SimpleFrame(AxesOrientation.XRIGHT_YOUT_ZDOWN)
  cbar = None
  if g==None:
    ipg = sf.addImagePanels(s1,s2,s3,f)
    if cmap!=None:
      ipg.setColorModel(cmap)
    if cmin!=None and cmax!=None:
      ipg.setClips(cmin,cmax)
    else:
      ipg.setClips(-2.0,2.0)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-2,2)
    if cmin!=None and cmax!=None:
      ipg.setClips2(cmin,cmax)
    if cmap==None:
      cmap = jetFill(0.8)
    ipg.setColorModel2(cmap)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMap2Listener(cbar)
    sf.world.addChild(ipg)
  if cbar:
    cbar.setWidthMinimum(120)
  if xyz:
    pg = PointGroup(0.2,xyz)
    ss = StateSet()
    cs = ColorState()
    cs.setColor(Color.YELLOW)
    ss.add(cs)
    pg.setStates(ss)
    #ss = StateSet()
    #ps = PointState()
    #ps.setSize(5.0)
    #ss.add(ps)
    #pg.setStates(ss)
    sf.world.addChild(pg)
  if cells:
    ss = StateSet()
    lms = LightModelState()
    lms.setTwoSide(True)
    ss.add(lms)
    ms = MaterialState()
    ms.setSpecular(Color.GRAY)
    ms.setShininess(100.0)
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE)
    ms.setEmissiveBack(Color(0.0,0.0,0.5))
    ss.add(ms)
    cmap = ColorMap(0.0,0.5,ColorMap.JET)
    xyz,uvw,rgb = FaultCell.getXyzUvwRgbForLikelihood(0.7,cmap,cells,False)
    qg = QuadGroup(xyz,uvw,rgb)
    qg.setStates(ss)
    sf.world.addChild(qg)
  if horizon:
    sd = SurfaceDisplay()
    ts = sd.horizonWithAmplitude([-0.5,0.5],horizon,f)
    tg = TriangleGroup(True,ts[0],ts[1])
    #tg = TriangleGroup(True,s3,s2,horizon)
    #tg.setColor(Color.CYAN)
    sf.world.addChild(tg)
  if skins:
    sg = Group()
    ss = StateSet()
    lms = LightModelState()
    lms.setLocalViewer(True)
    lms.setTwoSide(True)
    ss.add(lms)
    ms = MaterialState()
    ms.setSpecular(Color.GRAY)
    ms.setShininess(100.0)
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE)
    if not smax:
      ms.setEmissiveBack(Color(0.0,0.0,0.5))
    ss.add(ms)
    sg.setStates(ss)
    size = 2.0
    if links:
      size = 0.65 
      ls = LineState()
      ls.setWidth(1.5)
      ls.setSmooth(True)
      ss.add(ls)
    ct = 0
    for skin in skins:
      if smax>0.0: # show fault throws
        cmap = ColorMap(-smax,smax,ColorMap.JET)
        xyz,uvw,rgb = skin.getCellXyzUvwRgbForThrow(size,cmap,False)
      else: # show fault likelihood
        cmap = ColorMap(0.0,1.0,ColorMap.JET)
        xyz,uvw,rgb = skin.getCellXyzUvwRgbForLikelihood(size,cmap,False)
      qg = QuadGroup(xyz,uvw,rgb)
      qg.setStates(None)
      sg.addChild(qg)
      if curve or trace:
        cell = skin.getCellNearestCentroid()
        if curve:
          xyz = cell.getFaultCurveXyz()
          pg = PointGroup(0.5,xyz)
          sg.addChild(pg)
        if trace:
          xyz = cell.getFaultTraceXyz()
          pg = PointGroup(0.5,xyz)
          sg.addChild(pg)
      if links:
        if ct==0:
          r,g,b=0,0,0
        if ct==1:
          r,g,b=0,0,1
        if ct==2:
          r,g,b=0,1,1
        if ct==3:
          #r,g,b=0.627451,0.12549,0.941176
          r,g,b=1,1,1
        r,g,b=0,0,1
        xyz = skin.getCellLinksXyz()
        #rgb = skin.getCellLinksRgb(r,g,b,xyz)
        #lg = LineGroup(xyz,rgb)
        lg = LineGroup(xyz)
        sg.addChild(lg)
        #ct = ct+1
    sf.world.addChild(sg)
  #ipg.setSlices(450,530,393)
  ipg.setSlices(k1,596,n3)
  if cbar:
    sf.setSize(1037,900)
  else:
    sf.setSize(900,900)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  zscale = 0.8*max(n2*d2,n3*d3)/(n1*d1)
  #zscale = 1.5*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.6)
  #ov.setScale(2.5)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(0.0,-0.15,-0.01))
  ov.setAzimuthAndElevation(225.0,40.0)
  #ov.setAzimuthAndElevation(-55.0,35.0)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")


#############################################################################
run(main)
