"""
Demonstrate 3D seismic image processing for faults
Author: Xinming Wu, Colorado School of Mines
Version: 2016.01.22
"""

from utils import *
#setupForSubset("nathan")
setupForSubset("nathanSub5")
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

# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
minPhi,maxPhi = 0,360
minTheta,maxTheta = 65,85
sigmaPhi,sigmaTheta = 20,60

# These parameters control the construction of fault skins.
# See the class FaultSkinner for more information.
lowerLikelihood = 0.3
upperLikelihood = 0.7
minSkinSize = 500

# These parameters control the computation of fault dip slips.
# See the class FaultSlipper for more information.
minThrow = 0.0
maxThrow = 85.0

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
#pngDir = "../../../png/beg/hongliu/"
pngDir = "../../../png/beg/nathan/sub5/"
pngDir = None
plotOnly = False

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goSlopes()
  #goScan()
  #goThin()
  #goSkin()
  #goSkinTv()
  #goReskin()
  #goSmooth()
  #goSlip()
  #goUnfaultS()
  #goDisplay()
  #goFaultImages()
  #goFaultPoints()
  #getOceanBottom()
  #goSeisResample()
  #goRose()
  goStrikeRotation()
  #goMask()
  '''
  gx = readImage(gxfile)
  writeImage("gx70",gx[70])
  fp = readImage(fptfile)
  fps = copy(300,n2,n3,fp)
  writeImage("fps",fps)
  '''
def goMask():
  fp = readImage(fptfile)
  ob = readImage2D(n2,n3,"ob")
  hp = Helper()
  hp.mask(ob,fp)
  writeImage(fptfile,fp)
  plot3(fp,fp,cmin=0,cmax=180,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=20,png="fpt")

def goFaultPoints():
  fp = readImage(fptfile)
  rp = RosePlot()
  ob = readImage2D(n2,n3,"ob")
  ps = rp.faultPoints(ob,fp)
  print len(ps)
  print len(ps[0])
  writeImage("fps",ps)

def goSeisResample():
  hp = Helper()
  d3 = 12.5
  gx = readImage(gxfile)
  c3 = Sampling(n3,18.75,0.0)
  gi = zerofloat(n1,n2,round(n3*1.5))
  hp.resample(s1,s2,c3,d3,gx,gi)
  writeImage("gi",gi)
  #plot3(gx)
  #plot3(gx,cmin=-10000,cmax=10000)
def goStrikeRotation():
  '''
  gx = readImage(gxfile)
  fp = readImage(fptfile)
  '''
  #gx = readImage("gxs")
  #fp = readImage("fpp")
  gx = readImage("gxs2")
  fp = readImage("fps2")
  hpr = Helper()
  hpr.rotateX(299,fp)
  hpr.convert(fp)
  plot3(gx,fp,cmin=0,cmax=180,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=20,png="fpt")

def goRose():
  rp = RosePlot()
  ob = readImage2D(n2,n3,"ob")
  #fp = readImage2D(104068862,4,"fpp")
  fp = readImage2D(107757761,4,"fps")
  rp.rotateX(299,fp[3])
  rp.convert(fp[3])
  fc = rp.removeSignature(29,fp)
  #rp.rose(fp[3],36)
  #c2,c3=10,2
  c2,c3=15,3
  #c2,c3=20,4
  for tp  in range(0,480,20):
  #for tp  in range(200,400,200):
    bt = tp+20
    title1 = str(tp*5)
    title2 = str(bt*5)
    title =title1+"~"+title2+" m"
    pp = rp.applyForRosePlotsX(tp,bt,c2,c3,n2,n3,36,fc,ob)
    #pp = rp.applyForRosePlotsX(99,99,c2,c3,n2,n3,36,fc,ob)
    pp.addTitle(title)
    pf = PlotFrame(pp)
    #wx,wy = 1450,round((c3*1450)/c2)+50
    wx,wy = 1700,round((c3*1750)/c2)+50
    #wx,wy = 2100,round((c3*2100)/c2)+100
    pf.setSize(wx,wy)
    pf.setVisible(True)
    pf.paintToPng(720,6,pngDir+str(tp)+".png")

def getOceanBottom():
  hp = Helper()
  gx = readImage(gxfile)
  ob = hp.getOceanBottom(0.1,gx)
  writeImage("ob",ob)
  #ob = readImage2D(n2,n3,"ob")
  #plot3(gx,horizon=ob)
def goDisplay():
  '''
  gx = readImage(gxfile)
  fx = gain(gx)
  writeImage(fxfile,fx)
  '''
  fx = readImage(gxfile)
  fx = copy(n1,3675,550,0,1100,60,fx)
  #fx = slog(fx)
  fx = gain(fx)
  writeImage("fx",fx)
  #fx = div(fx,600000)
  fxmin = min(fx)
  fxmax = max(fx)
  print fxmin
  print fxmax
  #plot3(fx,cmin=-5,cmax=5)
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
  '''
  plot3(gx,p2, cmin=-1,cmax=1,cmap=bwrNotch(1.0),
        clab="Inline slope (sample/sample)",png="p2")
  plot3(gx,p3, cmin=-1,cmax=1,cmap=bwrNotch(1.0),
        clab="Crossline slope (sample/sample)",png="p3")
  plot3(gx,sub(1,ep),cmin=0,cmax=1,cmap=jetRamp(1.0),
        clab="Planarity")
  '''

def goScan():
  print "goScan ..."
  if not plotOnly:
    gx = readImage(gxfile)
    fs = FaultScanner(sigmaPhi,sigmaTheta)
    sig1,sig2,smooth=4.0,2.0,4.0
    fl,fp,ft = fs.scan(minPhi,maxPhi,minTheta,maxTheta,sig1,sig2,smooth,gx)
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
  plot3(gx,ft,cmin=60,cmax=85,cmap=jetFill(1.0),
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
  '''
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  plot3(gx,ftt,cmin=60,cmax=85,cmap=jetFillExceptMin(1.0),
        clab="Fault dip (degrees)",png="ftt")
  plot3(gx,fpt,cmin=0,cmax=360,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=45,png="fpt")
  '''

def goReskin():
  print "go skin..."
  #gx = readImage(gxfile)
  if not plotOnly:
    fl = readImage(flfile)
    fp = readImage(fpfile)
    ft = readImage(ftfile)
    fsk = FaultSkinner()
    fsk.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
    fsk.setMaxDeltaStrike(10)
    fsk.setMaxPlanarDistance(0.2)
    fsk.setMinSkinSize(minSkinSize)
    cells = fsk.findCells([fl,fp,ft])
    sks = fsk.findSkins(cells)
    print len(sks)
    print "fault skins load finish..."
    fcs = FaultSkin.getCells(sks)
    cells = []
    for ic in range(0,len(fcs),4):
      cells.append(fcs[ic])
    print len(cells)
    print "fault cells load finish..."
    fs = FaultScanner(sigmaPhi,sigmaTheta)
    sp = fs.makePhiSampling(minPhi,maxPhi)
    dp = sp.getDelta()
    fr = FaultReconstructor(n1,n2,n3,cells)
    skins = fr.reskin(minSkinSize,dp)
    writeSkins(fsktv,skins)
    fd = FaultDisplay()
    print "fault skins load finish..."
    fd = FaultDisplay()
    fd.getFlt(skins,fl)
    fd.getFpt(skins,fp)
    fd.getFtt(skins,ft)
    writeImage(fltfile,fl)
    writeImage(fptfile,fp)
    writeImage(fttfile,ft)
  else:
    skins = readSkins(fsktv)
  '''
  print len(skins)
  fd = FaultDisplay()
  cells = FaultSkin.getCells(skins)
  flt = fillfloat(-0.001,n1,n2,n3)
  fd.getFlImage(cells,flt)
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),clab="Fault likelihood",png="smt")
  plot3(gx,skins=skins,png="skinsTv")
  '''


def goSkinTv():
  print "go skin..."
  #gx = readImage(gxfile)
  if not plotOnly:
    fl = readImage(flfile)
    fp = readImage(fpfile)
    ft = readImage(ftfile)
    fsk = FaultSkinner()
    fsk.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
    fsk.setMaxDeltaStrike(10)
    fsk.setMaxPlanarDistance(0.2)
    fsk.setMinSkinSize(minSkinSize)
    cells = fsk.findCells([fl,fp,ft])
    sks = fsk.findSkins(cells)
    print len(sks)
    print "fault skins load finish..."
    fcs = FaultSkin.getCells(sks)
    cells = []
    for ic in range(0,len(fcs),6):
      cells.append(fcs[ic])
    print len(cells)
    print "fault cells load finish..."

    fs = FaultScanner(sigmaPhi,sigmaTheta)
    sp = fs.makePhiSampling(minPhi,maxPhi)
    st = fs.makeThetaSampling(minTheta,maxTheta)
    fsx = FaultSkinnerX()
    fsx.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
    fsx.setMinSkinSize(minSkinSize)
    fsx.setMaxPlanarDistance(0.2)

    fsx.resetCells(cells)
    fsx.setGaussWeights(sp,st)
    skins = fsx.findSkins(n1,n2,n3,cells)
    removeAllSkinFiles(fsktv)
    writeSkins(fsktv,skins)
    fd = FaultDisplay()
    print "fault skins load finish..."
    fd = FaultDisplay()
    fd.getFlt(skins,fl)
    fd.getFpt(skins,fp)
    fd.getFtt(skins,ft)
    writeImage(fltfile,fl)
    writeImage(fptfile,fp)
    writeImage(fttfile,ft)
  else:
    skins = readSkins(fsktv)
  '''
  print len(skins)
  fd = FaultDisplay()
  cells = FaultSkin.getCells(skins)
  flt = fillfloat(-0.001,n1,n2,n3)
  fd.getFlImage(cells,flt)
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),clab="Fault likelihood",png="smt")
  plot3(gx,skins=skins,png="skinsTv")
  '''

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
    print "total number of cells =",len(cells)
    print "total number of skins =",len(skins)
    print "number of cells in skins =",FaultSkin.countCells(skins)
    fd = FaultDisplay()
    print "fault skins load finish..."
    fpt = fillfloat(-0.001,n1,n2,n3)
    fd = FaultDisplay()
    fd.getFpt(skins,gx,fpt)
    writeImage("fpk",fpt)
    removeAllSkinFiles(fskbase)
    writeSkins(fskbase,skins)
  else:
    skins = readSkins(fskbase)
  '''
  fd = FaultDisplay()
  sk = fd.getLargeFaults(20000,skins)
  print len(sk)
  plot3(gx,skins=sk)
  '''

def goFaultImages():
  gx = readImage(gxfile)
  if not plotOnly:
    fd = FaultDisplay()
    #skins = readSkins(fsktv)
    skins = readSkins(fskbase)
    print "fault skins load finish..."
    #flt = fillfloat(-0.001,n1,n2,n3)
    fpt = fillfloat(-0.001,n1,n2,n3)
    #ftt = fillfloat(-0.001,n1,n2,n3)
    fd = FaultDisplay()
    #fd.getFlt(skins,gx,flt)
    fd.getFpt(skins,gx,fpt)
    #fd.getFtt(skins,gx,ftt)
    #writeImage(fltfile,flt)
    writeImage("fpk",fpt)
    #writeImage(fttfile,ftt)
  else:
    #flt = readImage(fltfile)
    #fpt = readImage(fptfile)
    fpt = readImage("fpk")
    #ftt = readImage(fttfile)
  '''
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  plot3(gx,ftt,cmin=65,cmax=85,cmap=jetFillExceptMin(1.0),
        clab="Fault dip (degrees)",png="ftt")
  plot3(gx,fpt,cmin=0,cmax=180,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=10,png="fpt")
  '''
def goSmooth():
  print "goSmooth ..."
  flstop = 0.1
  fsigma = 8.0
  fl = readImage(flfile)
  gx = readImage(gxfile)
  skins = readSkins(fskgood)
  flt = zerofloat(n1,n2,n3)
  fsx = FaultSkinnerX()
  fsx.getFl(skins,flt)
  p2,p3,ep = FaultScanner.slopes(8.0,1.0,1.0,5.0,gx)
  gsx = FaultScanner.smooth(flstop,fsigma,p2,p3,flt,gx)
  writeImage(p2file,p2)
  writeImage(p3file,p3)
  writeImage(epfile,ep)
  writeImage(gsxfile,gsx)
  '''
  plot3(gx,flt,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
        clab="Fault likelihood",png="fli")
  plot3(gsx,png="gsx")
  '''

def goSlip():
  print "goSlip ..."
  gx = readImage(gxfile)
  if not plotOnly:
    gsx = readImage(gsxfile)
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    skins = readSkins(fskgood)
    fsl = FaultSlipper(gsx,p2,p3)
    fsl.setOffset(3.0) # the default is 2.0 samples
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
    '''
    smark = -999.999
    s1,s2,s3 = fsl.getDipSlips(skins,smark)
    s1,s2,s3 = fsl.interpolateDipSlips([s1,s2,s3],smark)
    gw = fsl.unfault([s1,s2,s3],gx)
    writeImage(gwfile,gw)
    writeImage(fs1file,s1)
    writeImage(fs2file,s2)
    writeImage(fs3file,s3)
    '''
  else:
    gw = readImage(gwfile)
    #s1 = readImage(fs1file)
    #skins = readSkins(fslbase)
  '''
  plot3(gx,skins=skins,smax=10.0,png="skinss1")
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
    fw = readImage(fwsfile)
    gw = readImage(gwfile)
  plot3(gx,png="gxuf")
  plot3(fw,png="fwuf")
  plot3(gw,png="fwuf")
  '''
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
          links=False,curve=False,trace=False,png=None):
  n3 = len(f)
  n2 = len(f[0])
  n1 = len(f[0][0])
  s1,s2,s3=Sampling(n1),Sampling(n2),Sampling(n3)
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  sf = SimpleFrame(AxesOrientation.XRIGHT_YIN_ZDOWN)
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
    cmap = ColorMap(0.0,1.0,ColorMap.JET)
    xyz,uvw,rgb = FaultCell.getXyzUvwRgbForLikelihood(0.7,cmap,cells,False)
    qg = QuadGroup(xyz,uvw,rgb)
    qg.setStates(ss)
    sf.world.addChild(qg)
  if horizon:
    tg = TriangleGroup(True,s3,s2,horizon)
    tg.setColor(Color.CYAN)
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
        cmap = ColorMap(0.2,0.8,ColorMap.JET)
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
  ipg.setSlices(150,5,56)
  #ipg.setSlices(85,5,43)
  #ipg.setSlices(85,5,102)
  #ipg.setSlices(n1,0,n3) # use only for subset plots
  ipg.setSlices(n1,376,308)
  if cbar:
    sf.setSize(1037,700)
  else:
    sf.setSize(900,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  zscale = 0.5*max(n2*d2,n3*d3)/(n1*d1)
  #zscale = 1.5*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.5)
  #ov.setScale(2.5)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(0.0,-0.00,-0.05))
  ov.setAzimuthAndElevation(45.0,35.0)
  #ov.setAzimuthAndElevation(-55.0,35.0)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")


#############################################################################
run(main)
