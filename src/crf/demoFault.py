"""
Demonstrate 3D seismic image processing for faults
Author: Xinming Wu, Colorado School of Mines
Version: 2016.01.22
"""

from utils import *
#setupForSubset("nathan")
setupForSubset("nathanSub8")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count
# Names and descriptions of image files used below.
gxfile  = "gx" # input image (maybe after bilateral filtering)
gsxfile = "gsx"
p2file = "p2"
p3file = "p3"
gxpfile  = "gxp" # input image (maybe after bilateral filtering)
epfile  = "ep" # fault likelihood
eppfile  = "epp" # fault likelihood
flfile  = "fl" # fault likelihood
fpfile  = "fp" # fault strike (phi)
ftfile  = "ft" # fault dip (theta)
flvfile  = "flv" # fault likelihood
fpvfile  = "fpv" # fault strike (phi)
ftvfile  = "ftv" # fault dip (theta)
fltfile = "flt" # fault likelihood thinned
fptfile = "fpt" # fault strike thinned
fttfile = "ftt" # fault dip thinned
fltvfile = "fltv" # fault likelihood thinned
fptvfile = "fptv" # fault strike thinned
fptmfile = "fptm" # fault strike masked out 20~40 degrees
fttvfile = "fttv" # fault dip thinned
fskbase = "fsk" # fault skin (basename only)
fslbase = "fsl" # fault skin (basename only)
fskgood = "fsg" # fault skin (basename only)
fs1file = "fs1"
fw1file = "fw1"
fw2file = "fw2"
fw3file = "fw3"
gwfile = "gw"
fsktv = "fst" # fault skin (basename only)
fppfile = "fpp"
hl1file = "hl1"
hu1file = "hu1"
hm1file = "hm1"

# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
minPhi,maxPhi = 0,360
minTheta,maxTheta = 65,85
sigmaPhi,sigmaTheta = 20,50

# These parameters control the construction of fault skins.
# See the class FaultSkinner for more information.
lowerLikelihood = 0.6
upperLikelihood = 0.7
minSkinSize = 200

# These parameters control the computation of fault dip slips.
# See the class FaultSlipper for more information.
minThrow = 0.0
maxThrow = 20.0

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
#pngDir = "../../../png/beg/hongliu/"
pngDir = "../../../png/beg/nathan/sub8/"
pngDir = None
plotOnly = False

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goPlanar()
  #goMask()
  #goFaultScan()
  #goThin()
  #goSkinTv()
  #goSmooth()
  goSlip()
  #goFaultImages()
  #goSurfaces()
  #goFaultPoints()
  #getOceanBottom()
  #goSeisResample()
  #goHorizon()
  #goRosePlots()
  #goRosePlotsN()
  #goResetSurfaces()
  #goFaultsAndSurfs()
  #goFaultDensity()
  #goSetFaultImages()
  #goStrikeMask()


def goStrikeMask():
  gx = readImage(gxfile)
  fp = readImage(fptvfile)
  hp = Helper()
  hp.strikeMask(15,40,fp)
  writeImage(fptmfile,fp)
  plot3(gx,fp,cmin=0,cmax=180,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=10,png="fpt")
def goFaultsAndSurfs():
  #fpt = readImage(fptvfile)
  hu1 = readImage2D(n2,n3,hu1file)
  hm1 = readImage2D(n2,n3,hm1file)
  hl1 = readImage2D(n2,n3,hl1file)
  hpr = Helper()
  spm = readImage2D(888165,3,"hzs/M1")
  hpr.mergeU1AndTop(round(s2.getFirst()),spm,hm1)
  writeImage(hm1file,hm1)
  gx = readImage(gxfile)
  flt = readImage(fltvfile)
  hpr.horizonToImage(div(hu1,5),flt)
  hpr.horizonToImage(div(hm1,5),flt)
  hpr.horizonToImage(div(hl1,5),flt)
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  '''
  plot3(gx,fpt,cmin=0,cmax=180,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=10,png="fpt")
  '''

def goFaultDensity():
  stfile = hm1file
  sbfile = hl1file
  st = readImage2D(n2,n3,stfile)
  sb = readImage2D(n2,n3,sbfile)
  if not plotOnly:
    fp = readImage(fptvfile)
    st = div(st,5)
    sb = div(sb,5)
    hp = Helper()
    fd = hp.faultDensity(st,sb,fp)
    writeImage("fd"+stfile+sbfile,fd)
  else:
    fd = readImage2D(n2,n3,"fd"+stfile+sbfile)
    print min(fd)
    print max(fd)
    rgf = RecursiveGaussianFilterP(20)
    rgf.apply00(fd,fd)
  #gx = readImage(gxfile)
  #plot3(gx,horizon=sb,fd=fd)
  hp = Helper()
  fd = hp.transpose(fd)
  c3 = Sampling(n3,1.0/1.5,s3.getFirst())
  plot2(c3,s2,fd,cmin=0.01,cmax=0.15,png="fd"+stfile+sbfile)


def goRosePlots():
  tpfile = hu1file
  btfile = hm1file
  tp = readImage2D(n2,n3,tpfile)
  bt = readImage2D(n2,n3,btfile)
  tp = div(tp,5)
  bt = div(bt,5)
  fpp = readImage2D(65230747,4,fppfile)
  c2,c3=20,4
  rp = RosePlot()
  title = tpfile+'~'+btfile
  #pp = rp.applyForRosePlotsX(npm,tp,bt,c2,c3,n2,n3,36,fc,ob)
  pp = rp.applyForRosePlotsN(64,tp,bt,c2,c3,n2,n3,36,fpp)
  #pp = rp.applyForRosePlotsN(99,99,c2,c3,n2,n3,18,fc,ob)
  pp.addTitle(title)
  pf = PlotFrame(pp)
  #wx,wy = 1450,round((c3*1450)/c2)+50
  #wx,wy = 1700,round((c3*1750)/c2)+50
  wx,wy = 2100,round((c3*2100)/c2)+100
  wx = round(wx*1.5)
  wy = round(wy*1.3)
  pf.setSize(wx,wy)
  pf.setVisible(True)
  pf.paintToPng(720,6,pngDir+title+".png")

def goRosePlotsN():
  hu = readImage2D(n2,n3,hu1file)
  hm = readImage2D(n2,n3,hm1file)
  hl = readImage2D(n2,n3,hl1file)
  fp = readImage2D(65230747,4,fppfile)
  hu = div(hu,5)
  hm = div(hm,5)
  hl = div(hl,5)
  tpfile = hm1file
  btfile = hl1file
  tp = hm
  bt = hl
  c2,c3=20,4
  rp = RosePlot()
  #npm = rp.findMaxSamples(c2,c3,hu,hm,hl,fp)
  npm = 300000
  title = tpfile+'~'+btfile
  pp = rp.applyForRosePlotsX(64,npm,tp,bt,c2,c3,36,fp)
  pp.addTitle(title)
  pf = PlotFrame(pp)
  wx,wy = 2100,round((c3*2100)/c2)+100
  wx = round(wx*1.5)
  wy = round(wy*1.3)
  pf.setSize(wx,wy)
  pf.setVisible(True)
  pf.paintToPng(720,6,pngDir+title+"N"+".png")

def goPlanar():
  gx = readImage(gxfile)
  if not plotOnly:
    lof = LocalOrientFilter(8,2)
    ets = lof.applyForTensors(gx)
    ets.setEigenvalues(1.0,0.01,0.1)
    fer = FaultEnhancer(sigmaPhi,sigmaTheta)
    ep = fer.applyForPlanar(20,ets,gx)
    writeImage(epfile,ep)
    print min(ep)
    print max(ep)
  else:
    ep = readImage(epfile)
  #plot3(gx,cmin=-3,cmax=3)
  plot3(ep,cmin=0.1,cmax=0.9)

def goMask():
  f1 = s1.getFirst()
  d1 = s1.getDelta()
  ep = readImage(epfile)
  plot3(ep,cmin=0.1,cmax=0.9)
  tp = readImage2D(n2,n3,"ob")
  bt = readImage2D(n2,n3,"hl1")
  bt = sub(bt,f1)
  bt = div(bt,d1)
  hp = Helper()
  hp.padValues(tp,bt,ep)
  writeImage(eppfile,ep)
  plot3(ep,cmin=0.1,cmax=0.9)

def goFaultScan():
  ep = readImage(eppfile)
  ep = clip(0.0,1.0,ep)
  if not plotOnly:
    fe = FaultEnhancer(sigmaPhi,sigmaTheta)
    flpt = fe.scan(minPhi,maxPhi,minTheta,maxTheta,ep)
    writeImage(flfile,flpt[0])
    writeImage(fpfile,flpt[1])
    writeImage(ftfile,flpt[2])
  else:
    fl = readImage(flfile)
  #gx = readImage(gxfile)
  #plot3(gx,flpt[0],cmin=0.01,cmax=1,cmap=jetRamp(1.0))
  '''
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  flt,fpt,ftt=fe.thin(flpt)
  writeImage(fltfile,flt)
  writeImage(fptfile,ftt)
  writeImage(fttfile,ftt)
  plot3(gx)
  plot3(ep,cmin=0.1,cmax=0.9)
  plot3(gx,flt,cmin=0.2,cmax=1.0,cmap=jetFillExceptMin(1.0))
  '''

def goSurfaces():
  gx = readImage(gxfile)
  fl = readImage(fltvfile)
  hp = Helper()
  fns = ["sm1","su1"]
  for fni in fns:
    ndfs = readImage2D(3,2,fni+"ndfs")
    ny = round(ndfs[0][0])
    nx = round(ndfs[1][0])
    sf = readImage2D(ny,nx,fni)
    b2 = round(ndfs[0][2]-s2.getFirst())
    b3 = round(ndfs[1][2]-s3.getFirst())
    hp.horizonToImage(s1,b2,b3,sf,fl)
  sl1 = readImage2D(n2,n3,hl1file)
  hp.horizonToImage(s1,0,0,sl1,fl)
  plot3(gx,fl,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")

def goFaultPoints():
  fp = readImage(fptmfile)
  rp = RosePlot()
  ps = rp.faultPoints(fp)
  print len(ps)
  print len(ps[0])
  writeImage(fppfile,ps)

def getOceanBottom():
  hp = Helper()
  gx = readImage(gxfile)
  ob = hp.getOceanBottom(0.1,gx)
  writeImage("ob",ob)
  #ob = readImage2D(n2,n3,"ob")
  #plot3(gx,horizon=ob)

def goScan():
  print "goScan ..."
  if not plotOnly:
    gx = readImage(gxpfile)
    fs = FaultScanner(sigmaPhi,sigmaTheta)
    sig1,sig2,smooth=16.0,1.0,4.0
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
  plot3(gx,fl,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
      clab="Fault likelihood",png="fl")
  plot3(gx,fp,cmin=0,cmax=360,cmap=hueFill(1.0),
      clab="Fault strike (degrees)",cint=45,png="fp")
  plot3(gx,ft,cmin=60,cmax=85,cmap=jetFill(1.0),
      clab="Fault dip (degrees)",png="ft")

def goThin():
  print "goThin ..."
  gx = readImage(gxfile)
  if not plotOnly:
    fl = readImage(flfile)
    fp = readImage(fpfile)
    ft = readImage(ftfile)
    flt,fpt,ftt = FaultEnhancer.thin([fl,fp,ft])
    writeImage(fltfile,flt)
    writeImage(fptfile,fpt)
    writeImage(fttfile,ftt)
  else:
    flt = readImage(fltfile)
    fpt = readImage(fptfile)
    ftt = readImage(fttfile)
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  plot3(gx,ftt,cmin=60,cmax=85,cmap=jetFillExceptMin(1.0),
        clab="Fault dip (degrees)",png="ftt")
  plot3(gx,fpt,cmin=0,cmax=180,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=45,png="fpt")

def goSkinTv():
  print "go skin..."
  gx = readImage(gxfile)
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
    st = fs.makeThetaSampling(minTheta,maxTheta)
    fsx = FaultSkinnerX()
    fsx.setGrowLikelihoods(0.3,0.5)
    fsx.setMinSkinSize(1000)
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
    writeImage("fltv",fl)
    writeImage("fptv",fp)
    writeImage("fttv",ft)
  else:
    fl = readImage("fltv")
    fp = readImage("fptv")
    ft = readImage("fttv")
  plot3(gx,fp,cmin=0,cmax=180,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=10,png="fpt")

  plot3(gx,fl,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  plot3(gx,ft,cmin=60,cmax=85,cmap=jetFillExceptMin(1.0),
        clab="Fault dip (degrees)",png="ftt")

def goSmooth():
  print "goSmooth ..."
  flstop = 0.01
  fsigma = 8.0
  gx = readImage(gxfile)
  if not plotOnly:
    skins = readSkins(fsktv)
    flt = zerofloat(n1,n2,n3)
    fsx = FaultSkinnerX()
    fsx.getFl(100,skins,flt)
    p2,p3,ep = FaultScanner.slopes(8.0,2.0,2.0,5.0,gx)
    gsx = FaultScanner.smooth(flstop,fsigma,p2,p3,flt,gx)
    writeImage(p2file,p2)
    writeImage(p3file,p3)
    writeImage(gsxfile,gsx)
  else:
    gsx = readImage(gsxfile)
  plot3(gsx,png="gsx")

def goSlip():
  print "goSlip ..."
  #gx = readImage(gsxfile)
  if not plotOnly:
    gsx = readImage(gsxfile)
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    skins = readSkins(fsktv)
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
    s1,s2,s3 = fsl.getDipSlips(skins,smark)
    writeImage(fs1file,s1)
    '''
    w1,w2,w3 = fsl.interpolateDipSlips([s1,s2,s3],smark)
    gw = fsl.unfault([w1,w2,w3],gx)
    writeImage(fw1file,w1)
    writeImage(fw2file,w2)
    writeImage(fw3file,w3)
    writeImage(gwfile,gw)
    '''
  else:
    s1 = readImage(fs1file)
    #w1 = readImage(fw1file)
    #w2 = readImage(fw2file)
    #w3 = readImage(fw3file)
    #gw = readImage(gwfile)

  plot3(gsx,s1,cmin=0,cmax=20.0,cmap=jetFillExceptMin(1.0),
        clab="Fault throw (samples)",png="gxs1")
  '''
  plot3(gx,s1,cmin=0.0,cmax=30.0,cmap=jetFill(0.3),
        clab="Vertical shift (samples)",png="gxs1i")
  plot3(gx,s2,cmin=-2.0,cmax=2.0,cmap=jetFill(0.3),
        clab="Inline shift (samples)",png="gxs2i")
  plot3(gx,s3,cmin=-1.0,cmax=1.0,cmap=jetFill(0.3),
        clab="Crossline shift (samples)",png="gxs3i")
  plot3(gx)
  plot3(gw,clab="Amplitude",png="gw")
  '''


def goFaultImages():
  gx = readImage(gxfile)
  if not plotOnly:
    fl = fillfloat(-0.01,n1,n2,n3)
    fp = fillfloat(-0.01,n1,n2,n3)
    ft = fillfloat(-0.01,n1,n2,n3)
    skins = readSkins(fsktv)
    fsx = FaultSkinnerX()
    fsx.getFlpt(3000,skins,fl,fp,ft)
    hp = Helper()
    hp.rotateX(26,fp)
    hp.rotate(90,fp)
    hp.convert(fp)
    writeImage(fltvfile,fl)
    writeImage(fptvfile,fp)
    writeImage(fttvfile,ft)
  else:
    #fl = readImage(fltvfile)
    #ft = readImage(fttvfile)
    fp = readImage(fptvfile)
  '''
  plot3(gx,fl,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  plot3(gx,ft,cmin=65,cmax=85,cmap=jetFillExceptMin(1.0),
        clab="Fault dip (degrees)",png="ftt")
  '''
  plot3(gx,fp,cmin=0,cmax=180,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=10,png="fpt")

def goResetSurfaces():
  gx = readImage(gxfile)
  fl = readImage(fltvfile)
  fns = ["m1","u1"]
  hp = Helper()
  for fni in fns:
    ndfs = readImage2D(3,2,"s"+fni+"ndfs")
    ny = round(ndfs[0][0])
    nx = round(ndfs[1][0])
    sf = readImage2D(ny,nx,"s"+fni)
    b2 = round(ndfs[0][2]-s2.getFirst())
    b3 = round(ndfs[1][2]-s3.getFirst())
    hz = zerofloat(n2,n3)
    for ix in range(nx):
      for iy in range(ny):
        hz[ix+b3][iy+b2] = sf[ix][iy]
    writeImage("h"+fni,hz)
    hp.horizonToImage(s1,0,0,hz,fl)
  plot3(gx,fl,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")

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

def plot3(f,g=None,cmin=-2,cmax=2,cmap=None,clab=None,cint=None,
          sx=None,sy=None,horizon=None,fd=None,cells=None,skins=None,png=None):
  n3 = len(f)
  n2 = len(f[0])
  n1 = len(f[0][0])
  #s1,s2,s3=Sampling(n1),Sampling(n2),Sampling(n3)
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
  if horizon and not fd:
    hp = Helper()
    print "fvalues"
    print min(f)
    print max(f)
    ts = hp.horizonWithAmplitude(s1,s2,s3,s1,sy,sx,[cmin,cmax],horizon,f)
    tg = TriangleGroup(True,ts[0],ts[1])
    sf.world.addChild(tg)
  if horizon and fd:
    hp = Helper()
    ts = hp.horizonWithFaultDensity(s1,s2,s3,[0.0,0.15],horizon,fd)
    tg = TriangleGroup(True,ts[0],ts[1])
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
    ms.setEmissiveBack(Color(0.0,0.0,0.5))
    ss.add(ms)
    sg.setStates(ss)
    size = 2.0
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


def plot2(s1,s2,f,cmap=None,cmin=None,cmax=None,cint=None,clab=None,png=None): 
  f1 = s1.getFirst()
  f2 = s2.getFirst()
  d1 = s1.getDelta()
  d2 = s2.getDelta()
  n1 = s1.getCount()
  orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT;
  panel = PlotPanel(1,1,orientation)#,PlotPanel.AxesPlacement.NONE)
  #panel.setVInterval(0.1)
  #panel.setHInterval(1.0)
  panel.setHLabel("Inline")
  panel.setVLabel("Crossline")
  pxv = panel.addPixels(0,0,s1,s2,f);
  pxv.setColorModel(ColorMap.JET)
  #pxv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmin and cmax:
    pxv.setClips(cmin,cmax)
  cb = panel.addColorBar();
  if cint:
    cb.setInterval(cint)
  if clab:
    cb.setLabel(clab)
  panel.setColorBarWidthMinimum(50)
  moc = panel.getMosaic();
  frame = PlotFrame(panel);
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  #frame.setTitle("normal vectors")
  frame.setVisible(True);
  wx,wy = 2100,2100/5+100+50
  frame.setSize(wx,wy) 
  frame.setFontSize(13)
  if pngDir and png:
    frame.paintToPng(720,3.333,pngDir+png+".png")


#############################################################################
run(main)
