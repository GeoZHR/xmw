"""
Demonstrate 3D seismic image processing for faults
Author: Xinming Wu, Colorado School of Mines
Version: 2016.01.22
"""

from utils import *
#setupForSubset("nathan")
setupForSubset("nathanTest")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count
# Names and descriptions of image files used below.
gxfile  = "gx" # input image (maybe after bilateral filtering)
gxpfile  = "gxp" # input image (maybe after bilateral filtering)
epfile  = "ep" # fault likelihood
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
fttvfile = "fttv" # fault dip thinned
fskbase = "fsk" # fault skin (basename only)
fslbase = "fsl" # fault skin (basename only)
fskgood = "fsg" # fault skin (basename only)
fsktv = "fst" # fault skin (basename only)

# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
minPhi,maxPhi = 0,360
minTheta,maxTheta = 65,85
sigmaPhi,sigmaTheta = 20,50

# These parameters control the construction of fault skins.
# See the class FaultSkinner for more information.
lowerLikelihood = 0.3
upperLikelihood = 0.5
minSkinSize = 200

# These parameters control the computation of fault dip slips.
# See the class FaultSlipper for more information.
minThrow = 0.0
maxThrow = 85.0

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
#pngDir = "../../../png/beg/hongliu/"
pngDir = "../../../png/beg/nathan/sub8/"
pngDir = None
plotOnly = False

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goMask()
  #goPlanar()
  #goFaultScan()
  #goThin()
  #goSkin()
  goSkinTv()
  #goFaultImages()
  #goSurfaces()
  #goFaultPoints()
  #getOceanBottom()
  #goSeisResample()
  #goHorizon()


def goMask():
  f1 = s1.getFirst()
  d1 = s1.getDelta()
  gx = readImage(gxfile)
  tp = readImage2D(n2,n3,"ob")
  bt = readImage2D(n2,n3,"hl1")
  bt = sub(bt,f1)
  bt = div(bt,d1)
  hp = Helper()
  hp.padValues(tp,bt,gx)
  plot3(gx)
  writeImage(gxpfile,gx)

def goPlanar():
  gx = readImage(gxfile)
  lof = LocalOrientFilter(8,2)
  ets = lof.applyForTensors(gx)
  ets.setEigenvalues(1.0,0.01,0.1)
  fer = FaultEnhancer(sigmaPhi,sigmaTheta)
  ep = fer.applyForPlanar(20,ets,gx)
  writeImage("ep",ep)
  print min(ep)
  print max(ep)
  #plot3(gx,cmin=-3,cmax=3)
  #plot3(ep,cmin=0.1,cmax=0.9)

def goFaultScan():
  ep = readImage("ep")
  ep = clip(0.0,1.0,ep)
  gx = readImage(gxfile)
  if not plotOnly:
    fe = FaultEnhancer(sigmaPhi,sigmaTheta)
    flpt = fe.scan(minPhi,maxPhi,minTheta,maxTheta,ep)
    writeImage(flfile,flpt[0])
    writeImage(fpfile,flpt[1])
    writeImage(ftfile,flpt[2])
  else:
    fl = readImage(flfile)
  plot3(gx)
  plot3(ep,cmin=0.1,cmax=0.9)
  plot3(gx,fl,cmin=0.2,cmax=1.0,cmap=jetFillExceptMin(1.0))


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
  plot3(gx,fpt,cmin=0,cmax=360,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=45,png="fpt")

def goSurfaces():
  fn = "sm1"
  gx = readImage(gxfile)
  if not plotOnly:
    nl1 = 1646796
    nu1 = 1330218
    nm1 =  888165
    sp = readImage2D(nl1,3,"hzs/L1")
    #sp = readImage2D(nu1,3,"hzs/U1")
    #sp = readImage2D(nm1,3,"hzs/M1")
    print min(sp[0])
    print max(sp[0])
    print min(sp[2])
    print max(sp[2])
    hp = Helper()
    ndfs = zerofloat(3,2)
    sf = hp.surfaceResample(s2,s3,1.5,sp,ndfs)
    sy = Sampling(round(ndfs[0][0]),ndfs[0][1],ndfs[0][2])
    sx = Sampling(round(ndfs[1][0]),ndfs[1][1],ndfs[1][2])
    writeImage(fn,sf)
    writeImage(fn+"ndfs",ndfs)
  else:
    ndfs = readImage2D(3,2,fn+"ndfs")
    ny = round(ndfs[0][0])
    nx = round(ndfs[1][0])
    sf = readImage2D(ny,nx,fn)
    sy = Sampling(round(ndfs[0][0]),ndfs[0][1],ndfs[0][2])
    sx = Sampling(round(ndfs[1][0]),ndfs[1][1],ndfs[1][2])
  plot3(gx,cmin=-3,cmax=3,sx=sx,sy=sy,horizon=sf)

def getOceanBottom():
  hp = Helper()
  gx = readImage(gxfile)
  ob = hp.getOceanBottom(0.1,gx)
  writeImage("ob",ob)
  #ob = readImage2D(n2,n3,"ob")
  #plot3(gx,horizon=ob)


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
    for ic in range(0,len(fcs),2):
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
    fsx.setGrowLikelihoods(0.3,upperLikelihood)
    #fsx.setMinSkinSize(minSkinSize)
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
    fp = zerofloat(n1,n2,n3)
    skins = readSkins(fsktv)
    fsx = FaultSkinnerX()
    fsx.getFp(4000,skins,fp)

    #fl = readImage("fltv")
    #fp = readImage("fptv")
    #ft = readImage("fttv")
  plot3(gx,fp,cmin=0,cmax=180,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=10,png="fpt")

  '''
  plot3(gx,fl,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  plot3(gx,ft,cmin=60,cmax=85,cmap=jetFillExceptMin(1.0),
        clab="Fault dip (degrees)",png="ftt")
  plot3(gx,fp,cmin=0,cmax=180,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=45,png="fpt")
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
    fd.getFlt(skins,fl)
    #fd.getFpt(skins,fp)
    #fd.getFtt(skins,ft)
    writeImage("flk",fl)
    removeAllSkinFiles(fskbase)
    writeSkins(fskbase,skins)
  else:
    #skins = readSkins(fskbase)
    fl  = readImage("flk")
  flt = readImage(fltfile)
  plot3(gx,fl,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")

def goFaultImages():
  gx = readImage(gxfile)
  if not plotOnly:
    fl = fillfloat(-0.001,n1,n2,n3)
    fp = fillfloat(-0.001,n1,n2,n3)
    ft = fillfloat(-0.001,n1,n2,n3)
    skins = readSkins(fsktv)
    fp = zerofloat(n1,n2,n3)
    skins = readSkins(fsktv)
    fsx = FaultSkinnerX()
    fsx.getFlpt(4000,skins,fl,fp,ft)
    writeImage(fltvfile,fl)
    writeImage(fptvfile,fp)
    writeImage(fttvfile,ft)
  else:
    fl = readImage(fltvfile)
    fp = readImage(fptvfile)
    ft = readImage(fttvfile)
  plot3(gx,fl,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  plot3(gx,ft,cmin=65,cmax=85,cmap=jetFillExceptMin(1.0),
        clab="Fault dip (degrees)",png="ftt")
  plot3(gx,fp,cmin=0,cmax=180,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=10,png="fpt")

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
    ts = hp.horizonWithFaultDensity(n1-2,[0.0,0.15],horizon,fd)
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


#############################################################################
run(main)
