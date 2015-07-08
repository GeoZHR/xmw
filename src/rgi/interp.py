"""
Interpolates scattered data, such as data from well logs.
"""
from tputils import *

#setupForSubset("subz_51_4_1400")
#setupForSubset("subz_401_4_400")
setupForSubset("subz_401_4_600")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.getCount(),s2.getCount(),s3.getCount()
method = "b" # blended
logSet = "d" # deep logs only
logType = "v"; logLabel = "Velocity (km/s)"; vmin,vmax = 2.4,5.6
logType = "d"; logLabel = "Density (g/cc)"; vmin,vmax = 2.0,2.8
#logType = "p"; logLabel = "Porosity"; vmin,vmax = 0.0,0.4
#logType = "g"; logLabel = "Gamma ray (API units)"; vmin,vmax = 0.0,200.0
smin,smax = -5.5,5.5
#smooth = 50 # half-width of smoothing filter for logs
smooth = 0 # no smoothing

sfile = "tpsz" # seismic image
efile = "tpet" # eigen-tensors (structure tensors)
esfile = "tpets" # eigen-tensors scaled by semblances
s1file = "tps1" # semblance w,uv
s2file = "tps2" # semblance vw,u
s3file = "tps3" # semblance uvw,
gfile = "tpg"+logType # simple gridding with null for unknown samples
pfile = "tpp"+logType+method # values of nearest known samples
qfile = "tpq"+logType+method # output of blended gridder
tfile = "tpt"+logType+method # times to nearest known samples

gxfile = "tpsz" # seismic image
gsxfile = "gsx" # image after lsf with fault likelihoods
epfile  = "ep" # eigenvalue-derived planarity
p2file  = "p2" # inline slopes
p3file  = "p3" # crossline slopes
flfile  = "fl" # fault likelihood
fpfile  = "fp" # fault strike (phi)
ftfile  = "ft" # fault dip (theta)
fltfile = "flt" # fault likelihood thinned
fptfile = "fpt" # fault strike thinned
fttfile = "ftt" # fault dip thinned
fs1file = "fs1" # fault slip (1st component)
fs2file = "fs2" # fault slip (2nd component)
fs3file = "fs3" # fault slip (3rd component)
fskbase = "fsk" # fault skin (basename only)
fslbase = "fsl" # fault skin (basename only)
fskgood = "fsg" # fault skin (basename only)
fqfile = "fq"+logType # output of blended gridder

horizons = ["CrowMountainCRMT"]
"""
horizons = [
  "CrowMountainCRMT",
  "TensleepASand",
  "TensleepBbaseC1Dolo"]
"""

pngDir = "../../../png/figi/"

#k1,k2,k3 = 228,170,74 # 2D displays
#k1,k2,k3 = 228,170,106 # 2D displays
#k1,k2,k3 = 228,170,96 # 2D displays
k1,k2,k3 = 366,15,96 # 3D displays


# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
minPhi,maxPhi = 0,360
minTheta,maxTheta = 80,85
sigmaPhi,sigmaTheta = 5,20

# These parameters control the construction of fault skins.
# See the class FaultSkinner for more information.
lowerLikelihood = 0.2
upperLikelihood = 0.5
minSkinSize = 3500

# These parameters control the computation of fault dip slips.
# See the class FaultSlipper for more information.
minThrow =  0.0
maxThrow = 20.0

plotOnly = False

def main(args):
  #goSlopes()
  #goScan()
  #goThin()
  #goSkin()
  #goReSkin()
  #goSmooth()
  #goSlip()
  #goInterp()
  #goInterpO()
  #goFigures()
  #goImpedance()
  goTest()

def goSlopes():
  print "goSlopes ..."
  gx = readImage(gxfile)
  sigma1,sigma2,sigma3,pmax = 16.0,1.0,1.0,5.0
  p2,p3,ep = FaultScanner.slopes(sigma1,sigma2,sigma3,pmax,gx)
  zm = ZeroMask(0.1,1,1,1,gx)
  zero,tiny=0.0,0.01
  zm.setValue(zero,p2)
  zm.setValue(zero,p3)
  zm.setValue(tiny,ep)
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
    gx = FaultScanner.taper(10,0,0,gx)
    fs = FaultScanner(sigmaPhi,sigmaTheta)
    fl,fp,ft = fs.scan(minPhi,maxPhi,minTheta,maxTheta,p2,p3,gx)
    zm = ZeroMask(0.1,1,1,1,gx)
    zm.setValue(0.0,fl)
    zm.setValue(0.0,fp)
    zm.setValue(0.0,ft)
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
  plot3(gx,convertDips(ft),cmin=25,cmax=55,cmap=jetFill(1.0),
      clab="Fault dip (degrees)",png="ft")

def goThin():
  print "goThin ..."
  gx = readImage(gxfile)
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  flt,fpt,ftt = FaultScanner.thin([fl,fp,ft])
  writeImage(fltfile,flt)
  writeImage(fptfile,fpt)
  writeImage(fttfile,ftt)
  plot3(gx,clab="Amplitude",png="gx")
  plot3(gx,fl,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
        clab="Fault likelihood",png="fl")
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  plot3(gx,fpt,cmin=0,cmax=360,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=45,png="fpt")
  plot3(gx,convertDips(ftt),cmin=15,cmax=55,cmap=jetFillExceptMin(1.0),
        clab="Fault dip (degrees)",png="ftt")

def goSkin():
  print "goSkin ..."
  gx = readImage(gxfile)
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  fs = FaultSkinner()
  fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  fs.setMaxDeltaStrike(10)
  fs.setMaxPlanarDistance(0.1)
  fs.setMinSkinSize(minSkinSize)
  cells = fs.findCells([fl,fp,ft])
  skins = fs.findSkins(cells)
  for skin in skins:
    skin.smoothCellNormals(4)
  print "total number of cells =",len(cells)
  print "total number of skins =",len(skins)
  print "number of cells in skins =",FaultSkin.countCells(skins)
  removeAllSkinFiles(fskbase)
  writeSkins(fskbase,skins)
  plot3(gx,cells=cells,png="cells")
  plot3(gx,skins=skins)

def goReSkin():
  print "goReSkin ..."
  useOldCells = True
  gx = readImage(gxfile)
  if not plotOnly:
    fl = readImage(flfile)
    sk = readSkins(fskbase)
    fsx = FaultSkinnerX()
    fsx.setParameters(20,10,0.5)
    fsx.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
    fsx.setMinSkinSize(3000)
    fsx.setMaxPlanarDistance(0.5)
    fsx.setSkinning(useOldCells)
    cells = FaultSkin.getCells(sk)
    fsx.resetCells(cells)
    skins = fsx.findSkinsXX(cells,fl)
    removeAllSkinFiles(fskgood)
    writeSkins(fskgood,skins)
  skins = readSkins(fskgood)
  for skin in skins:
    skin.smoothCellNormals(4)
  plot3(gx,skins=skins,png="skins")
  plot3(gx,skins=skins,links=True,png="links")
  for iskin,skin in enumerate(skins):
    plot3(gx,skins=[skin],links=True,png="skin"+str(iskin))

def goSmooth():
  print "goSmooth ..."
  flstop = 0.1
  fsigma = 8.0
  fl = readImage(flfile)
  gx = readImage(gxfile)
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  skins = readSkins(fskgood)
  flt = zerofloat(n1,n2,n3)
  fsx = FaultSkinnerX()
  fsx.getFl(skins,flt)
  p2,p3,ep = FaultScanner.slopes(8.0,2.0,2.0,5.0,gx)
  gsx = FaultScanner.smooth(flstop,fsigma,p2,p3,flt,gx)
  writeImage(p2file,p2)
  writeImage(p3file,p3)
  writeImage(gsxfile,gsx)
  plot3(gx)
  plot3(gsx)


def goSlip():
  print "goSlip ..."
  gx = readImage(gxfile)
  gsx = readImage(gsxfile)
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  skins = readSkins(fskgood)
  fsl = FaultSlipper(gsx,p2,p3)
  fsl.setOffset(2.0) # the default is 2.0 samples
  #fsl.setErrorPower(2.0)
  fsl.setZeroSlope(False) # True only if we want to show the error
  fsl.computeDipSlips(skins,minThrow,maxThrow)
  print "  dip slips computed, now reskinning ..."
  print "  number of skins before =",len(skins),
  fsk = FaultSkinner() # as in goSkin
  fsk.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  fsk.setMinSkinSize(minSkinSize)
  fsk.setMaxDeltaThrow(1.5)
  fsk.setMinMaxThrow(minThrow,maxThrow)
  skins = fsk.reskin(skins)

  print ", after =",len(skins)
  removeAllSkinFiles(fslbase)
  writeSkins(fslbase,skins)
  smark = -999.999
  s1,s2,s3 = fsl.getDipSlips(skins,smark)
  writeImage(fs1file,s1)
  writeImage(fs2file,s2)
  writeImage(fs3file,s3)
  print min(s1)
  print max(s1)
  print min(s2)
  print max(s2)
  print min(s3)
  print max(s3)
  ss = zerofloat(n1,n2,n3)
  FaultSkin.setValueOnFaults(-100,skins,ss)
  t1,t2,t3 = fsl.getDipSlipsX(skins,smark)
  plot3(add(gx,ss),t1,cmin=0,cmax=20.0,cmap=jetFillExceptMin(1.0),
        clab="Fault throw (samples)",png="gxss")
  plot3(gx,skins=skins,smax=20.0,png="skinss1")


def goInterp():
  if not plotOnly:
    gx = readImage(gxfile)

    '''
    lof = LocalOrientFilter(4.0,1.0,1.0)
    et = lof.applyForTensors(gx)
    et.setEigenvalues(0.001,1.0,1.0)
    '''
    ws = fillfloat(1.0,n1,n2,n3)
    mk = fillfloat(0.0,n1,n2,n3)
    skins = readSkins(fslbase)
    fsc = FaultSlipConstraints(skins)
    sp = fsc.screenPoints(ws)
    et = getEigenTensors()
    p = readImage(gfile)
    gp = GetPoints()
    wp = gp.getWeights(p)
    fx = gp.getValues(p)
    [x1,x2,x3] = gp.getCoordinates(p)
    figi = FastImageGuidedInterp(fx,x1,x2,x3)
    figi.setTensors(et)
    figi.setIters(400)
    figi.setSmoothings(10.0,10.0)
    fq = figi.apply(sp,wp)
    writeImage(fqfile,fq)
  else:
    gx = readImage(gxfile)
    fq = readImage(fqfile)
  display(gx,fq,vmin,vmax,logType,png=logType+"New")
  plot3(gx,fq,cmin=vmin,cmax=vmax,cmap=jetFill(0.3),
        clab=logLabel,png=logType)


def goInterpO():
  global k1,k2,k3
  k1,k2,k3 = 366,15,96
  #gridBlendedP()
  #gridBlendedQ()
  s = readImage(sfile); print "s min =",min(s)," max =",max(s)
  display1(s,True,cmin=vmin,cmax=vmax,png=logType+"Wells")
  #display1(s,False)
  #display1(s,False,["CrowMountainCRMT","TensleepASand"])
  #display1(s,True,["CrowMountainCRMT","TensleepASand"])
  #p = readImage(pfile); print "p min =",min(p)," max =",max(p)
  q = readImage(qfile); print "q min =",min(q)," max =",max(q)
  #t = readImage(tfile); print "t min =",min(t)," max =",max(t)
  #display(s,p,vmin,vmax,logType)
  display(s,q,vmin,vmax,logType,png=logType+"Old")
  #display(s,t,0.0,100.0,logType)
  #display(s,q,vmin,vmax,logType,["CrowMountainCRMT"])
  #display(s,q,vmin,vmax,logType,["TensleepASand"])

def goFigures():
  global k1,k2,k3
  k1,k2,k3 = 228,170,96 # intersect low-velocity layers
  s = readImage(sfile); print "s min =",min(s)," max =",max(s)
  #p = readImage(pfile); print "p min =",min(p)," max =",max(p)
  #q = readImage(qfile); print "q min =",min(q)," max =",max(q)
  #display3(s,None,0.0,0.0,"tpsz")
  #display3(s,p,vmin,vmax,"tppvb")
  #display3(s,q,vmin,vmax,"tpqvb")
  p = readImage("ig6/tppvb"); print "p min =",min(p)," max =",max(p)
  display3(s,p,vmin,vmax,"tppvb")
  p = readImage("tppvo09"); print "p min =",min(p)," max =",max(p)
  display3(s,p,vmin,vmax,"tppvb")

def gridBlendedP():
  e = getEigenTensors()
  bi = BlendedGridder3(e)
  p = readImage(gfile)
  t = bi.gridNearest(0.0,p)
  writeImage(pfile,p)
  writeImage(tfile,t)

def gridBlendedQ():
  e = getEigenTensors()
  bg = BlendedGridder3(e)
  bg.setSmoothness(1.0)
  p = readImage(pfile)
  t = readImage(tfile)
  t = clip(0.0,50.0,t)
  q = copy(p)
  bg.gridBlended(t,p,q)
  writeImage(qfile,q)

def getEigenTensors():
  e = readTensors(esfile)
  return e

def goImpedance():
  global k1,k2,k3,logLabel
  k1,k2,k3 = 366,15,96
  logType = None
  logLabel = "Impedance (g/cc x km/s)"
  d = readImage("tpqdb"); print "d min =",min(d)," max =",max(d)
  v = readImage("tpqvb"); print "v min =",min(v)," max =",max(v)
  s = readImage("tpsz"); print "s min =",min(s)," max =",max(s)
  i = mul(d,v)
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply2XX(mul(0.5,log(i)),i)
  #imin,imax = min(i),max(i)
  imin,imax = -0.05,0.05
  display2S(s,i,imin,imax,logType)

def display2S(s,g,cmin,cmax,logType,horizons=[]):
  world = World()
  ipg = addImageToWorld(world,s)
  ipg.setClips(smin,smax)
  ipg.setSlices(k1,k2,k3)
  ipg = addImageToWorld(world,g)
  ipg.setClips(cmin,cmax)
  ipg.setSlices(k1,k2,k3)
  if logType:
    addLogsToWorld(world,logSet,logType,cmin,cmax,smooth=smooth)
  for horizon in horizons:
    addHorizonToWorld(world,horizon)
  frame = makeFrame(world)

def display(s,g,cmin,cmax,logType,horizons=[],png=None):
  world = World()
  ipg = addImage2ToWorld(world,s,g)
  ipg.setClips1(smin,smax)
  ipg.setClips2(cmin,cmax)
  ipg.setSlices(k1,k2,k3)
  '''
  if logType:
    addLogsToWorld(world,logSet,logType,cmin,cmax,smooth=smooth)
  '''
  for horizon in horizons:
    addHorizonToWorld(world,horizon)
  frame = makeFrame(world)
  cbar = addColorBar(frame,logLabel)
  cbar.setWidthMinimum(120)
  ipg.addColorMap2Listener(cbar)
  if png and pngDir:
    frame.paintToFile(pngDir+png+".png")

def display1(s,wells=True,horizons=[],cmin=0,cmax=0,png=None):
  world = World()
  ipg = addImageToWorld(world,s)
  ipg.setClips(smin,smax)
  ipg.setSlices(k1,k2,k3)
  frame = makeFrame(world)
  if wells:
    cbar = addColorBar(frame,logLabel)
    addLogsToWorld(world,logSet,logType,cmin,cmax,cbar,smooth=smooth)
  for horizon in horizons:
    addHorizonToWorld(world,horizon)
  if png and pngDir:
    frame.paintToFile(pngDir+png+".png")

def addColorBar(frame,label):
  cbar = ColorBar(logLabel)
  cbar.setFont(cbar.getFont().deriveFont(36.0))
  frame.add(cbar,BorderLayout.EAST)
  cbar.setWidthMinimum
  #frame.viewCanvas.setBackground(frame.getBackground())
  return cbar

def display2(s,g=None,cmin=0,cmax=0):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.addColorBar()
  sp.getPlotPanel().setColorBarWidthMinimum(80)
  pv = sp.addPixels(s)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if g!=None:
    pv = sp.addPixels(g)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setColorModel(ColorMap.getJet(0.3))
    if cmin!=cmax:
      pv.setClips(cmin,cmax)

def display3(s,g=None,cmin=0,cmax=0,png=None):
  pp = PlotPanelPixels3(
    PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT,
    PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM,
    s1,s2,s3,s)
  pp.setSlices(k1,k2,k3)
  pp.setLabel1("Depth (km)")
  pp.setLabel2("Crossline (km)")
  pp.setLabel3("Inline (km)")
  pp.setClips(smin,smax)
  if g:
    pp.setLineColor(Color.BLACK)
    cb = pp.addColorBar(logLabel)
    cb.setInterval(1.0)
  else:
    pp.setLineColor(Color.YELLOW)
    cb = pp.addColorBar("Amplitude")
    cb.setInterval(5.0)
  pp.setColorBarWidthMinimum(100)
  pp.setInterval1(0.5)
  pp.setInterval2(2.0)
  pp.setInterval3(2.0)
  pp.mosaic.setHeightElastic(0,100)
  pp.mosaic.setHeightElastic(1,200)
  if g:
    pv12 = PixelsView(s1,s2,slice12(k3,g))
    pv12.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv13 = PixelsView(s1,s3,slice13(k2,g))
    pv13.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv23 = PixelsView(s2,s3,slice23(k1,g))
    pv23.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
    for pv in [pv12,pv13,pv23]:
      pv.setColorModel(ColorMap.getJet(0.5))
      if cmin!=cmax:
        pv.setClips(cmin,cmax)
    pp.pixelsView12.tile.addTiledView(pv12)
    pp.pixelsView13.tile.addTiledView(pv13)
    pp.pixelsView23.tile.addTiledView(pv23)
  pf = PlotFrame(pp)
  pf.setFontSizeForSlide(1.0,1.0)
  pf.setSize(996,815)
  pf.setVisible(True)
  if png and pngDir:
    pf.paintToPng(300,6,pngDir+png+".png")

#############################################################################

def like(x):
  n3 = len(x)
  n2 = len(x[0])
  n1 = len(x[0][0])
  return zerofloat(n1,n2,n3)


def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(10.0)
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

def wellTeapot(np,k2,k3):
  n1,n2,n3 =401,357,161
  image = zerofloat(n1,n2,n3)
  ais = ArrayInputStream("../../../data/seis/rgi/tpgd.dat")
  ais.readFloats(image)
  ais.close()
  gp = GetPoints()
  return gp.getXf(np,k2,k3,image)



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
          xyz=None,cells=None,skins=None,smax=0.0,
          links=False,curve=False,trace=False,png=None):
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  world = World()
  sf = SimpleFrame(world)
  cbar = None
  if g==None:
    ipg = sf.addImagePanels(s1,s2,s3,f)
    if cmap!=None:
      ipg.setColorModel(cmap)
    if cmin!=None and cmax!=None:
      ipg.setClips(cmin,cmax)
    else:
      ipg.setClips(smin,smax)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setColorModel1(ColorMap.getGray())
    ipg.setClips1(smin,smax)
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
    xyz,uvw,rgb = FaultCell.getXyzUvwRgbForLikelihood(0.5,cmap,cells,False)
    qg = QuadGroup(xyz,uvw,rgb)
    qg.setStates(ss)
    sf.world.addChild(qg)
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
      size = 0.6 
    ct = 0
    for skin in skins:
      if smax>0.0: # show fault throws
        cmap = ColorMap(-10,smax,ColorMap.JET)
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
          r,g,b=0,0,1
        if ct==1:
          r,g,b=0,1,0
        if ct==2:
          r,g,b=0,1,1
        if ct==3:
          r,g,b=0.627451,0.12549,0.941176
        xyz = skin.getCellLinksXyz()
        rgb = skin.getCellLinksRgb(r,g,b,xyz)
        #lg = LineGroup(xyz,rgb)
        lg = LineGroup(xyz)
        sg.addChild(lg)
        ct = ct+1
    sf.world.addChild(sg)
  ipg.setSlices(498,296,110)
  ipg.setSlices(k1,k2,k3)
  #ipg.setSlices(n1-1,0,n3-1)
  if cbar:
    sf.setSize(837,700)
    sf.setSize(1250,900)
  else:
    sf.setSize(700,700)
    sf.setSize(1250,900)

  view = sf.getOrbitView()
  zscale = 0.75*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.3)
  #view.setAzimuth(75.0)
  #view.setAzimuth(-75.0)
  view.setAzimuth(-65.0)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  sf.viewCanvas.setBackground(sf.getBackground())
  sf.setSize(1250,900)

  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")

#############################################################################
run(main)
