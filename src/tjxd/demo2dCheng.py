"""
Demonstrate velocity estimation
Author: Xinming Wu, Colorado School of Mines
Version: 2016.07.20
"""
from utils2 import *
setupForSubset("fake")
s1,s2,s3= getSamplings()
n1,n2= s1.count,s2.count

# Names and descriptions of image files used below.
gxfile  = "image" # migrated image 
gffile  = "gf"
smfile  = "sm" # migrated image 
pkfile  = "pk" # picked velocity 
epfile  = "ep" # eigenvalue-derived planarity
p2file  = "p2" # inline slopes
flfile  = "fl" # fault likelihood
ftfile  = "ft" # fault dip (theta)
fltfile = "flt" # fault likelihood thinned
fttfile = "ftt" # fault dip thinned

# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
minTheta,maxTheta = 10,45
sigmaTheta = 20

# These parameters control the construction of fault skins.
# See the class FaultSkinner for more information.
lowerLikelihood = 0.1
upperLikelihood = 0.50
minSize = 10

minThrow = -35
maxThrow =  0

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
plotOnly = False
pngDir = None
pngDir = "../../../png/tjxd/fake/"

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goFault()
  #goScan()
  #goFaultThrow()
  goInterp()
  goMigVel()
def goMigVel():
  gx = readImage2D(n1,n2,gxfile)
  gx = gain(gx)
  vx = readImage2D(401,801,"vel_mode")
  c1 = Sampling(401)
  c2 = Sampling(801)
  plot2x(gx,s1,s2,label="Amplitude",png="seis")
  plot2x(vx,c1,c2,g=vx,gmin=2500,gmax=6000,label="Velocity",png="migVel")
def goFault():
  p11 = [  0, 10, 42, 45,125]
  p12 = [360,354,300,270,125]
  p21 = [  5, 57, 76,115,139]
  p22 = [450,341,266,190,131]
  gx = readImage2D(n1,n2,gxfile)
  gx = gain(gx)
  hpr = Helper()
  fl = hpr.fault(n1,n2,[p11,p21],[p12,p22])
  writeImage(flfile,fl)
  plot2(s1,s2,gx,g=fl,cmin=0.20,cmax=1,cmap=jetRamp(1.0),
      label="Fault likelihood",png="fl")
def goInterp():
  vx = readImage2D(n1,n2,"vel")
  gx = readImage2D(n1,n2,gxfile)
  gx = gain(gx)
  x1,x2,fx = getWells(50,750)
  et = getTensors()
  bg = BlendedGridder2(et,fx,x1,x2)
  bg.setSmoothness(50)
  vi = bg.grid(s1,s2)
  plot2x(gx,s1,s2,f=fx,x1=x1,x2=x2,gmin=2500,gmax=6000,label="Velocity",png="seisAndLogs")
  plot2x(gx,s1,s2,g=vi,f=fx,x1=x1,x2=x2,gmin=2500,gmax=6000,label="Velocity",png="interpVel")
  plot2x(gx,s1,s2,g=vx,gmin=2500,gmax=6000,label="Velocity",png="trueVel")
  au = zerofloat(n1,n2)
  av = zerofloat(n1,n2)
  et.getEigenvalues(au,av)
  au = clip(0.003,1,au)
  et.setEigenvalues(au,av)
  e1 = Sampling(19,10,1)
  e2 = Sampling(21,40,1)
  plotTensors(gx,s1,s2,d=et,dscale=3.5,ne=20,e1=e1,e2=e2,cmin=0,cmax=0,png="tensors")

def getWells(k1,k2):
  vx = readImage2D(n1,n2,"vel")
  x1 = zerofloat(n1*2)
  x2 = zerofloat(n1*2)
  fx = zerofloat(n1*2)
  k = 0
  for i1 in range(n1):
    x1[k] = i1
    x2[k] = k1
    fx[k] = vx[k1][i1]
    k = k+1
    x1[k] = i1
    x2[k] = k2
    fx[k] = vx[k2][i1]
    k = k+1
  return x1,x2,fx
  

def getTensors():
  gx = readImage2D(n1,n2,gxfile)
  fl = readImage2D(n1,n2,flfile)
  gx = gain(gx)
  fl = sub(1,fl)
  av = clip(0.00001,1.0,fl)
  au = fillfloat(0.0001,n1,n2)
  lof = LocalOrientFilter(3,2)
  ets = lof.applyForTensors(gx)
  ets.setEigenvalues(au,av)
  return ets

def goScan():
  print "goScan ..."
  gx = readImage2D(n1,n2,gxfile)
  p11 = [ 10, 42, 45,125]
  p12 = [354,300,266,125]
  p21 = [  2, 57, 79,115,129]
  p22 = [463,346,260,192,152]
  gx = readImage2D(n1,n2,gxfile)
  gx = gain(gx)
  hpr = Helper()
  fxs = hpr.fault(n1,n2,[p11,p21],[p12,p22])
  if not plotOnly:
    gx = FaultScanner2.taper(10,0,gx)
    fs = FaultScanner2(sigmaTheta)
    fl,ft = fs.scan(minTheta,maxTheta,fxs)
    print "fl min =",min(fl)," max =",max(fl)
    print "ft min =",min(ft)," max =",max(ft)
    writeImage(flfile,fl)
    writeImage(ftfile,ft)
  else:
    fl = readImage2D(n1,n2,flfile)
    ft = readImage2D(n1,n2,ftfile)
  fl = sub(fl,min(fl))
  fl = div(fl,max(fl))
  plot2(s1,s2,gx,g=fl,cmin=0.10,cmax=1,cmap=jetRamp(1.0),
      label="Fault likelihood",png="fl")
  '''
  plot2(s1,s2,gx,g=abs(ft),cmin=minTheta,cmax=maxTheta,cmap=jetFill(1.0),
      label="Fault dip (degrees)",png="ft")
  '''
def goThin():
  print "goThin ..."
  gx = readImage2D(n1,n2,gxfile)
  gx = gain(gx)
  if not plotOnly:
    fl = readImage2D(n1,n2,flfile)
    ft = readImage2D(n1,n2,ftfile)
    fs = FaultScanner2(sigmaTheta)
    flt,ftt = fs.thin([fl,ft])
    writeImage(fltfile,flt)
    writeImage(fttfile,ftt)
  else:
    flt = readImage2D(n1,n2,fltfile)
    ftt = readImage2D(n1,n2,fttfile)
  plot2(s1,s2,gx)
  plot2(s1,s2,gx,g=flt,cmin=0.2,cmax=1,cmap=jetFillExceptMin(1.0))
  plot2(s1,s2,gx,g=abs(ftt),cmin=minTheta,cmax=maxTheta,cmap=jetFill(1.0),
      label="Fault dip (degrees)",png="ft")

def goSmooth():
  gx = readImage2D(n1,n2,gxfile)
  gx = gain(gx)
  fl = readImage2D(n1,n2,flfile)
  fl = sub(fl,min(fl))
  fl = div(fl,max(fl))
  ft = readImage2D(n1,n2,ftfile)
  fc = FaultCurver()
  fc.setMinCurveSize(minSize)
  fc.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  ps = fc.findPoints([fl,ft])
  cc = fc.findCurves(ps)
  ftt = zerofloat(n1,n2)
  flt = zerofloat(n1,n2)
  FaultCurve.getFlsImage(cc,flt)
  wp = sub(1,flt);
  wp = pow(wp,10)
  gs = zerofloat(n1,n2)
  lof = LocalOrientFilterP(8,2)
  ets = lof.applyForTensors(gx)
  ets.setEigenvalues(0.001,1.0)
  lsf = LocalSmoothingFilter()
  lsf.apply(ets,400,wp,gx,gs)
  plot2(s1,s2,gx)
  plot2(s1,s2,gs)
  plot2(s1,s2,gx,g=flt,cmin=0.1,cmax=1,cmap=jetFillExceptMin(1.0))
  return gs
def goFaultCurve():
  gx = readImage2D(n1,n2,gxfile)
  gx = gain(gx)
  fl = readImage2D(n1,n2,flfile)
  ft = readImage2D(n1,n2,ftfile)
  fc = FaultCurver()
  fc.setMinCurveSize(minSize)
  fc.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  ps = fc.findPoints([fl,ft])
  print len(ps)
  cs = fc.findCurves(ps)
  ftt = zerofloat(n1,n2)
  flt = fillfloat(-1000,n1,n2)
  FaultCurve.getFtImage(cs,ftt)
  FaultCurve.getFlImage(cs,flt)
  writeImage(fltfile,flt)
  plot2(s1,s2,gx,g=flt,cmin=0.1,cmax=1,cmap=jetFillExceptMin(1.0))

def goFaultThrow():
  gs = readImage2D(n1,n2,"vel")
  gx = readImage2D(n1,n2,gxfile)
  gx = gain(gx)
  fl = readImage2D(n1,n2,flfile)
  ft = readImage2D(n1,n2,ftfile)
  fc = FaultCurver()
  fc.setMinCurveSize(minSize)
  fc.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  ps = fc.findPoints([fl,ft])
  cc = fc.findCurves(ps)
  cs = [cc[1],cc[0]]
  ftt = zerofloat(n1,n2)
  flt = zerofloat(n1,n2)
  FaultCurve.getFlsImage(cc,flt)
  plot2(s1,s2,gx,g=flt,cmin=0.1,cmax=1,cmap=jetFillExceptMin(1.0))
  wp = sub(1,flt);
  wp = pow(wp,10)
  lsp = LocalSlopeFinder(8,2,5)
  el = zerofloat(n1,n2)
  p2 = zerofloat(n1,n2)
  lsp.findSlopes(gs,p2,el)
  #gs = goSmooth()
  plot2(s1,s2,gs,cmin=2500,cmax=5500)
  fcr = FaultCorrelater(gs,p2)
  fcr.setZeroSlope(False) # True only if we want to show the error
  fcr.setOffset(2)
  fcr.computeThrow([cs[0]],minThrow,maxThrow)
  fcr.setOffset(4)
  fcr.computeThrow([cs[1]],minThrow,maxThrow)
  fst = zerofloat(n1,n2)
  FaultCurve.getFsImage(cs,fst)
  fs1 = fillfloat(-1000,n1,n2)
  fs2 = fillfloat(-1000,n1,n2)
  FaultCurve.getFsImage(cs,fs1,fs2)
  writeImage("seis",gx)
  writeImage("slip1",fs1)
  writeImage("slip2",fs2)
  print min(fst)
  print max(fst)
  plot2(s1,s2,gx)
  smark = -999.999
  gw = copy(gx)
  for ci in cs:
    p1,p2 = fcr.getDipSlips(n1,n2,[ci],smark)
    p1,p2 = fcr.interpolateDipSlips([p1,p2],smark)
    gw = fcr.unfault([p1,p2],gw)
  writeImage("gw",gw)
  plot2(s1,s2,gw,label="Amplitude",png="gw")
  return cs

def goTensors():
  gx = readImage2D(n1,n2,gxfile)
  fl = readImage2D(n1,n2,fltfile)
  lof = LocalOrientFilterP(4,2)
  et = lof.applyForTensors(gx)
  eu = zerofloat(n1,n2)
  ev = zerofloat(n1,n2)
  et.getEigenvalues(eu,ev)
  eu=clip(0.008,max(eu),eu)
  ev=clip(0.008,max(ev),ev)
  et.setEigenvalues(eu,ev)
  et.invertStructure(1.0,1.0)
  et.getEigenvalues(eu,ev)
  print min(eu)
  print max(eu)
  #plotTensors(gx,s1,s2,d=et,dscale=20,mk=mk,cmin=-2,cmax=2,png="tensors")
  plot2(s1,s2,gx,cmin=-2,cmax=2,label="Amplitdue",png="seis")
  plotTensors(gx,s1,s2,d=et,dscale=1,ne=25,cmin=-2,cmax=2,png="tensors")

def goUnfault():
  cc = goFaultThrow()
  cc = [cc[0]]
  fl = zerofloat(n1,n2)
  fw = zerofloat(n1,n2)
  wp = fillfloat(1.0,n1,n2)
  gx = readImage2D(n1,n2,gxfile)
  gx = gain(gx)
  lof = LocalOrientFilterP(8,2)
  et = lof.applyForTensors(gx)
  et.setEigenvalues(0.001,1.0)
  slp = FaultSlipConstraints2(cc)
  sp = slp.screenPointsX(wp)
  uf = UnfaultS2(10,10)
  uf.setIters(100)
  uf.setTensors(et)
  np =  len(sp[0][0])
  scale = (n1*n2/np)
  mul(sp[3][0],scale,sp[3][0])
  [t1,t2] = uf.findShifts(sp,wp)
  #[t1,t2] = uf.convertShifts(40,[t1,t2])
  uf.applyShifts([t1,t2],gx,fw)
  plot2(s1,s2,gx)
  plot2(s1,s2,fw)
  plot2(s1,s2,wp)


def goVelocity():
  cc = goFaultThrow()
  fl = zerofloat(n1,n2)
  wp = fillfloat(1.0,n1,n2)
  slp = FaultSlipConstraints2(cc)
  sp = slp.screenPointsX(wp)
  #sp = None
  FaultCurve.getFlsImage(cc,fl)
  gx = readImage2D(n1,n2,gxfile)
  sm = readImage2D(n1,n2,smfile)
  pk = readImage2D(n1,n2,pkfile)
  lof = LocalOrientFilterP(4,2)
  et = lof.applyForTensors(gx)
  eu = zerofloat(n1,n2)
  ev = zerofloat(n1,n2)
  et.getEigenvalues(eu,ev)
  eu=clip(0.005,max(eu),eu)
  ev=clip(0.005,max(ev),ev)
  et.setEigenvalues(eu,ev)
  et.invertStructure(1.0,1.0)
  sm = clip(0.0,1.0,sm)
  sm = pow(sm,2)
  sm = mul(sm,wp)
  ve = VelocityEstimator(10,10)
  ve.setIterations(0.00001,400)
  ve.setTensors(et)
  ve.setSmoothness(200.0)
  ps = pow(pk,2)
  vi = ve.applyForVelocity(sp,sm,wp,ps)
  vi = sqrt(vi)
  vp = ve.predictVelocity(vi)
  plot2(s1,s2,wp)
  clab = "Velocity (km/s)"
  clab1 = "Picked migration velocity (km/s)"
  clab2 = "Predicted migration velocity (km/s)"
  clab3 = "Interval velocity (km/s)"
  plot2(s1,s2,gx,g=pk,cmin=1.4,cmax=2.3,cmap=jetFill(1.0),label=clab1,png="pk")
  plot2(s1,s2,gx,g=vp,cmin=1.4,cmax=2.3,cmap=jetFill(1.0),label=clab2,png="vp")
  plot2(s1,s2,gx,g=vi,cmin=1.6,cmax=2.7,cmap=jetFill(1.0),label=clab3,png="vi")
  plot2(s1,s2,gx,g=vi,cmin=1.6,cmax=2.7,cmap=jetFill(0.6),label=clab3,png="seisvi")

def goShapping():
  dix = readImage2D(n1,n2,"dix")
  shp = readImage2D(n1,n2,"shp")
  shpp = readImage2D(n1,n2,"shpp")
  gx = readImage2D(n1,n2,gxfile)
  clab1 = "Picked migration velocity (km/s)"
  clab2 = "Predicted migration velocity (km/s)"
  clab3 = "Interval velocity (km/s)"
  plot2(s1,s2,gx,g=dix,cmin=1.6,cmax=2.7,cmap=jetFill(1.0),label=clab3,png="dix")
  plot2(s1,s2,gx,g=shp,cmin=1.6,cmax=2.7,cmap=jetFill(1.0),label=clab3,png="shp")
  plot2(s1,s2,gx,g=shp,cmin=1.6,cmax=2.7,cmap=jetFill(0.6),label=clab3,png="seishp")
  plot2(s1,s2,gx,g=shpp,cmin=1.4,cmax=2.3,cmap=jetFill(1.0),label=clab2,png="shpp")
def like(x):
  n2 = len(x)
  n1 = len(x[0])
  return zerofloat(n1,n2)



def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(80.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2)
  div(x,sqrt(g),y)
  return y

#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET
backgroundColor = Color(0xfd,0xfe,0xff) # easy to make transparent
def jetFill(alpha):
  return ColorMap.setAlpha(ColorMap.JET,alpha)

def jetFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.JET,a)
def bwrNotch(alpha):
  a = zerofloat(256)
  for i in range(len(a)):
    if i<128:
      a[i] = alpha*(128.0-i)/128.0
    else:
      a[i] = alpha*(i-127.0)/128.0
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,a)


def bwrFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,a)

def jetRamp(alpha):
  return ColorMap.setAlpha(ColorMap.JET,rampfloat(0.0,alpha/256,256))

def bwrRamp(alpha):
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,rampfloat(0.0,alpha/256,256))

def grayRamp(alpha):
  return ColorMap.setAlpha(ColorMap.GRAY,rampfloat(0.0,alpha/256,256))

def plotTensors(g,s1,s2,d=None,dscale=1,ne=20,mk=None,e1=None,e2=None,cmin=0,cmax=0,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setBackground(backgroundColor)
  sp.setHLabel("Lateral position (sample)")
  sp.setVLabel("Depth (sample)")
  sp.setHInterval(100)
  sp.setVInterval(50)
  sp.setFontSize(24)
  #sp.setFontSizeForPrint(8,240)
  #sp.setFontSizeForSlide(1.0,0.9)
  sp.setSize(793,500)
  pv = sp.addPixels(s1,s2,g)
  pv.setColorModel(ColorMap.GRAY)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setClips(-1,1)
  if d:
    tv = TensorsView(s1,s2,d)
    tv.setOrientation(TensorsView.Orientation.X1DOWN_X2RIGHT)
    tv.setLineColor(Color.YELLOW)
    tv.setLineWidth(2)
    if(e1 and e2):
      tv.setEllipsesDisplayed(e1,e2)
    else:
      tv.setEllipsesDisplayed(ne)
    tv.setScale(dscale)
    tile = sp.plotPanel.getTile(0,0)
    tile.addTiledView(tv)
  sp.setFontSize(18)
  if pngDir and png:
    sp.paintToPng(720,3.3,pngDir+png+".png")
    #sp.paintToPng(720,3.3,pngDir+png+".png")


def plot2x(s,s1,s2,g=None,f=None,x1=None,x2=None,gmin=None,gmax=None,label=None,png=None,et=None):
  n2 = len(s)
  n1 = len(s[0])
  s1,s2=Sampling(n1),Sampling(n2)
  panel = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.LEFT_TOP)
  panel.setVLimits(0,n1-1)
  panel.setHInterval(100)
  panel.setVInterval(50)
  panel.setHLabel("Lateral position (sample)")
  panel.setVLabel("Depth (sample)")
  #panel.setHInterval(100.0)
  #panel.setVInterval(100.0)
  #panel.setHLabel("Pixel")
  #panel.setVLabel("Pixel")
  if label:
    panel.addColorBar(label)
  else:
    panel.addColorBar()
  panel.setColorBarWidthMinimum(80)
  pv = panel.addPixels(s1,s2,s)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setColorModel(ColorMap.GRAY)
  pv.setClips(-1.0,1.0)
  if g:
    alpha = 1.0
  if g:
    pv = panel.addPixels(s1,s2,g)
    pv.setInterpolation(PixelsView.Interpolation.LINEAR)
    pv.setColorModel(ColorMap.getJet(alpha))
    if gmin:
      pv.setClips(gmin,gmax)
  if x1 and x2:
    cmap = ColorMap(gmin,gmax,ColorMap.JET)
    fs,x1s,x2s = makePointSets(cmap,f,x1,x2)
    for i in range(len(fs)):
      #color = cmap.getColor((fs[i][0]-min(fs))/max(fs))
      color = cmap.getColor((fs[i][0]))
      #color = Color(color.red,color.green,color.blue)
      pv = panel.addPoints(x1s[i],x2s[i])
      pv.setLineStyle(PointsView.Line.NONE)
      pv.setMarkStyle(PointsView.Mark.FILLED_SQUARE)
      pv.setMarkSize(4)
      pv.setMarkColor(color)
  frame = PlotFrame(panel)
  #frame.setFontSizeForPrint(8,240)
  #frame.setSize(1240,774)
  #frame.setFontSizeForSlide(1.0,0.8)
  frame.setSize(900,500)
  frame.setVisible(True)
  frame.setFontSize(18)
  if png and pngDir:
    frame.paintToPng(720,3.3,pngDir+png+".png")
  return frame
def makePointSets(cmap,f,x1,x2):
  sets = {}
  for i in range(len(f)):
    if f[i] in sets:
      points = sets[f[i]]
      points[0].append(f[i])
      points[1].append(x1[i])
      points[2].append(x2[i])
    else:
      points = [[f[i]],[x1[i]],[x2[i]]] # lists of f, x1, x2
      sets[f[i]] = points
  ns = len(sets)
  fs = zerofloat(1,ns)
  x1s = zerofloat(1,ns)
  x2s = zerofloat(1,ns)
  il = 0
  for points in sets:
    fl = sets[points][0]
    x1l = sets[points][1]
    x2l = sets[points][2]
    nl = len(fl)
    fs[il] = zerofloat(nl)
    x1s[il] = zerofloat(nl)
    x2s[il] = zerofloat(nl)
    copy(fl,fs[il])
    copy(x1l,x1s[il])
    copy(x2l,x2s[il])
    il += 1
  return fs,x1s,x2s

def plot2(s1,s2,f,g=None,cmin=None,cmax=None,cmap1=ColorMap.GRAY,cmap=None,label=None,png=None):
  n2 = len(f)
  n1 = len(f[0])
  f1,f2 = s1.getFirst(),s2.getFirst()
  d1,d2 = s1.getDelta(),s2.getDelta()
  panel = panel2Teapot()
  panel.setHInterval(5.0)
  panel.setVInterval(5.0)
  panel.setHLabel("Lateral position (km)")
  panel.setVLabel("Time (s)")
  #panel.setHInterval(100.0)
  #panel.setVInterval(100.0)
  #panel.setHLabel("Pixel")
  #panel.setVLabel("Pixel")
  if label:
    panel.addColorBar(label)
  else:
    panel.addColorBar()
  panel.setColorBarWidthMinimum(80)
  pv = panel.addPixels(s1,s2,f)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setColorModel(ColorMap.GRAY)
  pv.setClips(-2,2)
  if g:
    pv = panel.addPixels(s1,s2,g)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setColorModel(cmap)
    if label:
      panel.addColorBar(label)
    else:
      panel.addColorBar()
  if cmin and cmax:
    pv.setClips(cmin,cmax)
  frame2Teapot(panel,png)
def panel2Teapot():
  panel = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT)#,PlotPanel.AxesPlacement.NONE)
  return panel
def frame2Teapot(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  #frame.setFontSizeForPrint(8,240)
  #frame.setSize(1240,774)
  #frame.setFontSizeForSlide(1.0,0.9)
  frame.setFontSize(12)
  frame.setSize(n2,n1*2)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(400,3.2,pngDir+png+".png")
  return frame


#############################################################################
run(main)
