"""
Demonstrate structure-oriented processing
Author: Xinming Wu, Colorado School of Mines
Version: 2016.09.16
"""
from utils import *
setupForSubset("shapping")
s1,s2,s3= getSamplings()
n1,n2= s1.count,s2.count

# Names and descriptions of image files used below.
gxfile  = "gx" # migrated image 
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
minTheta,maxTheta = 75,85
sigmaTheta = 40

# These parameters control the construction of fault skins.
# See the class FaultSkinner for more information.
lowerLikelihood = 0.75
upperLikelihood = 0.85
minSize = 100

minThrow = 0.0
maxThrow = 20.0

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
pngDir = None
plotOnly = False

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  goTensors()
  goGaussianSmoothing()
  goIsotropicDiffusion()
  goAnisotropicDiffusion()
  #goVelocityInversion()

def goTensors():
  print "goTensors ..."
  gx = readImage2D(n1,n2,gxfile)
  fl = readImage2D(n1,n2,fltfile)
  el = zerofloat(n1,n2)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  lof = LocalOrientFilterP(8,1)
  et = lof.applyForTensors(gx)
  lof.applyForNormalLinear(gx,u1,u2,el)
  eu = zerofloat(n1,n2)
  ev = zerofloat(n1,n2)
  et.getEigenvalues(eu,ev)
  eu=clip(0.005,max(eu),eu)
  ev=clip(0.005,max(ev),ev)
  et.setEigenvalues(eu,ev)
  et.invertStructure(1.0,1.0)
  et.getEigenvalues(eu,ev)
  #plotTensors(gx,s1,s2,d=et,dscale=20,mk=mk,cmin=-2,cmax=2,png="tensors")
  plot2(s1,s2,gx)
  plotTensors(gx,s1,s2,d=et,dscale=1,ne=30,cmin=-2,cmax=2,png="tensors")
  plot2(s1,s2,el,cmin=0.1,cmax=1.0)

def goGaussianSmoothing():
  print "goGaussianSmoothing ..."
  gx = readImage2D(n1,n2,gxfile)
  gs = zerofloat(n1,n2)
  rgf = RecursiveGaussianFilter(10)
  rgf.apply00(gx,gs)
  plot2(s1,s2,gx)
  plot2(s1,s2,gs)
def goIsotropicDiffusion():
  print "goIsotropicDiffusion ..."
  gx = readImage2D(n1,n2,gxfile)
  gs = zerofloat(n1,n2)
  lof = LocalOrientFilterP(4,2)
  et = lof.applyForTensors(gx)
  et.setEigenvalues(1.0,1.0)
  lsf = LocalSmoothingFilter()
  lsf.apply(et,400,gx,gs)
  plot2(s1,s2,gx)
  plot2(s1,s2,gs)

def goAnisotropicDiffusion():
  print "goAnisotropicDiffusion ..."
  gx = readImage2D(n1,n2,gxfile)
  g1 = zerofloat(n1,n2)
  g2 = zerofloat(n1,n2)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lof = LocalOrientFilterP(4,2)
  lof.applyForNormalLinear(gx,u1,u2,el)
  el = pow(el,8)
  et = lof.applyForTensors(gx)
  et.setEigenvalues(0.001,1.0)
  lsf = LocalSmoothingFilter()
  lsf.apply(et,10,gx,g1)
  lsf.apply(et,10,el,gx,g2)
  plot2(s1,s2,gx)
  plot2(s1,s2,g1)
  plot2(s1,s2,g2)

def goVelocityInversion():
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
  #plot2(s1,s2,wp)
  clab = "Velocity (km/s)"
  clab1 = "Picked migration velocity (km/s)"
  clab2 = "Predicted migration velocity (km/s)"
  clab3 = "Interval velocity (km/s)"
  plot2(s1,s2,gx,g=pk,cmin=1.4,cmax=2.3,cmap=jetFill(1.0),label=clab1,png="pk")
  plot2(s1,s2,gx,g=vp,cmin=1.4,cmax=2.3,cmap=jetFill(1.0),label=clab2,png="vp")
  plot2(s1,s2,gx,g=vi,cmin=1.6,cmax=2.7,cmap=jetFill(1.0),label=clab3,png="vi")
  plot2(s1,s2,gx,g=vi,cmin=1.6,cmax=2.7,cmap=jetFill(0.6),label=clab3,png="seisvi")

def goFaultThrow():
  gx = readImage2D(n1,n2,gxfile)
  fl = readImage2D(n1,n2,flfile)
  ft = readImage2D(n1,n2,ftfile)
  for i2 in range(n2):
    for i1 in range(95):
      fl[i2][i1] = 0.01
  fc = FaultCurver()
  fc.setMinCurveSize(minSize)
  fc.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  ps = fc.findPoints([fl,ft])
  cs = fc.findCurves(ps)
  ftt = zerofloat(n1,n2)
  flt = zerofloat(n1,n2)
  cc = [cs[1],cs[3],cs[7],cs[9],cs[15],cs[16],cs[17]]
  #cc = [cs[1],cs[3],cs[7],cs[9],cs[11],cs[14],cs[15],cs[16],cs[17]]
  FaultCurve.getFlsImage(cc,flt)
  lof = LocalOrientFilterP(4,2)
  ets = lof.applyForTensors(gx)
  ets.setEigenvalues(0.001,1.0)
  lsf = LocalSmoothingFilter()
  wp = sub(1,flt);
  wp = pow(wp,10)
  gs = zerofloat(n1,n2)
  lsf.apply(ets,10,wp,gx,gs)
  lsp = LocalSlopeFinder(8,2,5)
  el = zerofloat(n1,n2)
  p2 = zerofloat(n1,n2)
  lsp.findSlopes(gx,p2,el)
  fcr = FaultCorrelater(gs,p2)
  fcr.setOffset(2)
  fcr.setZeroSlope(False) # True only if we want to show the error
  fcr.computeThrow(cc,minThrow,maxThrow)
  fst = zerofloat(n1,n2)
  FaultCurve.getFsImage(cc,fst)
  '''
  print min(fst)
  print max(fst)
  plot2(s1,s2,gx)
  plot2(s1,s2,gs)
  plot2(s1,s2,gx,g=flt,cmin=0.2,cmax=1,cmap=jetFillExceptMin(1.0))
  plot2(s1,s2,gx,g=mul(fst,4),cmin=minThrow,cmax=maxThrow,
        cmap=jetFillExceptMin(1.0),label="Fault throw (ms)",png="fst")
  smark = -999.999
  p1,p2 = fcr.getDipSlips(n1,n2,cc,smark)
  p1,p2 = fcr.interpolateDipSlips([p1,p2],smark)
  gw = fcr.unfault([p1,p2],gx)
  plot2(s1,s2,gw,label="Amplitude",png="gw")
  '''
  return cc


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

def plotTensors(g,s1,s2,d=None,dscale=1,ne=20,mk=None,cmin=0,cmax=0,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setBackground(backgroundColor)
  sp.setHLabel("Lateral position (km)")
  sp.setVLabel("Time (s)")

  sp.setHInterval(1.0)
  sp.setVInterval(1.0)
  sp.setFontSize(24)
  #sp.setFontSizeForPrint(8,240)
  #sp.setFontSizeForSlide(1.0,0.9)
  sp.setSize(523,800)
  pv = sp.addPixels(s1,s2,g)
  pv.setColorModel(ColorMap.GRAY)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  else:
    pv.setPercentiles(1,99)
  if d:
    tv = TensorsView(s1,s2,d)
    tv.setOrientation(TensorsView.Orientation.X1DOWN_X2RIGHT)
    tv.setLineColor(Color.YELLOW)
    tv.setLineWidth(2)
    if(mk):
      tv.setEllipsesDisplayed(mk)
    else:
      tv.setEllipsesDisplayed(ne)
    tv.setScale(dscale)
    tile = sp.plotPanel.getTile(0,0)
    tile.addTiledView(tv)
  if pngDir and png:
    sp.paintToPng(360,3.3,pngDir+png+".png")
    #sp.paintToPng(720,3.3,pngDir+png+".png")

def plot2(s1,s2,f,g=None,cmin=None,cmax=None,cmap=None,label=None,png=None):
  n2 = len(f)
  n1 = len(f[0])
  f1,f2 = s1.getFirst(),s2.getFirst()
  d1,d2 = s1.getDelta(),s2.getDelta()
  panel = panel2Teapot()
  panel.setHInterval(1.0)
  panel.setVInterval(1.0)
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
  frame.setFontSize(24)
  frame.setSize(550+80,800)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(400,3.2,pngDir+png+".png")
  return frame


#############################################################################
run(main)
