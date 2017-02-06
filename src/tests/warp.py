#############################################################################
# Dynamic warping for 1D sequences
import sys
from java.awt import *
from java.awt.image import *
from java.io import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from swt import *

#############################################################################
pngDir = "./png/geop"
pngDir = None

nrms = 0.0
strainMax = 1.0

def main(args):
  goWarpForSeanOne()
  #goWarpForSeanTwo()
def goWarpForSeanOne():
  global seismicDir
  seismicDir = "../../../data/seis/sean/"
  fxfile = "tp_synthetic"
  gxfile = "tp_seistrace"
  umin,umax = -0.10,0.20
  rmin,rmax = -0.02,0.02
  dmin = 0.1
  n1 = 590
  ns = 313
  ds = 0.002
  fs = 0.104434
  s1 = Sampling(n1,0.002,0.0)
  sf = Sampling(ns,0.002,0.104434)
  gx = readImageM(n1,gxfile)
  fx = readImageM(ns,fxfile)
  #f = mul(4000,fx)
  #g = gx
  f = normalize(fx)
  g = normalize(gx)
  dw = DynamicWarpingRD(umin,umax,sf)
  dw.setStrainLimits(rmin,rmax)
  dw.setSmoothness(dmin)
  ss = dw.samplingS
  e = dw.computeErrors(sf,f,s1,g)
  u = dw.findShifts(e)
  st,r = invertShifts(sf,u)
  nt,dt,ft = st.count,st.delta,st.first
  print ft
  s = sub(rampfloat(ft,dt,nt),r)
  si = SincInterpolator()
  hx = zerofloat(nt)
  fx = mul(4000,fx)
  si.interpolate(ns,ds,fs,fx,nt,s,hx)
  plot1s(s1,sf,gx,fx)
  plot1s(s1,st,gx,hx)
  plot1(sf,u,umin,umax)
  print len(hx)
  writeImageM("tp_shifts",u)
  writeImageM("tp_warped",hx)

def goWarpForSeanTwo():
  global seismicDir
  seismicDir = "../../../data/seis/sean/"
  fxfile = "trace"
  gxfile = "synthetic"
  n1 = 423
  ml = 30
  sm = 0.25
  s1 = Sampling(n1)
  f = readImageM(n1,fxfile)
  g = readImageM(n1,gxfile)
  dw = DynamicWarping(-ml,ml)
  dw.setErrorExtrapolation(DynamicWarping.ErrorExtrapolation.NEAREST)
  dw.setStrainMax(sm)
  dw.setShiftSmoothing(2)
  u = dw.findShifts(f,g)
  h = dw.applyShifts(u,g)
  plot1(s1,u,0,ml)
  plot1s(s1,s1,f,g)
  plot1s(s1,s1,f,h)
  writeImageM("shifts",u)
  writeImageM("warped",h)

def invertShifts(ss,u):
  """ 
  Given a uniformly sampled u(s) such that s+u(s) increases monotonically, 
  computes a uniformly sampled r(t) = u(t-r(t)) such that t-r(t) increases 
  monotonically. Uses the following 3-step algorithm:
  (1) computes t(s) = s+u(s)
  (2) computes s(t) = by inverse linear interpolation of t(s)
  (3) computes r(t) = t-s(t)
  Returns the sampling of time st and the sequence of sampled r(t).
  """
  ns,ds,fs,ls = ss.count,ss.delta,ss.first,ss.last
  dt = ds # make sampling intervals equal
  ft = fs+u[0] # ft = time t of first sample
  lt = ls+u[ns-1] # lt = time of last sample
  ft = int(ft/dt)*dt # force ft to be a multiple of interval dt
  nt = 1+int((lt-ft)/dt) # number of t samples
  st = Sampling(nt,dt,ft) # sampling of t
  t = add(rampfloat(fs,ds,ns),u) # t(s) = s+u(s)
  s = zerofloat(nt)
  ii = InverseInterpolator(ss,st)
  ii.invert(t,s) # both t(s) and s(t) increase monotonically
  r = sub(rampfloat(ft,dt,nt),s) # r(t) = t-s(t)
  return st,r


def goFigures():
  global nrms,strainMax
  ml = 33
  strainMax = 1.0
  for nrms in [0.0,0.5]:
    f,g,s = makeSequences()
    dw = DynamicWarping(-ml,ml)
    dw.setErrorExtrapolation(DynamicWarping.ErrorExtrapolation.NEAREST)
    dw.setStrainMax(strainMax)
    e = dw.computeErrors(f,g)
    d = dw.accumulateForward(e)
    u = dw.backtrackReverse(d,e)
    plotfg(f,g,png="fg")
    #plot2c(etran(e),s,u,perc=93,png="cesu")
    plotc(etran(e),None,None,perc=95,png="cesu")
    plotc(dtran(d),s,u,perc=97,png="cdsu")
    if nrms>0:
      ef = dw.accumulateForward(e)
      er = dw.accumulateReverse(e)
      es = sub(add(ef,er),e) # avoid normalization after smoothing
      #es = dw.smoothErrors(e) # this method normalizes after smoothing
      plotc(dtran(ef),None,None,perc=97,png="cef")
      plotc(dtran(er),None,None,perc=97,png="cer")
      plotc(dtran(es),None,None,perc=97,png="ces")
  for strainMax in [1.0,0.5,0.2]:
    dw.setStrainMax(strainMax)
    e = dw.computeErrors(f,g)
    d = dw.accumulateForward(e)
    u = dw.backtrackReverse(d,e)
    plotc(dtran(d),s,u,perc=97,png="cdsu")
    ea = dw.smoothErrors(e)
    da = dw.accumulateForward(ea)
    ua = dw.backtrackReverse(da,ea)
    plotc(dtran(da),s,ua,perc=97,png="cdasu")
    #print " esumu =",dw.sumErrors(e,u), "  esums =",dw.sumErrors(e,s);
    #print "easumu =",dw.sumErrors(ea,ua)," easums =",dw.sumErrors(ea,s);

def gain(sig,x):
  n1 = len(x)
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(sig)
  ref.apply1(g,g)
  y = zerofloat(n1)
  div(x,sqrt(g),y)
  return y

def smooth(u):
  v = copy(u)
  rgf = RecursiveGaussianFilter(4)
  rgf.apply0(u,v)
  return v

def normalize(e):
  emin = min(e)
  emax = max(e)
  return mul(sub(e,emin),1.0/(emax-emin))

def etran(e):
  #return transpose(pow(e,0.25))
  return transpose(e)

def dtran(d):
  return transpose(d)

def makeSequences():
  n = 500
  fpeak = 0.125
  shift = 2.0/fpeak
  #w = Warp1Function.constant(shift,n)
  w = WarpFunction1.sinusoid(shift,n)
  #f = makeCosine(fpeak,n)
  f = makeRandomEvents(n,seed=seed); 
  g = w.warp(f)
  f = addRickerWavelet(fpeak,f)
  g = addRickerWavelet(fpeak,g)
  f = addNoise(nrms,fpeak,f,seed=10*seed+1)
  g = addNoise(nrms,fpeak,g,seed=10*seed+2)
  s = zerofloat(n)
  for i in range(n):
    s[i] = w.ux(i)
  return f,g,s

def makeCosine(freq,n):
  return cos(mul(2.0*PI*freq,rampfloat(0.0,1.0,n)))

def makeRandomEvents(n,seed=0):
  if seed!=0:
    r = Random(seed)
  else:
    r = Random()
  return pow(mul(2.0,sub(randfloat(r,n),0.5)),15.0)

def addRickerWavelet(fpeak,f):
  n = len(f)
  ih = int(3.0/fpeak)
  nh = 1+2*ih
  h = zerofloat(nh)
  for jh in range(nh):
    h[jh] = ricker(fpeak,jh-ih)
  g = zerofloat(n)
  Conv.conv(nh,-ih,h,n,0,f,n,0,g)
  return g

def ricker(fpeak,time):
  x = PI*fpeak*time
  return (1.0-2.0*x*x)*exp(-x*x)

def addNoise(nrms,fpeak,f,seed=0):
  n = len(f)
  if seed!=0:
    r = Random(seed)
  else:
    r = Random()
  nrms *= max(abs(f))
  g = mul(2.0,sub(randfloat(r,n),0.5))
  g = addRickerWavelet(fpeak,g)
  #rgf = RecursiveGaussianFilter(3.0)
  #rgf.apply1(g,g)
  frms = sqrt(sum(mul(f,f))/n)
  grms = sqrt(sum(mul(g,g))/n)
  g = mul(g,nrms*frms/grms)
  return add(f,g)
 

# read/write images
def readImageM(n1,basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".rsf@"
  image = zerofloat(n1)
  ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  ais.readFloats(image)
  ais.close()
  return image

def writeImageM(name,image):
  """ 
  Writes an image to a file with specified name.
  name: base name of image file; e.g., "tpgp"
  image: the image
  """
  fileName = seismicDir+name+".rsf@"
  aos = ArrayOutputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  aos.writeFloats(image)
  aos.close()
  return image

#############################################################################
# plotting

backgroundColor = Color.WHITE

def plot1(s1,y1,hmin,hmax,hlabel="Shifts",vlabel="Time (sample)",png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv1 = sp.addPoints(s1,y1)
  pv1.setLineColor(Color.BLACK)
  sp.setHLimits(hmin,hmax)
  sp.setSize(300,800)
  sp.setHLabel(hlabel)
  sp.setVLabel(vlabel)
  if png and pngDir:
    sp.paintToPng(300,7.0,pngDir+png+".png")


def plot1s(s1,s2,y1,y2,hlabel="Amplitude",vlabel="Time (sample)",png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv1 = sp.addPoints(s1,y1)
  pv1.setLineColor(Color.BLACK)
  pv2 = sp.addPoints(s2,y2)
  pv2.setLineColor(Color.RED)
  sp.setHLimits(-8,8)
  sp.setSize(300,800)
  sp.setHLabel(hlabel)
  sp.setVLabel(vlabel)
  if png and pngDir:
    sp.paintToPng(300,7.0,pngDir+png+".png")

def plotfg(f,g,png=None):
  n = len(f)
  panel = PlotPanel(2,1,PlotPanel.Orientation.X1RIGHT_X2UP)
  panel.mosaic.setHeightElastic(0,25)
  panel.mosaic.setHeightElastic(1,25)
  panel.setHLimits(0,0,n-1)
  fv = panel.addPoints(0,0,f)
  gv = panel.addPoints(1,0,g)
  fv.setLineWidth(2)
  gv.setLineWidth(2)
  panel.setVLimits(0,-4.0,4.0)
  panel.setVLimits(1,-4.0,4.0)
  #panel.setHLabel("sample index i")
  #panel.setVLabel(0,"f")
  #panel.setVLabel(1,"g")
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(backgroundColor)
  frame.setFontSizeForPrint(8,240)
  frame.setSize(1000,650)
  frame.setVisible(True)
  if png and pngDir:
    png += "n"+str(int(10*nrms))
    frame.paintToPng(720,3.33333,pngDir+"/"+png+".png")

def plotc(c,s=None,u=None,cmin=0.0,cmax=0.0,perc=None,png=None):
  n,nlag = len(c[0]),len(c)
  s1 = Sampling(n,1.0,0.0)
  slag = Sampling(nlag,1.0,-(nlag-1)/2)
  panel = PlotPanel(1,1,PlotPanel.Orientation.X1RIGHT_X2UP)
  panel.setHLimits(0,0,n-1)
  panel.setVLimits(0,slag.first,slag.last)
  cv = panel.addPixels(0,0,s1,slag,c)
  cv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if perc:
    cv.setPercentiles(100-perc,perc)
  elif cmin<cmax:
    cv.setClips(cmin,cmax)
  if s:
    sv = panel.addPoints(0,0,s)
    sv.setLineColor(Color.WHITE)
    sv.setLineStyle(PointsView.Line.DOT)
    sv.setLineWidth(3)
  if u:
    uv = panel.addPoints(0,0,u)
    uv.setLineColor(Color.WHITE)
    uv.setLineWidth(3)
  panel.addColorBar()
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(backgroundColor)
  frame.setFontSizeForPrint(8,240)
  frame.setSize(1000,470)
  frame.setVisible(True)
  if png and pngDir:
    png += "n"+str(int(10*nrms))
    png += "s"+str(int(10*strainMax))
    frame.paintToPng(720,3.33333,pngDir+"/"+png+".png")

def plot2c(c,s,u,clip=None,perc=None,png=None):
  n,nlag = len(c[0]),len(c)
  s1 = Sampling(n,1.0,0.0)
  slag = Sampling(nlag,1.0,-(nlag-1)/2)
  panel = PlotPanel(2,1,PlotPanel.Orientation.X1RIGHT_X2UP)
  panel.setHLimits(0,0,n-1)
  panel.setVLimits(0,slag.first,slag.last)
  panel.setVLimits(1,slag.first,slag.last)
  cv0 = panel.addPixels(0,0,s1,slag,c)
  cv1 = panel.addPixels(1,0,s1,slag,c)
  cv0.setInterpolation(PixelsView.Interpolation.NEAREST)
  cv1.setInterpolation(PixelsView.Interpolation.NEAREST)
  cv0.setColorModel(ColorMap.getGray(0.0,0.8))
  cv1.setColorModel(ColorMap.getGray(0.0,0.8))
  if perc:
    cv0.setPercentiles(0,perc)
    cv1.setPercentiles(0,perc)
  elif clip:
    cv0.setClips(0.0,clip)
    cv1.setClips(0.0,clip)
  if s:
    sv = panel.addPoints(1,0,s)
    sv.setLineColor(Color.WHITE)
    sv.setLineStyle(PointsView.Line.DOT)
    sv.setLineWidth(3)
  if u:
    uv = panel.addPoints(1,0,u)
    uv.setLineColor(Color.WHITE)
    uv.setLineWidth(3)
  panel.addColorBar()
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(backgroundColor)
  frame.setFontSizeForPrint(8,240)
  frame.setSize(1000,850)
  frame.setVisible(True)
  if png and pngDir:
    png += "n"+str(int(10*nrms))
    png += "s"+str(int(10*strainMax))
    frame.paintToPng(720,3.33333,pngDir+"/"+png+".png")

def plot3c(c,s,u,cmin=0.0,cmax=0.0,png=None):
  print "c0: min =",min(c[0])," max =",max(c[0])
  print "c1: min =",min(c[1])," max =",max(c[1])
  print "c2: min =",min(c[2])," max =",max(c[2])
  n,nlag = len(c[0][0]),len(c[0])
  s1 = Sampling(n,1.0,0.0)
  slag = Sampling(nlag,1.0,-(nlag-1)/2)
  panel = PlotPanel(3,1,PlotPanel.Orientation.X1RIGHT_X2UP)
  panel.setHLimits(0,0,n-1)
  panel.setVLimits(0,slag.first,slag.last)
  panel.setVLimits(1,slag.first,slag.last)
  panel.setVLimits(2,slag.first,slag.last)
  cv0 = panel.addPixels(0,0,s1,slag,c[0])
  cv1 = panel.addPixels(1,0,s1,slag,c[1])
  cv2 = panel.addPixels(2,0,s1,slag,c[2])
  cv0.setInterpolation(PixelsView.Interpolation.NEAREST)
  cv1.setInterpolation(PixelsView.Interpolation.NEAREST)
  cv2.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmin<cmax:
    cv0.setClips(cmin,cmax)
    cv1.setClips(cmin,cmax)
    cv2.setClips(cmin,cmax)
  if s:
    sv = panel.addPoints(0,0,s[0])
    sv.setLineColor(Color.WHITE)
    sv.setLineStyle(PointsView.Line.DOT)
    sv.setLineWidth(3)
    sv = panel.addPoints(1,0,s[1])
    sv.setLineColor(Color.WHITE)
    sv.setLineStyle(PointsView.Line.DOT)
    sv.setLineWidth(3)
    sv = panel.addPoints(2,0,s[2])
    sv.setLineColor(Color.WHITE)
    sv.setLineStyle(PointsView.Line.DOT)
    sv.setLineWidth(3)
  if u:
    uv = panel.addPoints(0,0,u[0])
    uv.setLineColor(Color.WHITE)
    uv.setLineWidth(3)
    uv = panel.addPoints(1,0,u[1])
    uv.setLineColor(Color.WHITE)
    uv.setLineWidth(3)
    uv = panel.addPoints(2,0,u[2])
    uv.setLineColor(Color.WHITE)
    uv.setLineWidth(3)
  panel.addColorBar()
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(backgroundColor)
  frame.setFontSizeForPrint(8,240)
  frame.setSize(1000,1050)
  frame.setVisible(True)
  if png and pngDir:
    png += "n"+str(int(10*nrms))
    #png += "s"+str(int(10*strainMax))
    frame.paintToPng(720,3.33333,pngDir+"/"+png+".png")

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
