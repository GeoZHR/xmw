#############################################################################
# Dynamic warping for 1D sequences
"""
|f[i]-g[i+l]|

fg, nrms=0.0
c without s,u
c with s,u
fg, nrms=0.5
c without s,u
c with s,u
    sm=1
    sm=1/2
    sm=1/5
"""

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

from lsa import *
from pik import *

#############################################################################

pngDir = "../../../png/test/wlm/"
seismicDir = "../../../data/seis/tests/"
fxfile = "cbound"

n1,n2=400,10
s1 = Sampling(n1)
s2 = Sampling(n2)
seed = 99 
seed = 954 
seed = 2127 
seed = abs(Random().nextInt()/1000000)
seed = 877 
print "seed =",seed

nrms = 0.5
strainMax = 1.0

def main(args):
  goZhiguang()
  #d = goDynamicWarping()
  #d = goWarpWithSimilarity()
  #goPikWithError()
  #f,g,s = makeSequences()
  #plotfgX(f,d,png="fdx")
  #plotfgX(f,g,png="fgx")
  #l = goLocalSimilarity()
  #plotfgs(f,g,g,g,png="fgg")
  #plotfgs(f,g,d,l,png="fgs")
def goZhiguang():
  ml = 35
  strainMax = 0.2
  #f,g,s = makeSequences()
  fx = readImage2D(400,10,fxfile)
  f,g=fx[0],fx[5]
  dw = DynamicWarping(-ml,0)
  dw.setErrorExtrapolation(DynamicWarping.ErrorExtrapolation.NEAREST)
  dw.setStrainMax(strainMax)
  e = dw.computeErrors(f,g)
  d = dw.accumulateForward(e)
  u = dw.backtrackReverse(d,e)
  #plotfg(f,g,png="fg")
  #plot2c(etran(e),s,u,perc=93,png="cesu")
  #plotc(etran(e),None,None,perc=98,png="cesu")
  #plotc(dtran(d),s,u,perc=97,png="cdsu")
  ef = dw.accumulateForward(e)
  er = dw.accumulateReverse(e)
  es = sub(add(ef,er),e) # avoid normalization after smoothing
  #es = dw.smoothErrors(e) # this method normalizes after smoothing
  '''
  plotc(dtran(ef),None,None,cbar="Accumulated error",perc=97,png="cef")
  plotc(dtran(er),None,None,cbar="Accumulated error",perc=97,png="cer")
  plotc(dtran(es),None,None,cbar="Accumulated error",perc=97,png="ces")
  '''
  for strainMax in [0.4]:
    dw.setStrainMax(strainMax)
    e = dw.computeErrors(f,g)
    d = dw.accumulateForward(e)
    u = dw.backtrackReverse(d,e)
    #plotc(dtran(d),s,u,cbar="Accumulated error",perc=97,png="cdsu")
    ea = dw.smoothErrors(e)
    da = dw.accumulateForward(ea)
    ua = dw.backtrackReverse(da,ea)
    #plotc(dtran(da),s,ua,cbar="Accumulated error",perc=97,png="cdasu")
    #edge = RecursiveExponentialFilter.Edges.OUTPUT_ZERO_VALUE
    #ref = RecursiveExponentialFilter(10)
    #ref.setEdges(edge)
    #ref.apply(ua,ua)
    #plotc(dtran(d),s,ua,cbar="Accumulated error",perc=97,png="cdsus")
  fx[1] = fx[5]
  gx = copy(fx)
  gx[1] = dw.applyShifts(ua,fx[5])
  plot2(s1,s2,fx)
  plot2(s1,s2,gx)

  '''

  ml = 33
  sd = 1.0
  fx = readImage2D(400,10,fxfile)
  f = fx[0]
  g = fx[5]
  dw = DynamicWarping(-ml,ml)
  e = dw.computeErrors(f,g)
  sa = transpose(sub(1,e))
  sa = sub(sa,min(sa))
  sa = div(sa,max(sa))
  #plotc(sa,None,None,cbar="1-error",perc=99, png="lsdw")
  opp = OptimalPathPicker(3,0.5)
  wht = opp.applyForWeight(sa)
  n1,n2=len(sa[0]),len(sa)
  tms1 = zerofloat(n1,n2)
  tms2 = zerofloat(n1,n2)
  pik1 = opp.forwardPick(33,wht,tms1)
  piks = opp.applyForPath(33,0,sa)
  pik1 = sub(pik1,33)
  piks = sub(piks,33)
  #pik2 = opp.backwardPick(round(pik1[n2-1]),wht,tms2)
  #x2 = rampfloat(0,1,n2)
  tt = tms1
  plotc(tt,None,None,cbar="Travel time",perc=95,png="timedw")
  plotc(tt,None,None,cbar="Travel time",contour=True,perc=95,png="timeCdw")
  plotc(tt,pik1,pik1,cbar="Travel time",contour=True,perc=95,png="timePikdw")
  plotc(tt,piks,piks,cbar="Travel time",contour=True,perc=95,png="timePiksdw")
  fx[1] = fx[5]
  gx = copy(fx)
  gx[1] = dw.applyShifts(piks,fx[5])
  plot2(s1,s2,fx)
  plot2(s1,s2,gx)
  '''

def goWarpWithSimilarity():
  global nrms,strainMax
  ml = 33
  strainMax = 0.2
  f,g,s = makeSequences()
  dw = DynamicWarping(-ml,ml)
  dw.setErrorExtrapolation(DynamicWarping.ErrorExtrapolation.NEAREST)
  dw.setStrainMax(strainMax)
  ls = LocalSimilarity(-ml,ml,1)
  ls.setSmoothings(7,0)
  sa = ls.apply(g,f)
  sa = clip(0.001,1,sa)
  #sa = pow(sa,0.2)
  sa = sub(sa,min(sa))
  sa = div(sa,max(sa))
  e = sub(1,sa)
  e = transpose(e)
  d = dw.accumulateForward(e)
  u = dw.backtrackReverse(d,e)
  plotfg(f,g,png="fg")
  #plot2c(etran(e),s,u,perc=93,png="cesu")
  plotc(etran(e),None,None,cmin=0.1,cmax=0.9,cbar="1-similarity",png="cesusm")
  #plotc(dtran(d),s,u,perc=97,png="cdsu")
  ef = dw.accumulateForward(e)
  er = dw.accumulateReverse(e)
  es = sub(add(ef,er),e) # avoid normalization after smoothing
  #es = dw.smoothErrors(e) # this method normalizes after smoothing
  plotc(dtran(ef),None,None,cbar="Accumulated similarity",perc=97,png="cefsm")
  for strainMax in [0.2]:
    dw.setStrainMax(strainMax)
    ls = LocalSimilarity(-ml,ml,1)
    ls.setSmoothings(7,0)
    sa = ls.apply(g,f)
    sa = sub(sa,min(sa))
    sa = div(sa,max(sa))
    e = sub(1,sa)
    e = transpose(e)
    d = dw.accumulateForward(e)
    u = dw.backtrackReverse(d,e)
    plotc(dtran(d),s,u,cbar="Accumulated similarity",perc=97,png="cdsusm")
    ea = dw.smoothErrors(e)
    da = dw.accumulateForward(ea)
    ua = dw.backtrackReverse(da,ea)
    plotc(dtran(da),s,ua,cbar="Accumulated similarity",perc=97,png="cdasusm")
    edge = RecursiveExponentialFilter.Edges.OUTPUT_ZERO_VALUE
    ref = RecursiveExponentialFilter(10)
    ref.setEdges(edge)
    ref.apply(ua,ua)
    plotc(dtran(d),s,ua,cbar="Accumulated similarity",perc=97,png="cdsussm")

def goDynamicWarping():
  global nrms,strainMax
  ml = 33
  strainMax = 0.2
  f,g,s = makeSequences()
  dw = DynamicWarping(-ml,ml)
  dw.setErrorExtrapolation(DynamicWarping.ErrorExtrapolation.NEAREST)
  dw.setStrainMax(strainMax)
  e = dw.computeErrors(f,g)
  d = dw.accumulateForward(e)
  u = dw.backtrackReverse(d,e)
  plotfg(f,g,png="fg")
  #plot2c(etran(e),s,u,perc=93,png="cesu")
  plotc(etran(e),None,None,perc=98,png="cesu")
  #plotc(dtran(d),s,u,perc=97,png="cdsu")
  ef = dw.accumulateForward(e)
  er = dw.accumulateReverse(e)
  es = sub(add(ef,er),e) # avoid normalization after smoothing
  #es = dw.smoothErrors(e) # this method normalizes after smoothing
  plotc(dtran(ef),None,None,cbar="Accumulated error",perc=97,png="cef")
  plotc(dtran(er),None,None,cbar="Accumulated error",perc=97,png="cer")
  plotc(dtran(es),None,None,cbar="Accumulated error",perc=97,png="ces")
  for strainMax in [0.2]:
    dw.setStrainMax(strainMax)
    e = dw.computeErrors(f,g)
    d = dw.accumulateForward(e)
    u = dw.backtrackReverse(d,e)
    plotc(dtran(d),s,u,cbar="Accumulated error",perc=97,png="cdsu")
    ea = dw.smoothErrors(e)
    da = dw.accumulateForward(ea)
    ua = dw.backtrackReverse(da,ea)
    plotc(dtran(da),s,ua,cbar="Accumulated error",perc=97,png="cdasu")
    edge = RecursiveExponentialFilter.Edges.OUTPUT_ZERO_VALUE
    ref = RecursiveExponentialFilter(10)
    ref.setEdges(edge)
    ref.apply(ua,ua)
    plotc(dtran(d),s,ua,cbar="Accumulated error",perc=97,png="cdsus")
  return dw.applyShifts(ua,g)
    #print " esumu =",dw.sumErrors(e,u), "  esums =",dw.sumErrors(e,s);
    #print "easumu =",dw.sumErrors(ea,ua)," easums =",dw.sumErrors(ea,s);

def goPikWithError():
  ml = 33
  sd = 1.0
  f,g,s = makeSequences()
  dw = DynamicWarping(-ml,ml)
  e = dw.computeErrors(f,g)
  sa = transpose(sub(1,e))
  sa = sub(sa,min(sa))
  sa = div(sa,max(sa))
  plotc(sa,None,None,cbar="1-error",perc=99, png="lsdw")
  opp = OptimalPathPicker(3,1.5)
  wht = opp.applyForWeight(sa)
  n1,n2=len(sa[0]),len(sa)
  tms1 = zerofloat(n1,n2)
  tms2 = zerofloat(n1,n2)
  pik1 = opp.forwardPick(33,wht,tms1)
  piks = opp.applyForPath(33,10,sa)
  pik1 = sub(pik1,33)
  #piks = sub(piks,33)
  #pik2 = opp.backwardPick(round(pik1[n2-1]),wht,tms2)
  #x2 = rampfloat(0,1,n2)
  tt = tms1
  plotc(tt,None,None,cbar="Travel time",perc=95,png="timedw")
  plotc(tt,None,None,cbar="Travel time",contour=True,perc=95,png="timeCdw")
  plotc(tt,s,pik1,cbar="Travel time",contour=True,perc=95,png="timePikdw")
  plotc(tt,s,piks,cbar="Travel time",contour=True,perc=95,png="timePiksdw")
  '''
  gx = div(gx,max(gx))
  plot(gx,cmin=-1,cmax=1.0,w1=250,w2=1000,clab="Amplitude",png="seis")
  plot(fx,cmin=0.1,cmax=1.0,w1=250,w2=1000,clab="Amplitude",png="env")
  plot(tt,cmap=ColorMap.JET,cmin=0.1,cmax=max(tt),contour=True,w1=250,w2=1000,
       clab="Travel time",png="time")
  plot(tt,cmap=ColorMap.JET,cmin=0.1,cmax=max(tt),contour=True,
          xp=[pik2,x2],w1=250,w2=1000,clab="Travel time",png="timePik")
  plot(fx,cmin=0.1,cmax=1.0,xp=[pik2,x2],
        w1=250,w2=1000,clab="Amplitude",png="seisPik")
  '''
  return dw.applyShifts(piks,g)

def goLocalSimilarity():
  ml = 33
  sd = 1.0
  f,g,s = makeSequences()
  ls = LocalSimilarity(-ml,ml,sd)
  ls.setSmoothings(6,0)
  sa = ls.apply(g,f)
  dw = DynamicWarping(-ml,ml)
  '''
  dw.setErrorExtrapolation(DynamicWarping.ErrorExtrapolation.NEAREST)
  dw.setStrainMax(strainMax)
  e = dw.computeErrors(f,g)
  sa = transpose(sub(1,e))
  '''

  sa = sub(sa,min(sa))
  sa = div(sa,max(sa))
  plotc(sa,None,None,cbar="Local similarity",perc=99, png="ls")
  opp = OptimalPathPicker(3,0.5)
  wht = opp.applyForWeight(sa)
  n1,n2=len(sa[0]),len(sa)
  tms1 = zerofloat(n1,n2)
  tms2 = zerofloat(n1,n2)
  pik1 = opp.forwardPick(33,wht,tms1)
  piks = opp.applyForPath(33,10,sa)
  pik1 = sub(pik1,33)
  #piks = sub(piks,33)
  #pik2 = opp.backwardPick(round(pik1[n2-1]),wht,tms2)
  #x2 = rampfloat(0,1,n2)
  tt = tms1
  plotc(tt,None,None,cbar="Travel time",perc=95,png="time")
  plotc(tt,None,None,cbar="Travel time",contour=True,perc=95,png="timeC")
  plotc(tt,s,pik1,cbar="Travel time",contour=True,perc=95,png="timePik")
  plotc(tt,s,piks,cbar="Travel time",contour=True,perc=95,png="timePiks")
  '''
  gx = div(gx,max(gx))
  plot(gx,cmin=-1,cmax=1.0,w1=250,w2=1000,clab="Amplitude",png="seis")
  plot(fx,cmin=0.1,cmax=1.0,w1=250,w2=1000,clab="Amplitude",png="env")
  plot(tt,cmap=ColorMap.JET,cmin=0.1,cmax=max(tt),contour=True,w1=250,w2=1000,
       clab="Travel time",png="time")
  plot(tt,cmap=ColorMap.JET,cmin=0.1,cmax=max(tt),contour=True,
          xp=[pik2,x2],w1=250,w2=1000,clab="Travel time",png="timePik")
  plot(fx,cmin=0.1,cmax=1.0,xp=[pik2,x2],
        w1=250,w2=1000,clab="Amplitude",png="seisPik")
  '''
  return dw.applyShifts(piks,g)


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
def readImage2D(n1,n2,name):
  """ 
  Reads an image from a file with specified name.
  name: base name of image file; e.g., "tpsz"
  """
  fileName = seismicDir+name+".rsf@"
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  ais.readFloats(image)
  ais.close()
  return image
 
#############################################################################

# plotting

backgroundColor = Color.WHITE

def plotfgs(f,g,d,l,png=None):
  n = len(f)
  panel = PlotPanel(4,1,PlotPanel.Orientation.X1RIGHT_X2UP)
  panel.mosaic.setHeightElastic(0,25)
  panel.mosaic.setHeightElastic(1,25)
  panel.mosaic.setHeightElastic(2,25)
  panel.mosaic.setHeightElastic(3,25)
  panel.setHLimits(0,0,n-1)
  fv = panel.addPoints(0,0,f)
  gv = panel.addPoints(1,0,g)
  dv = panel.addPoints(2,0,d)
  lv = panel.addPoints(3,0,l)
  fv.setLineWidth(2)
  gv.setLineWidth(2)
  dv.setLineWidth(2)
  lv.setLineWidth(2)
  panel.setVLabel(0,"Amplitude")
  panel.setVLabel(1,"Amplitude")
  panel.setVLabel(2,"Amplitude")
  panel.setVLabel(3,"Amplitude")
  panel.setHLabel("Sample index")

  panel.setVLimits(0,-1.3,1.3)
  panel.setVLimits(1,-1.3,1.3)
  panel.setVLimits(2,-1.3,1.3)
  panel.setVLimits(3,-1.3,1.3)

  #panel.setHLabel("sample index i")
  #panel.setVLabel(0,"f")
  #panel.setVLabel(1,"g")
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(backgroundColor)
  frame.setFontSize(16)
  frame.setSize(1000,650)
  frame.setVisible(True)
  if png and pngDir:
    png += "n"+str(int(10*nrms))
    frame.paintToPng(720,3.33333,pngDir+"/"+png+".png")

def plotfgX(f,g,png=None):
  n = len(f)
  panel = PlotPanel(1,1,PlotPanel.Orientation.X1RIGHT_X2UP)
  panel.mosaic.setHeightElastic(0,25)
  panel.setHLimits(0,0,n-1)
  fv = panel.addPoints(0,0,f)
  gv = panel.addPoints(0,0,g)
  fv.setLineWidth(2)
  gv.setLineWidth(2)
  fv.setLineColor(Color.BLACK)
  gv.setLineColor(Color.RED)
  panel.setVLabel("Amplitude")
  panel.setHLabel("Sample index")
  panel.setVLimits(0,-1.3,1.3)
  #panel.setHLabel("sample index i")
  #panel.setVLabel(0,"f")
  #panel.setVLabel(1,"g")
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(backgroundColor)
  frame.setFontSize(24)
  frame.setSize(915,250)
  frame.setVisible(True)
  if png and pngDir:
    png += "n"+str(int(10*nrms))
    frame.paintToPng(720,3.33333,pngDir+"/"+png+".png")

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
  panel.setVLabel(0,"Amplitude")
  panel.setVLabel(1,"Amplitude")
  panel.setHLabel("Sample index")

  panel.setVLimits(0,-1.3,1.3)
  panel.setVLimits(1,-1.3,1.3)
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

def plotc(c,s=None,u=None,contour=False,cbar="Error",
    cmin=0.0,cmax=0.0,perc=None,png=None):
  n,nlag = len(c[0]),len(c)
  s1 = Sampling(n,1.0,0.0)
  slag = Sampling(nlag,1.0,-(nlag-1)/2)
  panel = PlotPanel(1,1,PlotPanel.Orientation.X1RIGHT_X2UP)
  panel.setVLabel("Shift lag")
  panel.setHLabel("Sample index")
  panel.setHLimits(0,0,n-1)
  panel.setVLimits(0,slag.first,slag.last)
  cv = panel.addPixels(0,0,s1,slag,c)
  cv.setColorModel(ColorMap.JET)
  cv.setInterpolation(PixelsView.Interpolation.NEAREST)

  if perc:
    cv.setPercentiles(100-perc,perc)
  if contour:
    cv = panel.addContours(s1,slag,c)
    cv.setContours(120)
    #cv.setColorModel(ColorMap.JET)
    cv.setLineWidth(3.0)
    cv.setLineColor(Color.BLACK)

  elif cmin<cmax:
    cv.setClips(cmin,cmax)
  if s:
    sv = panel.addPoints(0,0,s)
    sv.setLineColor(Color.WHITE)
    sv.setLineStyle(PointsView.Line.DOT)
    sv.setLineWidth(4)
  if u:
    uv = panel.addPoints(0,0,u)
    uv.setLineColor(Color.MAGENTA)
    uv.setLineWidth(4)
  panel.addColorBar(cbar)
  panel.setColorBarWidthMinimum(85)
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(backgroundColor)
  #frame.setFontSize(24)
  frame.setFontSize(24)
  frame.setSize(1000,350)
  frame.setVisible(True)
  if png and pngDir:
    png += "n"+str(int(10*nrms))
    png += "s"+str(int(10*strainMax))
    frame.paintToPng(720,3.33333,pngDir+"/"+png+".png")

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
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
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
  frame.setSize(n2/2,n1)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(400,3.2,pngDir+png+".png")
  return frame

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
