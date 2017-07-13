#############################################################################
"""
Demo of dynamic warping for automatic picking
Author: Xinming Wu, University of Texas at Austin
Version: 2016.06.01
"""


from utils import * 
setupForSubset("edge")
setupForSubset("fault2d")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count
f1,f2,f3 = s1.getFirst(),s2.getFirst(),s3.getFirst()
d1,d2,d3 = s1.getDelta(),s2.getDelta(),s3.getDelta()
#############################################################################
gxfile = "gx238" # input semblance image
gxfile = "fxnwc" # input semblance image
gxfile = "246009" # edge
gxfile = "388067" # edge
gxfile = "teste" # edge
gxfile = "f3d75s" # input semblance image
gxfile = "ep56" # input semblance image
gxfile = "clyde200" # input semblance image
elfile = "el" # picked path using Sergey's method
smfile = "sm"
pngDir = getPngDir()
pngDir = None
plotOnly = False

def main(args):
  #goSemblance()
  #goFaultPik()
  #goFaultPik1()
  #goTimeMarker()
  #goPolar()
  #goEdgeEnhance()
  goEdgeTime()
def goPolar():
  sm = readImage(smfile)
  sm = sub(1,sm)
  '''
  gx = readImage(gxfile)
  ss = copy(100,100,120,100,sm)
  gs = copy(100,100,120,100,gx)
  plot(gs,ss,cmin=0.1,cmax=0.5,cmap=jetFillExceptMin(0.6),cint=0.2)
  k = 1.2
  n2 = 100
  for i2 in range(n2):
    i1 = round(k*i2)
    if i1<100:
      ss[i2][i1] = 2
  rgf = RecursiveGaussianFilterP(1)
  rgf.apply00(ss,ss)
  '''
  pc = PolarCoordinates(555,125,100,0,360)
  sc,gx = pc.applyTransform(sm)
  plot(sm)
  plot(sc)
  plot(gx)

def goTimeMarker():
  gx = readImage(gxfile)
  gx = gain(gx)
  fe = FaultEnhance(4,1.0)
  st1 = zerofloat(n1,n2)
  seeds=fe.findSeeds(5,1,0.05,sm,st1)
  lof = LocalOrientFilter(4,2)
  ets = lof.applyForTensors(gx)
  ets.setEigenvalues(sm,sm)
  k1 = [108]
  k2 = [511]
  f = [1]
  t = fillfloat(1,n1,n2)
  p = zerofloat(n1,n2)
  for i in range(len(k1)):
    p[k2[i]][k1[i]] = 1
    t[k2[i]][k1[i]] = 0
  bd = BlendedGridder2(ets,f,k1,k2)
  bd.gridNearest(t,p)
  u = zerofloat(n1)
  for i1 in range(n1):
    tm = t[0][i1]
    for i2 in range(n2):
      if(tm>t[i2][i1]):
        tm = t[i2][i1]
        u[i1] = i2
  cmin = -1
  cmax =  1
  mp1 = ColorMap.GRAY
  mp2 = ColorMap.JET
  plot2(s1,s2,sm,u=u,vint=10,hint=20,cmin=0.0,cmax=0.5,cmap=mp1)
  plot(gx,t,cmin=0.0,cmax=max(t),cmap=jetFillExceptMin(0.6),cint=0.2)

def goSemblance():
  gx = readImage(gxfile)
  lof = LocalOrientFilter(2,2)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lof.applyForNormalLinear(gx,u1,u2,el)
  writeImage(smfile,el)

def goFaultPik1():
  gx = readImage(gxfile)
  gx = gain(gx)
  sm = readImage(smfile)
  plot(sub(1,sm))
  sm = pow(sm,8)
  sm = sub(1,sm)
  fe = FaultEnhance(4,0.2)
  st1 = zerofloat(n1,n2)
  seeds=fe.findSeeds(1,1,0.3,sm,st1)
  se = fe.enhanceInPolarSpace1(60,60,10,seeds,sm)
  se = sub(se,min(se))
  se = div(se,max(se))
  plot(st1)
  plot(se,cmin=0.01,cmax=0.6)
  cmin = -1
  cmax =  1
  mp1 = ColorMap.GRAY
  mp2 = ColorMap.JET
  plot(gx,sm,cmin=0.3,cmax=1.0,cmap=jetFillExceptMin(0.6),cint=0.2)
  plot(gx,se,cmin=0.1,cmax=0.3,cmap=jetFillExceptMin(0.6),cint=0.2)
  '''
  plot(gx,stt,cmin=0.1,cmax=0.3,cmap=jetFillExceptMin(0.6),cint=0.2)
  plot2(s1,s2,et,u=pik1,vint=20,hint=20,cmin=0,cmax=0.5,cmap=mp2)
  plot2(s1,s2,gx,u=pik1,vint=20,hint=20,cmin=cmin,cmax=cmax,cmap=mp1)
  plot2(s1,s2,gx,u=pik2,vint=20,hint=20,cmin=cmin,cmax=cmax,cmap=mp1)
  '''
  #plot2(s1,s2,gx,u=pik3,vint=20,hint=20,cmin=cmin,cmax=cmax,cmap=mp)
  
def goEdgeEnhance():
  gx = readImage(gxfile)
  '''
  g1 = zerofloat(n1,n2)
  g2 = zerofloat(n1,n2)
  plot(gx,cmap=ColorMap.GRAY)
  rgf = RecursiveGaussianFilterP(1)
  rgf.apply1X(gx,g1)
  rgf.applyX1(gx,g2)
  g1 = mul(g1,g1)
  g2 = mul(g2,g2)
  gs = add(g1,g2)
  gs = sqrt(gs)
  gs = sub(gs,min(gs))
  sm = div(gs,max(gs))
  plot(gs,cmap=ColorMap.GRAY)
  '''
  gx = sub(gx,min(gx))
  sm = div(gx,max(gx))
  sm = sub(1,sm)
  fe = FaultEnhance(10,2)
  st1 = zerofloat(n1,n2)
  seeds=fe.findSeedsX(2,2,0.1,sm,st1)
  plot(sm)
  plot(st1)
  #se = fe.applyForEnhanceX(40,40,8,seeds,sm)
  se = fe.enhanceInPolarSpace(30,30,8,seeds,st1)
  se = sub(se,min(se))
  se = div(se,max(se))
  plot(se,cmin=0.01,cmax=0.5)
  '''
  cmin = -1
  cmax =  1
  mp1 = ColorMap.GRAY
  mp2 = ColorMap.JET
  plot(gx,sm,cmin=0.3,cmax=1.0,cmap=jetFillExceptMin(0.6),cint=0.2)
  plot(gx,se,cmin=0.1,cmax=0.3,cmap=jetFillExceptMin(0.6),cint=0.2)
  plot(gx,stt,cmin=0.1,cmax=0.3,cmap=jetFillExceptMin(0.6),cint=0.2)
  plot2(s1,s2,et,u=pik1,vint=20,hint=20,cmin=0,cmax=0.5,cmap=mp2)
  plot2(s1,s2,gx,u=pik1,vint=20,hint=20,cmin=cmin,cmax=cmax,cmap=mp1)
  plot2(s1,s2,gx,u=pik2,vint=20,hint=20,cmin=cmin,cmax=cmax,cmap=mp1)
  '''
  #plot2(s1,s2,gx,u=pik3,vint=20,hint=20,cmin=cmin,cmax=cmax,cmap=mp)

def goEdgeTime():
  gx = readImage(smfile)
  gx = pow(gx,4)
  gx = sub(gx,min(gx))
  sm = div(gx,max(gx))
  sm = sub(1,sm)
  '''
  g1 = zerofloat(n1,n2)
  g2 = zerofloat(n1,n2)
  rgf = RecursiveGaussianFilterP(1)
  rgf.apply10(gx,g1)
  rgf.apply01(gx,g2)
  gs = sqrt(add(mul(g1,g1),mul(g2,g2)))
  gs = sub(gs,min(gs))
  sm = div(gs,max(gs))
  '''
  lof = LocalOrientFilter(4,4)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lof.applyForNormalLinear(gx,u1,u2,el)
  ee = EdgeEnhance()
  st = ee.findRidges(0.01,sm)
  sd = ee.pickSeeds(4,0.3,st)
  fd = ee.seedsToImage(sd,st)
  #ps1,ps2 = ee.applyEnhance(650,520,100,8,u1,u2,st1)
  ss = zerofloat(n1,n2)
  rgf = RecursiveGaussianFilterP(1)
  rgf.apply00(st,ss)
  ss = sub(ss,min(ss))
  ss = div(ss,max(ss))
  se = ee.applyEnhance(60,10,sd,u1,u2,ss)
  se = sub(se,min(se))
  se = div(se,max(se))
  #ps = ee.backTrack(350,550,600,450,t)
  cmin = -1
  cmax =  1
  mp1 = ColorMap.GRAY
  mp2 = ColorMap.JET
  st = sub(1,st)
  se = sub(1,se)
  plot(gx)
  plot(fd)
  plot(st,cmin=0.7,cmax=1.0)
  plot(se,cmin=0.7,cmax=1.0)

def goFaultPik():
  #gx = readImage(gxfile)
  #gx = gain(gx)
  sm = readImage(gxfile)
  plot(sub(1,sm))
  sm = pow(sm,8)
  sm = sub(1,sm)
  fe = FaultEnhance(4,0.5)
  st1 = zerofloat(n1,n2)
  seeds=fe.findSeedsX(1,1,0.3,sm,st1)
  #se = fe.applyForEnhanceX(20,20,8,seeds,sm)
  se = fe.enhanceInPolarSpace(40,40,10,seeds,sm)
  se = sub(se,min(se))
  se = div(se,max(se))
  plot(st1)
  plot(se,cmin=0.01,cmax=0.5)
  '''
  cmin = -1
  cmax =  1
  mp1 = ColorMap.GRAY
  mp2 = ColorMap.JET
  plot(gx,sm,cmin=0.3,cmax=1.0,cmap=jetFillExceptMin(0.6),cint=0.2)
  plot(gx,se,cmin=0.1,cmax=0.3,cmap=jetFillExceptMin(0.6),cint=0.2)
  plot(gx,stt,cmin=0.1,cmax=0.3,cmap=jetFillExceptMin(0.6),cint=0.2)
  plot2(s1,s2,et,u=pik1,vint=20,hint=20,cmin=0,cmax=0.5,cmap=mp2)
  plot2(s1,s2,gx,u=pik1,vint=20,hint=20,cmin=cmin,cmax=cmax,cmap=mp1)
  plot2(s1,s2,gx,u=pik2,vint=20,hint=20,cmin=cmin,cmax=cmax,cmap=mp1)
  '''
  #plot2(s1,s2,gx,u=pik3,vint=20,hint=20,cmin=cmin,cmax=cmax,cmap=mp)

def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(100.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2)
  div(x,sqrt(g),y)
  return y

def smooth(sig,u):
  v = copy(u)
  rgf = RecursiveGaussianFilterP(sig)
  rgf.apply0(u,v)
  return v

def smooth2(sig1,sig2,u):
  v = copy(u)
  rgf1 = RecursiveGaussianFilterP(sig1)
  rgf2 = RecursiveGaussianFilterP(sig2)
  rgf1.apply0X(u,v)
  rgf2.applyX0(v,v)
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

#############################################################################
# plotting
def jetFill(alpha):
  return ColorMap.setAlpha(ColorMap.JET,alpha)
def jetFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.JET,a)
def jetRamp(alpha):
  return ColorMap.setAlpha(ColorMap.JET,rampfloat(0.0,alpha/256,256))
def bwrRamp(alpha):
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,rampfloat(0.0,alpha/256,256))
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


backgroundColor = Color.WHITE

def plot1(f):
  n1 = len(f)
  s1 = Sampling(n1)
  sp = SimplePlot()
  pv = sp.addPoints(s1,f)

def plot(f,g=None,ps=None,t=None,cmap=None,cmin=None,cmax=None,cint=None,
        label=None,neareast=False,png=None): 
  orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT;
  n1,n2=len(f[0]),len(f)
  s1,s2=Sampling(n1),Sampling(n2)
  panel = PlotPanel(1,1,orientation)#,PlotPanel.AxesPlacement.NONE)
  panel.setVInterval(50)
  panel.setHInterval(50)
  #panel.setHLabel("Inline (traces)")
  #panel.setVLabel("Time (samples)")
  pxv = panel.addPixels(0,0,s1,s2,f);
  pxv.setColorModel(ColorMap.GRAY)
  #pxv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if g:
    pxv.setClips(-1,1)
  else:
    if cmin and cmax:
      pxv.setClips(cmin,cmax)
  if g:
    pv = panel.addPixels(s1,s2,g)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setColorModel(cmap)
    if cmin and cmax:
      pv.setClips(cmin,cmax)
  if ps:
    uv = panel.addPoints(0,0,ps[0],ps[1])
    uv.setLineColor(Color.YELLOW)
    uv.setLineWidth(2)
  moc = panel.getMosaic();
  frame = PlotFrame(panel);
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  #frame.setTitle("normal vectors")
  frame.setVisible(True);
  #frame.setSize(1400,700)
  frame.setSize(n2,n1)
  frame.setFontSize(24)
  if pngDir and png:
    frame.paintToPng(720,3.333,pngDir+png+".png")

def plot2(s1,s2,c,u=None,us=None,ss=None,cps=None,css=None,vint=1,hint=1,
          cmin=0.0,cmax=0.0,cmap=ColorMap.JET,title=None,perc=None,png=None):
  panel = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT)
          #PlotPanel.AxesPlacement.NONE)
  panel.setHLimits(0,s2.first,s2.last)
  panel.setVLimits(0,s1.first,s1.last)
  panel.setVInterval(0,vint)
  panel.setHInterval(0,hint)
  if title:
    panel.addTitle(title)
  cv = panel.addPixels(0,0,s1,s2,c)
  cv.setInterpolation(PixelsView.Interpolation.LINEAR)
  cv.setColorModel(cmap)
  if perc:
    cv.setPercentiles(100-perc,perc)
  elif cmin<cmax:
    cv.setClips(cmin,cmax)
  if u:
    uv = panel.addPoints(0,0,s1,u)
    uv.setLineColor(Color.RED)
    uv.setLineWidth(2)
  if us:
    colors = [Color.RED,Color.GREEN,Color.BLUE]
    for k in range(len(us)):
      sk = ss[k]
      uk = copy(ss[k].getCount(),0,us[k])
      uv = panel.addPoints(0,0,sk,uk)
      uv.setLineColor(colors[k])
      #uv.setLineColor(Color.WHITE)
      uv.setLineStyle(PointsView.Line.DASH)
      uv.setLineWidth(3)
  if cps and css:
    colors = [Color.RED,Color.GREEN,Color.BLUE]
    for k in range(len(cps)):
      pv = panel.addPoints(0,0,cps[k],css[k])
      pv.setMarkColor(colors[k])
      #uv.setLineColor(Color.WHITE)
      pv.setLineStyle(PointsView.Line.NONE)
      pv.setMarkStyle(PointsView.Mark.FILLED_SQUARE)
      pv.setMarkSize(8)
  #panel.addColorBar()
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(backgroundColor)
  #frame.setFontSizeForPrint(8,240)
  #frame.setSize(470,1000)
  frame.setFontSize(24)
  frame.setSize(1400,700)
  frame.setVisible(True)
  if png and pngDir:
    png += "n"+str(int(10*nrms))
    png += "s"+str(int(10*strainMax))
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

def plot3(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          xyz=None,cells=None,skins=None,fbs=None,surf=None,smax=0.0,
          links=False,curve=False,trace=False,png=None):
  n1 = len(f[0][0])
  n2 = len(f[0])
  n3 = len(f)
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
      ipg.setClips(-3.0,3.0)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-3.0,3.0)
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
  if surf:
    tg = TriangleGroup(True, s3, s2, surf)
    states = StateSet()
    cs = ColorState()
    cs.setColor(Color.CYAN)
    states.add(cs)
    lms = LightModelState()
    lms.setTwoSide(True)
    states.add(lms)
    ms = MaterialState()
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE)
    ms.setSpecular(Color.WHITE)
    ms.setShininess(100.0)
    states.add(ms)
    tg.setStates(states);
    sf.world.addChild(tg)
  ipg.setSlices(232,63,0)
  if cbar:
    sf.setSize(987,700)
  else:
    sf.setSize(850,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  zscale = 0.3*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.5)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(0.0,0.05,-0.08))
  ov.setAzimuthAndElevation(25,45.0)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(137,1,pngDir+png+"cbar.png")
#############################################################################
# Run the function main on the Swing thread
import sys
class _RunMain(Runnable):
  def __init__(self,main):
    self.main = main
  def run(self):
    self.main(sys.argv)
def run(main):
  SwingUtilities.invokeLater(_RunMain(main)) 
run(main)
