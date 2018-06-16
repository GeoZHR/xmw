#############################################################################
"""
Demo of dynamic warping for automatic picking
Author: Xinming Wu, University of Texas at Austin
Version: 2016.06.01
"""


from utils import * 
setupForSubset("path")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count
f1,f2,f3 = s1.getFirst(),s2.getFirst(),s3.getFirst()
d1,d2,d3 = s1.getDelta(),s2.getDelta(),s3.getDelta()
#############################################################################
gxfile = "gx1"
elfile = "el"
fpfile = "fp"
fefile = "fe"
ftfile = "ft"
flfile = "fl"
ptfile = "pt"
fvfile = "fv"
fvtfile = "fvt"
pngDir = getPngDir()
plotOnly = False

def main(args):
  goLinearity()
  goPath()

def goLinearity():
  gx = readImage(gxfile)
  el = zerofloat(n1,n2)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  lof = LocalOrientFilter(2,2)
  lof.applyForNormalLinear(gx,u1,u2,el)
  writeImage(elfile,el)
  el = pow(el,4)
  plot(gx,label="Amplitude",png="gx")
  plot(gx,sub(1,el),cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),label="1-linearity",png="el")

def goPath():
  gx = readImage(gxfile)
  el = readImage(elfile)
  el = pow(el,4)
  op = OptimalPath(50,77)
  gt = op.transpose(gx)
  et = op.transpose(el)
  plott(gt,sub(1,et),cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
          label="1-linearity",png="et")
  c1,c2=52,102
  op.setStrainMax(0.4)
  op.setPathSmoothing(2)
  op.setAttributeSmoothing(1)
  op.setControlPoint(c1,c2,et)
  ec = copy(et)
  p1 = op.findPath(c1,c2,et)
  p2 = rampfloat(0,1,n1)
  ps = [p1,p2]
  plot(gx,sub(1,el),cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
          label="1-linearity",png="el")
  plott(gt,sub(1,ec),cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
          label="1-linearity",png="ec")
  plot(gx,sub(1,el),ps=[p2,p1],cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
          label="1-linearity",png="elp")
  plot(gx,ps=[p2,p1],cmin=-2,cmax=2.0,label="Amplitude", png="gxp")
  plott(gt,sub(1,et),ps=ps,cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
          label="1-linearity",png="eap")

  plott(gt,sub(1,et),cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
          label="1-linearity",png="ea")

  plott(gt,sub(1,ec),ps=ps,cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
          label="1-linearity",png="etp")
  ef = zerofloat(n1,n2)
  es = zerofloat(n1,n2)
  ps = zerofloat(n1)
  pf = zerofloat(n1)
  for i1 in range(n1):
    i2 = round(p1[i1])
    ef[i2  ][i1] = 1-el[i2][i1]
    ef[i2+1][i1] = 1-el[i2][i1]
    pf[i1] = 1-el[i2][i1]
  rgf = RecursiveGaussianFilterP(20)
  rgf.apply0(pf,ps)
  for i1 in range(n1):
    i2 = round(p1[i1])
    es[i2  ][i1] = ps[i1]
    es[i2+1][i1] = ps[i1]
  plot(gx,ef,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),nearest=True,
          label="1-linearity",png="ef")
  plot(gx,es,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),nearest=True,
          label="1-linearity",png="es")

  ef = zerofloat(n2,n1)
  eb = zerofloat(n2,n1)

  op.accumulateForwardX(sub(1,ec),ef)
  print min(ec)
  print max(ec)
  op.accumulateBackwardX(sub(1,ec),eb)
  es = add(ef,eb)
  es = sub(es,sub(1,ec))
  op.normalizeErrors(es)
  op.setControlPoint(c1,c2,es)
  plott(gt,ef,cmin=0.25,cmax=100,cmap=jetRamp(1.0),
          label="1-linearity",png="eaf")
  plott(gt,es,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
          label="1-linearity",png="eaf")


def gain(x):
  n2 = len(x)
  n1 = len(x[0])
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(20.0)
  ref.apply(g,g)
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

def plot(f,g=None,ps=None,t=None,cmap=None,cmin=None,cmax=None,cint=None,
        label=None,nearest=False,png=None): 
  orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT;
  n1,n2=len(f[0]),len(f)
  s1,s2=Sampling(n1),Sampling(n2)
  panel = PlotPanel(1,1,orientation)#,PlotPanel.AxesPlacement.NONE)
  panel.setVInterval(20)
  panel.setHInterval(20)
  panel.setVLimits(0,n1-1)
  panel.setHLabel("Inline (traces)")
  panel.setVLabel("Depth (samples)")
  pxv = panel.addPixels(0,0,s1,s2,f);
  pxv.setColorModel(ColorMap.GRAY)
  pxv.setInterpolation(PixelsView.Interpolation.LINEAR)
  if g:
    pxv.setClips(-2,2)
  else:
    if cmin and cmax:
      pxv.setClips(cmin,cmax)
  if g:
    pv = panel.addPixels(s1,s2,g)
    if nearest:
      pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    else:
      pv.setInterpolation(PixelsView.Interpolation.LINEAR)
    pv.setColorModel(cmap)
    if cmin and cmax:
      pv.setClips(cmin,cmax)
  if ps:
    uv = panel.addPoints(0,0,ps[0],ps[1])
    uv.setLineColor(Color.CYAN)
    uv.setLineWidth(8)
  if label:
    cbar = panel.addColorBar(label)
    cbar.setInterval(0.1)
  panel.setColorBarWidthMinimum(70)
  moc = panel.getMosaic();
  frame = PlotFrame(panel);
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  #frame.setTitle("normal vectors")
  frame.setVisible(True);
  #frame.setSize(1400,700)
  frame.setSize(round(n2*5)+70,round(n1*5))
  frame.setFontSize(20)
  if pngDir and png:
    frame.paintToPng(720,3.333,pngDir+png+".png")
def plott(f,g=None,ps=None,t=None,cmap=None,cmin=None,cmax=None,cint=None,
        label=None,neareast=False,png=None): 
  orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT;
  n1,n2=len(f[0]),len(f)
  s1,s2=Sampling(n1),Sampling(n2)
  panel = PlotPanel(1,1,orientation)#,PlotPanel.AxesPlacement.NONE)
  panel.setVInterval(20)
  panel.setHInterval(20)
  panel.setHLimits(0,0,n2-1)
  panel.setHLabel("Depth (samples)")
  panel.setVLabel("Inline (traces)")
  pxv = panel.addPixels(0,0,s1,s2,f);
  pxv.setColorModel(ColorMap.GRAY)
  pxv.setInterpolation(PixelsView.Interpolation.LINEAR)
  if g:
    pxv.setClips(-2,2)
  else:
    if cmin and cmax:
      pxv.setClips(cmin,cmax)
  if g:
    pv = panel.addPixels(s1,s2,g)
    if neareast:
      pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    else:
      pv.setInterpolation(PixelsView.Interpolation.LINEAR)
    pv.setColorModel(cmap)
    if cmin and cmax:
      pv.setClips(cmin,cmax)
  if ps:
    uv = panel.addPoints(0,0,ps[0],ps[1])
    uv.setLineColor(Color.CYAN)
    uv.setLineWidth(5)
  if label:
    cbar = panel.addColorBar(label)
    cbar.setInterval(0.1)
  panel.setColorBarWidthMinimum(70)
  moc = panel.getMosaic();
  frame = PlotFrame(panel);
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  #frame.setTitle("normal vectors")
  frame.setVisible(True);
  #frame.setSize(1400,700)
  frame.setSize(round(n2*5)+70,round(n1*5))
  frame.setFontSize(20)
  if pngDir and png:
    frame.paintToPng(720,3.333,pngDir+png+".png")

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
