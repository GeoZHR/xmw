from utils2 import *


setupForSubset("seamSub2d")
pngDir = getPngDir()
s1,s2 = getSamplings()
n1,n2 = s1.count,s2.count

fxfile = "gx108" # input seismic

def main(args):
  pickWithEnvelope()
  pickWithSlatLikelihood()
  pickWithMultipleAttributes()

def pickWithEnvelope():
  gx = readImage(fxfile)
  gx = gain(gx)
  sp = SaltPicker2()
  pa = sp.applyForInsAmp(gx)
  w1 = n1
  w2 = n2
  pa = normalize(pa)
  c1 = [315,28,120,326,548,315]
  c2 = [ 20,50,380,215,150, 20]
  xu = sp.initialBoundary(1,c1,c2)
  xt = copy(xu)
  nr,dr = 65,1
  bs = sp.refine(nr,1,10,2,xu,pa)
  xu = sp.pickNext(10,1,5,2,xu[0],xu[1],pa)
  xu = sp.regridBoundary(1,xu[0],xu[1])
  f1,f2=s1.first,s2.first
  d1,d2=s1.delta,s2.delta
  np = len(xu[0])
  for ip in range(np):
    xu[0][ip] = xu[0][ip]*d1+f1
    xu[1][ip] = xu[1][ip]*d2+f2
  for ip in range(len(c1)):
    c1[ip] = c1[ip]*d1+f1
    c2[ip] = c2[ip]*d2+f2
  plot(gx,cmin=-2,cmax=2,w1=w1,w2=w2,png="seis")
  plot(gx,cmin=-2,cmax=2,pp=[c1,c2],w1=w1,w2=w2,png="initial")
  plot(pa,cmin=0.2,cmax=0.6,w1=w1,w2=w2,png="env")
  plot(gx,cmin=-2,cmax=2,xp=[xu[0],xu[1]],w1=w1,w2=w2,png="paPick")

def pickWithSlatLikelihood():
  c1 = [315,28,120,326,548,315]
  c2 = [ 20,50,380,215,150, 20]
  gx = readImage(fxfile)
  sp = SaltPicker2()
  w1 = n1
  w2 = n2
  sl = getSaltLikelihood(gx)
  xu = sp.initialBoundary(1,c1,c2)
  nr,dr = 65,1
  bs = sp.refine(nr,1,10,2,xu,sl)
  #xu = sp.pickNext(15,1,5,2,xu[0],xu[1],sl)
  xu = sp.regridBoundary(1,xu[0],xu[1])
  f1,f2=s1.first,s2.first
  d1,d2=s1.delta,s2.delta
  np = len(xu[0])
  for ip in range(np):
    xu[0][ip] = xu[0][ip]*d1+f1
    xu[1][ip] = xu[1][ip]*d2+f2
  for ip in range(len(c1)):
    c1[ip] = c1[ip]*d1+f1
    c2[ip] = c2[ip]*d2+f2
  plot(sl,cmin=0.2,cmax=0.6,w1=w1,w2=w2,png="saltLike")
  plot(gx,cmin=-2,cmax=2,xp=[xu[0],xu[1]],w1=w1,w2=w2,png="slPik")

def pickWithMultipleAttributes():
  gx = readImage(fxfile)
  sp = SaltPicker2()
  pa = sp.applyForInsAmp(gx)
  w1 = n1
  w2 = n2
  pa = normalize(pa)
  sl = getSaltLikelihood(gx)
  ps = div(add(mul(pa,2),sl),3.0)
  ps = normalize(ps)
  c1 = [315,28,120,326,548,315]
  c2 = [ 20,50,380,215,150, 20]
  xu = sp.initialBoundary(1,c1,c2)
  nr,dr = 65,1
  bs = sp.refine(nr,1,10,2,xu,ps)
  xu = sp.pickNext(15,1,5,2,xu[0],xu[1],ps)
  xu = sp.regridBoundary(1,xu[0],xu[1])
  f1,f2=s1.first,s2.first
  d1,d2=s1.delta,s2.delta
  np = len(xu[0])
  for ip in range(np):
    xu[0][ip] = xu[0][ip]*d1+f1
    xu[1][ip] = xu[1][ip]*d2+f2
  for ip in range(len(c1)):
    c1[ip] = c1[ip]*d1+f1
    c2[ip] = c2[ip]*d2+f2
  plot(ps,cmin=0.2,cmax=0.6,w1=w1,w2=w2,png="ps")
  plot(gx,cmin=-2,cmax=2,xp=[xu[0],xu[1]],w1=w1,w2=w2,png="psPik")
  '''
  opp = OptimalPathPicker(10,2)
  ft = opp.applyTransform(bs)
  m2,m1 = len(bs),len(bs[0])
  wht = opp.applyForWeight(ft)
  tms1 = zerofloat(m2,m1)
  tms2 = zerofloat(m2,m1)
  pik1 = opp.forwardPick(nr,wht,tms1)
  pik2 = opp.backwardPick(round(pik1[m2-1]),wht,tms2)
  x2 = rampfloat(0,1,m2)
  plot(bs,cmin=0.05,cmax=0.9,w1=180,w2=1200,png="bandMap")
  plot(bs,cmin=0.05,cmax=0.9,xp=[pik2,x2],w1=180,w2=1200,png="bandPik")
  '''
def getSaltLikelihood(gx):
  lof = LocalOrientFilter(4,4)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  e1 = zerofloat(n1,n2)
  e2 = zerofloat(n1,n2)
  lof.applyForNormalLinear(gx,u1,u2,el)
  for i2 in range(n2):
    for i1 in range(n1):
      if (isNaN(el[i2][i1])):
        el[i2][i1] = 1
  rgf = RecursiveGaussianFilterP(4)
  rgf.apply10(el,e1)
  rgf.apply01(el,e2)
  e1 = mul(e1,e1)
  e2 = mul(e2,e2)
  es = sqrt(add(e1,e2))
  es = normalize(es)
  return es


def isNaN(num):
    return num != num

def gain(x):
  g = mul(x,x) 
  n2 = len(x)
  n1 = len(x[0])
  ref = RecursiveExponentialFilter(50.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2)
  div(x,sqrt(g),y)
  return y

def normalize(x):
  x = sub(x,min(x))
  return div(x,max(x))

def smooth(sig1,sig2,sigs,x):
  n2 = len(x)
  n1 = len(x[0])
  lof = LocalOrientFilter(sig1,sig2)
  ets = lof.applyForTensors(x)
  ets.setEigenvalues(0.001,1.0)
  lsf = LocalSmoothingFilter()
  g = zerofloat(n1,n2)
  lsf.apply(ets,sigs,x,g)
  return g

#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET

def plot(f,xp=None,pp=None,xs=None,xu=None,nr=50,phi=None,v1=None,v2=None,
        cmin=None,cmax=None,w1=None,w2=None,clab=None,png=None): 
  orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT;
  panel = PlotPanel(1,1,orientation);
  #panel.setVInterval(0.2)
  l1,l2=s1.last,s2.last
  f1,f2=s1.first,s2.first
  panel.setHLabel("Crossline (km)")
  panel.setVLabel("Depth (km)")
  panel.setHLimits(0,f2,l2)
  panel.setVLimits(0,f1,l1)
  pxv = panel.addPixels(0,0,s1,s2,f);
  #pxv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pxv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pxv.setColorModel(ColorMap.GRAY)
  if cmin and cmax:
    pxv.setClips(cmin,cmax)
  else:
    pxv.setClips(min(f),max(f))
  #panel.setTitle("normal vectors")
  if phi:
    cv = panel.addContours(phi)
    cv.setContours([0])
    cv.setLineColor(Color.RED)
    cv.setLineWidth(1.0)
  if xp:
    ptv = panel.addPoints(0,0,xp[0],xp[1])
    ptv.setLineColor(Color.YELLOW)
    ptv.setLineWidth(3.0)
  if xu:
    np = len(xu[0])
    for ip in range(np):
      x1c = xu[0][ip]
      x2c = xu[1][ip]
      u1c = xu[2][ip]
      u2c = xu[3][ip]
      x1m = x1c-u1c*nr
      x2m = x2c-u2c*nr
      x1p = x1c+u1c*nr
      x2p = x2c+u2c*nr
      x1s = [x1m,x1c,x1p]
      x2s = [x2m,x2c,x2p]
      ptv2 = panel.addPoints(x1s,x2s)
      ptv2.setLineColor(Color.RED)
      ptv2.setLineWidth(1.0)
    if xp:
      ptv = panel.addPoints(0,0,xp[0],xp[1])
      ptv.setLineColor(Color.YELLOW)
      ptv.setLineWidth(3.0)
  if pp:
    ptvl = panel.addPoints(0,0,pp[0],pp[1])
    ptvl.setLineColor(Color.YELLOW)
    ptvl.setLineWidth(3.0)
    ptvp = panel.addPoints(0,0,pp[0],pp[1])
    ptvp.setLineStyle(PointsView.Line.NONE)
    ptvp.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE)
    #ptvp.setMarkStyle(PointsView.Mark.CROSS)
    ptvp.setMarkColor(Color.RED)
    ptvp.setMarkSize(12.0)
    ptvp.setLineWidth(3.0)
  if xs:
    for ip in range(len(xs)):
      ptv = panel.addPoints(xs[ip][0],xs[ip][1])
      ptv.setLineStyle(PointsView.Line.NONE)
      ptv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
      ptv.setMarkColor(Color.RED)
      ptv.setMarkSize(2.0)
  if(clab):
    cb = panel.addColorBar();
    cb.setLabel(clab)
  panel.setColorBarWidthMinimum(130)
  moc = panel.getMosaic();
  frame = PlotFrame(panel);
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setVisible(True);
  if w1 and w2:
    frame.setSize(w2,w1)
  else:
    frame.setSize(round(n2*0.8),round(n1*0.8))
  #frame.setSize(1190,760)
  frame.setFontSize(16)
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
