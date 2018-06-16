#############################################################################
"""
Demo of dynamic warping for automatic picking
Author: Xinming Wu, University of Texas at Austin
Version: 2016.06.01
"""


from utils import * 
setupForSubset("fake2d")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count
f1,f2,f3 = s1.getFirst(),s2.getFirst(),s3.getFirst()
d1,d2,d3 = s1.getDelta(),s2.getDelta(),s3.getDelta()
#############################################################################
gxfile = "gx2d02" # input semblance image
p2file = "p2" # inline slopes
p3file = "p3" # crossline slopes
epfile = "ep"  # planarity
c2file = "c2" # c2 coherence or semblance-based coherence
c2sfile = "c2s" # c2 coherence or semblance-based coherence
c3file = "c3" # c3 coherence or covariance-based coherence
smsfile = "sms" # structure-oriented semblance
epsfile = "eps" # structure-oriented planarity
effile = "ef"  # 1-planarity
fefile = "fe"  # 1-planarity
flfile = "fl"  # fault likelihood
fpfile = "fp"  # fault strike;
ftfile = "ft"  # fault dip;
fvfile = "fv"  # fault dip;
vpfile = "vp"  # fault dip;
vtfile = "vt"  # fault dip;
fetfile = "fet" # fault likelihood thinned
fptfile = "fpt" # fault strike thinned
fttfile = "ftt" # fault dip thinned
fskfile = "skin"

plotOnly = False
plotOnly = True
pngDir = None
pngDir = getPngDir()
# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
minTheta,maxTheta = 65,80
minPhi,maxPhi = 0,360
sigmaPhi,sigmaTheta=3,6

def main(args):
  #goSlopes()
  #goC2()
  #goC2s()
  #goC3()
  #goLinearity()
  #goStructureOrientedSemblance()
  #goStructureOrientedLinearity()
  #goFaultLikelihoodScan()
  #goFaultLikelihood()
  #goOptimalPathVoting()
  showPathes()
  #goCnnf()
  #goFaultOrientScan()
  #goSurfaceVoting()
  #goFaultSurfaces()
  #goFaultLikelihood()
  #goFaultSkins()
  #goVoters()

def goCnnf():
  gx = readImage2D(gxfile)
  fp = readImage2D("gxCnnf")
  plot(gx,fp,cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
        clab="CNN fault probability",png="cnnf")


# a simplified C2 coherence method
# (Marfurt, Kirlin, Farmer and Bahorich, 1998)
def goC2():
  sig1,sig2 = 2,1 
  gx = readImage2D(gxfile)
  ref1 = RecursiveExponentialFilter(sig1)
  ref2 = RecursiveExponentialFilter(sig2)
  ref1.setEdges(RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE)
  ref2.setEdges(RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE)
  gn = zerofloat(n1,n2)
  gd = zerofloat(n1,n2)
  # compute the numerator of coherence
  ref2.apply2(gx,gn)
  gn = mul(gn,gn)
  ref1.apply1(gn,gn)
  # compute the denominator of coherence
  gd = mul(gx,gx)
  ref2.apply2(gd,gd)
  ref1.apply1(gd,gd)
  c2 = div(gn,gd)
  c2 = pow(c2,8)
  #gx = copy(81,81,10,10,gx)
  #c2 = copy(81,81,10,10,c2)
  plot(gx,clab="Amplitude",png="gx")

  plot(gx,sub(1,c2),cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
        clab="1-coherence",png="c2")
def goC2s():
  sig1,sig2 = 4,1 
  gx = readImage2D(gxfile)
  lof = LocalOrientFilter(8,2)
  ets = lof.applyForTensors(gx)
  ets.setEigenvalues(0.0001,1.0)
  lsf = LocalSmoothingFilter()
  ref1 = RecursiveExponentialFilter(sig1)
  ref1.setEdges(RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE)
  gn = zerofloat(n1,n2)
  gd = zerofloat(n1,n2)
  # compute the numerator of coherence
  lsf.apply(ets,sig2,gx,gn)
  gn = mul(gn,gn)
  ref1.apply1(gn,gn)
  # compute the denominator of coherence
  lsf.apply(ets,sig2,mul(gx,gx),gd)
  ref1.apply1(gd,gd)
  c2s = div(gn,gd)
  c2s = pow(c2s,8)
  #gx = copy(81,81,10,10,gx)
  #c2s = copy(81,81,10,10,c2s)
  plot(gx,sub(1,c2s),cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
        clab="1-coherence",png="c2s")
# covariance-matrix-based semblance
# (Gersztenkorn and Marfurt, 1999)
def goC3():
  gx = readImage2D(gxfile)
  p2 = zerofloat(n1,n2)
  lsf = LocalSlopeFinder(8,2,5)
  lsf.findSlopes(gx,p2)
  cv = Covariance()
  em,es = cv.covarianceEigen(10,p2,gx)
  sm = div(em,es)
  sm = pow(sm,20)
  #gx = copy(81,81,10,10,gx)
  #sm = copy(81,81,10,10,sm)
  plot(gx,sub(1,sm),cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
        clab="1-coherence",png="sm")
def goLinearity():
  gx = readImage2D(gxfile)
  lof = LocalOrientFilter(3,1,1)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lof.applyForNormalLinear(gx,u1,u2,el)
  el = pow(el,8)
  #gx = copy(81,81,10,10,gx)
  #el = copy(81,81,10,10,el)
  plot(gx,sub(1,el),cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
      clab="1-linearity",png="el")

def goStructureOrientedSemblance():
  sigma1,sigma2 = 8,2
  gx = readImage2D(gxfile)
  lof = LocalOrientFilter(sigma1,sigma2)
  ets = lof.applyForTensors(gx)
  lsf = LocalSemblanceFilter(2,16)
  sms = lsf.semblance(LocalSemblanceFilter.Direction2.V,ets,gx)
  sms = pow(sms,8)
  #gx = copy(81,81,10,10,gx)
  #sms = copy(81,81,10,10,sms)
  plot(gx,sub(1,sms),cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
      clab="1-semblance",png="sms") 

def goStructureOrientedLinearity():
  sigma1,sigma2 = 8,2
  gx = readImage2D(gxfile)
  el = zerofloat(n1,n2)
  lof = LocalOrientFilter(sigma1,sigma2)
  ets = lof.applyForTensors(gx)
  sta = StructureTensorAttribute(ets,24)
  sta.setEigenvalues(1.0,0.001)
  sta.applyForLinear(gx,el)
  el = pow(el,8)
  #gx = copy(81,81,10,10,gx)
  #el = copy(81,81,10,10,el)
  plot(gx,sub(1,el),cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
      clab="1-planarity",png="els") 

def goFaultLikelihoodScan():
  gx = readImage2D(gxfile)
  gx = FaultScanner2.taper(10,0,gx)
  fs = FaultScanner2(20)
  sig1,sig2,smooth=16.0,2.0,4.0
  p2,ep = fs.slopes(sig1,sig2,5,gx)
  fls = fs.scanOverThetas(65,85,p2,gx)
  k=0
  for fli in fls:
    plot(gx,fli,cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
      clab="Falt likelihod",png="fl"+str(k)) 
    k = k+1


def goFaultLikelihood():
  gx = readImage2D(gxfile)
  gx = FaultScanner2.taper(10,0,gx)
  fs = FaultScanner2(15)
  sig1,sig2,smooth=16.0,2.0,4.0
  p2,ep = fs.slopes(sig1,sig2,5,gx)
  fl,ft = fs.scan(65,85,p2,gx)
  flt,ftt=fs.thin([fl,ft])
  #gx = copy(81,81,10,10,gx)
  #fl = copy(81,81,10,10,fl)
  #ft = copy(81,81,10,10,ft)
  plot(gx,fl,cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
      clab="Falt likelihod",png="fl") 
  plot(gx,flt,cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
      clab="Falt likelihod",png="flt") 
  plot(gx,ft,cmin=-85,cmax=85,cmap=jetFill(1.0),
      clab="Fault dip (degrees)",png="ft")

def goOptimalPathVoting():
  gx = readImage2D(gxfile)
  lof = LocalOrientFilter(3,1)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lof.applyForNormalLinear(gx,u1,u2,el)
  fos = FaultOrientScanner2(6)
  fe,fp = fos.scanDip(65,85,el)
  ef = pow(el,4)
  ef = sub(1,ef)
  ft,pt = fos.thin([fe,fp])
  osv = OptimalPathVoter(20,40)
  osv.setStrainMax(0.4)
  osv.setPathSmoothing(6)
  fv,w1,w2 = osv.applyVoting(3,0.4,ft,pt)
  #gx = copy(91,81,10,10,gx)
  #fv = copy(91,81,10,10,fv)
  plot(gx,fv,cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
      clab="Voting score",png="fv") 

def showPathes():
  gx = readImage2D(gxfile)
  lof = LocalOrientFilter(3,1)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lof.applyForNormalLinear(gx,u1,u2,el)
  fos = FaultOrientScanner2(6)
  fe,fp = fos.scanDip(65,85,el)
  ft,pt = fos.thin([fe,fp])
  ftc = copy(ft)
  ptc = copy(pt)
  osv = OptimalPathVoter(20,40)
  osv.setStrainMax(0.4)
  osv.setPathSmoothing(6)
  fes,xs = osv.applyVotingAndPathes(3,0.4,ft,pt)
  ft,pt = fos.thin([fe,fp])
  sds = osv.pickSeeds(3,0.4,ftc,ptc)
  sps = osv.seedToPoint(sds)
  el = pow(el,8)
  ep = sub(1,el)
  fas = osv.accumulate(fes)
  jk = 0
  p1s = []
  p2s = []
  ns = len(fas)
  for ik in range(len(fes)):
    fmax = max(fes[ik])
    if fmax>0:
      jk = jk+1
      p1s.append(sps[0][ik])
      p2s.append(sps[1][ik])
      '''
      plot(fes[ik],fmin=0,fmax=fmax,fmap=jetRamp(1.0),ps=[[sps[0][ik]],[sps[1][ik]]],
         title="Path "+str(jk), clab="Voting score",png="fe"+str(jk)) 
      plot(fas[ik],fmin=0,fmax=max(fas[ik]),fmap=jetRamp(1.0),
         title="Accumulated voting scores",clab="Voting score",png="fa"+str(jk)) 
      '''

  plot(ep,fmin=0.25,fmax=1.0,fmap=jetRamp(1.0),title="Input attribute map", 
          clab="1-linearity",png="el") 
  plot(zerofloat(n1,n2),fmin=0,fmax=1,fmap=jetRamp(1.0),ps=[p1s,p2s],
       title="Seed points", clab="Voting score",png="seeds") 
  plot(fas[ns-1],fmin=0,fmax=max(fas),fmap=jetRamp(1.0),
       clab="Voting score",png="fas") 


def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(10.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2,n3)
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

def plot(f,g=None,ps=None,t=None,fmap=ColorMap.GRAY,
        fmin=-2,fmax=2,
        cmap=None,cmin=None,cmax=None,cint=None,
        clab=None,title=" ",nearest=False,png=None): 
  orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT;
  n1,n2=len(f[0]),len(f)
  s1,s2=Sampling(n1),Sampling(n2)
  panel = PlotPanel(1,1,orientation,PlotPanel.AxesPlacement.NONE)
  if title:
    panel.addTitle(title)
  '''
  panel = PlotPanel(1,1,orientation)#,PlotPanel.AxesPlacement.NONE)
  panel.setHLabel("Inline (traces)")
  panel.setVLabel("Depth (samples)")
  panel.setVInterval(10)
  panel.setHInterval(10)
  '''

  panel.setHLimits(0,n2-1)
  #panel.setVLimits(0,n1-1)
  pxv = panel.addPixels(0,0,s1,s2,f);
  pxv.setColorModel(fmap)
  pxv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pxv.setClips(fmin,fmax)
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
    uv.setMarkStyle(PointsView.Mark.PLUS)
    uv.setMarkColor(Color.BLACK)
    uv.setMarkSize(12)
    uv.setLineWidth(4)
    uv.setLineStyle(PointsView.Line.NONE)
  if clab:
    panel.addColorBar(clab)
  panel.setColorBarWidthMinimum(65)
  moc = panel.getMosaic();
  frame = PlotFrame(panel);
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  #frame.setTitle("normal vectors")
  frame.setVisible(True);
  #frame.setSize(1400,700)
  frame.setSize(round(n2*4)+75,round(n1*4))
  frame.setFontSize(14)
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
