from utils3 import *
setupForSubset("part")

pngDir = getPngDir()
pngDir = None
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count
f1,f2,f3 = s1.getFirst(),s2.getFirst(),s3.getFirst()
d1,d2,d3 = s1.getDelta(),s2.getDelta(),s3.getDelta()
l1,l2,l3 = s1.getLast(),s2.getLast(),s3.getLast()
k3,k2,k1 = 69,168,1200; azimuth=240; elevation=20 # for 3D views
fmin,fmax = -5.5,5.5

gxfile = "gx" # input seismic image
epfile = "ep" # planarity
p2file = "p2" # inline slope
p3file = "p3" # crossline slope
ggfile = "gg" # flattened seismic image
x1file = "x1" # horizon volume with slopes
c1file = "c1" # corrected horizon volume with correlations
gtfile = "gt" # relative geologic time
gifile = "gi" # interpolated volume
cpfile = "cp" # most positive curvature
cnfile = "cn" # most negative curvature

#log information
c2s = [670916,660065,659116,664891,661082,677668,
       676790,677017,682230,676987,671153,672222]	

c3s = [4530909,4532272,4533012,4532232,4530130,4533492,
       4531236,4533302,4533286,4532326,4529709,4530982]

wns =["YC1","YJ1-3","YJ1-5","YJ1-9X","YJ1X","YJ2-3","YJ2-7X",
      "YJ2-9","YJ2X","YJ3","YJ3-2","YJ3-3H"]
logDir = "../../../data/seis/tjxd/3d/logs/las/"
dvtDir = "../../../data/seis/tjxd/3d/logs/dvt/"
logType = "velocity"
wmin,wmax=1.8,6.0

plotOnly = True
plotOnly = False

def main(args):
  #goSubset()
  #goDisplay()
  #goFlatten()
  #goCorrection(4,2,1,20,0.5,0.5)
  #goCorrection(4,2,2,5,0.25,0.5)
  #goRgtInterp()
  goCurvature()

def goSubset():
  gx = readImage3D(gxfile)
  gs = copy(1401,n2,n3,100,0,0,gx)
  writeImage("gxSub",gs)
def goDisplay():
  gx = readImage3D(gxfile)
  gx = gain(gx)
  writeImage(gxfile,gx)
  samples=getLogSamples("velocity")
  plot3(gx,samples=samples)
  #plot3(gx)

def goSlopes():
  gx = readImage3D(gxfile)
  p2 = copy(gx)
  p3 = copy(gx)
  ep = copy(gx)
  lsf = LocalSlopeFinder(8.0,2.0)
  lsf.findSlopes(gx,p2,p3,ep);
  '''
  zm = ZeroMask(0.1,4,2,2,gx)
  zero = 0.00;
  tiny = 0.01;
  zm.setValue(zero,p2);
  zm.setValue(zero,p3);
  zm.setValue(tiny,ep);
  '''
  writeImage(p2file,p2)
  writeImage(p3file,p3)
  writeImage(epfile,ep)

def goFlatten():
  gx = readImage3D(gxfile)
  print min(gx)
  print max(gx)
  c1,c2,c3=Sampling(n1),Sampling(n2),Sampling(n3)
  if not plotOnly:
    p2 = zerofloat(n1,n2,n3)
    p3 = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    lsf = LocalSlopeFinder(8.0,4.0)
    lsf.findSlopes(gx,p2,p3,ep);
    ep = pow(ep,10)
    fl = Flattener3Dw()
    fl.setWeight1(0.2)
    fl.setIterations(0.01,300)
    fm = fl.getMappingsFromSlopes(c1,c2,c3,p2,p3,ep)
    gg = fm.flatten(gx)
    x1 = fm.x1
    writeImage(ggfile,gg)
    writeImage(x1file,x1)
  else:
    gg = readImage3D(ggfile)
  plot3(gx,clab="Amplitude",png="seis")
  plot3(gg,clab="Amplitude",png="flat")

def goCorrection(sig1,sig2,k,smax,strain,w1):
  if(k==1):
    gx = readImage3D(ggfile)
  else:
    gx = readImage3D("gc"+str(k))
  print min(gx)
  print max(gx)
  c1,c2,c3=Sampling(n1),Sampling(n2),Sampling(n3)
  if not plotOnly:
    p2 = zerofloat(n1,n2,n3)
    p3 = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    lsf = LocalSlopeFinder(sig1,sig2)
    lsf.findSlopes(gx,p2,p3,ep);
    wp = zerofloat(n2,n3)
    for i3 in range(n3):
      for i2 in range(n2):
        wp[i3][i2] = sum(ep[i3][i2])/n1
    gc = GlobalCorrelationFinder(-smax,smax)
    gc.setStrainMax(strain)
    ks = gc.getTraceIndexes(8,8,120,0.2,wp)
    ts = gc.findCorrelations(ks,gx)
    ep = pow(ep,10)
    fl = Flattener3Dw()
    fl.setWeight1(w1)
    fl.setIterations(0.01,200)
    fm = fl.getMappingsFromSlopesAndCorrelations(c1,c2,c3,0.001,p2,p3,ep,ks,ts)
    gg = fm.flatten(gx)
    xc1 = fm.x1
    if(k==1):
      x1 = readImage3D(x1file)
    else:
      x1 = readImage3D("xc"+str(k))
    fl.updateHorizonVolume(x1,xc1)
    writeImage("gc"+str(k+1),gg)
    writeImage("xc"+str(k+1),xc1)
  else:
    gg = readImage3D("gc"+str(k+1))
    xc = readImage3D("xc"+str(k+1))
  plot3(gx,clab="Amplitude (before)",png="seis")
  plot3(gg,clab="Amplitude (after)",png="flat")
  plot3(gx,g=xc,cmap=jetFill(1),cmin=1,cmax=n1)
 
def goCurvature():
  gx = readImage3D(gxfile)
  if not plotOnly:
    x1 = readImage3D("xc3")
    cv = Curvature()
    cn,cp = cv.horizonVolumetricCurvature(2,x1)
    writeImage(cnfile,cn)
    writeImage(cpfile,cp)
  else:
    cn = readImage3D(cnfile)
    cp = readImage3D(cpfile)
  plot3(gx,g=cn,cmap=bwrFill(1.0),cmin=-0.001,cmax=0.001)
  plot3(gx,g=cp,cmap=bwrFill(1.0),cmin=-0.001,cmax=0.001)

def goRgtInterp():
  if not plotOnly:
    x1 = readImage3D("xc3")
    fl = Flattener3Dw()
    gt = fl.rgtFromHorizonVolume(Sampling(n1),x1)
    ri = RgtInterpolator(s1,s2,s3,0.001)
    fx,x1,x2,x3=getLogSamples(logType)
    samples=fx,x1,x2,x3
    gi = ri.apply(fx,x1,x2,x3,gt)
    writeImage(gifile,gi)
  else:
    fx,x1,x2,x3=getLogSamples(logType)
    samples=fx,x1,x2,x3
    gi = readImage3D(gifile)
  print max(gi)
  print min(gi)
  gx = readImage3D(gxfile)
  #gi = copy(n1-400,n2,n3,0,0,0,gi)
  #gx = copy(n1-400,n2,n3,0,0,0,gx)
  #c1 = Sampling(n1-400,d1,f1)
  k3 = 72
  plot3(gx,s1=s1,k1=816,k3=k3,clab="Amplitude",png="seisSub")
  plot3(gx,s1=s1,samples=samples,k1=816,k3=k3,png="seisWellSub")
  plot3(gx,g=gi,s1=s1,samples=samples,k1=816,k3=k3,
        cmin=wmin,cmax=wmax,cint=0.5,cmap=jetFill(0.9),clab="Velocity",png="interp")
  
def getLogSamples(curve):
  k1,k2,k3,fk=[],[],[],[]
  wldata = WellLog.Data(logDir, dvtDir)
  logs = wldata.getLogsWith(curve)
  wi = 0
  for log in logs:
    print log.name
    #if(log.name=="YJ2-7X"):
    #  continue
    log.smooth(10)
    c2,c3 = getLogLocation(log.name)
    fxs = log.getSamplesX(curve)
    if(fxs!=None):
      fx,x1 = fxs[0],fxs[1]
      if(fx!=None and x1!=None):
        np = len(x1)
        x2 = fillfloat(c2,np)
        x3 = fillfloat(c3,np)
        k1.append(x1)
        k2.append(x2)
        k3.append(x3)
        fk.append(fx)
  return fk,k1,k2,k3

def getLogLocation(name):
  wi = 0
  for wni in wns:
    if(wni==name):
      return c2s[wi],c3s[wi]
    wi = wi+1

def gain(x):
    g = mul(x,x) 
    ref = RecursiveExponentialFilter(100.0)
    ref.apply1(g,g)
    y = zerofloat(n1,n2,n3)
    div(x,sqrt(g),y)
    return y
#############################################################################
# graphics
jet = ColorMap.JET
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
    cbar.setInterval(0.5)
  cbar.setFont(Font("Arial",Font.PLAIN,24)) # size by experimenting
  cbar.setWidthMinimum
  cbar.setBackground(Color.WHITE)
  frame.add(cbar,BorderLayout.EAST)
  return cbar

def convertDips(ft):
  return FaultScanner.convertDips(0.2,ft) # 5:1 vertical exaggeration

def makePointGroup(f,x1,x2,x3,cmin,cmax,cbar):
  n = len(x1)
  xyz = zerofloat(3*n)
  copy(n,0,1,x3,0,3,xyz)
  copy(n,0,1,x2,1,3,xyz)
  copy(n,0,1,x1,2,3,xyz)
  rgb = None
  if cmin<cmax:
    cmap = ColorMap(cmin,cmax,ColorMap.getJet(1.0))
    if cbar:
      cmap.addListener(cbar)
    rgb = cmap.getRgbFloats(f)
  pg = PointGroup(xyz,rgb)
  ps = PointState()
  ps.setSize(10)
  ps.setSmooth(False)
  ss = StateSet()
  ss.add(ps)
  pg.setStates(ss)
  return pg

def plot3(f,g=None,s1=s1,s2=s2,s3=s3,k1=1100,k2=0,k3=22,au=60,cmin=-4,cmax=4,
          cmap=None,clab=None,cint=None,curve=False,hs=None,samples=None,png=None):
  n3 = len(f)
  n2 = len(f[0])
  n1 = len(f[0][0])
  l1,l2,l3 = s1.last,s2.last,s3.last
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
      ipg.setClips(-4.0,4.0)
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
      cbar = addColorBar(sf,clab,0.1)
      ipg.addColorMap2Listener(cbar)
    sf.world.addChild(ipg)
  if cbar:
    cbar.setWidthMinimum(85)
  if samples:
    fx,x1,x2,x3 = samples
    vmin,vmax,vmap= wmin,wmax,ColorMap.JET
    nw = len(fx)
    for iw in range(nw):
      pgi = makePointGroup(fx[iw],x1[iw],x2[iw],x3[iw],vmin,vmax,None)
      sf.world.addChild(pgi)
  if hs:
    for hi in hs:
      if not curve:
        tg = TriangleGroup(True,hi[0],hi[1])
        sf.world.addChild(tg)
      else:
        lg = LineGroup(hi[0],hi[1])
        ss = StateSet()
        lg.setStates(ss)
        ls = LineState()
        ls.setWidth(4)
        ls.setSmooth(False)
        ss.add(ls)
        sf.world.addChild(lg)
  ipg.setSlices(k1,k2,k3)
  if cbar:
    sf.setSize(1051,850)
  else:
    sf.setSize(950,850)
  view = sf.getOrbitView()
  zscale = 0.8*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.1)
  view.setAzimuth(au)
  view.setElevation(25)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  view.setTranslate(Vector3(0.06,-0.00,0.0))
  sf.viewCanvas.setBackground(sf.getBackground())
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")
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
