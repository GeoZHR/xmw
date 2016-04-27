"""
Demo of surface reconstruction from fault cells/oriented points
Author: Xinming Wu and Dave Hale, Colorado School of Mines
Version: 2015.02.09
"""

from utils import *
setupForSubset("fake")
#setupForSubset("tp")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

# Names and descriptions of image files used below.
gxfile  = "gx" # input image (maybe after bilateral filtering)
rxfile  = "rx" # reflectivity image
rnfile  = "rn" # reflectivity image
pxfile  = "px" # impedance image
gsxfile = "gsx" # image after lsf with fault likelihoods
epfile  = "ep" # eigenvalue-derived planarity
p2file  = "p2" # inline slopes
p3file  = "p3" # crossline slopes
p2kfile = "p2k" # inline slopes (known)
p3kfile = "p3k" # crossline slopes (known)
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
u1file = "u1"
gffile = "gf"
gufile = "gu"
gsfile = "gs"
gcfile = "gc"
ghfile = "gh"
grfile = "gr"
dwfile = "dw"
grifile = "gri"

# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
minPhi,maxPhi = 0,145
minTheta,maxTheta = 65,85
sigmaPhi,sigmaTheta = 8,30

# These parameters control the construction of fault skins.
# See the class FaultSkinner for more information.
lowerLikelihood = 0.5
upperLikelihood = 0.7
minSkinSize = 3000

# These parameters control the computation of fault dip slips.
# See the class FaultSlipper for more information.
minThrow = 0.01
maxThrow = 15.0

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
plotOnly = False
pngDir = "../../../png/aii/fake/"
pngDir = None

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goFakeData()
  #goImpedance2()
  goImpedance3()

def goImpedance3():
  gx = readImage(gxfile)
  px = readImage(pxfile)
  rx = readImage(rnfile)
  smooth = 0.5 # for noisy data
  ep = fillfloat(1.0,n1,n2,n3)
  u1 = fillfloat(1.0,n1,n2,n3)
  u2 = fillfloat(1.0,n1,n2,n3)
  u3 = fillfloat(1.0,n1,n2,n3)
  lof = LocalOrientFilter(8.0,2.0)
  et3 = lof.applyForTensors(gx)
  lof.applyForNormalPlanar(gx,u1,u2,u3,ep)
  wp = pow(ep,6.0)
  et3.setEigenvalues(0.000001,1.0,1.0)
  ai3 = AcousticImpedanceInv3(6.0,6.0)
  ai3.setIterations(0.001,200)
  ai3.setTensors(et3)
  ai3.setSmoothness(smooth)
  pt = zerofloat(n1,n2,n3)
  k1,k2,k3,fp = getLogs(px)
  ai3.setInitial(pt,k1,k2,k3,fp)
  pi3 = ai3.applyForImpedance(copy(pt),rx,wp,k1,k2,k3,fp)
  samples = fp,k1,k2,k3
  print min(px)
  print max(px)
  plot3(gx,clab="Amplitude",png="seismic")
  plot3(rx,cmin=min(rx),cmax=max(rx),clab="Impedance",png="refx")
  plot3(gx,px,cmin=min(px),cmax=max(px),clab="Impedance",png="pTrue")
  plot3(gx,pt,cmin=min(px),cmax=max(px),clab="Impedance",
        samples=samples,png="pInitial")
  plot3(gx,pi3,cmin=min(px),cmax=max(px),clab="Impedance",
        samples=samples,png="pRecover01")

def goImpedance2():
  gx = readImage(gxfile)
  px = readImage(pxfile)
  rx = readImage(rnfile)
  #rx = readImage(rxfile)
  #smooth = 0.5 # for clean data
  rx = readImage(rnfile)
  smooth = 0.9 # for noisy data
  r3,p3,g3 = rx[50],px[50],gx[50]
  k1,k2,fp = setWells(p3)
  ep = fillfloat(1.0,n1,n2)
  u1 = fillfloat(1.0,n1,n2)
  u2 = fillfloat(1.0,n1,n2)
  lof = LocalOrientFilter(8.0,2.0)
  et2 = lof.applyForTensors(g3)
  lof.applyForNormalLinear(r3,u1,u2,ep)
  wp = pow(ep,6.0)
  et2.setEigenvalues(0.000001,1.0)
  ai2 = AcousticImpedanceInv2(6.0,6.0)
  ai2.setIterations(0.001,800)
  ai2.setTensors(et2)
  ai2.setSmoothness(smooth)
  pt = zerofloat(n1,n2)
  ai2.setInitial(u1,pt,k1,k2,fp)
  pi2 = ai2.applyForImpedance(r3,wp,k1,k2,fp)
  plot2(s1,s2,r3,clab="Reflectivity",cmin=-0.05,cmax=0.05,png="rx")
  plot2(s1,s2,p3,cmap=ColorMap.JET,clab="True impedance",
        cmin=min(p3),cmax=max(p3),png="p3")
  plot2(s1,s2,pt,cmap=ColorMap.JET,clab="Initial impedance",
        cmin=min(p3),cmax=max(p3),interp=False,png="pt")
  plot2(s1,s2,pi2,cmap=ColorMap.JET,clab="Recovered impedance",
        cmin=min(p3),cmax=max(p3),png="pi")

def getLogs(px):
  k2u = [81,56,20,8,96,50]
  k3u = [90,93,65,8,20,30]
  np = len(k2u)
  k1,k2,k3,fp = [],[],[],[]
  for ip in range(np):
    i2 = round(k2u[ip])
    i3 = round(k3u[ip])
    for i1 in range(n1):
      k1.append(i1)
      k2.append(i2)
      k3.append(i3)
      fp.append(px[i3][i2][i1])
  return k1,k2,k3,fp

def setWells3D(p3):
  k1 = zerofloat(n1*2)
  k2 = zerofloat(n1*2)
  k3 = zerofloat(n1*2)
  fp = zerofloat(n1*2)
  k = 0
  for i1 in range(n1):
    k1[k] = i1
    k2[k] = 15
    fp[k] = p3[15][i1]
    k = k+1
  for i1 in range(n1):
    k1[k] = i1
    k2[k] = 61
    fp[k] = p3[61][i1]
    k = k+1
  return k1,k2,fp

def setWells(p3):
  k1 = zerofloat(n1*2)
  k2 = zerofloat(n1*2)
  fp = zerofloat(n1*2)
  k = 0
  for i1 in range(n1):
    k1[k] = i1
    k2[k] = 15
    fp[k] = p3[15][i1]
    k = k+1
  for i1 in range(n1):
    k1[k] = i1
    k2[k] = 61
    fp[k] = p3[61][i1]
    k = k+1
  return k1,k2,fp

def goFakeData():
  #sequence = 'A' # 1 episode of faulting only
  sequence = 'OA' # 1 episode of folding, followed by one episode of faulting
  #sequence = 'OOOOOAAAAA' # 5 episodes of folding, then 5 of faulting
  #sequence = 'OAOAOAOAOA' # 5 interleaved episodes of folding and faulting
  nplanar = 0 # number of planar faults
  conjugate = False # if True, two large planar faults will intersect
  conical = False # if True, may want to set nplanar to 0 (or not!)
  impedance = False # if True, data = impedance model
  wavelet = True # if False, no wavelet will be used
  noise = 0.6 # (rms noise)/(rms signal) ratio
  gx,p2,p3 = FakeData.seismicAndSlopes3d2014A(
      sequence,nplanar,conjugate,conical,impedance,wavelet,noise)
  writeImage(gxfile,gx)
  writeImage(p2kfile,p2)
  writeImage(p3kfile,p3)
  print "gx min =",min(gx)," max =",max(gx)
  print "p2 min =",min(p2)," max =",max(p2)
  print "p3 min =",min(p3)," max =",max(p3)
  gmin,gmax,gmap = -3.0,3.0,ColorMap.GRAY
  if impedance:
    gmin,gmax,gmap = 0.0,1.4,ColorMap.JET
  plot3(gx,cmin=gmin,cmax=gmax,cmap=gmap,clab="Amplitude",png="gx")

def goSlopes():
  print "goSlopes ..."
  gx = readImage(gffile)
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

def like(x):
  n2 = len(x)
  n1 = len(x[0])
  return zerofloat(n1,n2)

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

def makePointGroup(f,x1,x2,x3,cmin,cmax,cbar):
  n = len(x1)
  xyz = zerofloat(3*n)
  copy(n,0,1,x3,0,3,xyz)
  copy(n,0,1,x2,1,3,xyz)
  copy(n,0,1,x1,2,3,xyz)
  rgb = None
  if cmin<cmax:
    cmap = ColorMap(cmin,cmax,ColorMap.getJet(0.3))
    if cbar:
      cmap.addListener(cbar)
    rgb = cmap.getRgbFloats(f)
  pg = PointGroup(xyz,rgb)
  ps = PointState()
  ps.setSize(4)
  ps.setSmooth(False)
  ss = StateSet()
  ss.add(ps)
  pg.setStates(ss)
  return pg

def plot1s(s1,ys,rs=None,vmin=None,vmax=None,color=Color.RED,
  hlabel="Seismic traces",vlabel="time (ms)",png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sf = 1.0
  yf = sf
  sp.setVLimits(0,n1)
  if vmin and vmax:
    sp.setVLimits(vmin,vmax)
  sp.setHLimits(0,len(ys))
  for il,y in enumerate(ys):
    ya = sum(y)/len(y)
    y = sub(y,ya)
    y = add(y,yf)
    pv = sp.addPoints(s1,y)
    pv.setLineColor(color)
    yf = yf+sf
  rf = sf
  if rs:
    for il,r in enumerate(rs):
      ra = sum(r)/len(r)
      r = sub(r,ra)
      r = add(r,rf)
      pv = sp.addPoints(s1,r)
      pv.setLineColor(Color.BLACK)
      rf = rf+sf
  sp.setSize(600,500)
  sp.setHLabel(hlabel)
  sp.setVLabel(vlabel)
  if png and pngDir:
    sp.paintToPng(300,7.0,pngDir+png+".png")

def plot2(s1,s2,x,cmap=ColorMap.GRAY,clab=None,cmin=0,cmax=0,
         title=None,interp=True,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  if title:
    sp.setTitle(title)
  sp.addColorBar(clab)
  sp.setSize(680,600)
  sp.plotPanel.setColorBarWidthMinimum(100)
  pv = sp.addPixels(s1,s2,x)
  pv.setColorModel(cmap)
  if interp:
    pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  else:
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if pngDir and png:
    sp.paintToPng(300,3.333,pngDir+png+".png")


def plot3X(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          slices=None, samples=None,png=None):
  n3 = len(f)
  n2 = len(f[0])
  n1 = len(f[0][0])
  s1,s2,s3=Sampling(n1),Sampling(n2),Sampling(n3)
  l1,l2,l3 = s1.last,s2.last,s3.last
  f1,f2,f3 = s1.first,s2.first,s3.first
  d1,d2,d3=s1.getDelta(),s2.getDelta(),s3.getDelta()
  sf = SimpleFrame(AxesOrientation.XRIGHT_YOUT_ZDOWN)
  cbar = None
  if g==None:
    ipg = sf.addImagePanels(s1,s2,s3,f)
    if cmap!=None:
      ipg.setColorModel(cmap)
    if cmin!=None and cmax!=None:
      ipg.setClips(cmin,cmax)
    else:
      #ipg.setClips(-2.0,2.0)
      ipg.setClips(-2.0,1.5) # use for subset plots
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-2.0,1.5)
    if cmin!=None and cmax!=None:
      ipg.setClips2(cmin,cmax)
    if cmap==None:
      cmap = jetFill(1.0)
    ipg.setColorModel2(cmap)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMap2Listener(cbar)
    sf.world.addChild(ipg)
  if cbar:
    cbar.setWidthMinimum(120)
  #ipg.setSlices(109,138,31)
  ipg.setSlices(2850,850,75) # for logs only
  if samples:
    fx,x1,x2,x3 = samples
    vmin,vmax,vmap= min(fx),max(fx)/2,ColorMap.JET
    pg = makePointGroup(fx,x1,x2,x3,vmin,vmax,None)
    sf.world.addChild(pg)
  if cbar:
    sf.setSize(887,700)
  else:
    sf.setSize(750,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  zscale = 0.5*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(3.5)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(-0.04,0.00,0.05))
  ov.setAzimuthAndElevation(130,35.0)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")

def plot3(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          xyz=None,cells=None,skins=None,smax=0.0,slices=None,
          links=False,curve=False,trace=False,htgs=None,
          uncs=None,samples=None,png=None):
  n3 = len(f)
  n2 = len(f[0])
  n1 = len(f[0][0])
  s1,s2,s3=Sampling(n1),Sampling(n2),Sampling(n3)
  d1,d2,d3=s1.getDelta(),s2.getDelta(),s3.getDelta()
  sf = SimpleFrame(AxesOrientation.XRIGHT_YOUT_ZDOWN)
  cbar = None
  if g==None:
    ipg = sf.addImagePanels(s1,s2,s3,f)
    if cmap!=None:
      ipg.setColorModel(cmap)
    if cmin!=None and cmax!=None:
      ipg.setClips(cmin,cmax)
    else:
      #ipg.setClips(-2.0,2.0)
      ipg.setClips(-2.0,1.5) # use for subset plots
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-2.0,1.5)
    if cmin!=None and cmax!=None:
      ipg.setClips2(cmin,cmax)
    if cmap==None:
      cmap = jetFill(1.0)
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
    xyz,uvw,rgb = FaultCell.getXyzUvwRgbForLikelihood(0.7,cmap,cells,False)
    qg = QuadGroup(xyz,uvw,rgb)
    qg.setStates(ss)
    sf.world.addChild(qg)
  if htgs:
    for htg in htgs:
      sf.world.addChild(htg)
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
      size = 0.65 
      ls = LineState()
      ls.setWidth(1.5)
      ls.setSmooth(True)
      ss.add(ls)
    ct = 0
    for skin in skins:
      if smax>0.0: # show fault throws
        cmap = ColorMap(0.0,smax,ColorMap.JET)
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
          r,g,b=0,0,0
        if ct==1:
          r,g,b=0,0,1
        if ct==2:
          r,g,b=0,1,1
        if ct==3:
          #r,g,b=0.627451,0.12549,0.941176
          r,g,b=1,1,1
        r,g,b=0,0,1
        xyz = skin.getCellLinksXyz()
        rgb = skin.getCellLinksRgb(r,g,b,xyz)
        lg = LineGroup(xyz,rgb)
        #lg = LineGroup(xyz)
        sg.addChild(lg)
        #ct = ct+1
    sf.world.addChild(sg)
  #ipg.setSlices(109,138,31)
  ipg.setSlices(109,138,7) # for logs only
  if uncs:
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
    ss.add(ms)
    sg.setStates(ss)
    us = UncSurfer()
    ul=readImage(ulfile)
    #ul = div(exp(ul),exp(1.0))
    for unc in uncs:
      [xyz,rgb]=us.buildTrigs(n1,s3,s2,-0.1,unc,ul)
      #[xyz,rgb]=us.buildTrigs(n1,s3,s2,0.01,unc,ul)
      tg  = TriangleGroup(True,xyz,rgb)
      sg.addChild(tg)
    sf.world.addChild(sg)
  if samples:
    fx,x1,x2,x3 = samples
    vmin,vmax,vmap= min(fx),max(fx),ColorMap.JET
    pg = makePointGroup(fx,x1,x2,x3,vmin,vmax,None)
    sf.world.addChild(pg)
  if cbar:
    sf.setSize(887,700)
  else:
    sf.setSize(750,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.48*sqrt(n1*n1+n2*n2+n3*n3)
  zscale = 0.80*max(n2*d2,n3*d3)/(n1*d1)
  ov = sf.getOrbitView()
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setWorldSphere(BoundingSphere(0.5*n1,0.4*n2,0.4*n3,radius))
  ov.setAzimuthAndElevation(120.0,25.0)
  #ov.setTranslate(Vector3(0.02,0.16,-0.27))
  ov.setTranslate(Vector3(0.02,0.23,-0.27))
  ov.setScale(1.25)
  # for subset plots
  #ov.setWorldSphere(BoundingSphere(0.5*n1,0.5*n2,0.5*n3,radius))
  #ov.setAzimuthAndElevation(-40.0,25.0)
  #ov.setTranslate(Vector3(0.0241,-0.0400,0.0103))
  #ov.setScale(1.3) #use only for subset plots
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")

#############################################################################
run(main)
