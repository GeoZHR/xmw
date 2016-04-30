"""
Demo of surface reconstruction from fault cells/oriented points
Author: Xinming Wu and Dave Hale, Colorado School of Mines
Version: 2015.02.09
"""

from utils import *
setupForSubset("f3dSub")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

# Names and descriptions of image files used below.
gxfile  = "gx" # input image (maybe after bilateral filtering)
rxfile  = "rx" # shallower reflectivity image
rxffile  = "rxf" # deeper reflectivity image
pxfile  = "px" # impedance image
gsxfile = "gsx" # image after lsf with fault likelihoods
epfile  = "ep" # eigenvalue-derived planarity
p2file  = "p2" # inline slopes
p3file  = "p3" # crossline slopes
p2kfile = "p2k" # inline slopes (known)
p3kfile = "p3k" # crossline slopes (known)
flfile  = "fli" # fault likelihood
fpfile  = "fpi" # fault strike (phi)
ftfile  = "fti" # fault dip (theta)
fltfile = "flit" # fault likelihood thinned
fptfile = "fpit" # fault strike thinned
fttfile = "ftit" # fault dip thinned
fskbase = "fsk"

# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
minPhi,maxPhi = 0,360
minTheta,maxTheta = 70,80
sigmaPhi,sigmaTheta = 10,20

# These parameters control the construction of fault skins.
# See the class FaultSkinner for more information.
lowerLikelihood = 0.3
upperLikelihood = 0.5
minSkinSize = 1000

# These parameters control the computation of fault dip slips.
# See the class FaultSlipper for more information.
minThrow = 0.0
maxThrow = 20.0


# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
plotOnly = False
pngDir = None
pngDir = "../../../png/aii/f3d/sub/"

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goLogs()
  #goTie()
  #goLinearity()
  #goImpedance()
  #goInitial() # display only
  #goSeisTracesAtWells()
  goFaults()

def goFaults():
  print "goSkin ..."
  gx = readImage(gxfile)
  if not plotOnly:
    fl = readImage(flfile)
    fp = readImage(fpfile)
    ft = readImage(ftfile)
    fs = FaultSkinner()
    fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
    fs.setMaxDeltaStrike(10)
    fs.setMaxPlanarDistance(0.2)
    fs.setMinSkinSize(minSkinSize)
    cells = fs.findCells([fl,fp,ft])
    flt = like(gx)
    FaultSkin.getLikelihoods(cells,flt)
    writeImage(fltfile,flt)
  else:
    flt = readImage(fltfile)
  '''
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="fls")
  '''

def goInitial():
  wps,frs=goTie()
  wpm = goWellSeisFit(wps,frs)
  x2 = [ 33, 84]
  x3 = [259,141]
  k1,k2,k3,fp = [],[],[],[]
  for i2 in range(2):
    m1 = len(wpm[i2])
    for i1 in range(1350,m1):
      k2.append(x2[i2])
      k3.append(x3[i2])
      k1.append(i1-1350)
      fp.append(exp(wpm[i2][i1]*2))
  '''
  ai3 = AcousticImpedanceInv3(8.0,8.0)
  pt = zerofloat(n1,n2,n3)
  ai3.setInitial(pt,k1,k2,k3,fp)
  plot3X(pt,pt,cmin=3500,cmax=5500,clab="Impedance",
        samples=samples,png="initial")
  '''
  samples = fp,k1,k2,k3
  rx = readImage(rxfile)
  gx = readImage(gxfile)
  gx = gain(gx)
  writeImage(gxfile,gx)
  plot3X(rx,cmin=-0.15,cmax=0.15,clab="Reflectivity",
        samples=samples,png="ref")
  print min(gx)
  print max(gx)
  plot3X(gx,clab="Amplitude", samples=samples,png="seis")

def goReflectivity():
  #rx = readImage(rxfile)
  gx = readImage(gxfile)
  gx = div(gx,10000)
  #plot3X(rx,cmin=-0.15,cmax=0.15,clab="Reflectivity", png="ref")
  plot3X(gx,cmin=-1.5,cmax=1.5,clab="Amplitude", png="seis")

def goImpSub():
  print "goImpedance..."
  m2 = 2
  smooth = 0.2
  wpm = readImage2D(n1,m2,'wpm')
  x2 = [ 33, 84]
  x3 = [159, 41]
  k1,k2,k3,fp = [],[],[],[]
  for i2 in range(m2):
    for i1 in range(n1):
      k1.append(i1)
      k2.append(x2[i2])
      k3.append(x3[i2])
      fp.append(wpm[i2][i1])
  if not plotOnly:
    rx = readImage(rxfile)
    rx = div(rx,2.5)
    gx = readImage(gxfile)
    #ep = readImage(epfile)
    #wp = pow(ep,0.0)
    wp = fillfloat(1.0,n1,n2,n3)
    lof = LocalOrientFilter(8.0,2.0)
    et3 = lof.applyForTensors(gx,True)
    print "tensor computation done..."
    et3.setEigenvalues(0.000001,1.0,1.0)
    ai3 = AcousticImpedanceInv3(8.0,8.0)
    ai3.setIterations(0.001,20)
    ai3.setTensors(et3)
    ai3.setSmoothness(smooth)
    pt = zerofloat(n1,n2,n3)
    ai3.setInitial(pt,k1,k2,k3,fp)
    px = ai3.applyForImpedance(pt,rx,wp,k1,k2,k3,fp)
    writeImage(pxfile,px)
  else:
    px = readImage(pxfile)
  samples = fp,k1,k2,k3
  '''
  plot3(gx,clab="Amplitude",png="seismic")
  plot3(rx,cmin=min(rx),cmax=max(rx),clab="Impedance",png="refx")
  plot3(gx,px,cmin=min(px),cmax=max(px),clab="Impedance",png="pTrue")
  plot3(gx,pt,cmin=min(px),cmax=max(px),clab="Impedance",
        samples=samples,png="pInitial")
  '''
  rx = readImage(rxfile)
  plot3X(rx,cmin=min(rx)/10,cmax=max(rx)/10,clab="Impedance",
        samples=samples)

  plot3X(px,px,cmin=min(fp),cmax=max(fp),clab="Impedance",
        samples=samples)

def goSeisTracesAtWells():
  rx = readImage(rxffile)
  frs = zerofloat(n1,2)
  x2s = [ 33, 84]
  x3s = [259,141]
  for k in range(2):
    x2 = x2s[k]
    x3 = x3s[k]
    frs[k] = rx[x3][x2]
  writeImage("frs",frs)

def goTie():
  wpt = readImage2D(2873,4,'logs')
  rxs = readImage2D(2897,2,'frs')
  wp1 = copy(2121,0,wpt[0])
  wp2 = copy(2473,0,wpt[3])
  wps = [wp1,wp2]
  frs = []
  wrs = []
  wrc = []
  wpc = []
  k = 0
  ss = []
  rgf = RecursiveGaussianFilterP(1.0)
  for wpi in wps:
    m1 = len(wpi)
    wri = zerofloat(m1)
    rgf.apply0(wpi,wpi)
    for i1 in range(m1-1):
      wri[i1] = 0.5*(log(wpi[i1+1])-log(wpi[i1]))
      wpi[i1] = 0.5*log(wpi[i1])
    wpi[m1-1] = 0.5*log(wpi[m1-1])
    s1 = Sampling(m1)
    ss.append(s1)
    fri = copy(m1,0,rxs[k])
    frs.append(fri)
    wrs.append(wri)
    dwk = DynamicWarpingK(8,-160,160,s1)
    dwk.setStrainLimits(-0.2,0.2)
    dwk.setSmoothness(4)
    rs = dwk.findShifts(s1,fri,s1,wri)
    wrci = dwk.applyShifts(s1,wri,rs)
    wpci = dwk.applyShifts(s1,wpi,rs)
    wpc.append(wpci)
    wrc.append(wrci)
    k=k+1
  plot1s(ss,wrs,rs=frs,color=Color.BLUE)
  plot1s(ss,wrc,rs=frs)
  return wpc,frs

def goWellSeisFit(wps,frs):
  wrm = []
  wpm = []
  wpc = []
  frc = []
  ss = []
  for k in range(2):
    fis = FitLogImpWithSeisRef(6.0)
    fis.setSeisBalance(0.8)
    wpk = wps[k]
    frk = frs[k]
    if (k==0):
      wpk = copy(1959,0,wps[k])
      frk = copy(1959,0,frs[k])
    wpc.append(wpk)
    frc.append(frk)
    wpmk = fis.fitImpedance(wpk,frk)
    m1 = len(wpmk)
    wrmk = zerofloat(m1)
    for i1 in range(m1-1):
      wrmk[i1] = wpmk[i1+1]-wpmk[i1]
    wrm.append(wrmk)
    wpm.append(wpmk)
    ss.append(Sampling(m1))
  plot1s(ss,wrm,rs=frc,color=Color.MAGENTA)
  plot1s(ss,wpm,rs=wpc,color=Color.MAGENTA)
  return wpm


def goLinearity():
  print "goLinearity..."
  gx = readImage(gxfile)
  ep = fillfloat(1.0,n1,n2,n3)
  u1 = fillfloat(1.0,n1,n2,n3)
  u2 = fillfloat(1.0,n1,n2,n3)
  u3 = fillfloat(1.0,n1,n2,n3)
  lof = LocalOrientFilter(8.0,2.0)
  lof.applyForNormalPlanar(gx,u1,u2,u3,ep)
  writeImage(epfile,ep)

def goImpedance():
  print "goImpedance..."
  if not plotOnly:
    m2 = 4
    smooth = 0.2
    wpm = readImage2D(n1,m2,'wpm')
    x2 = [ 33,545,704, 84]
    x3 = [259,619,339,141]
    k1,k2,k3,fp = [],[],[],[]
    for i2 in range(m2):
      for i1 in range(n1):
        k1.append(i1)
        k2.append(x2[i2])
        k3.append(x3[i2])
        fp.append(wpm[i2][i1])
    rx = readImage(rxfile)
    rx = div(rx,2.5)
    gx = readImage(gxfile)
    wp = fillfloat(1.0,n1,n2,n3)#readImage(epfile)
    lof = LocalOrientFilter(8.0,2.0)
    et3 = lof.applyForTensors(gx,True)
    print "tensor computation done..."
    et3.setEigenvalues(0.000001,1.0,1.0)
    ai3 = AcousticImpedanceInv3(8.0,8.0)
    ai3.setIterations(0.001,max(n2,n3))
    ai3.setTensors(et3)
    ai3.setSmoothness(smooth)
    pt = zerofloat(n1,n2,n3)
    ai3.setInitial(pt,k1,k2,k3,fp)
    px = ai3.applyForImpedance(pt,rx,wp,k1,k2,k3,fp)
    px = mul(px,2.0)
    px = exp(px)
    writeImage(pxfile,px)
  else:
    m2 = 4
    px = readImage(pxfile)
    wpm = readImage2D(n1,m2,'wpm')
    x2 = [ 33,545,704, 84]
    x3 = [259,619,339,141]
    k1,k2,k3,fp = [],[],[],[]
    for i2 in range(m2):
      for i1 in range(n1):
        k1.append(i1)
        k2.append(x2[i2])
        k3.append(x3[i2])
        fp.append(exp(wpm[i2][i1]*2))
  '''
  plot3(gx,clab="Amplitude",png="seismic")
  plot3(gx,px,cmin=min(px),cmax=max(px),clab="Impedance",png="pTrue")
  plot3(gx,pt,cmin=min(px),cmax=max(px),clab="Impedance",
        samples=samples,png="pInitial")
  '''
  samples = fp,k1,k2,k3
  plot3X(px,px,cmin=3000,cmax=5400,clab="Impedance",
        samples=samples,png="pRecover02")

def goLogs():
  rx = readImage(rxfile)
  print min(rx)
  print max(rx)
  wp,k1,k2,k3,wps,wrs = getF3dLogs()
  r1 = rx[259][ 33]
  r2 = rx[141][ 84]
  frs = [r1,r2]
  #plot1s(s1,wrs,rs=frs)
  samples=wp,k1,k2,k3
  plot3X(rx,cmin=min(rx)/10,cmax=max(rx)/10,samples=samples)

def getF3dLogs():
  m2 = 4 # number of logs
  m1 = 2873 # number of samples for each log
  x2 = [ 33, 84]
  x3 = [259,141]
  lgt = readImage2D(m1,m2,'logs')
  lgs = zerofloat(m1,2)
  lgs[0] = lgt[0]
  lgs[1] = lgt[3]
  k1 = []
  k2 = []
  k3 = []
  fp = []
  wrs = zerofloat(m1,2)
  wps = zerofloat(m1,2)
  for i2 in range(2):
    for i1 in range(m1-1):
      if(lgs[i2][i1]!=-999.25):
        k1.append(i1+0)
        k2.append(x2[i2])
        k3.append(x3[i2])
        fp.append(0.5*log(lgs[i2][i1]))
        if(lgs[i2][i1+1]!=-999.25):
          wps[i2][i1] = 0.5*log(lgs[i2][i1])
          wrs[i2][i1] = 0.5*(log(lgs[i2][i1+1])-log(lgs[i2][i1]))
  return fp,k1,k2,k3,wps,wrs

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
  ref = RecursiveExponentialFilter(80.0)
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

def plot1s(ss,ys,rs=None,vmin=None,vmax=None,color=Color.RED,
  hlabel="Seismic traces",vlabel="time (ms)",png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sf = 1.0
  yf = sf
  sp.setVLimits(0,2500)
  if vmin and vmax:
    sp.setVLimits(vmin,vmax)
  sp.setHLimits(0,len(ys)+1)
  for il,y in enumerate(ys):
    ya = sum(y)/len(y)
    y = sub(y,ya)
    y = add(y,yf)
    pv = sp.addPoints(ss[il],y)
    pv.setLineColor(color)
    yf = yf+sf
  rf = sf
  if rs:
    for il,r in enumerate(rs):
      ra = sum(r)/len(r)
      r = sub(r,ra)
      r = add(r,rf)
      pv = sp.addPoints(ss[il],r)
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
    cbar.setWidthMinimum(140)
  #ipg.setSlices(109,138,31)
  ipg.setSlices(1396,20,128) # for logs only
  if samples:
    fx,x1,x2,x3 = samples
    #vmin,vmax,vmap= min(fx),max(fx),ColorMap.JET
    vmin,vmax,vmap= 3000,6500,ColorMap.JET
    pg = makePointGroup(fx,x1,x2,x3,vmin,vmax,None)
    sf.world.addChild(pg)
  if cbar:
    sf.setSize(907,700)
  else:
    sf.setSize(750,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  zscale = 0.5*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(3.6)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(-0.06,0.005,0.015))
  ov.setAzimuthAndElevation(50,36.0)
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
  ipg.setSlices(1396,20,128) # for logs only
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
    sf.setSize(907,700)
  else:
    sf.setSize(750,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  zscale = 0.5*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(3.6)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(-0.06,0.005,0.015))
  ov.setAzimuthAndElevation(50,36.0)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")

#############################################################################
run(main)
