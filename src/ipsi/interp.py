"""
3D seismic image processing for faults
Author: Xinming Wu, Colorado School of Mines
Version: 2015.05.07
"""

from utils import *
setupForSubset("ufs")
#setupForSubset("unc")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

# Names and descriptions of image files used below.
gxfile = "gx" # input image
rgtfile = "rgt" #relative geologic time image
denfile = "tpgd" # array of density log data
tfile = "time" # array of density log data
pfile = "nearest" # array of density log data
qfile = "blended" # array of density log data
sx1file = "sx1"
sx2file = "sx2"
sx3file = "sx3"

# Directory for saved png images. If None, png images will not be saved.
pngDir = None
pngDir = "../../../png/ipsi/"
plotOnly = False

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  goInterp()
def goInterp():
  if not plotOnly:

    gx = readImage(gxfile)
    rgt = readImage(rgtfile)
    den = readImage(denfile)
    sx1 = readImage(sx1file)
    sx2 = readImage(sx2file)
    sx3 = readImage(sx3file)
    rgtx = zerofloat(n1,n2,n3)
    uf = UnfaultS(4.0,2.0)
    uf.applyShiftsX([sx1,sx2,sx3],rgt,rgtx)

    gp = GetPoints()
    x1,x2,x3,fx=wellTeapot(n1,150,150)
    samples = fx,x1,x2,x3
    rgi = RgtInterp3(fx,x1,x2,x3)
    rgi.setRgt(rgtx)
    rgi.setScales(0.01,1.00)
    t,p,q = rgi.grid(s1,s2,s3)
    '''
    tx = zerofloat(n1,n2,n3)
    px = zerofloat(n1,n2,n3)
    qx = zerofloat(n1,n2,n3)
    uf.applyShiftsX([sx1,sx2,sx3],t,tx)
    uf.applyShiftsX([sx1,sx2,sx3],p,px)
    uf.applyShiftsX([sx1,sx2,sx3],q,qx)
    '''
    writeImage(tfile,t)
    writeImage(pfile,p)
    writeImage(qfile,q)
  else:
    t = readImage(tfile)
    p = readImage(pfile)
    q = readImage(qfile)
  gx = gain(gx)
  k3 = str(249)
  plot3(gx,samples=samples,png="log"+k3)
  plot3(gx,rgtx,cmin=10.0,cmax=n1,cmap=jetFill(1.0),
        clab="Relative geologic time (samples)",png="rgt"+k3)
  plot3(gx,t,cmin=min(t),cmax=100,cmap=jetFill(0.5),
    clab="Time (samples)",png="time"+k3)
  plot3(gx,p,cmin=2.35,cmax=2.6,cmap=jetFill(0.5),
    clab="Density (g/cc)",samples=samples,png="nearest"+k3)
  plot3(gx,q,cmin=min(q),cmax=max(q),cmap=jetFill(0.5),
    clab="Blended values",png="blended"+k3)


def smoothF(x):
  fsigma = 4.0
  flstop = 0.9
  flt = fillfloat(0.0,n1,n2,n3)
  sigma1,sigma2,sigma3,pmax = 8.0,1.0,1.0,1.0
  p2,p3,ep = FaultScanner.slopes(sigma1,sigma2,sigma3,pmax,x)
  return FaultScanner.smooth(flstop,fsigma,p2,p3,flt,x)

def like(x):
  n3 = len(x)
  n2 = len(x[0])
  n1 = len(x[0][0])
  return zerofloat(n1,n2,n3)


def gain2(x,sigma):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(sigma)
  ref.apply1(g,g)
  #y = like(x)
  div(x,sqrt(g),x)
  return x

def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(20.0)
  ref.apply1(g,g)
  y = like(x)
  div(x,sqrt(g),y)
  return y

def array(x1,x2,x3=None,x4=None):
  if x3 and x4:
    return jarray.array([x1,x2,x3,x4],Class.forName('[[[F'))
  elif x3:
    return jarray.array([x1,x2,x3],Class.forName('[[[F'))
  else:
    return jarray.array([x1,x2],Class.forName('[[[F'))




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

def rgbFromHeight(h,r,g,b):
  n1 = len(h[0])
  n2 = len(h)
  ht = zerofloat(n1*n2)
  mp = ColorMap(-max(h),-min(h),ColorMap.JET)
  i = 0
  for i1 in range(n1):
    for i2 in range(n2):
      ht[i] = -h[i2][i1]
      i=i+1
  htRGB = mp.getRgbFloats(ht)
  i = 0
  for i1 in range(n1):
    for i2 in range(n2):
      r[i2][i1] = htRGB[i  ] 
      g[i2][i1] = htRGB[i+1] 
      b[i2][i1] = htRGB[i+2] 
      i = i+3

def makePointGroup(f,x1,x2,x3,cmin,cmax,cbar):
  n = len(x1)
  xyz = zerofloat(3*n)
  copy(n,0,1,x3,0,3,xyz)
  copy(n,0,1,x2,1,3,xyz)
  copy(n,0,1,x1,2,3,xyz)
  rgb = None
  if cmin<cmax:
    cmap = ColorMap(cmin,cmax,ColorMap.JET)
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

def plot3(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          slices=None, samples=None, png=None):
  n1,n2,n3 = s1.count,s2.count,s3.count
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
      cmap = jetFill(0.8)
    ipg.setColorModel2(cmap)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMap2Listener(cbar)
    sf.world.addChild(ipg)
  if cbar:
    cbar.setWidthMinimum(120)
  if samples:
    fx,x1,x2,x3=samples
    pg = makePointGroup(fx,x1,x2,x3,2.35,2.6,None)
    sf.world.addChild(pg)
  ipg.setSlices(117,40,220)
  ipg.setSlices(108,25,249)
  #ipg.setSlices(108,25,150)
  if cbar:
    sf.setSize(987,720)
  else:
    sf.setSize(850,720)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  zscale = 0.5*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.3)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(0.0,0.15,0.05))
  ov.setAzimuthAndElevation(-56.0,35.0)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")


#############################################################################
run(main)
