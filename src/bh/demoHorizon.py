"""
Demo of generating avo models
Author: Xinming Wu, University of Texas at Austin
Version: 2017.05.05
"""

from utils import *
setupForSubset("bh")
#setupForSubset("tp")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

# Names and descriptions of image files used below.
gxfile = "gx"
gdfile = "rou"
gsfile = "gxs"
dsfile = "gds"
gcfile = "gxc"
p2file = "p2"
p3file = "p3"
epfile = "ep"
flfile = "fl"
fpfile = "fp"
ftfile = "ft"
fltfile = "flt"
fptfile = "fpt"
fttfile = "ftt"
fskbase = "fsk"
gwfile = "gw"
fwfile = "fw"
gsxfile = "gsc"
sf1file = "sf1"

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
pngDir = None
pngDir = "../../../png/bh/"
plotOnly = True
def main(args):
  #goSlopes()
  goHorizonOne()
  #goScan()
  #goThin()
  #goSkin()
  #goSmooth()
  #goSlip()
  #goTopBottomHorizons()

def goDen():
  gx = readImage3D(n1,n2,n3,gxfile)
  gd = readImage3D(n1,n2,n3,gdfile)
  for i3 in range(n3):
    gt = copy(gd[i3])
    for i2 in range(n2):
      gd[i3][i2] = gt[n2-i2-1]
  ds = copy(550,n2,n3,780,0,0,gd)
  gs = copy(550,n2,n3,780,0,0,gx)
  writeImage("gds",ds)
  writeImage("gxs",gs)
  plot3(gs)
  plot3(gs,ds,cmin=2.2,cmax=2.7,cmap=jetRamp(0.4))

def goSlopes():
  print "goSlopes ..."
  gx = readImage(gsfile)
  p2 = zerofloat(n1,n2,n3)
  p3 = zerofloat(n1,n2,n3)
  ep = zerofloat(n1,n2,n3)
  sigma1,sigma2,sigma3,pmax = 8.0,2.0,2.0,5.0
  lsf = LocalSlopeFinder(sigma1,sigma2,sigma3,pmax)
  lsf.findSlopes(gx,p2,p3,ep)
  ep = pow(ep,6)
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

def goHorizonOne():
  k1 = [173,183,185,187,248,291,242,249,250,256,251,263,285,299,281,262,297, 313,340,348, 284,363,249,260,255,269,262]
  k2 = [134,242,355,218,261,348,446,608,699,853,908,777,623,618,803,853,873,1039,799,672,1031,630,759,714,809,640,661]
  k3 = [236,433,270,536,160, 61,295,327,327,432,306,324,125, 85,119,188,103, 103, 29, 24, 171,0.0,284,233,350,192,195]
  gx = readImage(gsfile)
  dx = readImage(dsfile)
  gh = GlobalHorizon3()
  dxc = copy(dx)
  if not plotOnly:
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    ep = readImage(epfile)
    ep = pow(ep,3)
    ep = sub(ep,min(ep))
    ep = div(ep,max(ep))
    lmt = n1-1
    gh.setWeights(0)
    gh.setSmoothings(8,8)
    gh.setCG(0.01,100)
    gh.setExternalIterations(20)
    sf1 = gh.surfaceInitialization(n2,n3,lmt,k1,k2,k3)
    sf1 = gh.surfaceUpdateFromSlopesAndCorrelations(5,20,gx,ep,p2,p3,k2,k3,sf1)
    writeImage(sf1file,sf1) 
  else:
    sf1 = readImage2D(n2,n3,sf1file)
  plot3(gx,cmap=gray,clab="Amplitude",png="seis")
  mp = ColorMap(-1,1,gray)
  r1,g1,b1 = gh.amplitudeRgb(mp,gx,sf1) 
  surf1 = [sf1,r1,g1,b1]
  plot3(gx,hz=surf1,cmap=gray,png="sfSeis")
  hp = Helper();
  hp.horizonToImage(4,sf1,dx)
<<<<<<< HEAD
  plot3(gx,g=dx,ks=[k1,k2,k3],cmin=2.2,cmax=2.7,cmap=jetRamp(0.6),png="surf1m")
  mp = ColorMap(2.2,2.7,jet)
  r2,g2,b2 = gh.amplitudeRgb(mp,dxc,sf1) 
  surf2 = [sf1,r2,g2,b2]
  plot3(gx,dxc,hz=surf2,cmin=2.2,cmax=2.7,cmap=jetRamp(0.4))
=======
  plot3(gx,g=dxc,cmin=2.2,cmax=2.7,cmap=jetRamp(1.0),clab="Density",png="den")
  mp = ColorMap(2.2,2.7,jet)
  r2,g2,b2 = gh.amplitudeRgb(mp,dxc,sf1) 
  surf2 = [sf1,r2,g2,b2]
  plot3(gx,dxc,hz=surf2,cmin=2.2,cmax=2.7,cmap=jetRamp(1.0),clab="Density",png="sfDen")
>>>>>>> bfeba1eb132e7d1054c11a88f8b5e8e06986dd85

def goTopBottomHorizons():
  gs = readImage3D(n1,n2,n3,"gxs")
  ds = readImage3D(n1,n2,n3,"gds")
  hp = Helper()
  hs = hp.getTopBottomHorizons(ds)
  gc = hp.getSection(hs,gs)
  writeImage(gcfile,gc)
  plot3(gs)
  plot3(gs,ds,cmin=2.2,cmax=2.7,cmap=jetRamp(1.0))
  plot3(gc,ds,cmin=2.2,cmax=2.7,cmap=jetRamp(1.0))


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
  ref = RecursiveExponentialFilter(100.0)
  ref.apply1(g,g)
  y = like(x)
  div(x,sqrt(g),y)
  return y

#############################################################################
# graphics
gray = ColorMap.GRAY
jet = ColorMap.JET
bwr = ColorMap.BLUE_WHITE_RED
rwb = ColorMap.RED_WHITE_BLUE
#############################################################################
# plotting
backgroundColor = Color.WHITE
cjet = ColorMap.JET
alpha = fillfloat(1.0,256); alpha[0] = 0.0
ajet = ColorMap.setAlpha(cjet,alpha)
def plotLogs(sz,sw,fx,wh=500,wv=800,cmin=None,cmax=None,
    hint=100,vlab="Depth (samples)",cbar=None,png=None):
  fz = sz.first
  dz = sz.delta
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(wh,wv)
  #sp.setHInterval(hint)
  sp.setVLabel(vlab)
  sp.setHLabel("Lateral position (samples)")
  sp.addColorBar(cbar)
  sp.plotPanel.setColorBarWidthMinimum(85)
  #sp.setHLimits(sw.first-sw.delta/2,sw.last)
  #sp.setVLimits(sz.first,sz.last)
  pv = sp.addPixels(sz,sw,fx)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cjet)
  n1 = len(fx[0])
  sp.setFontSize(20)
  if cmin and cmax:
    pv.setClips(cmin,cmax)
  if png and pngDir:
    sp.paintToPng(1080,3.33333,pngDir+png+".png")


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

def addTensorsInImage(ip,et,esize):
  tp = TensorsPanel(s1,s2,s3,et)
  tp.setEllipsoidSize(esize)
  ip.getFrame().addChild(tp)
  return tp

def addColorBar(frame,clab=None,cint=None):
  cbar = ColorBar(clab)
  if cint:
    cbar.setInterval(cint)
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
    cmap = ColorMap(cmin,cmax,ColorMap.getJet(0.3))
    if cbar:
      cmap.addListener(cbar)
    rgb = cmap.getRgbFloats(f)
  pg = PointGroup(xyz,rgb)
  ps = PointState()
  ps.setSize(8)
  ps.setSmooth(False)
  ss = StateSet()
  ss.add(ps)
  pg.setStates(ss)
  return pg

def setPointGroup(k1,k2,k3,size):
  np  = len(k1)
  xyz = zerofloat(np*3)
  rgb = zerofloat(np*3)
  ki = 0
  for i in range(np):
    xyz[ki  ] = k3[i]
    xyz[ki+1] = k2[i]
    xyz[ki+2] = k1[i]
    rgb[ki  ]  = 0#1/225 
    rgb[ki+1]  = 1#225/225 
    rgb[ki+2]  = 0#1/225 
    ki = ki+3
  pg = PointGroup(size,xyz,rgb);
  states = StateSet();
  cs = ColorState();
  cs.setColor(Color.GREEN);
  states.add(cs);
  lms = LightModelState();
  lms.setTwoSide(True);
  states.add(lms);
  ms = MaterialState();
  ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE);
  ms.setShininess(100.0);
  states.add(ms);
  pg.setStates(states);
  return pg;

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

def plot2(f,x1,x2,s,s1,s2,g=None,gmin=None,gmax=None,
                label=None,png=None,et=None):
  n2 = len(s)
  n1 = len(s[0])
  panel = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.LEFT_TOP)
  panel.setVLimits(0,500)
  panel.setHInterval(100)
  panel.setVInterval(100)
  panel.setHLabel("Lateral position (sample)")
  panel.setVLabel("Time (sample)")
  #panel.setHInterval(100.0)
  #panel.setVInterval(100.0)
  #panel.setHLabel("Pixel")
  #panel.setVLabel("Pixel")
  if label:
    panel.addColorBar(label)
  else:
    panel.addColorBar()
  panel.setColorBarWidthMinimum(90)
  pv = panel.addPixels(s1,s2,s)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setColorModel(ColorMap.GRAY)
  pv.setClips(-0.1,0.1)
  if g:
    alpha = 1.0
  else:
    g = zerofloat(s1.count,s2.count)
    alpha = 0.0
  pv = panel.addPixels(s1,s2,g)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setColorModel(ColorMap.getJet(alpha))
  if label and label[0]=="T":
    pv.setClips(0.0,1000.0)
    pv.setClips(min(g),max(g))
  else:
    if gmin:
      pv.setClips(gmin,gmax)
    else:
      pv.setClips(2.0,2.8)
  cmap = ColorMap(gmin,gmax,ColorMap.JET)
  fs,x1s,x2s = makePointSets(cmap,f,x1,x2)
  for i in range(len(fs)):
    #color = cmap.getColor((fs[i][0]-min(fs))/max(fs))
    color = cmap.getColor((fs[i][0]))
    #color = Color(color.red,color.green,color.blue)
    pv = panel.addPoints(x1s[i],x2s[i])
    pv.setLineStyle(PointsView.Line.NONE)
    pv.setMarkStyle(PointsView.Mark.FILLED_SQUARE)
    pv.setMarkSize(4)
    pv.setMarkColor(color)
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  #frame.setFontSizeForPrint(8,240)
  #frame.setSize(1240,774)
  #frame.setFontSizeForSlide(1.0,0.8)
  frame.setSize(1290,777)
  frame.setSize(900,777)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(400,3.2,pngDir+png+".png")
  return frame
def makePointSets(cmap,f,x1,x2):
  sets = {}
  for i in range(len(f)):
    if f[i] in sets:
      points = sets[f[i]]
      points[0].append(f[i])
      points[1].append(x1[i])
      points[2].append(x2[i])
    else:
      points = [[f[i]],[x1[i]],[x2[i]]] # lists of f, x1, x2
      sets[f[i]] = points
  ns = len(sets)
  fs = zerofloat(1,ns)
  x1s = zerofloat(1,ns)
  x2s = zerofloat(1,ns)
  il = 0
  for points in sets:
    fl = sets[points][0]
    x1l = sets[points][1]
    x2l = sets[points][2]
    nl = len(fl)
    fs[il] = zerofloat(nl)
    x1s[il] = zerofloat(nl)
    x2s[il] = zerofloat(nl)
    copy(fl,fs[il])
    copy(x1l,x1s[il])
    copy(x2l,x2s[il])
    il += 1
  return fs,x1s,x2s

def plot3(f,g=None,hz=None,ks=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          slices=None,png=None):
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
      ipg.setClips(-1,1) # use for subset plots
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-1,1)
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
    cbar.setWidthMinimum(100)
  if hz:
    tg=TriangleGroup(True,s3,s2,hz[0],hz[1],hz[2],hz[3])
    sf.world.addChild(tg)
  if ks:
    pg = setPointGroup(ks[0],ks[1],ks[2],10)
    sf.world.addChild(pg)
  ipg.setSlices(n1-1,85,n3-1)
  if cbar:
    sf.setSize(967,720)
  else:
    sf.setSize(850,720)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  zscale = 0.5*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.6)
  ov.setScale(1.4)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(-0.08,0.05,0.05))
  ov.setAzimuthAndElevation(-50.0,50.0)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")

#############################################################################
run(main)
