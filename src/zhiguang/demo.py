"""
Demo of surface reconstruction from fault cells/oriented points
Author: Xinming Wu and Dave Hale, Colorado School of Mines
Version: 2015.02.09
"""

from fakeutils import *
setupForSubset("zhiguang")
#setupForSubset("tp")
#setupForSubset("f3d")
s1,s2,s3 = getSamplings()
n1,n2 = s1.count,s2.count

# Names and descriptions of image files used below.
g1file = "data1" # input image (maybe after bilateral filtering)
g2file = "data2" # input image (maybe after bilateral filtering)
gwfile = "gw"
ssfile = "ss"

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
pngDir = None
plotOnly = False
pngDir = "../../../png/zhiguang/"

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goWarp()
  goTest()
def goTest():
  g = zerofloat(n1,n2)
  g1 = zerofloat(n1,n2)
  for i1 in range(20,1000,1):
    g[100][i1] = 1
    g[101][i1] = 1
    g[102][i1] = 1
  for i1 in range(1010,2000,1):
    g[100][i1] = 1
    g[101][i1] = 1
    g[102][i1] = 1
  rgf = RecursiveGaussianFilter(10)
  rgf.applyX1(g,g1)
  plot2(s1,s2,g,clab="g",cmin=0,cmax=1)
  plot2(s1,s2,g1,clab="g1",cmin=min(g1),cmax=max(g1))
def goWarp():
  g1 = readImage2(g1file)
  g2 = readImage2(g2file)
  g1 = clip(-0.0005,0.0005,g1)
  g2 = clip(-0.0005,0.0005,g2)
  g1 = mul(g1,10000)
  g2 = mul(g2,10000)
  smin,smax = -5.0,70.0
  r1min,r1max = -0.8,0.8
  r2min,r2max = -0.8,0.8
  dwk = DynamicWarpingR(smin,smax,s1,s2)
  dwk.setStrainLimits(r1min,r1max,r2min,r2max)
  dwk.setSmoothness(8,8)
  ss = dwk.findShifts(s1,g1,s1,g2)
  gw = dwk.applyShifts(s1,g2,ss)
  writeImage(gwfile,gw)
  writeImage(ssfile,ss)
  plot2(s1,s2,g1,clab="g1",cmin=-2,cmax=2,png="g1")
  plot2(s1,s2,g2,clab="g2",cmin=-2,cmax=2,png="g2")
  plot2(s1,s2,gw,clab="gw",cmin=-2,cmax=2,png="gw")
  plot2(s1,s2,ss,cmap=ColorMap.JET,
   clab="Shifts (samples)",cmin=min(ss),cmax=max(ss),png="ss")

def like2(x):
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
  y = like2(x)
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

def plot1s(s1,ss,ys,rs=None,vmin=None,vmax=None,color=Color.RED,
  hlabel="Seismic traces",vlabel="time (ms)",png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sf = 10.0
  yf = sf
  sp.setVLimits(0,n1)
  if vmin and vmax:
    sp.setVLimits(vmin,vmax)
  sp.setHLimits(5,5+len(ys)*10)
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
      pv = sp.addPoints(s1,r)
      pv.setLineColor(Color.BLACK)
      rf = rf+sf
  sp.setSize(600,500)
  sp.setHLabel(hlabel)
  sp.setVLabel(vlabel)
  if png and pngDir:
    sp.paintToPng(300,7.0,pngDir+png+".png")

def plot2(s1,s2,x,cmap=ColorMap.GRAY,clab=None,cmin=0,cmax=0,
         title=None,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  if title:
    sp.setTitle(title)
  sp.addColorBar(clab)
  sp.setSize(680,600)
  sp.plotPanel.setColorBarWidthMinimum(80)
  pv = sp.addPixels(s1,s2,x)
  pv.setColorModel(cmap)
  #pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if pngDir and png:
    sp.paintToPng(300,3.333,pngDir+png+".png")


def plot3(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          xyz=None,cells=None,skins=None,fbs=None,smax=0.0,
          links=False,curve=False,trace=False,png=None):
  n1 = len(f[0][0])
  n2 = len(f[0])
  n3 = len(f)
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
    xyz,uvw,rgb = FaultCell.getXyzUvwRgbForLikelihood(0.5,cmap,cells,False)
    qg = QuadGroup(xyz,uvw,rgb)
    qg.setStates(ss)
    sf.world.addChild(qg)
  if fbs:
    mc = MarchingCubes(s1,s2,s3,fbs)
    ct = mc.getContour(0.0)
    tg = TriangleGroup(ct.i,ct.x,ct.u)
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
  if skins:
    sg = Group()
    ss = StateSet()
    lms = LightModelState()
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
      size = 0.5 
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
        xyz = skin.getCellLinksXyz()
        lg = LineGroup(xyz)
        sg.addChild(lg)
    sf.world.addChild(sg)
  #ipg.setSlices(198,0,89)
  ipg.setSlices(198,0,58)
  if cbar:
    sf.setSize(837,600)
  else:
    sf.setSize(700,600)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  ov.setWorldSphere(BoundingSphere(0.5*n1,0.5*n2,0.5*n3,radius))
  ov.setAzimuthAndElevation(-55.0,25.0)
  ov.setTranslate(Vector3(0.03,0.33,0.15))
  ov.setScale(1.4)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(137,1,pngDir+png+"cbar.png")

def plot3p(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          xyz=None,cells=None,skins=None,fbs=None,smax=0.0,
          links=False,curve=False,trace=False,png=None):
  n1 = len(f[0][0])
  n2 = len(f[0])
  n3 = len(f)
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
    xyz,uvw,rgb = FaultCell.getXyzUvwRgbForLikelihood(0.5,cmap,cells,False)
    qg = QuadGroup(xyz,uvw,rgb)
    qg.setStates(ss)
    sf.world.addChild(qg)
  if fbs:
    mc = MarchingCubes(s1,s2,s3,fbs)
    ct = mc.getContour(0.0)
    tg = TriangleGroup(ct.i,ct.x,ct.u)
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
  if skins:
    sg = Group()
    ss = StateSet()
    lms = LightModelState()
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
      size = 0.5 
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
        xyz = skin.getCellLinksXyz()
        lg = LineGroup(xyz)
        sg.addChild(lg)
    sf.world.addChild(sg)
  #ipg.setSlices(198,0,89)
  ipg.setSlices(189,0,0)
  if cbar:
    sf.setSize(937,600)
  else:
    sf.setSize(800,600)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  ov.setAxesScale(1.0,1.0,1.3)
  ov.setWorldSphere(BoundingSphere(0.5*n1,0.5*n2,0.5*n3,radius))
  ov.setAzimuthAndElevation(50,35.0)
  ov.setTranslate(Vector3(-0.2,-0.3,0.3))
  ov.setScale(1.3)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(137,1,pngDir+png+"cbar.png")

#############################################################################
run(main)
