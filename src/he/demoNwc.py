from utils import *

#pngDir = None
setupForSubset("nwc")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

n1,n2,n3 = s1.count,s2.count,s3.count
d1,d2,d3 = s1.delta,s2.delta,s3.delta
#k1,k2,k3 = 88,60,160; azimuth=285; elevation=11 # for 3D view of all horizons
k1,k2,k3 = n1-1,950,540; azimuth=240; elevation=25 # for 3D view of strips
fmin,fmax = -5.5,5.5
k1f,k2f,k3f = 65,406,114
k1f,k2f,k3f = 48,406,114
k1f,k2f,k3f = 48,406,0
gmin,gmax,gint,glab = -2.0,2.0,0.5,"Amplitude"
background = Color.WHITE

pngDir = "../../../png/he/"
pngDir = None
gxfile = "gx"
p2file = "p2"
p3file = "p3"
epfile = "ep"
fskfile = "fsr"
ssfile = "surfs"
srfile = "surfsr"
smfile = "surfm"
mrfile = "surfmr"

plotOnly = False

def main(args):
  #slopes()
  #oneControlPoint()
  goSurfaceRefine()
  #multipleControlPoints()
def slopes():
  gx = readImage(gxfile)
  p2 = zerofloat(n1,n2,n3)
  p3 = zerofloat(n1,n2,n3)
  ep = zerofloat(n1,n2,n3)
  #good to use sigma2>sigma1 for extracting sequence boundaries
  sigma1,sigma2=4.0,2.0
  lsf = LocalSlopeFinder(sigma1,sigma2,sigma2,5) 
  lsf.findSlopes(gx,p2,p3,ep);
  writeImage(p2file,p2)
  writeImage(p3file,p3)
  writeImage(epfile,ep)
  plot3(gx)
  plot3(gx,p2,cmin=-2,cmax=2,cmap=jetFill(1.0))
  plot3(gx,p3,cmin=-2,cmax=2,cmap=jetFill(1.0))
  plot3(gx,ep,cmin=0.1,cmax=0.9,cmap=jetRamp(1.0))

def oneControlPoint():
  k11 = [453.0]
  k12 = [349.0]
  k13 = [147.0]
  if not plotOnly:
    gx  = readImage(gxfile)
    p2  = readImage(p2file)
    p3  = readImage(p3file)
    ep  = readImage(epfile)
    ep = pow(ep,6.0) 
    horizonExtraction(gx,p2,p3,ep,k11,k12,k13,ssfile)
  pg = setPointGroup(k11,k12,k13,9.0)
  displayHorizon(pg,ssfile)
  
def multipleControlPoints():
  k11 = [100, 43, 35, 91, 39, 38, 82, 76, 47, 76, 
          86, 57, 39, 37, 35,106, 58,101, 39,  6]
  k12 = [335,706,832,624,945,920,620,620,650,640,
          635,519,875,821,950,370,556,365,768,940]
  k13 = [433,200,495,  0,353,  9, 95,165,286,120, 
          22,547, 26,150,168,280,500,380,200,530]
  gx  = readImage(gxfile)
  p2  = readImage(p2file)
  p3  = readImage(p3file)
  ep  = readImage(epfile)
  wp = pow(ep,10.0) 
  horizonExtraction(gx,p2,p3,wp,k11,k12,k13,smfile)
  pg = setPointGroup(k11,k12,k13,10.0)
  displayHorizon(pg,smfile)


def goSurfaceRefine():
  k11 = [453.0]
  k12 = [349.0]
  k13 = [147.0]
  d1 = .2
  c1 = 60
  m1 = c1*2+1
  if not plotOnly:
    sk = readSkins(fskfile)
    gx = readImage(gxfile)
    sf = readImage2(srfile)
    srd = SurfaceRefinerDp()
    srd.setStrainMax(0.25,0.25)
    srd.setErrorSmoothing(3)
    #fm,em = srd.getErrorMatrix(m1,d1,sk,gx,sf)
    em = srd.getErrorMatrix(m1,d1,gx,sf)
    sfr = zerofloat(n2,n3)
    #srd.findSurface(fm,em,sfr)
    srd.findSurface(em,sfr)
    sfr = mul(sfr,d1)
    sfr = sub(sfr,c1*d1)
    sfr = add(sf,sfr)
    writeImage(srfile,sfr) 
  pg1 = setPointGroup(k11,k12,k13,9.0)
  pg2 = setPointGroup(k11,k12,k13,9.0)
  displayHorizon(pg1,ssfile)
  displayHorizon(pg2,srfile)

def horizonExtraction(f,p2,p3,wp,k11,k12,k13,filename):
  lmt = n1-1
  se = SurfaceExtractorC()
  #k11=se.refineConstraints(k11,k12,k13,f)
  se.setWeights(0.0)
  se.setSmoothings(4.0,4.0)
  se.setCG(0.01,100)
  surf = se.surfaceInitialization(n2,n3,lmt,k11,k12,k13)
  se.surfaceUpdateFromSlopes(wp,p2,p3,k11,k12,k13,surf)
  writeImage(filename,surf) 

def display(filename):
  f = readImage(filename)
  world = World()
  ipg = addImageToWorld(world,f,cmap=gray)
  ipg.setSlices(k1,k2,k3)
  makeFrame(world)

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

def rgbFromAmplitude(f,h,r,g,b):
  amp = zerofloat(n2*n3)
  si = SincInterpolator()
  #si.setUniform(n1,1,0,n2,1,0,n3,1,0,f)
  i = 0
  for i3 in range(n3):
    for i2 in range(n2):
      amp[i] = si.interpolate(n1,1,0,n2,1,0,n3,1,0,f,h[i3][i2],i2,i3)
      i = i+1
  aMin,aMax = -2.5,2.5
  mp = ColorMap(aMin,aMax,ColorMap.RED_WHITE_BLUE)
  ampRGB = mp.getRgbFloats(amp)
  i = 0
  for i3 in range(n3):
    for i2 in range(n2):
      r[i3][i2] = ampRGB[i  ] 
      g[i3][i2] = ampRGB[i+1] 
      b[i3][i2] = ampRGB[i+2] 
      i = i+3

def displayHorizon(pg,filename):
  f = readImage(gxfile)
  f = gain(f)
  h = readImage2(filename) 
  r = zerofloat(n2,n3)
  g = zerofloat(n2,n3)
  b = zerofloat(n2,n3)
  #rgbFromHeight(h,r,g,b)
  rgbFromAmplitude(f,h,r,g,b)
  world = World()
  print min(f)
  print max(f)
  ipg = addImageToWorld(world,f,cmap=rwb,cmin=-2.5,cmax=2.5)
  ipg.setSlices(k1,k2,k3)
  tg  = TriangleGroup(True,s3,s2,add(h,0.0),r,g,b)
  world.addChild(pg)
  world.addChild(tg)
  makeFrame(world,png=filename)
def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(20.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2,n3)
  div(x,sqrt(g),y)
  return y

#############################################################################
def slice12(k3,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n1,n2)
  SimpleFloat3(f).get12(n1,n2,0,0,k3,s)
  return s

def slice13(k2,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n1,n3)
  SimpleFloat3(f).get13(n1,n3,0,k2,0,s)
  return s

def slice23(k1,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n2,n3)
  SimpleFloat3(f).get23(n2,n3,k1,0,0,s)
  return s


#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET
bwr = ColorMap.BLUE_WHITE_RED
rwb = ColorMap.RED_WHITE_BLUE

def addImageToWorld(world,image,cmap=gray,cmin=0,cmax=0):
  ipg = ImagePanelGroup(s1,s2,s3,image)
  ipg.setColorModel(cmap)
  if cmin<cmax:
    ipg.setClips(cmin,cmax)
  world.addChild(ipg)
  return ipg

def addImage2ToWorld(world,image1,image2):
  ipg = ImagePanelGroup2(s1,s2,s3,image1,image2)
  ipg.setColorModel1(ColorMap.getGray())
  ipg.setColorModel2(ColorMap.getJet())
  #ipg.setColorModel2(ColorMap.getHue(0.0,20.0,0.3))
  world.addChild(ipg)
  return ipg

def makeFrame(world,png=None):
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  #lightPosition=[-0.18,-0.4,0.8,0.0] # good for horizons 1 and 2
  lightPosition=[0.2,0.0,0.30,0.0] # good for horizons 1 and 2
  #lightPosition=[0.,0.,1.0,0.0] #default position
  frame = SimpleFrame(world)
  view = frame.getOrbitView()
  #view.setLightPosition(lightPosition)
  zscale = 0.5*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.1)
  view.setAzimuth(azimuth)
  view.setElevation(elevation)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  #frame.viewCanvas.setBackground(frame.getBackground())
  frame.setSize(1000,900)
  frame.setVisible(True)
  if png:
    frame.paintToFile(png+".png")
  return frame
"""
  frame = SimpleFrame(world)
  view = frame.getOrbitView()
  zscale = 0.75*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.3)
  #view.setAzimuth(75.0)
  #view.setAzimuth(-75.0)
  view.setAzimuth(-65.0)
"""
 
def display2(k3,s,png=None):
  pp = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.LEFT_BOTTOM)
  pp.setHInterval(10.0)
  pp.setVInterval(10.0)
  pp.setHLabel("Crossline (km)")
  pp.setVLabel("Inline (km)")
  pv = pp.addPixels(s1,s2,slice12(k3,s))
  pv.setClips(fmin,fmax)
  pf = PlotFrame(pp)
  pf.setFontSizeForPrint(6.0,480)
  pf.setSize(926,510)
  pf.setVisible(True)
  if png and pngDir:
    pf.paintToPng(300,6,pngDir+png+".png")
 
def display3(s,c=None,clabel="",cmin=0,cmax=0,png=None):
  pp = PlotPanelPixels3(
    PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT,
    PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM,
    s1,s2,s3,s)
  pp.setSlices(k1,k2,k3)
  pp.setLabel1("Time (s)")
  pp.setLabel2("Crossline (km)")
  pp.setLabel3("Inline (km)")
  pp.setClips(fmin,fmax)
  if c:
    cb = pp.addColorBar(clabel)
    #cb.setInterval(1.0)
    pp.setColorBarWidthMinimum(140)
    pp.setLineColor(Color.BLACK)
  else:
    pp.setLineColor(Color.YELLOW)
    #cb = pp.addColorBar("Amplitude")
    #cb.setInterval(5.0)
  pp.setInterval1(0.5)
  pp.setInterval2(2.0)
  pp.setInterval3(2.0)
  pp.mosaic.setHeightElastic(0,100)
  pp.mosaic.setHeightElastic(1,200)
  if c:
    pv12 = PixelsView(s1,s2,slice12(k3,c))
    pv12.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv13 = PixelsView(s1,s3,slice13(k2,c))
    pv13.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv23 = PixelsView(s2,s3,slice23(k1,c))
    pv23.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
    for pv in [pv12,pv13,pv23]:
      pv.setColorModel(ColorMap.getJet(0.5))
      if cmin!=cmax:
        pv.setClips(cmin,cmax)
    pp.pixelsView12.tile.addTiledView(pv12)
    pp.pixelsView13.tile.addTiledView(pv13)
    pp.pixelsView23.tile.addTiledView(pv23)
  pf = PlotFrame(pp)
  pf.setFontSizeForSlide(1.0,1.0)
  if c:
    pf.setSize(1036,814)
  else:
    pf.setSize(859,814)
  pf.setVisible(True)
  if png and pngDir:
    pf.paintToPng(300,6,pngDir+png+".png")

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


def plot3(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          horizon=None,skins=None,slices=None, k1=n1/2,png=None):
  n3 = len(f)
  n2 = len(f[0])
  n1 = len(f[0][0])
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
      ipg.setClips(-2.0,2.0)
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
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMap2Listener(cbar)
    sf.world.addChild(ipg)
  if cbar:
    cbar.setWidthMinimum(120)
  if horizon:
    sd = SurfaceDisplay()
    ts = sd.horizonWithAmplitude([-0.5,0.5],horizon,f)
    tg = TriangleGroup(True,ts[0],ts[1])
    #tg = TriangleGroup(True,s3,s2,horizon)
    #tg.setColor(Color.CYAN)
    sf.world.addChild(tg)
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
        cmap = ColorMap(-smax,smax,ColorMap.JET)
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
        #rgb = skin.getCellLinksRgb(r,g,b,xyz)
        #lg = LineGroup(xyz,rgb)
        lg = LineGroup(xyz)
        sg.addChild(lg)
        #ct = ct+1
    sf.world.addChild(sg)
  #ipg.setSlices(450,530,393)
  ipg.setSlices(k1,596,n3)
  if cbar:
    sf.setSize(1037,900)
  else:
    sf.setSize(900,900)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  zscale = 0.8*max(n2*d2,n3*d3)/(n1*d1)
  #zscale = 1.5*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.6)
  #ov.setScale(2.5)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(0.0,-0.15,-0.01))
  ov.setAzimuthAndElevation(225.0,40.0)
  #ov.setAzimuthAndElevation(-55.0,35.0)
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
