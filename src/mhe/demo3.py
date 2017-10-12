from utils import *

#pngDir = None
setupForSubset("f3d")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

n1,n2,n3 = s1.count,s2.count,s3.count
d1,d2,d3 = s1.delta,s2.delta,s3.delta
#k1,k2,k3 = 88,60,160; azimuth=285; elevation=11 # for 3D view of all horizons
k1,k2,k3 = 154,950,540; azimuth=240; elevation=25 # for 3D view of strips
fmin,fmax = -5.5,5.5
k1f,k2f,k3f = 65,406,114
k1f,k2f,k3f = 48,406,114
k1f,k2f,k3f = 48,406,0
gmin,gmax,gint,glab = -2.0,2.0,0.5,"Amplitude"
background = Color.WHITE

pngDir = None
pngDir = "../../../png/mhe/f3d/"
gxfile = "gx"
p2file = "p2"
p3file = "p3"
epfile = "ep"
mh1file = "mh1"
mh2file = "mh2"
sh1file = "sh1"
sh2file = "sh2"
plotOnly = True

def main(args):
  #slopes()
  goHorizonOne()
  #goHorizonTwo()

def slopes():
  gx = readImage(gxfile)
  p2 = zerofloat(n1,n2,n3)
  p3 = zerofloat(n1,n2,n3)
  ep = zerofloat(n1,n2,n3)
  sigma1,sigma2=8.0,2.0
  lsf = LocalSlopeFinder(sigma1,sigma2,sigma2,5) 
  lsf.findSlopes(gx,p2,p3,ep);
  writeImage(p2file,p2)
  writeImage(p3file,p3)
  writeImage(epfile,ep)
  for g in [p2,p3,ep]:
    world = World()
    addImage2ToWorld(world,gx,g)
    makeFrame(world)

def goHorizonOne():
  k1 = [144,172]
  k2 = [667,340]
  k3 = [217,302]
  gx  = readImage(gxfile)
  se = GlobalHorizon3()
  if not plotOnly:
    p2  = readImage(p2file)
    p3  = readImage(p3file)
    ep  = readImage(epfile)
    wp = pow(ep,4.0) 
    lmt = n1-1
    se.setWeights(0.0)
    se.setSmoothings(8.0,8.0)
    se.setCG(0.01,100)
    se.setExternalIterations(20)
    surf = se.surfaceInitialization(n2,n3,lmt,k1,k2,k3)
    surf2 = copy(surf)
    surf1 = se.surfaceUpdateFromSlopesAndCorrelations(5,20,gx,wp,p2,p3,k2,k3,surf)
    surf2 = se.surfaceUpdateFromSlopesAndCorrelations(wp,p2,p3,k2,k3,surf2)
    writeImage(mh1file,surf1) 
    writeImage(sh1file,surf2) 
  else:
    surf1 = readImage2(mh1file)
    surf2 = readImage2(sh1file)
  mp = ColorMap(-2.5,2.5,rwb)
  r1,g1,b1 = se.amplitudeRgb(mp,gx,surf1) 
  r2,g2,b2 = se.amplitudeRgb(mp,gx,surf2) 
  '''
  cmin = min(surf1)+5
  cmax = max(surf1)
  mp = ColorMap(cmin,cmax,jet)
  r1,g1,b1 = se.heightRgb(mp,surf1) 
  r2,g2,b2 = se.heightRgb(mp,surf2) 
  '''
  sf1 = [surf1,r1,g1,b1]
  sf2 = [surf2,r2,g2,b2]
  plot3(gx,cmap=rwb,clab="Amplitude",png="seis")
  plot3(gx,hz=sf1,ks=[k1,k2,k3],cmap=rwb,png="surf1m")
  plot3(gx,hz=sf2,ks=[k1,k2,k3],cmap=rwb,png="surf1s")
def goHorizonTwo():
  k1 = [90 ]
  k2 = [413]
  k3 = [312]
  gx  = readImage(gxfile)
  se = GlobalHorizon3()
  if not plotOnly:
    p2  = readImage(p2file)
    p3  = readImage(p3file)
    ep  = readImage(epfile)
    wp = pow(ep,8.0) 
    lmt = n1-1
    se.setWeights(0.0)
    se.setSmoothings(8.0,8.0)
    se.setCG(0.01,100)
    se.setExternalIterations(20)
    surf = se.surfaceInitialization(n2,n3,lmt,k1,k2,k3)
    surf2 = copy(surf)
    surf1 = se.surfaceUpdateFromSlopesAndCorrelations(5,20,gx,wp,p2,p3,k2,k3,surf)
    surf2 = se.surfaceUpdateFromSlopesAndCorrelations(wp,p2,p3,k2,k3,surf2)
    writeImage(mh2file,surf1) 
    writeImage(sh2file,surf2) 
  else:
    surf1 = readImage2(mh2file)
    surf2 = readImage2(sh2file)
  mp = ColorMap(-2.5,2.5,rwb)
  r1,g1,b1 = se.amplitudeRgb(mp,gx,surf1) 
  r2,g2,b2 = se.amplitudeRgb(mp,gx,surf2) 
  '''
  cmin = min(surf1)+5
  cmax = max(surf1)
  mp = ColorMap(cmin,cmax,jet)
  r1,g1,b1 = se.heightRgb(mp,surf1) 
  r2,g2,b2 = se.heightRgb(mp,surf2) 
  '''
  sf1 = [surf1,r1,g1,b1]
  sf2 = [surf2,r2,g2,b2]
  plot3(gx,hz=sf1,ks=[k1,k2,k3],cmap=rwb,png="surfm2")
  plot3(gx,hz=sf2,ks=[k1,k2,k3],cmap=rwb,png="surfs2")
 
def display(filename):
  f = readImage(filename)
  world = World()
  ipg = addImageToWorld(world,f,cmap=gray)
  ipg.setSlices(k1,k2,k3)
  makeFrame(world)


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

def addColorBar(frame,clab=None,cint=None):
  cbar = ColorBar(clab)
  if cint:
    cbar.setInterval(cint)
  cbar.setFont(Font("Arial",Font.PLAIN,32)) # size by experimenting
  cbar.setWidthMinimum
  cbar.setBackground(Color.WHITE)
  frame.add(cbar,BorderLayout.EAST)
  return cbar

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
      ipg.setClips(-2.5,2.5) # use for subset plots
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-2.5,2.5)
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
    pg = setPointGroup(ks[0],ks[1],ks[2],18)
    sf.world.addChild(pg)
  ipg.setSlices(n1-1,n2-1,0)
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
  ov.setTranslate(Vector3(-0.08,0.00,0.05))
  ov.setTranslate(Vector3( 0.02,0.00,0.05))
  ov.setAzimuthAndElevation(130.0,45.0)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")
 

def plot3f(g,a=None,amin=None,amax=None,amap=None,alab=None,aint=None,
           png=1):
  pp = PlotPanelPixels3(
    PlotPanelPixels3.Orientation.X1DOWN_X3RIGHT,
    PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM,
    s1,s2,s3,g)
  pp.setSlices(k1f,k2f,k3f)
  pp.setLabel1("Time (s)")
  pp.setLabel2("Inline (km)")
  pp.setLabel3("Crossline (km)")
  pp.mosaic.setHeightElastic(0,180)
  pp.mosaic.setHeightElastic(1, 70)
  pp.setClips(gmin,gmax)
  pp.setColorModel(gray)
  if a:
    pp.setLineColor(Color.WHITE)
    cb = pp.addColorBar(alab)
    if aint:
      cb.setInterval(aint)
  else:
    pp.setLineColor(Color.YELLOW)
    cb = pp.addColorBar("Amplitude")
    cb.setInterval(0.5)
  pp.setInterval1(0.1)
  pp.setInterval2(2.0)
  pp.setInterval3(2.0)
  if a:
    pv12 = PixelsView(s1,s2,slice12(k3,a))
    pv12.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv12.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv13 = PixelsView(s1,s3,slice13(k2,a))
    pv13.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv13.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv23 = PixelsView(s2,s3,slice23(k1f,a))
    pv23.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
    pv23.setInterpolation(PixelsView.Interpolation.NEAREST)
    for pv in [pv12,pv13,pv23]:
      pv.setColorModel(ColorMap.getJet(0.5))
      if amin!=amax:
        pv.setClips(amin,amax)
    pp.pixelsView12.tile.addTiledView(pv12)
    pp.pixelsView13.tile.addTiledView(pv13)
    pp.pixelsView23.tile.addTiledView(pv23)
  pf = PlotFrame(pp)
  pf.setBackground(background)
  pp.setColorBarWidthMinimum(170)
  #pf.setFontSizeForSlide(1.0,0.8)
  pf.setSize(1200,700)
  pf.setVisible(True)
  if png and pngDir:
    png = pngDir+"f"+str(k1f)
    pf.paintToPng(360,7.0,png+".png")

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

def plotFrame(s1,s2,f,h,i3t):
  orient = PlotPanel.Orientation.X1DOWN_X2RIGHT
  panel  = PlotPanel(2,1,orient)
  pxv    = panel.addPixels(0,0,s1,s2,f)
  pxv.setColorModel(ColorMap.GRAY)
  pxv    = panel.addPixels(1,0,s1,s2,f)
  pxv.setColorModel(ColorMap.GRAY)
  ptv1 = panel.addPoints(0,0,h[0],h[2])
  ptv2 = panel.addPoints(0,0,h[1],h[2])
  ptv1.setStyle("b-")
  ptv2.setStyle("g-")
  panel.setTitle("section "+i3t)
  panel.addColorBar()
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setVisible(True)
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
