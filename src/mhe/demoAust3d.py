from utils import *

#pngDir = None
setupForSubset("aust")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

n1,n2,n3 = s1.count,s2.count,s3.count
d1,d2,d3 = s1.delta,s2.delta,s3.delta
#k1,k2,k3 = 88,60,160; azimuth=285; elevation=11 # for 3D view of all horizons
k1,k2,k3 = 354,800,200; azimuth=240; elevation=25 # for 3D view of strips
fmin,fmax = -5.5,5.5
k1f,k2f,k3f = 65,406,114
k1f,k2f,k3f = 48,406,114
k1f,k2f,k3f = 250,640,385
gmin,gmax,gint,glab = -2.0,2.0,0.5,"Amplitude"
background = Color.WHITE

pngDir = None
pngDir = "../../../png/mhe/aust/"
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
  #goSlopes()
  goHorizonOne()

def goSlopes():
  print "go slopes..."
  gx = readImage(gxfile)
  if not plotOnly:
    p2 = zerofloat(n1,n2,n3)
    p3 = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    sigma1,sigma2=6.0,2.0
    lsf = LocalSlopeFinder(sigma1,sigma2,sigma2,5) 
    lsf.findSlopes(gx,p2,p3,ep);
    writeImage(p2file,p2)
    writeImage(p3file,p3)
    writeImage(epfile,ep)
  else:
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    ep = readImage(epfile)
  plot3(gx)
  plot3(ep)
  plot3(p2)
  plot3(p3)

def goHorizonOne():
  k1 = [298]
  k2 = [459]
  k3 = [668]
  gx = readImage(gxfile)
  #gx = gain(gx)
  #writeImage(gxfile,gx)
  se = GlobalHorizon3()
  print "go horizonOne..."
  if not plotOnly:
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    ep = readImage(epfile)
    wp = pow(ep,6.0) 
    lmt = n1-1
    se.setWeights(0.0)
    se.setSmoothings(8.0,8.0)
    se.setCG(0.01,150)
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
  cmin = -2.0
  cmax =  2.0
  cmap = gray
  mp = ColorMap(-2.0,3.0,cmap)
  r1,g1,b1 = se.amplitudeRgb(mp,gx,surf1) 
  r2,g2,b2 = se.amplitudeRgb(mp,gx,surf2) 
  #cmin = min(surf1)+5
  #cmax = max(surf1)
  #mp = ColorMap(cmin,cmax,jet)
  #r1,g1,b1 = se.heightRgb(mp,surf1) 
  #r2,g2,b2 = se.heightRgb(mp,surf2) 
  sf1 = [surf1,r1,g1,b1]
  sf2 = [surf2,r2,g2,b2]
  plot3(gx,cmap=cmap,cmin=cmin,cmax=cmax,png="seis")
  plot3(gx,hz=sf1,ks=[k1,k2,k3],cmap=cmap,cmin=cmin,cmax=cmax,png="surf")
  #plot3(gx,hz=sf2,ks=[k1,k2,k3],cmap=cmap,cmin=cmin,cmax=cmax)
  plot3p(s1,s2,s3,gx)
  plot3p(s1,s2,s3,gx,hv=[surf1])
def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(60.0)
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
      ipg.setClips(-2.0,2.0) # use for subset plots
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-2.0,2.0)
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
  if hz:
    tg=TriangleGroup(True,s3,s2,hz[0],hz[1],hz[2],hz[3])
    sf.world.addChild(tg)
  if ks:
    pg = setPointGroup(ks[0],ks[1],ks[2],18)
    sf.world.addChild(pg)
  ipg.setSlices(498,759,0)
  if cbar:
    sf.setSize(987,720)
  else:
    sf.setSize(850,620)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  zscale = 0.4*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.6)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(-0.00,0.05,0.05))
  ov.setAzimuthAndElevation(160.0,40.0)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")
 

def plot3p(s1,s2,s3,f,hv=None,k1=None,k2=None,k3=None,cmap=ColorMap.GRAY,
        cmin=-2,cmax=2,clab=None,cint=0.1,png=None):
  width,height,cbwm = 1200,800,200
  n1,n2,n3 = s1.count,s2.count,s3.count
  orient = PlotPanelPixels3.Orientation.X1DOWN_X3RIGHT;
  axespl = PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM
  panel = PlotPanelPixels3(orient,axespl,s1,s2,s3,f)
  #panel.mosaic.setWidthElastic(0,100)
  #panel.mosaic.setWidthElastic(1,75)
  panel.mosaic.setHeightElastic(0,200)
  #panel.mosaic.setHeightElastic(1,100)
  panel.setSlice23(k1f)
  panel.setSlice13(k2f)
  panel.setSlice12(k3f)
  #panel.setSlice103(70)
  panel.setClips(cmin,cmax)
  if clab:
    cbar = panel.addColorBar(clab)
    cbar.setInterval(cint)
  panel.setColorBarWidthMinimum(50)
  panel.setLabel1("Samples")
  panel.setLabel2("Inline (traces)")
  panel.setLabel3("Crossline (traces)")
  panel.setInterval2(100)
  panel.setInterval3(100)
  panel.setColorModel(ColorMap.GRAY)
  panel.setLineColor(Color.WHITE)
  panel.setHLimits(0,s2.first,s2.last)
  panel.setVLimits(1,s1.first,s1.last)
  panel.setHLimits(1,s2.first,s2.last)
  if hv:
    nh = len(hv)
    hd = HorizonDisplay()
    cv12 = hd.slice12(k3f,s2,hv)
    cv13 = hd.slice13(k2f,s3,hv)
    cv23 = hd.slice23X(k1f,s2,s3,hv)
    mp = ColorMap(0,nh,ColorMap.JET)
    print nh
    for ih in range(nh):
      pv12 = PointsView(cv12[ih][1],cv12[ih][0])
      pv13 = PointsView(cv13[ih][1],cv13[ih][0])
      pv12.setLineWidth(2.0)
      pv13.setLineWidth(2.0)
      #pv12.setLineColor(mp.getColor(ih))
      #pv13.setLineColor(mp.getColor(ih))
      pv12.setLineColor(Color.YELLOW)
      pv13.setLineColor(Color.YELLOW)
      panel.pixelsView12.tile.addTiledView(pv12)
      panel.pixelsView13.tile.addTiledView(pv13)
      nc = len(cv23[ih][0])
      for ic in range(nc):
        pv23 = PointsView(cv23[ih][1][ic],cv23[ih][0][ic])
        pv23.setLineWidth(2.0)
        #pv23.setLineColor(mp.getColor(ih))
        pv23.setLineColor(Color.YELLOW)
        panel.pixelsView23.tile.addTiledView(pv23)
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(Color(0xfd,0xfe,0xff)) # easy to make transparent
  frame.setFontSize(16)#ForSlide(1.0,0.8)
  frame.setSize(width,height)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(1080,4,pngDir+"/"+png+".png")


def plot3f(g,hv=None,amin=None,amax=None,amap=None,alab=None,aint=None,
           png=1):
  pp = PlotPanelPixels3(
    PlotPanelPixels3.Orientation.X1DOWN_X3RIGHT,
    PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM,
    s1,s2,s3,g)
  pp.setSlice23(k1f)
  pp.setSlice13(k2f)
  pp.setSlice12(k3f)
  pp.setLabel1("Time (samples)")
  pp.setLabel2("Inline (samples)")
  pp.setLabel3("Crossline (samples)")
  pp.mosaic.setHeightElastic(0,180)
  pp.mosaic.setHeightElastic(1, 70)
  pp.setClips(gmin,gmax)
  pp.setColorModel(gray)
  pp.setLineColor(Color.WHITE)
  cb = pp.addColorBar("Amplitude")
  cb.setInterval(0.5)
  pp.setInterval1(100.0)
  pp.setInterval2(100.0)
  pp.setInterval3(100.0)
  if hv:
    nh = len(hv)
    hd = HorizonDisplay()
    cv12 = hd.slice12(k3f,s2,hv)
    cv13 = hd.slice13(k2f,s3,hv)
    cv23 = hd.slice23X(k1f,s2,s3,hv)
    mp = ColorMap(0,nh,ColorMap.JET)
    print nh
    for ih in range(nh):
      pv12 = PointsView(cv12[ih][1],cv12[ih][0])
      pv13 = PointsView(cv13[ih][1],cv13[ih][0])
      pv12.setLineWidth(2.0)
      pv13.setLineWidth(2.0)
      pv12.setLineColor(Color.YELLOW)
      pv13.setLineColor(Color.YELLOW)
      #pv12.setLineColor(mp.getColor(ih))
      #pv13.setLineColor(mp.getColor(ih))
      pp.pixelsView12.tile.addTiledView(pv12)
      pp.pixelsView13.tile.addTiledView(pv13)
      nc = len(cv23[ih][0])
      for ic in range(nc):
        pv23 = PointsView(cv23[ih][1][ic],cv23[ih][0][ic])
        pv23.setLineWidth(2.0)
        pv23.setLineColor(Color.YELLOW)
        #pv23.setLineColor(mp.getColor(ih))
        pp.pixelsView23.tile.addTiledView(pv23)
  pf = PlotFrame(pp)
  pf.setBackground(background)
  pp.setColorBarWidthMinimum(65)
  pf.setFontSize(16)
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
