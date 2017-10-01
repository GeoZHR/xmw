#############################################################################
"""
Demo of dynamic warping for automatic picking
Author: Xinming Wu, University of Texas at Austin
Version: 2016.06.01
"""


from utils import * 
setupForSubset("fake")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count
f1,f2,f3 = s1.getFirst(),s2.getFirst(),s3.getFirst()
d1,d2,d3 = s1.getDelta(),s2.getDelta(),s3.getDelta()


#############################################################################
fxfile = "x" # 
gtfile = "gt" # 
flfile = "fl" # 
ftfile = "ft" # 
fltfile = "flt" # 
fttfile = "ftt" # 
p2file = "p2" # seismic slopes
p3file = "p3" # seismic slopes
epfile = "ep" # seismic slopes
uxfile = "ux"
usfile = "us"
#pngDir = getPngDir()
pngDir = None
plotOnly = False

def main(args):
  gx = readImage(fxfile)
  print min(gx)
  print max(gx)
  plot3(gx,cmin=-0.01,cmax=0.01)

def goPath():
  n2 = 5
  n3 = 5
  ws = fillfloat(0.1,n2,n3)
  ws[2][2] = 1
  dk = Dijkstra(ws)
  dt = zerofloat(n2*n3)
  pd = zeroint(n2*n3)
  dk.apply(0,dt,pd)
  dk.printPath(n3,pd,0,n2*n3-1)

def goSlopes():
  fx = readImage(fxfile)
  if not plotOnly:
    p2 = zerofloat(n1,n2,n3)
    p3 = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    lsf = LocalSlopeFinder(4.0,1.0,1.0,5.0)
    lsf.findSlopes(fx,p2,p3,ep)
    writeImage(p2file,p2)
    writeImage(p3file,p3)
    #writeImage(epfile,ep)
  else:
    ep = readImage(epfile)
  plot3(fx)
  ep = pow(ep,8)
  plot3(ep,cmin=0.1,cmax=0.8)
def goFlattenWithSlopes():
  fx = readImage(fxfile)
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  ep = readImage(epfile)
  ep = pow(ep,6.0)
  fl = Flattener3()
  fl.setIterations(0.01,100)
  fm = fl.getMappingsFromSlopes(s1,s2,s3,p2,p3,ep)
  gt = fm.flatten(fx)
  writeImage(gtfile,gt)
  gt = readImage(gtfile)
  plot3(fx)
  plot3(gt)

def goFlatten():
  fx = readImage(fxfile)
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  ep = readImage(epfile)
  ep = pow(ep,12)
  sub(ep,min(ep),ep)
  div(ep,max(ep),ep)
  rgf = RecursiveGaussianFilterP(1)
  #rgf.apply000(fx,fx)
  df = DynamicFlattener(-10,30)
  df.setErrorExponent(1)
  df.setStrainMax(0.5)
  df.setWindow(2,150)
  df.setGate(10)
  df.setErrorExponent(1)
  df.setShiftSmoothing(1)
  df.setErrorSmoothing(3)
  #ux = zerofloat(n1,n2,n3)
  #gx = df.flattenXL(0,0,ep,copy(fx),ux)
  #writeImage(uxfile,ux)
  ux = readImage(uxfile)
  ss = SmoothWithSlopes()
  gx = ss.flatten(ux,fx)
  us = ss.smooth(ep,p2,p3,ux)
  writeImage(usfile,us)
  us = readImage(usfile)
  gh = ss.flatten(us,fx)
  plot3(fx)
  plot3(gx)
  plot3(gh)
  gt = readImage(gtfile)
  plot3(gt)
  hs = ss.getHorizons(n1,1,0,us)
  hs = mul(hs,2)
  c1 = Sampling(n1,2,0)
  plot3p(c1,s2,s3,fx,k1=51,k2=445,k3=208,cmin=-1.5,cmax=1.5)
  plot3p(c1,s2,s3,gh,k1=46,k2=445,k3=208,cmin=-1.5,cmax=1.5)
  plot3p(c1,s2,s3,fx,hv=hs,k1=120,k2=445,k3=208,cmin=-1.5,cmax=1.5)
  plot3(fx,surf=div(hs[46],2))

  '''
  n2 = 1000
  n3 = 1000
  ws = fillfloat(1,n2,n3)
  dk = Dijkstra(ws)
  dt = zerofloat(n2*n3)
  pd = zeroint(n2*n3)
  dk.apply(0,dt,pd)
  dk.printPath(n3,pd,0,n2*n3-1)
  '''

def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(80.0)
  ref.apply1(g,g)
  y = like(x);
  div(x,sqrt(g),y)
  return y

def like(x):
  n3 = len(x)
  n2 = len(x[0])
  n1 = len(x[0][0])
  return zerofloat(n1,n2,n3)

def smooth(sig,u):
  v = copy(u)
  rgf = RecursiveGaussianFilterP(sig)
  rgf.apply0(u,v)
  return v


#############################################################################
# plotting
def jetFill(alpha):
  return ColorMap.setAlpha(ColorMap.JET,alpha)
def jetFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.JET,a)
def jetRamp(alpha):
  return ColorMap.setAlpha(ColorMap.JET,rampfloat(0.0,alpha/256,256))
def bwrRamp(alpha):
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,rampfloat(0.0,alpha/256,256))
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


backgroundColor = Color.WHITE

def plot3p(s1,s2,s3,f,g=None,hv=None,k1=None,k2=None,k3=None,cmap=ColorMap.GRAY,
        cmin=-1,cmax=1,clab=None,cint=0.1,png=None):
  width,height,cbwm = 800,1200,200
  n1,n2,n3 = s1.count,s2.count,s3.count
  print n1
  print n2
  print n3
  orient = PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT;
  axespl = PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM
  panel = PlotPanelPixels3(orient,axespl,s1,s2,s3,f)
  #panel.mosaic.setWidthElastic(0,100)
  #panel.mosaic.setWidthElastic(0,50)
  panel.mosaic.setHeightElastic(0,65)
  #panel.mosaic.setHeightElastic(1,100)
  panel.setSlice23(k1)
  panel.setSlice13(k2)
  panel.setSlice12(k3)
  #panel.setSlice103(70)
  panel.setClips(cmin,cmax)
  if clab:
    cbar = panel.addColorBar(clab)
    cbar.setInterval(cint)
  panel.setColorBarWidthMinimum(50)
  panel.setLabel1("Depth (samples)")
  panel.setLabel2("Inline (traces)")
  panel.setLabel3("Crossline (traces)")
  panel.setInterval2(50)
  panel.setInterval3(50)
  panel.setColorModel(ColorMap.GRAY)
  panel.setLineColor(Color.YELLOW)
  panel.setHLimits(0,s2.first,s2.last)
  panel.setVLimits(0,s3.first,s3.last)
  panel.setVLimits(1,s1.first,s1.last)
  panel.setHLimits(1,s2.first,s2.last)
  if g:
    pv12 = PixelsView(s1,s2,slice12(k3,g))
    pv12.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv13 = PixelsView(s1,s3,slice13(k2,g))
    pv13.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv23 = PixelsView(s2,s3,slice23(k1,g))
    pv23.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
    for pv in [pv12,pv13,pv23]:
      pv.setColorModel(cmap)
      if cmin!=cmax:
        pv.setClips(cmin,cmax)
    panel.pixelsView12.tile.addTiledView(pv12)
    panel.pixelsView13.tile.addTiledView(pv13)
    panel.pixelsView23.tile.addTiledView(pv23)
  if hv:
    nh = len(hv)
    hd = HorizonDisplay()
    cv12 = hd.slice12(k3,s2,hv)
    cv13 = hd.slice13(k2,s3,hv)
    cv23 = hd.slice23X(k1,s2,s3,div(hv,(float)(s1.getDelta())))
    mp = ColorMap(0,nh,ColorMap.JET)
    print nh
    k = 1
    for ih in range(0,nh-3,15):
      color = Color.MAGENTA
      if(k%4==0): color=Color.MAGENTA
      if(k%4==1): color=Color.GREEN
      if(k%4==2): color=Color.RED
      if(k%4==3): color=Color.BLUE
      k = k+1
      pv12 = PointsView(cv12[ih][1],cv12[ih][0])
      pv13 = PointsView(cv13[ih][1],cv13[ih][0])
      pv12.setLineWidth(3.0)
      pv13.setLineWidth(3.0)
      pv12.setLineColor(color)#mp.getColor(ih))
      pv13.setLineColor(color)#mp.getColor(ih))
      panel.pixelsView12.tile.addTiledView(pv12)
      panel.pixelsView13.tile.addTiledView(pv13)
      nc = len(cv23[ih][0])
      for ic in range(nc):
        pv23 = PointsView(cv23[ih][0][ic],cv23[ih][1][ic])
        pv23.setLineWidth(3.0)
        pv23.setLineColor(color)#mp.getColor(ih))
        panel.pixelsView23.tile.addTiledView(pv23)
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(Color(0xfd,0xfe,0xff)) # easy to make transparent
  frame.setFontSize(14)#ForSlide(1.0,0.8)
  frame.setSize(width,height)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(720,4,pngDir+"/"+png+".png")

def plot2f(s1,s2,f,g=None,cmin=None,cmax=None,cmap=None,label=None,png=None):
  n2 = len(f)
  n1 = len(f[0])
  f1,f2 = s1.getFirst(),s2.getFirst()
  d1,d2 = s1.getDelta(),s2.getDelta()
  panel = panel2Teapot()
  panel.setHInterval(1.0)
  panel.setVInterval(1.0)
  panel.setHLabel("Lateral position (km)")
  panel.setVLabel("Time (s)")
  #panel.setHInterval(100.0)
  #panel.setVInterval(100.0)
  #panel.setHLabel("Pixel")
  #panel.setVLabel("Pixel")
  if label:
    panel.addColorBar(label)
  else:
    panel.addColorBar()
  panel.setColorBarWidthMinimum(80)
  pv = panel.addPixels(s1,s2,f)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setColorModel(ColorMap.GRAY)
  pv.setClips(-2,2)
  if g:
    pv = panel.addPixels(s1,s2,g)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setColorModel(cmap)
    if label:
      panel.addColorBar(label)
    else:
      panel.addColorBar()
  if cmin and cmax:
    pv.setClips(cmin,cmax)
  frame2Teapot(panel,png)
def panel2Teapot():
  panel = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT)#,PlotPanel.AxesPlacement.NONE)
  return panel
def frame2Teapot(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  #frame.setFontSizeForPrint(8,240)
  #frame.setSize(1240,774)
  #frame.setFontSizeForSlide(1.0,0.9)
  frame.setFontSize(24)
  frame.setSize(450+80,700)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(400,3.2,pngDir+png+".png")
  return frame


def plot(f,g=None,t=None,cmap=None,cmin=None,cmax=None,cint=None,
        label=None,neareast=False,png=None): 
  orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT;
  panel = PlotPanel(1,1,orientation)#,PlotPanel.AxesPlacement.NONE)
  panel.setVInterval(50)
  panel.setHInterval(200)
  #panel.setHLabel("Inline (traces)")
  #panel.setVLabel("Time (samples)")
  pxv = panel.addPixels(0,0,s1,s2,f);
  pxv.setColorModel(ColorMap.GRAY)
  #pxv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if g:
    pxv.setClips(-1,1)
  else:
    if cmin and cmax:
      pxv.setClips(cmin,cmax)
  if g:
    pv = panel.addPixels(s1,s2,g)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setColorModel(cmap)
    if cmin and cmax:
      pv.setClips(cmin,cmax)
  moc = panel.getMosaic();
  frame = PlotFrame(panel);
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  #frame.setTitle("normal vectors")
  frame.setVisible(True);
  frame.setSize(1400,700)
  frame.setFontSize(24)
  if pngDir and png:
    frame.paintToPng(720,3.333,pngDir+png+".png")

def plot2(s1,s2,c,u=None,vint=30,hint=200,hw=None,vw=None,
          cmin=0.0,cmax=0.0,cmap=ColorMap.GRAY,color=Color.RED,
          title=None,perc=None,png=None):
  n2 = s2.getCount()
  n1 = s1.getCount()
  panel = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT)
          #PlotPanel.AxesPlacement.NONE)
  panel.setHLimits(0,s2.first,s2.last)
  panel.setVLimits(0,s1.first,s1.last)
  panel.setVInterval(0,vint)
  panel.setHInterval(0,hint)
  if title:
    panel.addTitle(title)
  cv = panel.addPixels(0,0,s1,s2,c)
  cv.setInterpolation(PixelsView.Interpolation.NEAREST)
  cv.setColorModel(cmap)
  if perc:
    cv.setPercentiles(100-perc,perc)
  elif cmin<cmax:
    cv.setClips(cmin,cmax)
  if u:
    nu = len(u)
    x2 = rampfloat(0,1,n2)
    cp = ColorMap(0,nu-1,ColorMap.PRISM)
    for iu in range(0,nu,10):
      uv = panel.addPoints(0,0,u[iu],x2)
      uv.setLineColor(cp.getColor(iu))
      uv.setLineWidth(2.5)
    uv = panel.addPoints(0,0,u[nu-1],x2)
    uv.setLineColor(cp.getColor(nu-1))
    uv.setLineWidth(2.5)
  #panel.addColorBar()
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(backgroundColor)
  #frame.setFontSizeForPrint(8,240)
  #frame.setSize(470,1000)
  frame.setFontSize(12)
  if hw and vw:
    frame.setSize(hw,vw)
  else:
    frame.setSize(round(n2*0.4),round(n1*2))
  frame.setVisible(True)
  if png and pngDir:
    png += "n"+str(int(10*nrms))
    png += "s"+str(int(10*strainMax))
    frame.paintToPng(720,3.33333,pngDir+"/"+png+".png")

def plotfg(f,g,png=None):
  n = len(f)
  panel = PlotPanel(2,1,PlotPanel.Orientation.X1RIGHT_X2UP)
  panel.mosaic.setHeightElastic(0,25)
  panel.mosaic.setHeightElastic(1,25)
  panel.setHLimits(0,0,n-1)
  fv = panel.addPoints(0,0,f)
  gv = panel.addPoints(1,0,g)
  fv.setLineWidth(2)
  gv.setLineWidth(2)
  panel.setVLimits(0,-1.3,1.3)
  panel.setVLimits(1,-1.3,1.3)
  #panel.setHLabel("sample index i")
  #panel.setVLabel(0,"f")
  #panel.setVLabel(1,"g")
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(backgroundColor)
  frame.setFontSizeForPrint(8,240)
  frame.setSize(1000,650)
  frame.setVisible(True)
  if png and pngDir:
    png += "n"+str(int(10*nrms))
    frame.paintToPng(720,3.33333,pngDir+"/"+png+".png")

def plotc(c,s=None,u=None,cmin=0.0,cmax=0.0,perc=None,png=None):
  n,nlag = len(c[0]),len(c)
  s1 = Sampling(n,1.0,0.0)
  slag = Sampling(nlag,1.0,-(nlag-1)/2)
  panel = PlotPanel(1,1,PlotPanel.Orientation.X1RIGHT_X2UP)
  panel.setHLimits(0,0,n-1)
  panel.setVLimits(0,slag.first,slag.last)
  cv = panel.addPixels(0,0,s1,slag,c)
  cv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if perc:
    cv.setPercentiles(100-perc,perc)
  elif cmin<cmax:
    cv.setClips(cmin,cmax)
  if s:
    sv = panel.addPoints(0,0,s)
    sv.setLineColor(Color.WHITE)
    sv.setLineStyle(PointsView.Line.DOT)
    sv.setLineWidth(3)
  if u:
    uv = panel.addPoints(0,0,u)
    uv.setLineColor(Color.WHITE)
    uv.setLineWidth(3)
  panel.addColorBar()
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(backgroundColor)
  frame.setFontSizeForPrint(8,240)
  frame.setSize(1000,470)
  frame.setVisible(True)
  if png and pngDir:
    png += "n"+str(int(10*nrms))
    png += "s"+str(int(10*strainMax))
    frame.paintToPng(720,3.33333,pngDir+"/"+png+".png")

def plot2c(c,s,u,clip=None,perc=None,png=None):
  n,nlag = len(c[0]),len(c)
  s1 = Sampling(n,1.0,0.0)
  slag = Sampling(nlag,1.0,-(nlag-1)/2)
  panel = PlotPanel(2,1,PlotPanel.Orientation.X1RIGHT_X2UP)
  panel.setHLimits(0,0,n-1)
  panel.setVLimits(0,slag.first,slag.last)
  panel.setVLimits(1,slag.first,slag.last)
  cv0 = panel.addPixels(0,0,s1,slag,c)
  cv1 = panel.addPixels(1,0,s1,slag,c)
  cv0.setInterpolation(PixelsView.Interpolation.NEAREST)
  cv1.setInterpolation(PixelsView.Interpolation.NEAREST)
  cv0.setColorModel(ColorMap.getGray(0.0,0.8))
  cv1.setColorModel(ColorMap.getGray(0.0,0.8))
  if perc:
    cv0.setPercentiles(0,perc)
    cv1.setPercentiles(0,perc)
  elif clip:
    cv0.setClips(0.0,clip)
    cv1.setClips(0.0,clip)
  if s:
    sv = panel.addPoints(1,0,s)
    sv.setLineColor(Color.WHITE)
    sv.setLineStyle(PointsView.Line.DOT)
    sv.setLineWidth(3)
  if u:
    uv = panel.addPoints(1,0,u)
    uv.setLineColor(Color.WHITE)
    uv.setLineWidth(3)
  panel.addColorBar()
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(backgroundColor)
  frame.setFontSizeForPrint(8,240)
  frame.setSize(1000,850)
  frame.setVisible(True)
  if png and pngDir:
    png += "n"+str(int(10*nrms))
    png += "s"+str(int(10*strainMax))
    frame.paintToPng(720,3.33333,pngDir+"/"+png+".png")

def plot3c(c,s,u,cmin=0.0,cmax=0.0,png=None):
  print "c0: min =",min(c[0])," max =",max(c[0])
  print "c1: min =",min(c[1])," max =",max(c[1])
  print "c2: min =",min(c[2])," max =",max(c[2])
  n,nlag = len(c[0][0]),len(c[0])
  s1 = Sampling(n,1.0,0.0)
  slag = Sampling(nlag,1.0,-(nlag-1)/2)
  panel = PlotPanel(3,1,PlotPanel.Orientation.X1RIGHT_X2UP)
  panel.setHLimits(0,0,n-1)
  panel.setVLimits(0,slag.first,slag.last)
  panel.setVLimits(1,slag.first,slag.last)
  panel.setVLimits(2,slag.first,slag.last)
  cv0 = panel.addPixels(0,0,s1,slag,c[0])
  cv1 = panel.addPixels(1,0,s1,slag,c[1])
  cv2 = panel.addPixels(2,0,s1,slag,c[2])
  cv0.setInterpolation(PixelsView.Interpolation.NEAREST)
  cv1.setInterpolation(PixelsView.Interpolation.NEAREST)
  cv2.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmin<cmax:
    cv0.setClips(cmin,cmax)
    cv1.setClips(cmin,cmax)
    cv2.setClips(cmin,cmax)
  if s:
    sv = panel.addPoints(0,0,s[0])
    sv.setLineColor(Color.WHITE)
    sv.setLineStyle(PointsView.Line.DOT)
    sv.setLineWidth(3)
    sv = panel.addPoints(1,0,s[1])
    sv.setLineColor(Color.WHITE)
    sv.setLineStyle(PointsView.Line.DOT)
    sv.setLineWidth(3)
    sv = panel.addPoints(2,0,s[2])
    sv.setLineColor(Color.WHITE)
    sv.setLineStyle(PointsView.Line.DOT)
    sv.setLineWidth(3)
  if u:
    uv = panel.addPoints(0,0,u[0])
    uv.setLineColor(Color.WHITE)
    uv.setLineWidth(3)
    uv = panel.addPoints(1,0,u[1])
    uv.setLineColor(Color.WHITE)
    uv.setLineWidth(3)
    uv = panel.addPoints(2,0,u[2])
    uv.setLineColor(Color.WHITE)
    uv.setLineWidth(3)
  panel.addColorBar()
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(backgroundColor)
  frame.setFontSizeForPrint(8,240)
  frame.setSize(1000,1050)
  frame.setVisible(True)
  if png and pngDir:
    png += "n"+str(int(10*nrms))
    #png += "s"+str(int(10*strainMax))
    frame.paintToPng(720,3.33333,pngDir+"/"+png+".png")

def plot3(f,g=None,k2=0,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          xyz=None,cells=None,skins=None,fbs=None,surf=None,smax=0.0,
          links=False,curve=False,trace=False,png=None):
  n1 = len(f[0][0])
  n2 = len(f[0])
  n3 = len(f)
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
  if surf:
    sd = SurfaceDisplay()
    xyz,rgb = sd.buildTrigs(s3,s2,surf)
    #tg = TriangleGroup(True, s3, s2, surf)
    tg = TriangleGroup(True,xyz,rgb)
    states = StateSet()
    cs = ColorState()
    #cs.setColor(Color.ORANGE)
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
  ipg.setSlices(329,445,208)
  if cbar:
    sf.setSize(987,700)
  else:
    sf.setSize(850,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  zscale = 0.9*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.4)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(-0.05,-0.05,0.02))
  ov.setAzimuthAndElevation(135,25.0)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(137,1,pngDir+png+"cbar.png")
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
