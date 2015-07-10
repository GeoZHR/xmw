"""
Demonstrate simultaneous multiple-well ties
Author: Xinming Wu, Colorado School of Mines
Version: 2015.06.03
"""

from utils import *
setupForSubset("subt")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count
d1,d2,d3 = s1.delta,s2.delta,s3.delta

# Names and descriptions of image files used below.
gxfile  = "gx" # input seismic image 
flfile  = "fl" # fault likelihood
fpfile  = "fp" # fault strike (phi)
ftfile  = "ft" # fault dip (theta)
gtfile  = "gt" # RGT volume
ghfile  = "gh" # horizon volume
gufile  = "gu" # flattened image 
p2file  = "p2" # inline slopes
p3file  = "p3" # crossline slopes
epfile  = "ep" # eigenvalue-derived planarity
wpfile  = "wp" # weight image for flattening

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
pngDir = None
pngDir = "../../../png/swt/"
plotOnly = True

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  goDisplay()
  #goFaultLikelihoods()
  #goSlopes()
  #goFlatten()
  #goFlattenC()
  #goFlattenD()

def goTest():
  gx = readImage(gxfile)
  #plot3(gx)
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  fs = FaultSkinner()
  cells = fs.findCells([fl,fp,ft])
  fci = FaultCellsToImage(cells)
  fli = fci.getLikelihoodThick(0.6,fl)
  plot3(gx,fl,cmin=0.2,cmax=1,cmap=jetRamp(1.0),
      clab="Fault likelihood",png="fl")
  plot3(gx,fli,cmin=0.2,cmax=1,cmap=jetRamp(1.0),
      clab="Fault likelihood",png="fl")

def goFaultLikelihoods():
  if not plotOnly:
    minPhi,maxPhi = -20,20
    minTheta,maxTheta = 80,85
    sigmaPhi,sigmaTheta = 8,30
    gx = readImage(gxfile)
    sig1,sig2,sig3,pmax = 16.0,1.0,1.0,5.0
    p2,p3,ep = FaultScanner.slopes(sig1,sig2,sig3,pmax,gx)
    gx = FaultScanner.taper(10,0,0,gx)
    fs = FaultScanner(sigmaPhi,sigmaTheta)
    fl,fp,ft = fs.scan(minPhi,maxPhi,minTheta,maxTheta,p2,p3,gx)
    print "fl min =",min(fl)," max =",max(fl)
    print "fp min =",min(fp)," max =",max(fp)
    print "ft min =",min(ft)," max =",max(ft)
    writeImage(flfile,fl)
    writeImage(fpfile,fp)
    writeImage(ftfile,ft)
  else:
    gx = readImage(gxfile)
    fl = readImage(flfile)
    fp = readImage(fpfile)
    ft = readImage(ftfile)
  plot3(gx,fl,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
      clab="Fault likelihood",png="fl")
  plot3(gx,fp,cmin=-20,cmax=20,cmap=jetFill(1.0),
      clab="Fault strike (degrees)",cint=5,png="fp")
  plot3(gx,ft,cmin=80,cmax=85,cmap=jetFill(1.0),
      clab="Fault dip (degrees)",png="ft")


def goSlopes():
  print "goSlopes ..."
  if not plotOnly:
    # set half-width of smoother for computing structure tensors
    sig1 = 4.0 #half-width in vertical direction
    sig2 = 1.0 #half-width in one literal direction
    sig3 = 1.0 #half-width in another literal direction
    pmax = 5.0 #maximum slope returned by this slope finder
    gx = readImage(gxfile)
    p2 = zerofloat(n1,n2,n3)
    p3 = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    lsf = LocalSlopeFinder(sig1,sig2,sig3,pmax)
    lsf.findSlopes(gx,p2,p3,ep);

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
    print "ep  min =",min(ep)," max =",max(ep)
  else:
    gx = readImage(gxfile)
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    ep = readImage(epfile)
  plot3(gx)
  plot3(gx,p2, cmin=-1,cmax=1,cmap=jetRamp(1.0),
      clab="Inline slope (sample/sample)",png="p2")
  plot3(gx,p3, cmin=-1,cmax=1,cmap=jetRamp(1.0),
      clab="Crossline slope (sample/sample)",png="p3")
  plot3(gx,pow(ep,4.0),cmin=0,cmax=1,cmap=jetRamp(1.0),
      clab="Planarity")

def goFlatten():
  print "goFlatten ..."
  if not plotOnly:
    gx = readImage(gxfile)
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    ep = readImage(epfile)
    p2 = mul(d1/d2,p2)
    p3 = mul(d1/d3,p3)
    ep = pow(ep,12.0)
    fl = Flattener3()
    #fl.setWeight1(0.05)
    fl.setIterations(0.01,200)
    fm = fl.getMappingsFromSlopes(s1,s2,s3,p2,p3,ep)
    gu = fm.flatten(gx) # flattened image
    gt = fm.u1 # rgt volume
    gh = fm.x1 # horizon volume
    writeImage(gufile,gu)
    writeImage(gtfile,gt)
    writeImage(ghfile,gh)
  else:
    gx = readImage(gxfile)
    gu = readImage(gufile)
    gt = readImage(gtfile)
    gh = readImage(ghfile)
  plot3(gx,png="seismic")
  plot3(gu,png="flattened")
  plot3(gx,gt,cmin=min(gt),cmax=max(gt),cmap=jetRamp(1.0),
        clab="Relative geologic time",png="rgt")
  '''
  ha = []
  hs = [800,750,700,650,600,550,500,450,400,350,300]
  hs = mul(hs,0.002)
  for ih, h in enumerate(hs):
    print h
    ha.append(h)
    plot3(gx,hs=ha,png="horizon"+str(ih))
  '''
def goFlattenC():
  print "goFlattenC ..."
  if not plotOnly:
    k11 = [351,358,328,337]
    k12 = [ 75,113,199, 19]
    k13 = [ 68, 68, 24, 29]
    k21 = [296,306,286,286,284,296,297]
    k22 = [ 89,102,163, 14,112,104, 76]
    k23 = [ 68, 68, 10, 13, 28, 52,  8]
    k31 = [269,275,256,260]
    k32 = [ 77,134,179, 21]
    k33 = [ 68, 68, 10, 10]
    k41 = [326,331,313,305]
    k42 = [ 67,127, 24,217]
    k43 = [ 68, 68, 15, 18]

    k51 = [567,559,587,601]
    k52 = [167,193, 51,214]
    k53 = [ 28, 28, 63, 60]

    gx = readImage(gxfile)
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    ep = readImage(epfile)
    fl = readImage(flfile)
    fp = readImage(fpfile)
    ft = readImage(ftfile)
    fs = FaultSkinner()
    cells = fs.findCells([fl,fp,ft])
    fci = FaultCellsToImage(cells)
    fli = fci.getLikelihoodThick(0.6,fl)
    wp = sub(1.0,fli)
    wp = mul(wp,ep)
    wp = pow(wp,16)
    #wh = pow(ep,8.0)
    #wp = pow(wp,20.0)
    #ep = pow(ep,20.0)
    #wp = mul(wp,ep)
    zm = ZeroMask(0.1,1,1,1,gx)
    zero,tiny=0.0,0.01
    zm.setValue(tiny,wp)

    sc = SetupConstraints()
    kk1 = sc.extend(k11,k12,k13,n2,n3,p2,p3,wp,gx)
    kk2 = sc.extend(k21,k22,k23,n2,n3,p2,p3,wp,gx)
    kk3 = sc.extend(k31,k32,k33,n2,n3,p2,p3,wp,gx)
    kk4 = sc.extend(k41,k42,k43,n2,n3,p2,p3,wp,gx)

    k1 = [kk1[0],kk2[0],kk3[0],kk4[0]]
    k2 = [kk1[1],kk2[1],kk3[1],kk4[1]]
    k3 = [kk1[2],kk2[2],kk3[2],kk4[2]]
    k4 = [kk1[3],kk2[3],kk3[3],kk4[3]]
    p2 = mul(d1/d2,p2)
    p3 = mul(d1/d3,p3)
    fl = Flattener3C()
    fl.setWeight1(0.08)
    fl.setIterations(0.01,200)
    fl.setSmoothings(8.0,8.0)
    fm = fl.getMappingsFromSlopes(s1,s2,s3,p2,p3,wp,k4,k1,k2,k3)
    gx = normalize(gx)
    gu = fm.flatten(gx) # flattened image
    gt = fm.u1 # rgt volume
    gh = fm.x1 # horizon volume
    writeImage(gufile,gu)
    writeImage(gtfile,gt)
    writeImage(ghfile,gh)
    writeImage(wpfile,wp)
  else:
    gx = readImage(gxfile)
    gu = readImage(gufile)
    gt = readImage(gtfile)
    gh = readImage(ghfile)
    wp = readImage(wpfile)
  '''
  plot3(gx,png="seismic")
  plot3(gu,png="flattenedC")
  plot3(gx,gt,cmin=min(gt),cmax=max(gt),cmap=jetRamp(1.0),
        clab="Relative geologic time",png="rgtC")
  plot3(gx,wp,cmin=0.0,cmax=1.0,cmap=jetRamp(1.0),
        clab="Weights",png="weights")
  '''

  logs = getLogs()
  gxs = zerofloat(n1,len(logs))
  gus = zerofloat(n1,len(logs))
  for il, log in enumerate(logs):
    model = SynSeis.getModel(log)
    i2 = s2.indexOfNearest(model.x2)
    i3 = s3.indexOfNearest(model.x3)
    gxs[il] = gx[i3][i2]
    gus[il] = gu[i3][i2]
  plot1(s1,gxs,vlabel="depth (km)",png="originalSeisTraces")
  plot1(s1,gus,vlabel="Relative geologic time",png="seisTracesFlattenC")


  '''
  ha = []
  hs = [800,750,700,650,600,550,500,450,400,350,300]
  hs = mul(hs,0.002)
  for ih, h in enumerate(hs):
    print h
    ha.append(h)
    plot3(gx,hs=ha,png="horizon"+str(ih))
  '''

def goFlattenD():
  gx = readImage(gxfile)
  logs = getLogs()
  gxs = zerofloat(n1,len(logs),1)
  for il, log in enumerate(logs):
    model = SynSeis.getModel(log)
    i2 = s2.indexOfNearest(model.x2)
    i3 = s3.indexOfNearest(model.x3)
    gxs[0][il] = gx[i3][i2]
  #maxShift = 50
  #errorPow = 2.0
  weight = 1.0
  maxShift = 40
  errorPow = 0.05

  wlw = WellLogWarping()
  wlw.setMaxShift(maxShift)
  wlw.setPowError([errorPow])
  s = wlw.findShifts([weight],gxs)
  gus = wlw.applyShifts(gxs[0],s)
  for i2 in range(len(gus)):
    for i1 in range(len(gus[0])):
      if(gus[i2][i1]<-10):
        gus[i2][i1] = 0.0
  plot1(s1,gxs[0],vlabel="depth (km)",png="originalSeisTraces")
  plot1(s1,gus,vlabel="Relative geologic time",png="seisTracesFlattenD")

def goDisplay():
  gx = readImage(gxfile)
  wl = getLogs()
  plot3(gx,logs=wl,curve="den",wmin=2,wmax=3,png="seisDen")
  plot3(gx,logs=wl,curve="vel",wmin=2,wmax=6,png="seisVel")

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

def normalize(x): 
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(100.0)
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

def convertDips(ft):
  return FaultScanner.convertDips(0.2,ft) # 5:1 vertical exaggeration

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

def wellGroup(logs,curve,cmin=0,cmax=0,cbar=None):
  print "number of logs =",len(logs)
  #s1 = Sampling(2762,0.002,0.000)
  #s2 = Sampling(357,0.025,0.000)
  #s3 = Sampling(161,0.025,0.000)
  fl,x1l,x2l,x3l = [],[],[],[]
  for log in logs:
    samples = log.getSamples(curve,s1,s2,s3)
    f,x1,x2,x3 = samples
    fl.append(f)
    x1l.append(x1)
    x2l.append(x2)
    x3l.append(x3)
  samples = fl,x1l,x2l,x3l
  lg = makeLogPoints(samples,cmin,cmax,cbar)
  return lg

def makeLogPoints(samples,cmin,cmax,cbar):
  lg = Group()
  fl,x1l,x2l,x3l = samples
  for i,f in enumerate(fl):
    f = fl[i]
    x1 = x1l[i]
    x2 = x2l[i]
    x3 = x3l[i]
    pg = makePointGroup(f,x1,x2,x3,cmin,cmax,cbar)
    lg.addChild(pg)
  return lg

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
          slices=None,surf=None,hs=None,logs=None,curve=None,
          wmin=0,wmax=0,png=None):
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
  if logs:
    wg = wellGroup(logs,curve,wmin,wmax)
    sf.world.addChild(wg)
  if hs:
    x1 = readImage(ghfile)
    u1 = readImage(gtfile)
    hfr = HorizonFromRgt(s1,s2,s3,x1,u1)
    for hi in hs:
      [xyz,rgb] = hfr.singleHorizon(hi)
      tg = TriangleGroup(True,xyz,rgb)
      sf.world.addChild(tg)
  if surf:
    tgs = Triangle()
    xyz = tgs.trianglesForSurface(surf,0,n1-1)
    tg  = TriangleGroup(True,xyz)
    sf.world.addChild(tg)
  ipg.setSlices(924,224,68)
  #ipg.setSlices(n1,0,n3) # use only for subset plots
  if cbar:
    sf.setSize(837,700)
  else:
    sf.setSize(700,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  ov = sf.getOrbitView()
  zscale = 0.8*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.1)
  ov.setAzimuthAndElevation(230,25)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(0.0,0.02,0.05))
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")


def plot1(s1,ys,hlabel="Seismic traces",vlabel="depth (km)",png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  yf = 10
  for y in ys:
    y = add(y,yf)
    pv = sp.addPoints(s1,y)
    pv.setLineColor(Color.BLACK)
    yf = yf+10
  sp.setHLimits(5.0,115)
  sp.setVLimits(0.15,1.5)
  sp.setSize(600,800)
  sp.setHLabel(hlabel)
  sp.setVLabel(vlabel)
  if png and pngDir:
    sp.paintToPng(300,7.0,pngDir+png+".png")

#############################################################################
run(main)
