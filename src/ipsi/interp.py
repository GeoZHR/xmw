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
pufile = "pu"
sx1file = "sx1"
sx2file = "sx2"
sx3file = "sx3"
fwsfile = "fws"
rgtfile = "rgt"
fgfile = "fg"
fslbase = "fsl" # fault skin fault estimating slips(basename only)
fskgood = "fsg" # fault skin fault estimating slips(basename only)
sx1file = "sx1"
sx2file = "sx2"
sx3file = "sx3"
ft1file = "ft1" # fault slip interpolated (1st component)
ft2file = "ft2" # fault slip interpolated (2nd component)
ft3file = "ft3" # fault slip interpolated (3rd component)
uncfile = "unc" # unconformity surface


# Directory for saved png images. If None, png images will not be saved.
pngDir = None
pngDir = "../../../png/ipsi/"
plotOnly = False

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goInterp()
  #goTest()
  #goLog()
  goUnfault()
def goTest():
  if not plotOnly:
    gx = readImage(gxfile)

    lof = LocalOrientFilter(4.0,1.0,1.0)
    et = lof.applyForTensors(gx)
    et.setEigenvalues(0.001,1.0,1.0)

    skins = readSkins(fslbase)
    fmark = fillshort(0,n1,n2,n3)
    fsc = FaultSlipConstraints(n1,n2,n3,skins)
    fps = fsc.faultMarks(fmark)
    x1,x2,x3,fx=wellTeapot(n1,150,150)
    pnull = -999
    p = fillfloat(pnull,n1,n2,n3)
    for k in range(len(x1)):
      k1 = round(x1[k])
      k2 = round(x2[k])
      k3 = round(x3[k])
      p[k3][k2][k1] = fx[k]
    #bgx = BlendedGridder3(et,fx,x1,x2,x3)
    bgx = BlendedGridder3X(et,fx,x1,x2,x3)
    bgx.setFaults(fmark,fps[0],fps[1])
    t = bgx.gridNearest(pnull,p)
    samples = fx,x1,x2,x3
    print max(t)
    print max(p)
    plot3(gx,t,cmin=min(t),cmax=10,cmap=jetFill(0.5),
      clab="Time (samples)")
    plot3(gx,p,cmin=2.35,cmax=2.6,cmap=jetFill(0.5),
      clab="Density (g/cc)",samples=samples)
def goUnfault():
  smark = -999.999
  gx = readImage(gxfile)
  gx = gain(gx)
  skins = readSkins(fslbase)
  mark = -1000000
  fss = fillfloat(mark,n1,n2,n3)
  FaultSkin.getThrowThick(mark,skins,fss)
  fss = mul(fss,4) #convert to ms
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(n1):
        fssi = fss[i3][i2][i1]
        if fssi<=0.0 and fssi>mark:
          fss[i3][i2][i1]=0.05
  rgt = readImage(rgtfile)
  sx1 = readImage(sx1file)
  sx2 = readImage(sx2file)
  sx3 = readImage(sx3file)
  rw1 = readImage(ft1file)
  rw2 = readImage(ft2file)
  rw3 = readImage(ft3file)
  uncs = readUncs(uncfile)
  rgtx = zerofloat(n1,n2,n3)
  uf = UnfaultS(4.0,2.0)
  uf.applyShiftsX([sx1,sx2,sx3],rgt,rgtx)
  hfr = HorizonExtraction(s1,s2,s3,None,rgtx)
  k3,k2=249,25
  for unc in uncs:
    uf.applyShiftsR([rw1,rw2,rw3],unc,unc)
  sub(uncs,1,uncs)
  sks = readSkins(fslbase)
  uls = hfr.horizonCurves(uncs,k2,k3,sks)
  plot3p(gx,fss,cmin=-0.05,cmax=25,cmap=jetFillExceptMin(1.0),
          clab="Fault throw (ms)",cint=8,curve=True,uncx=uls,png="throwUnc")

def goLog():
  gx = readImage(gxfile)
  gw = readImage(fwsfile)
  gu = readImage(fgfile)
  sx1 = readImage(sx1file)
  sx2 = readImage(sx2file)
  sx3 = readImage(sx3file)
  rgt = readImage(rgtfile)
  sks = readSkins(fskgood)
  x1,x2,x3,fx=wellTeapot(n1,150,150)
  log = zerofloat(len(x1),1,4)
  for k in range(len(x1)):
    log[0][0][k] = fx[k]
    log[1][0][k] = x1[k]
    log[2][0][k] = x2[k]
    log[3][0][k] = x3[k]
  gx = gain(gx)
  gw = gain(gw)
  gu = gain(gu)
  cp = ConvertPoints()
  ps = cp.setUnfaultCoord(log,sks,sx1,sx2,sx3)
  ps = cp.setFlattenedCoord(s1,s2,s3,rgt,ps)
  fw,w1,w2,w3=cp.getSamplesW(ps)
  fu,u1,u2,u3=cp.getSamplesU(ps)
  samplesW = fw,w1,w2,w3
  samplesX = fx,x1,x2,x3
  samplesU = fu,u1,u2,u3
  plot3(gx,samples=samplesX,png="gxw")
  plot3(gw,samples=samplesW,png="gww")
  plot3(gu,samples=samplesU,png="guw")

def goInterp():
  gx = readImage(gxfile)
  gu = readImage(fgfile)
  if not plotOnly:
    rgt = readImage(rgtfile)
    den = readImage(denfile)
    sx1 = readImage(sx1file)
    sx2 = readImage(sx2file)
    sx3 = readImage(sx3file)
    rgtx = zerofloat(n1,n2,n3)
    uf = UnfaultS(4.0,2.0)
    uf.applyShiftsX([sx1,sx2,sx3],rgt,rgtx)
    x1,x2,x3,fx=wellTeapot(n1,150,150)
    samples = fx,x1,x2,x3
    rgi = RgtInterp3(fx,x1,x2,x3)
    rgi.setRgt(rgtx)
    rgi.setScales(0.01,1.00)
    pu,t,p,q = rgi.grid(s1,s2,s3)
    fu,u1,u2,u3 = rgi.getPoints(s1)
    samplesU = fu,u1,u2,u3
    '''
    tx = zerofloat(n1,n2,n3)
    px = zerofloat(n1,n2,n3)
    qx = zerofloat(n1,n2,n3)
    uf.applyShiftsX([sx1,sx2,sx3],t,tx)
    uf.applyShiftsX([sx1,sx2,sx3],p,px)
    uf.applyShiftsX([sx1,sx2,sx3],q,qx)
    '''
    print len(pu[0][0])
    pu = copy(n1,n2,n3,pu)
    writeImage(pufile,pu)
    writeImage(tfile,t)
    writeImage(pfile,p)
    writeImage(qfile,q)
  else:
    pu = readImage(pufile)
    t = readImage(tfile)
    p = readImage(pfile)
    q = readImage(qfile)
  gx = gain(gx)
  gu = gain(gu)
  k3 = str(249)
  x1,x2,x3,fx=wellTeapot(n1,150,150)
  samples = fx,x1,x2,x3
  plot3(gx,samples=samples,png="log"+k3)
  '''
  plot3(gx,rgtx,cmin=10.0,cmax=n1,cmap=jetFill(1.0),
        clab="Relative geologic time (samples)",png="rgt"+k3)
  '''
  plot3(gu,pu,cmin=2.35,cmax=2.6,cmap=jetFill(0.5),
    clab="Density (g/cc)",samples=samplesU,png="time"+k3)
  plot3(gx,p,cmin=2.35,cmax=2.6,cmap=jetFill(0.5),
    clab="Density (g/cc)",samples=samples,png="nearest"+k3)
  plot3(gx,q,cmin=2.35,cmax=2.6,cmap=jetFill(0.5),
    clab="Density (g/cc)",png="blended"+k3)

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

def plot3p(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          xyz=None,cells=None,skins=None,smax=0.0,slices=None,
          links=False,curve=False,trace=False,hs=None,uncs=None,uncx=None,png=None):
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
    uc=readImage(ulfile)
    uc = gain2(uc,12)
    uc = sub(uc,min(uc))
    uc = div(uc,max(uc))
    ul = zerofloat(n1,n2,n3)
    copy(n1-4,n2,n3,0,0,0,uc,4,0,0,ul)
    for unc in uncs:
      [xyz,rgb]=us.buildTrigs(n1,s3,s2,-0.1,unc,ul)
      #[xyz,rgb]=us.buildTrigs(n1,s3,s2,0.01,unc,ul)
      tg  = TriangleGroup(True,xyz,rgb)
      sg.addChild(tg)
    sf.world.addChild(sg)
  if uncx:
    for unc in uncx:
      if not curve:
        tg = TriangleGroup(True,unc[0])
        tg.setColor(Color.MAGENTA)
        sf.world.addChild(tg)
      else:
        lg = LineGroup(unc[0],unc[1])
        ss = StateSet()
        lg.setStates(ss)
        ls = LineState()
        ls.setWidth(8)
        ls.setSmooth(False)
        ss.add(ls)
        sf.world.addChild(lg)
  if hs:
    for hi in hs:
      if not curve:
        tg = TriangleGroup(True,hi[0],hi[1])
        sf.world.addChild(tg)
      else:
        lg = LineGroup(hi[0],hi[1])
        ss = StateSet()
        lg.setStates(ss)
        ls = LineState()
        ls.setWidth(2)
        ls.setSmooth(False)
        ss.add(ls)
        sf.world.addChild(lg)

  if skins:
    sg = Group()
    ss = StateSet()
    lms = LightModelState()
    lms.setLocalViewer(True)
    #lms.setTwoSide(True)
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
      ls.setWidth(4.0)
      ls.setSmooth(True)
      ss.add(ls)
    ct = 0
    for skin in skins:
      if smax>0.0: # show fault throws
        cmap = ColorMap(-1.0,smax,ColorMap.JET)
        xyz,uvw,rgb = skin.getCellXyzUvwRgbForThrow(size,cmap,False)
      else: # show fault likelihood
        cmap = ColorMap(0.0,1.0,ColorMap.JET)
        xyz,uvw,rgb = skin.getCellXyzUvwRgbForLikelihood(size,cmap,False)
      qg = QuadGroup(xyz,uvw,rgb)
      qg.setStates(None)
      sg.addChild(qg)
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
  ipg.setSlices(117,40,220)
  ipg.setSlices(110,25,249)
  #ipg.setSlices(115,25,167)
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
  #ov.setAzimuthAndElevation(-56.0,40.0)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")

#############################################################################
run(main)
