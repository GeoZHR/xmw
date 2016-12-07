"""
Demo of salt likelihoods and salt boundaries
Author: Xinming Wu, Colorado School of Mines
Version: 2015.12.19
"""

from utils import *

setupForSubset("3d")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

# Names and descriptions of image files used below.
gxfile  = "gx" # input image 
slfile  = "sl" # salt likelihoods
epfile  = "ep" # eigenvalue-derived planarity
sffile  = "sf" # salt indicator function
mkfile  = "mk" # mask file
sfcfile  = "sfc" # salt indicator function with constraints

pngDir = None
pngDir = getPngDir()
pngDir = "../../../png/slt/3d/"

plotOnly = True

def main(args):
  #goSaltLike()  # compute salt likelihoods
  goSaltSurfer() # compute salt surfaces


def goSaltLike():
  gx = readImage(gxfile)
  if not plotOnly:
    u1 = zerofloat(n1,n2,n3)
    u2 = zerofloat(n1,n2,n3)
    u3 = zerofloat(n1,n2,n3)
    ss = SaltScanner()
    '''
    lof = LocalOrientFilterP(4,2)
    ets = lof.applyForTensors(gx)
    ets.setEigenvalues(0.01,1.0,1.0)
    ep = ss.applyForPlanar(200,ets,gx)
    '''
    ep = readImage(epfile)
    lof = LocalOrientFilterP(12,8)
    lof.applyForNormal(gx,u1,u2,u3)
    sl = ss.saltLikelihood(8,ep,u1,u2,u3)
    writeImage(epfile,ep)
    writeImage(slfile,sl)
  else:
    ep = readImage(epfile)
    sl = readImage(slfile)
  plot3(gx,clab="Amplitude",png="gx")
  plot3(gx,ep,cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),clab="Planarity",png="ep")
  plot3(gx,sl,cmin=0.25,cmax=0.9,cmap=jetRamp(1.0),clab="Salt likelihood",png="sl")

def goSaltSurfer():
  gx = readImage(gxfile)
  ep = readImage(epfile)
  sl = readImage(slfile)
  if not plotOnly:
    u1 = zerofloat(n1,n2,n3)
    u2 = zerofloat(n1,n2,n3)
    u3 = zerofloat(n1,n2,n3)
    g1 = zerofloat(n1,n2,n3)
    g2 = zerofloat(n1,n2,n3)
    g3 = zerofloat(n1,n2,n3)
    lof = LocalOrientFilterP(8,2)
    lof.applyForNormal(ep,u1,u2,u3)
    ss = SaltSurfer()
    fc = ss.findPoints(0.3,sl,u1,u2,u3)
    sf = readImage(sffile)
    plot3(gx,sl,cmin=0.25,cmax=0.9,cells=fc,cmap=jetRamp(1.0),png="points")
    #plot3(gx,sl,cmin=0.25,cmax=0.9,cells=fc,cmap=jetRamp(1.0),fbs=sf,png="pointsSf")
    mul(u1,sl,g1)
    mul(u2,sl,g2)
    mul(u3,sl,g3)
    sps = ScreenPoissonSurfer()
    sps.setSmoothings(20,20,20)
    mk = sps.getScreenMark(n1,n2,n3,fc)
    sf = sps.saltIndicator(mk,g1,g2,g3)
    writeImage(sffile,sf)
    writeImage(mkfile,mk)
  else:
    sf = readImage(sffile)
  print min(sf)
  print max(sf)
  plot3(gx,sl,cmin=0.25,cmax=0.9,cmap=jetRamp(1.0),fbs=sf,png="saltSl")
  plot3(sf,cmin=-max(sf)+1,cmax=max(sf)-1,cmap=ColorMap.BLUE_WHITE_RED,
        clab="Indicator function",png="sf")
  plot3(gx,sf,cmin=-max(sf)+1,cmax=max(sf)-1,cmap=bwrRamp(1.0),
        clab="Indicator function",png="sfSeis")
  plot3(gx,sf,cmin=-max(sf)+1,cmax=max(sf)-1,cmap=bwrRamp(1.0),fbs=sf,
        clab="Indicator function",png="saltSf")

def goSaltSurferC():
  mk = readImage(mkfile)
  gx = readImage(gxfile)
  ep = readImage(epfile)
  sl = readImage(slfile)
  if not plotOnly:
    u1 = zerofloat(n1,n2,n3)
    u2 = zerofloat(n1,n2,n3)
    u3 = zerofloat(n1,n2,n3)
    g1 = zerofloat(n1,n2,n3)
    g2 = zerofloat(n1,n2,n3)
    g3 = zerofloat(n1,n2,n3)
    lof = LocalOrientFilterP(8,2)
    lof.applyForNormal(gx,u1,u2,u3)
    mul(u1,sl,g1)
    mul(u2,sl,g2)
    mul(u3,sl,g3)
    sps = ScreenPoissonSurferC()
    sps.setSmoothings(20,20,20)
    mk = sps.getScreenMark(n1,n2,n3,fc)
    sf = sps.saltIndicator(mk,k1,k2,k3,g1,g2,g3)
    writeImage(sfcfile,sf)
  else:
    sf = readImage(sfcfile)
  print min(sf)
  print max(sf)
  plot3(gx,sl,cmin=0.1,cmax=0.8,cmap=jetRamp(1.0),fbs=sf,png="saltSl")
  plot3(gx,sf,cmin=-max(sf)+1,cmax=max(sf)-1,cmap=bwrRamp(1.0),
        clab="Indicator function",png="sf")
  plot3(gx,sf,cmin=-max(sf)+1,cmax=max(sf)-1,cmap=bwrRamp(1.0),fbs=sf,
        clab="Indicator function",png="saltSf")

def goSalt():
  gx = readImage(gxfile)
  ae = readImage("aet")#file)
  ce = readImage("cet")#file)
  ca = abs(sub(ce,ae))
  cd = sub(ce,0.5)
  ad = sub(ae,0.5)
  plot3(gx,ca,cmin=0.01,cmax=max(ca)-0.1,cmap=jetRamp(1.0))#,fbs=cd,png="pss")
  plot3(gx,ce,cmin=0.1,cmax=1.0,cmap=jetRamp(1.0),fbs=cd,png="pss")
  plot3(gx,ae,cmin=0.1,cmax=1.0,cmap=jetRamp(1.0),fbs=ad,png="pss")


def goPSS():
  print "point set surface method ..."
  gx = readImage(gxfile)
  sk = readSkins(fskbase)
  fc = FaultSkin.getCells(sk)
  sn = SaltNormal()
  xu = sn.applyForNormals(fc)
  ps = PointSetSurface()
  #sf = ps.findScalarField(n1,n2,n3,fc)
  sf = ps.findScalarField(n1,n2,n3,xu)
  print min(sf)
  print max(sf)
  plot3(sf,cmin=min(sf),cmax=max(sf))
  plot3(gx,cells=fc,png="cells")
  plot3(gx,sf,cmin=min(sf),cmax=max(sf),cells=fc,fbs=sf,cmap=jetRamp(1.0),
    clab="PointSetSurface",png="pss")


def goFakeData():
  #sequence = 'A' # 1 episode of faulting only
  sequence = 'OA' # 1 episode of folding, followed by one episode of faulting
  #sequence = 'OOOOOAAAAA' # 5 episodes of folding, then 5 of faulting
  #sequence = 'OAOAOAOAOA' # 5 interleaved episodes of folding and faulting
  nplanar = 0 # number of planar faults
  conjugate = False # if True, two large planar faults will intersect
  conical = True # if True, may want to set nplanar to 0 (or not!)
  impedance = False # if True, data = impedance model
  wavelet = True # if False, no wavelet will be used
  noise = 0.5 # (rms noise)/(rms signal) ratio
  gx,p2,p3 = FakeData.seismicAndSlopes3d2014A(
      sequence,nplanar,conjugate,conical,impedance,wavelet,noise)
  writeImage(gxfile,gx)
  writeImage(p2kfile,p2)
  writeImage(p3kfile,p3)
  print "gx min =",min(gx)," max =",max(gx)
  print "p2 min =",min(p2)," max =",max(p2)
  print "p3 min =",min(p3)," max =",max(p3)
  gmin,gmax,gmap = -3.0,3.0,ColorMap.GRAY
  if impedance:
    gmin,gmax,gmap = 0.0,1.4,ColorMap.JET
  plot3(gx,cmin=gmin,cmax=gmax,cmap=gmap,clab="Amplitude",png="gx")
  #plot3(gx,p2,cmap=bwrNotch(1.0))
  #plot3(gx,p3,cmap=bwrNotch(1.0))

def goSlopes():
  print "goSlopes ..."
  gx = readImage(gxfile)
  sigma1,sigma2,sigma3,pmax = 16.0,1.0,1.0,5.0
  p2,p3,ep = FaultScanner.slopes(sigma1,sigma2,sigma3,pmax,gx)
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

def goScan():
  print "goScan ..."
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  gx = readImage(gxfile)
  gx = FaultScanner.taper(10,0,0,gx)
  fs = FaultScanner(sigmaPhi,sigmaTheta)
  fl,fp,ft = fs.scan(minPhi,maxPhi,minTheta,maxTheta,p2,p3,gx)
  print "fl min =",min(fl)," max =",max(fl)
  print "fp min =",min(fp)," max =",max(fp)
  print "ft min =",min(ft)," max =",max(ft)
  writeImage(flfile,fl)
  writeImage(fpfile,fp)
  writeImage(ftfile,ft)
  plot3(gx,clab="Amplitude")
  plot3(gx,fl,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
        clab="Fault likelihood",png="fl")
  plot3(gx,fp,cmin=0,cmax=360,cmap=hueFill(1.0),
        clab="Fault strike (degrees)",cint=45,png="fp")
  plot3(gx,convertDips(ft),cmin=25,cmax=65,cmap=jetFill(1.0),
        clab="Fault dip (degrees)",png="ft")

def goThin():
  print "goThin ..."
  gx = readImage(gxfile)
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  flt,fpt,ftt = FaultScanner.thin([fl,fp,ft])
  writeImage(fltfile,flt)
  writeImage(fptfile,fpt)
  writeImage(fttfile,ftt)
  plot3(gx,clab="Amplitude")
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  plot3(gx,fpt,cmin=0,cmax=360,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=45,png="fpt")
  plot3(gx,convertDips(ftt),cmin=25,cmax=65,cmap=jetFillExceptMin(1.0),
        clab="Fault dip (degrees)",png="ftt")

def goStat():
  def plotStat(s,f,slabel=None):
    sp = SimplePlot.asPoints(s,f)
    sp.setVLimits(0.0,max(f))
    sp.setVLabel("Frequency")
    if slabel:
      sp.setHLabel(slabel)
  fl = readImage(fltfile)
  fp = readImage(fptfile)
  ft = readImage(fttfile)
  fs = FaultScanner(sigmaPhi,sigmaTheta)
  sp = fs.getPhiSampling(minPhi,maxPhi)
  st = fs.getThetaSampling(minTheta,maxTheta)
  pfl = fs.getFrequencies(sp,fp,fl); pfl[-1] = pfl[0] # 360 deg = 0 deg
  tfl = fs.getFrequencies(st,ft,fl)
  plotStat(sp,pfl,"Fault strike (degrees)")
  plotStat(st,tfl,"Fault dip (degrees)")

def goSkin():
  print "goSkin ..."
  gx = readImage(gxfile)
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  fs = FaultSkinner()
  fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  fs.setMinSkinSize(minSkinSize)
  fs.setMaxDeltaStrike(10)
  fs.setMaxPlanarDistance(0.2)
  cells = fs.findCells([fl,fp,ft])
  skins = fs.findSkins(cells)
  for skin in skins:
    skin.smoothCellNormals(4)
  print "total number of cells =",len(cells)
  print "total number of skins =",len(skins)
  print "number of cells in skins =",FaultSkin.countCells(skins)
  removeAllSkinFiles(fskbase)
  writeSkins(fskbase,skins)
  plot3(gx,cells=cells,png="cells")
  plot3(gx,skins=skins,png="skins")
  for iskin,skin in enumerate(skins):
    plot3(gx,skins=[skin],links=True,png="skin"+str(iskin))


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
          xyz=None,cells=None,skins=None,fbs=None,smax=0.0,
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
  if fbs:
    mc = MarchingCubes(s1,s2,s3,fbs)
    ct = mc.getContour(0.0)
    tg = TriangleGroup(ct.i,ct.x,ct.u)
    states = StateSet()
    cs = ColorState()
    cs.setColor(Color.MAGENTA)
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
  ipg.setSlices(232,29,585)
  if cbar:
    sf.setSize(987,700)
  else:
    sf.setSize(850,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  zscale = 0.5*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.5)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(0.0,0.05,-0.08))
  ov.setAzimuthAndElevation(-58.0,36.0)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(137,1,pngDir+png+"cbar.png")

#############################################################################
run(main)
