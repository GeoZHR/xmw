"""
Demonstrate simultaneous multiple-well ties
Author: Xinming Wu, Colorado School of Mines
Version: 2016.05.11
"""

from utils import *
setupForSubset("cfd2007")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count
d1,d2,d3 = s1.delta,s2.delta,s3.delta
f1,f2,f3 = s1.first,s2.first,s3.first

# Names and descriptions of image files used below.
sfile = "cfs" # input seismic image 
ssfile = "cfss" # smoothed seismic image 
#logType = "v"; logLabel = "Velocity (km/s)"; vmin,vmax = 2.4,5.0
logType = "d"; logLabel = "Density (g/cc)"; vmin,vmax = 2.2,2.8
gfile = "cfg"+logType # simple gridding with null for unknown samples
pfile = "cfp"+logType # values of nearest known samples
qfile = "cfq"+logType # output of blended gridder
tfile = "cft"+logType # times to nearest known samples
p2file = "p2"
p3file = "p3"
epfile = "ep"
gffile = "gf"
u1file = "u1"
fpfile = "fp"
flfile = "fl"
ftfile = "ft"
dwfile = "dw"
gufile = "gu"
fskbase = "fsk"
fskgood = "fsg"
fslbase = "fsl"

# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
minPhi,maxPhi = 180,300
minTheta,maxTheta = 75,85
sigmaPhi,sigmaTheta = 20,40

# These parameters control the construction of fault skins.
# See the class FaultSkinner for more information.
lowerLikelihood = 0.3
upperLikelihood = 0.7
minSkinSize = 20000

# These parameters control the computation of fault dip slips.
# See the class FaultSlipper for more information.
minThrow = 0.01
maxThrow = 15.0


# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
#pngDir = ".././../png/swt/print/"
pngDir = "../../../png/swt/slides/"
pngDir = None
plotOnly = True

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goSeisAndWells()
  #goSlopes()
  #goScan()
  #goSkin()
  #goReSkin()
  #gridNearest()
  gridBlendedP()
  #goInterp()

def goSeisAndWells():
  gx = readImage(sfile)
  x12,x13,w1s = getLog242()
  x22,x23,w2s = getLog281()
  mds=[]
  mds.append(SynSeis.getModel(x12,x13,w1s[0],w1s[1],w1s[2]))
  mds.append(SynSeis.getModel(x22,x23,w2s[0],w2s[1],w2s[2]))
  swt = SeismicWellTie()
  sps = swt.getSamples(s1,mds)
  plot3(gx,sps=sps[1],wmin=2.2,wmax=2.8,clab="Density (g/cc)",png="seisDen")
  plot3(gx,sps=sps[0],wmin=2.4,wmax=5.0,clab="Velocity (km/s)",png="seisVel")

def goSlopes():
  print "goSlopes ..."
  gx = readImage(sfile)
  sigma1,sigma2,sigma3,pmax = 4.0,2.0,2.0,5.0
  p2,p3,ep = FaultScanner.slopes(sigma1,sigma2,sigma3,pmax,gx)
  zm = ZeroMask(0.3,5,1,1,gx)
  zero,tiny=0.0,0.01
  zm.setValue(zero,p2)
  zm.setValue(zero,p3)
  zm.setValue(tiny,ep)
  writeImage(p2file,p2)
  writeImage(p3file,p3)
  print "p2  min =",min(p2)," max =",max(p2)
  print "p3  min =",min(p3)," max =",max(p3)
  plot3(gx,p2, cmin=-1,cmax=1,cmap=bwrNotch(1.0),
        clab="Inline slope (sample/sample)",png="p2")
  plot3(gx,p3, cmin=-1,cmax=1,cmap=bwrNotch(1.0),
        clab="Crossline slope (sample/sample)",png="p3")
  plot3(gx,sub(1,ep),cmin=0,cmax=1,cmap=jetRamp(1.0),
        clab="Planarity")

def goScan():
  print "goScan ..."
  gx = readImage(sfile)
  if not plotOnly:
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    gx = FaultScanner.taper(10,0,0,gx)
    fs = FaultScanner(sigmaPhi,sigmaTheta)
    fl,fp,ft = fs.scan(minPhi,maxPhi,minTheta,maxTheta,p2,p3,gx)
    zm = ZeroMask(0.3,5,1,1,gx)
    zero=0.0
    zm.setValue(zero,fl)
    zm.setValue(zero,fp)
    zm.setValue(zero,ft)
    print "fl min =",min(fl)," max =",max(fl)
    print "fp min =",min(fp)," max =",max(fp)
    print "ft min =",min(ft)," max =",max(ft)
    writeImage(flfile,fl)
    writeImage(fpfile,fp)
    writeImage(ftfile,ft)
  else:
    fl = readImage(flfile)
    fp = readImage(fpfile)
    ft = readImage(ftfile)
  plot3(gx,clab="Amplitude")
  plot3(gx,fl,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
        clab="Fault likelihood",png="fl")
  plot3(gx,fp,cmin=minPhi,cmax=maxPhi,cmap=jetFill(1.0),
        clab="Fault strike (degrees)",cint=45,png="fp")
  plot3(gx,ft,cmin=minTheta,cmax=maxTheta,cmap=jetFill(1.0),
        clab="Fault dip (degrees)",png="ft")

def goSkin():
  print "goSkin ..."
  gx = readImage(sfile)
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  fs = FaultSkinner()
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(690):
        fl[i3][i2][i1] = 0
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
  plot3F(gx,cells=cells,png="cells")
  plot3F(gx,skins=skins,png="skins")

def goReSkin():
  print "goReSkin ..."
  useOldCells = True
  gx = readImage(sfile)
  if not plotOnly:
    fl = readImage(flfile)
    sk = readSkins(fskbase)
    fsx = FaultSkinnerX()
    fsx.setParameters(10,10,0.2)
    fsx.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
    fsx.setMinSkinSize(minSkinSize)
    fsx.setMaxPlanarDistance(0.2)
    fsx.setSkinning(useOldCells)
    cells = FaultSkin.getCells(sk)
    fsx.resetCells(cells)
    skins = fsx.findSkinsXX(cells,fl)
    removeAllSkinFiles(fskgood)
    writeSkins(fskgood,skins)
  skins = readSkins(fskgood)
  for skin in skins:
    skin.smoothCellNormals(4)
  plot3F(gx,skins=skins,png="skinsNew")
  plot3F(gx,skins=skins,links=True,png="skinsNewLinks")

def goSmooth():
  print "goSmooth ..."
  flstop = 0.1
  fsigma = 8.0
  gx = readImage(sfile)
  skins = readSkins(fskgood)
  flt = zerofloat(n1,n2,n3)
  fsx = FaultSkinnerX()
  fsx.getFl(skins,flt)
  p2,p3,ep = FaultScanner.slopes(4.0,2.0,2.0,5.0,gx)
  gsx = FaultScanner.smooth(flstop,fsigma,p2,p3,flt,gx)
  writeImage(p2file,p2)
  writeImage(p3file,p3)
  writeImage(epfile,ep)
  writeImage(ssfile,gsx)
  plot3(gx,flt,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
        clab="Fault likelihood",png="fli")
  plot3(gsx,png="gsx")


def goInterp():
  gx = readImage(sfile)
  tm = TensorMaker()
  mk = tm.mask(0.3,5.0,1.0,1.0,gx)
  et = tm.applyForTensors(4.0,2.0,2.0,mk,gx)
  et.setEigenvalues(0.0001,1.0,1.0)
  k1,k2,k3,fx=getSamples()
  wp = fillfloat(1.0,n1,n2,n3)
  fp = FastInterp(6.0,6.0)
  fp.setTensors(et)
  fp.setIterations(0.001,500)
  px = fp.interpolate(wp,k1,k2,k3,fx)
  writeImage(qfile,px)
  mds=[]
  x12,x13,w1s = getLog242()
  x22,x23,w2s = getLog281()
  mds.append(SynSeis.getModel(x12,x13,w1s[0],w1s[1],w1s[2]))
  mds.append(SynSeis.getModel(x22,x23,w2s[0],w2s[1],w2s[2]))
  swt = SeismicWellTie()
  sps = swt.getSamples(s1,mds)
  plot3(gx,px,cmin=2.2,cmax=2.8,sps=sps[1],wmin=2.2,wmax=2.8,
        clab="Density (g/cc)",png="seisDen")

def gridBlendedP():
  tm = TensorMaker()
  gx = readImage(sfile)
  mk = tm.mask(0.3,5.0,1.0,1.0,gx)
  et = tm.applyForTensors(4.0,2.0,2.0,mk,gx)
  #fs = FaultSkinnerX()
  #sks = readSkins(fskgood)
  #fls = fillfloat(0.01,n1,n2,n3)
  #fs.getFls(sks,fls)
  #et.scale(fls) # scale structure tensors by fls
  #et.invertStructure(1.0,1.0,1.0) # invert and normalize
  et.setEigenvalues(0.001,1.0,1.0)
  bi = BlendedGridder3(et)
  p = readImage(gfile)
  t = bi.gridNearest(0.0,p)
  writeImage(pfile,p)
  writeImage(tfile,t)

def makeImageTensors(s):
  """ 
  Returns tensors for guiding along features in specified image.
  """
  sigma = 3
  n1,n2 = len(s[0]),len(s)
  lof = LocalOrientFilter(sigma)
  t = lof.applyForTensors(s) # structure tensors
  c = coherence(sigma,t,s) # structure-oriented coherence c
  c = clip(0.0,0.99,c) # c clipped to range [0,1)
  t.scale(sub(1.0,c)) # scale structure tensors by 1-c
  t.invertStructure(1.0,1.0) # invert and normalize
  return t


def getSamples():
  mds=[]
  x12,x13,w1s = getLog242()
  x22,x23,w2s = getLog281()
  mds.append(SynSeis.getModel(x12,x13,w1s[0],w1s[1],w1s[2]))
  mds.append(SynSeis.getModel(x22,x23,w2s[0],w2s[1],w2s[2]))
  swt = SeismicWellTie()
  sps = swt.getSamples(s1,mds)
  i12 = round((x12-f2)/d2)
  i13 = round((x13-f3)/d3)
  i22 = round((x22-f2)/d2)
  i23 = round((x23-f3)/d3)
  i2s = [i12,i22]
  i3s = [i13,i23]
  k1s,k2s,k3s,fxs=[],[],[],[]
  for il in range(2):
    i2 = i2s[il]
    i3 = i3s[il]
    w1 = sps[0][1][il]
    wv = sps[0][0][il]
    wd = sps[1][0][il]
    for k1 in range(len(w1)):
      i1 = round((w1[k1]-f1)/d1)
      k1s.append(i1)
      k2s.append(i2)
      k3s.append(i3)
      if logType=="d":
        fxs.append(wd[k1])
      else:
        fxs.add(wv[k1])
  return k1s,k2s,k3s,fxs


def gridNearest():
  mds=[]
  x12,x13,w1s = getLog242()
  x22,x23,w2s = getLog281()
  mds.append(SynSeis.getModel(x12,x13,w1s[0],w1s[1],w1s[2]))
  mds.append(SynSeis.getModel(x22,x23,w2s[0],w2s[1],w2s[2]))
  swt = SeismicWellTie()
  sps = swt.getSamples(s1,mds)
  i12 = round((x12-f2)/d2)
  i13 = round((x13-f3)/d3)
  i22 = round((x22-f2)/d2)
  i23 = round((x23-f3)/d3)
  i2s = [i12,i22]
  i3s = [i13,i23]
  gvs = zerofloat(n1,n2,n3)
  gds = zerofloat(n1,n2,n3)
  for il in range(2):
    i2 = i2s[il]
    i3 = i3s[il]
    w1 = sps[0][1][il]
    wv = sps[0][0][il]
    wd = sps[1][0][il]
    for k1 in range(len(w1)):
      i1 = round((w1[k1]-f1)/d1)
      gvs[i3][i2][i1] = wv[k1]
      gds[i3][i2][i1] = wd[k1]
  writeImage("cfgv",gvs)
  writeImage("cfgd",gds)

def like(x):
  n2 = len(x)
  n1 = len(x[0])
  return zerofloat(n1,n2)

def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(10.0)
  ref.apply1(g,g)
  y = like(x)
  div(x,sqrt(g),y)
  return y
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

def plot1(s1,ys,hlabel="Seismic traces",vlabel="depth (km)",png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  for y in ys:
    pv = sp.addPoints(s1,y)
    pv.setLineColor(Color.BLACK)
  #sp.setVLimits(0.1,1.1)
  sp.setSize(800,800)
  sp.setHLabel(hlabel)
  sp.setVLabel(vlabel)
  if png and pngDir:
    sp.paintToPng(300,7.0,pngDir+png+".png")

def plot1s(s1,ss,ys,rs=None,vmin=None,vmax=None,color=Color.RED,
  hlabel="Log index",vlabel="Time (s)",png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sf = 1.0
  yf = sf
  sp.setVLimits(0.1,1.0)
  if vmin and vmax:
    sp.setVLimits(vmin,vmax)
  sp.setHLimits(0.5,11.5)
  sp.setHInterval(2)
  for il,y in enumerate(ys):
    ya = sum(y)/len(y)
    y = sub(y,ya)
    y = div(y,10)
    y = add(y,yf)
    pv = sp.addPoints(ss[il],y)
    pv.setLineColor(color)
    pv.setLineWidth(1.5)
    yf = yf+sf
  rf = sf
  if rs:
    for il,r in enumerate(rs):
      ra = sum(r)/len(r)
      r = sub(r,ra)
      r = div(r,10)
      r = add(r,rf)
      pv = sp.addPoints(s1,r)
      pv.setLineColor(Color.BLACK)
      pv.setLineWidth(1.5)
      rf = rf+sf
  sp.setSize(600,650)
  sp.setHLabel(hlabel)
  sp.setVLabel(vlabel)
  #sp.setFontSize(20) #for print
  sp.setFontSize(30) #for slides
  sp.setVInterval(0.2)
  if png and pngDir:
    sp.paintToPng(300,7.0,pngDir+png+".png")


def plot2(w,sz,sl,wmin=0.0,wmax=0.0,vlabel="Time (s)",cbar=None,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(500,900)
  sp.setVLabel(vlabel)
  sp.setHLabel("Log index")
  sp.addColorBar(cbar)
  sp.plotPanel.setColorBarWidthMinimum(90)
  pv = sp.addPixels(sz,sl,w)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ColorMap.GRAY)
  pv.setClips(wmin,wmax)
  if png and pngDir:
    sp.paintToPng(300,7.0,pngDir+png+".png")

def plot3F(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
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


def plot3X(s1,f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          slices=None,surf=None,hs=None,logs=None,sps=None,curve=None,
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
    if wmin!=0 and wmax!=0:
      ipg.setClips(wmin,wmax)
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
    cbar.setWidthMinimum(120) # for slides
    #cbar.setWidthMinimum(80)
  if logs:
    wg = wellGroup(logs,curve,wmin,wmax)
    sf.world.addChild(wg)
  if sps:
    #samples = sps[0],sps[1],sps[2],sps[3]
    wg = makeLogPoints(sps,wmin,wmax,cbar)
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
    sf.setSize(700,700) # for slides
    #sf.setSize(740,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  ov = sf.getOrbitView()
  zscale = 0.8*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.1)
  ov.setAzimuthAndElevation(235,25)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(0.0,0.05,0.08))
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.setFont(Font("Arial", Font.PLAIN, 36)) #for slides
      #cbar.setFont(Font("Arial", Font.PLAIN, 24)) #for print
      cbar.setInterval(0.5)
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")

def plot3(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          slices=None,surf=None,hs=None,logs=None,sps=None,curve=None,
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
    if wmin!=0 and wmax!=0:
      ipg.setClips(wmin,wmax)
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
      cmap = jetFill(1.0)
    ipg.setColorModel2(cmap)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMap2Listener(cbar)
    sf.world.addChild(ipg)
  if cbar:
    cbar.setWidthMinimum(120) # for slides
    #cbar.setWidthMinimum(80)
  if logs:
    wg = wellGroup(logs,curve,wmin,wmax)
    sf.world.addChild(wg)
  if sps:
    #samples = sps[0],sps[1],sps[2],sps[3]
    wg = makeLogPoints(sps,wmin,wmax,cbar)
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
    sf.setSize(700,700) # for slides
    #sf.setSize(740,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  ov = sf.getOrbitView()
  zscale = 0.8*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.1)
  ov.setAzimuthAndElevation(235,25)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(0.0,0.05,0.08))
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.setFont(Font("Arial", Font.PLAIN, 36)) #for slides
      #cbar.setFont(Font("Arial", Font.PLAIN, 24)) #for print
      cbar.setInterval(0.5)
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")

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

def makePoint(f,x1,x2,x3,cmin,cmax,cbar):
  xyz = zerofloat(3)
  xyz[0],xyz[1],xyz[2]=x3,x2,x1
  rgb = None
  if cmin<cmax:
    cmap = ColorMap(cmin,cmax,ColorMap.JET)
    if cbar:
      cmap.addListener(cbar)
    rgb = cmap.getRgbFloats([f])
  pg = PointGroup(xyz,rgb)
  ps = PointState()
  ps.setSize(4)
  ps.setSmooth(False)
  ss = StateSet()
  ss.add(ps)
  pg.setStates(ss)
  return pg

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

#############################################################################
run(main)
