"""
3D seismic image processing for faults
Author: Xinming Wu, Colorado School of Mines
Version: 2015.05.07
"""
from utils import *
sys.setrecursionlimit(1500)
setupForSubset("sub1")
#setupForSubset("unc")
s1,s2,s3 = getSamplings()
#n1,n2,n3 = s1.count,s2.count,s3.count
n1,n2,n3 = 362,951,591
s1,s2,s3 = Sampling(n1),Sampling(n2),Sampling(n3)

# Names and descriptions of image files used below.
gxfile = "gx" # input image
gsxfile = "gsx" # image after lsf with fault likelihoods
epfile  = "ep" # eigenvalue-derived planarity
p2file  = "p2" # inline slopes
p3file  = "p3" # crossline slopes
flfile  = "fl" # fault likelihood
fpfile  = "fp" # fault strike (phi)
ftfile  = "ft" # fault dip (theta)
fltfile = "flt" # fault likelihood thinned
fptfile = "fpt" # fault strike thinned
fttfile = "ftt" # fault dip thinned
fskbase = "fsk" # fault skin (basename only)
fslbase = "fsl" # fault skin fault estimating slips(basename only)
fskgood = "fsg" # fault skin with interpolated cells (basename only)
ft1file = "ft1" # fault slip interpolated (1st component)
ft2file = "ft2" # fault slip interpolated (2nd component)
ft3file = "ft3" # fault slip interpolated (3rd component)
fwsfile = "fws" # image after unfaulting
ulfile  = "ul" # unconformity likelihood
ultfile = "ult" # thinned unconformity likelihood
uncfile = "unc" # unconformity surface
fgfile = "fg" #flattened image
rgtfile = "rgt" #relative geologic time image
sx1file = "sx1"
sx2file = "sx2"
sx3file = "sx3"
sffile = "sf"
hz1file = "hz1"
hz2file = "hz2"
unc1file = "unc1"
unc2file = "unc2"


# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
#minPhi,maxPhi = 0,150
minPhi,maxPhi = 0,360
minTheta,maxTheta = 70,80
sigmaPhi,sigmaTheta = 4,20

# These parameters control the construction of fault skins.
# See the class FaultSkinner for more information.
lowerLikelihood = 0.2
upperLikelihood = 0.5
minSkinSize = 2000

# These parameters control the computation of fault dip slips.
# See the class FaultSlipper for more information.
minThrow = 0.0
maxThrow = 20.0


# Directory for saved png images. If None, png images will not be saved.
#pngDir = None
pngDir = "../../../png/f3d/"
plotOnly = True

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goFaults()
  #goSalts()
  #goUncs()
  goSurfaces()
  #goDisplay()
  #goSlopes()
  #goScan()
  #goThin()
  #goThinImages()
  #goSkin()
  #goReSkin()
  #goSmooth()
  #goSlip()
  #goUnfault()
  #goUnfaultS()
  #goUncScan()
  #goUncConvert()
  #goFlatten()
  #goHorizons()
def goSurfaces():
  gx = readImage3D(362,951,591,gxfile)
  sf = readImage3D(242,611,591,sffile)
  gx = gain(gx)
  skins = readSkins(fskgood)
  sk1 = []
  sk2 = []
  ks = 0
  for ski in skins:
    if(ski.size()>1500):
      sk2.append(ski)
      if ks<(len(skins)-11):
        sk1.append(ski)
    ks = ks+1
  hz2 = readImage2D(n2,n3,hz2file)
  uc1 = readImage2D(n2,n3,unc1file)
  uc2 = readImage2D(n2,n3,unc2file)
  '''
  plot3(gx,png="seismic")
  plot3(gx,skins=sk1,png="faults")
  plot3(gx,skins=sk1,slt=sf,png="faultSalts")
  plot3(gx,skins=sk2,slt=sf,png="faultSaltsM")
  plot3(gx,skins=sk2,slt=sf,hs=[hz2],png="faultSaltsHorizonM")
  plot3(gx,skins=sk2,slt=sf,hs=[hz2],uncs=[uc2],png="faultSaltsHorizonUncM")
  '''
  plot3(gx,skins=sk2,slt=sf,hs=[hz2],uncs=[uc2])

def goSalts():
  gx = readImage3D(362,951,591,gxfile)
  sf = readImage3D(242,611,591,sffile)
  gx = gain(gx)
  plot3(gx,slt=sf,clab="skinsNew")

def goUncs():
  gx = readImage3D(362,951,591,gxfile)
  uc = readImage3D(362,951,591,uncfile)
  us = zerofloat(362,951,591)
  rgf = RecursiveGaussianFilter(2.0)
  rgf.apply000(uc,us)
  uh = UnconformityHelper()
  uh.setUncValues(30,172,us)
  sub(us,min(us),us)
  div(us,max(us),us)
  ut = zerofloat(362,951,591)
  uh.thin(0.1,us,ut)
  uncs = uh.surfer(n2,n3,0.1,100000,ut,us)
  uh.surfaceUpdate(2,2,gx,uncs)
  writeImage(unc1file,uncs[0])
  writeImage(unc2file,uncs[1])
  plot3(gx,uncs=uncs)
  plot3(gx,ut,cmin=0.1,cmax=1,cmap=jetFillExceptMin(1.0),
      clab="Thinned unconformity likelihood",png="unct")
  '''
  plot3(gx,uc,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
      clab="Unconformity likelihood",png="unc")
  '''

def goDisplay():
  gx = readImage(gxfile)
  '''
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  '''
  gx = copy(n1-20,n2,n3,20,0,0,gx)
  '''
  fl = copy(n1-20,n2,n3,20,0,0,fl)
  fp = copy(n1-20,n2,n3,20,0,0,fp)
  ft = copy(n1-20,n2,n3,20,0,0,ft)
  '''
  writeImage(gxfile,gx)
  '''
  writeImage(flfile,fl)
  writeImage(fpfile,fp)
  writeImage(ftfile,ft)
  '''
def goUnc():
  print "goUnc ..."
  gx = readImage(gxfile)
  uc = readImage(uncfile)
  fs = FaultSelect()
  fs.setUncValues(30,240,uc)
  sub(uc,min(uc),uc)
  div(uc,max(uc),uc)
  print min(uc)
  print max(uc)
  gx = gain(gx)
  plot3(gx,uc,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
      clab="Unconformity likelihood",png="unc")

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
  ''' 
  plot3(gx,p2, cmin=-1,cmax=1,cmap=bwrNotch(1.0),
        clab="Inline slope (sample/sample)",png="p2")
  plot3(gx,p3, cmin=-1,cmax=1,cmap=bwrNotch(1.0),
        clab="Crossline slope (sample/sample)",png="p3")
  plot3(gx,sub(1,ep),cmin=0,cmax=1,cmap=jetRamp(1.0),
        clab="Planarity")
  ''' 

def goScan():
  print "goScan ..."
  if not plotOnly:
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
  else:
    gx = readImage(gxfile)
    gx = gain(gx)
    fl = readImage(flfile)
    fp = readImage(fpfile)
    ft = readImage(ftfile)
  ''' 
  plot3(gx,fl,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
      clab="Fault likelihood",png="fl")
  plot3(gx,fp,cmin=0,cmax=360,cmap=hueFill(1.0),
      clab="Fault strike (degrees)",cint=45,png="fp")
  plot3(gx,convertDips(ft),cmin=35,cmax=50,cmap=jetFill(1.0),
      clab="Fault dip (degrees)",png="ft")
  ''' 

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
  gx = gain(gx)
  plot3(gx,clab="Amplitude",png="gx")
  plot3(gx,fl,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
        clab="Fault likelihood",png="fl")
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  plot3(gx,fpt,cmin=0,cmax=180,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=45,png="fpt")
  plot3(gx,convertDips(ftt),cmin=35,cmax=50,cmap=jetFillExceptMin(1.0),
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

def goThinImages():
  gx = readImage(gxfile)
  gx = gain(gx)
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  fs = FaultSkinner()
  fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  cells = fs.findCells([fl,fp,ft])
  flt = fillfloat(0.0,n1,n2,n3)
  fpt = fillfloat(0.0,n1,n2,n3)
  ftt = fillfloat(0.0,n1,n2,n3)
  FaultCell.getFlThick(0.0,cells,flt)
  FaultCell.getFpThick(0.0,cells,fpt)
  FaultCell.getFtThick(0.0,cells,ftt)
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  plot3(gx,fpt,cmin=0,cmax=360,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=45,png="fpt")
  plot3(gx,convertDips(ftt),cmin=35,cmax=50,cmap=jetFillExceptMin(1.0),
        clab="Fault dip (degrees)",png="ftt")

def goSkin():
  print "goSkin ..."
  gx = readImage(gxfile)
  gx = gain(gx)
  if not plotOnly:
    fl = readImage(flfile)
    fp = readImage(fpfile)
    ft = readImage(ftfile)
    fs = FaultSkinner()
    fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
    fs.setMaxDeltaStrike(10)
    fs.setMaxPlanarDistance(0.2)
    fs.setMinSkinSize(minSkinSize)
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
  skins = readSkins(fskbase)
  flt = like(gx)
  FaultSkin.getLikelihood(skins,flt)
  plot3(gx,skins=skins)
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="fls")

def goReSkin():
  useOldCells=True
  gx = readImage(gxfile)
  fl = readImage(flfile)
  sk = readSkins(fskbase)
  if not plotOnly:
    fsx = FaultSkinnerX()
    cells = fsx.resetCells(sk)
    fsx = FaultSkinnerX()
    fsx.setParameters(10.0,10.0,2.0)
    fsx.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
    fsx.setMinSkinSize(minSkinSize)
    fsx.setMaxPlanarDistance(0.2)
    fsx.setSkinning(useOldCells)
    sks = fsx.findSkinsXX(cells,fl)
    removeAllSkinFiles(fskgood)
    writeSkins(fskgood,sks)
  skins = readSkins(fskgood)
  '''
  for ski in skins:
    if(ski.size()>3000):
      sks.append(ski)
  #FaultSkin.getLikelihood(skins,flt)
  plot3(gx)
  '''
  gx = gain(gx)
  ks = [46,43,38,34,31,30,28,17,16,13,10,
   95,93,88,87,79,77,76,75,73,72,67,
   66,63,62,60,58,145,144,139,134,132,
   131,129,124,121,115,112,110,105,208,
   206,200,199,197,192,175,173,172,163,
   161,160,159,157,152,150]
  sks = []
  flt = zerofloat(n1,n2,n3)
  for k in ks:
    ski = skins[k]
    if(ski.size()>5000):
      sks.append(ski)
      FaultSkin.getLikelihoods([ski],flt)
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),png="flc")
  plot3(gx,skins=sks,clab="skinsNew")

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
  y = like(x)
  div(x,sqrt(g),y)
  return y

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

def readImage3D(n1,n2,n3,basename):
  """ 
  Reads an image from a file with specified basename
  """
  datDir = "../../../data/seis/f3d/"
  fileName = datDir+basename+".dat"
  image = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def readImage2D(n1,n2,basename):
  """ 
  Reads an image from a file with specified basename
  """
  datDir = "../../../data/seis/f3d/"
  fileName = datDir+basename+".dat"
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image



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

def plot3(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          xyz=None,cells=None,skins=None,smax=0.0,slices=None,
          links=False,curve=False,trace=False,slt=None,uv=None,
          hs=None,uncs=None,uncx=None,png=None):
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
  if uv:
    m3 = len(uv)
    m2 = len(uv[0])
    m1 = len(uv[0][0])
    c1,c2,c3=Sampling(m1),Sampling(m2),Sampling(m3)
    mc = MarchingCubes(c1,c2,c3,uv)
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
  if slt:
    m3 = len(slt)
    m2 = len(slt[0])
    m1 = len(slt[0][0])
    c1,c2,c3=Sampling(m1),Sampling(m2),Sampling(m3)
    mc = MarchingCubes(c1,c2,c3,slt)
    ct = mc.getContour(0.0)
    xyz = ct.x
    dh = DisplayHelper()
    dh.resetCellPositions(0,340,120,xyz)
    tg = TriangleGroup(ct.i,xyz,ct.u)
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
    uh = UnconformityHelper()
    uc=readImage3D(n1,n2,n3,uncfile)
    rgf = RecursiveGaussianFilter(2.0)
    rgf.apply000(uc,uc)
    uh = UnconformityHelper()
    uh.setUncValues(30,172,uc)
    sub(uc,min(uc),uc)
    div(uc,max(uc),uc)
    #copy(n1-4,n2,n3,0,0,0,uc,4,0,0,ul)
    print s2.getCount()
    print s3.getCount()
    for unc in uncs:
      [xyz,rgb]=uh.buildTrigs(n1,s3,s2,-0.1,unc,uc)
      #[xyz,rgb]=uh.buildTrigs(362,s3,s2,0.01,unc,uc)
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
        ls.setWidth(6)
        ls.setSmooth(False)
        ss.add(ls)
        sf.world.addChild(lg)
  if hs:
    for hi in hs:
      tg = TriangleGroup(True,s3,s2,hi)
      tg.setColor(Color.YELLOW)
      sf.world.addChild(tg)
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
    size = 3.0
    if links:
      size = 0.65 
      ls = LineState()
      ls.setWidth(4.0)
      ls.setSmooth(True)
      ss.add(ls)
    ct = 0
    dh = DisplayHelper()
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
  ipg.setSlices(n1-2,25,n3)
  if cbar:
    sf.setSize(1037,700)
  else:
    sf.setSize(1010,750)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  zscale = 0.65*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.4)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(-0.08,0.15,0.05))
  #ov.setAzimuthAndElevation(-60.0,35.0)
  ov.setAzimuthAndElevation(-65.0,22.0)
  #ov.setAzimuthAndElevation(-56.0,40.0)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(1080,1,pngDir+png+"cbar.png")


#############################################################################
run(main)
