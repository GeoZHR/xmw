"""
Demonstrate 3D seismic image processing for faults and horizons
Author: Dave Hale, Colorado School of Mines
Version: 2014.07.17
"""

from fakeutils import *
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

# Names and descriptions of image files used below.
gxfile  = "gx" # input image (maybe after bilateral filtering)
gsxfile = "gsx" # image after lsf with fault likelihoods
epfile  = "ep" # eigenvalue-derived planarity
p2file  = "p2" # inline slopes
p3file  = "p3" # crossline slopes
p2kfile = "p2k" # inline slopes (known)
p3kfile = "p3k" # crossline slopes (known)
semfile = "sem"
flfile  = "fl" # fault likelihood
fpfile  = "fp" # fault strike (phi)
ftfile  = "ft" # fault dip (theta)
fltfile = "flt" # fault likelihood thinned
fptfile = "fpt" # fault strike thinned
fttfile = "ftt" # fault dip thinned
flefile  = "fle" # fault likelihood
fpefile  = "fpe" # fault strike (phi)
ftefile  = "fte" # fault dip (theta)
fletfile = "flet" # fault likelihood thinned
fpetfile = "fpet" # fault strike thinned
ftetfile = "ftet" # fault dip thinned
flvfile  = "flv" # fault likelihood
fpvfile  = "fpv" # fault strike (phi)
ftvfile  = "ftv" # fault dip (theta)
fs1file = "fs1" # fault slip (1st component)
fs2file = "fs2" # fault slip (2nd component)
fs3file = "fs3" # fault slip (3rd component)
fskbaseF = "fskf" # fault skin (basename only)
fskbaseE = "fske" # fault skin (basename only)
fsktv = "fst"
smfile = "sm"
cmfile = "cm"
semfile = "sem"


# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
minPhi,maxPhi = 0,360
minTheta,maxTheta = 65,85
sigmaPhi,sigmaTheta = 6,20

# These parameters control the construction of fault skins.
# See the class FaultSkinner for more information.
lowerLikelihood = 0.3
upperLikelihood = 0.5
minSkinSize = 2000

# These parameters control the computation of fault dip slips.
# See the class FaultSlipper for more information.
minThrow = -15.0
maxThrow =  15.0

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
plotOnly = True
pngDir = "../../../png/mef/fake/"

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goFakeData()
  #goScan()
  goSkinF()
  #goSemblance()
  #goOrientScan()
  #goThinE()
  #goSkinE()
  #goTv()
  #goSkinTv()
  #goTest()
def goTest():
  gx = readImage(gxfile)
  skins = readSkins(fskbaseF)
  plot3(gx,skins=skins)
  fcs = FaultSkin.getCells(skins)
  cells = []
  for ic in range(0,len(fcs),2):
    cells.append(fcs[ic])
  fs = FaultScanner(sigmaPhi,sigmaTheta)
  sp = fs.makePhiSampling(minPhi,maxPhi)
  dp = sp.getDelta()
  fr = FaultReconstructor(n1,n2,n3,cells)
  sks = fr.reskin(minSkinSize,dp)
  plot3(gx,skins=sks)

def goFakeData():
  #sequence = 'A' # 1 episode of faulting only
  #sequence = 'OA' # 1 episode of folding, followed by one episode of faulting
  sequence = 'OOOOOAAAAA' # 5 episodes of folding, then 5 of faulting
  #sequence = 'OAOAOAOAOA' # 5 interleaved episodes of folding and faulting
  nplanar = 4 # number of planar faults
  conjugate = True # if True, two large planar faults will intersect
  conical = False # if True, may want to set nplanar to 0 (or not!)
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
    gmin,gax,gmap = 0.0,1.4,ColorMap.JET
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
  gx = readImage(gxfile)
  if not plotOnly:
    p2 = readImage(p2file)
    p3 = readImage(p3file)
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
    fl = readImage(flfile)
    fp = readImage(fpfile)
    ft = readImage(ftfile)
  plot3(gx,fl,k3=56,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
      clab="Fault likelihood",png="fl")
  plot3(gx,fp,k3=56,cmin=0,cmax=360,cmap=hueFill(1.0),
      clab="Fault strike (degrees)",cint=45,png="fp")
  plot3(gx,convertDips(ft),k3=56,cmin=15,cmax=55,cmap=jetFill(1.0),
      clab="Fault dip (degrees)",png="ft")

def goSemblance():
  print "go semblance ..."
  gx = readImage(gxfile)
  if not plotOnly:
    lof = LocalOrientFilterP(16.0,1.0,1.0)
    ets = lof.applyForTensors(gx)
    lsf = LocalSemblanceFilter(3,4)
    sem = lsf.semblance(LocalSemblanceFilter.Direction3.VW,ets,gx)
    writeImage(semfile,sem)
  else:
    sem = readImage(semfile)
  sem=sub(1,sem)
  sem = pow(sem,0.4)
  sub(sem,min(sem),sem)
  div(sem,max(sem),sem)
  plot3(gx,sem,k3=56,cmin=0.2,cmax=1.0,cmap=jetRamp(1.0),clab="Semblance")

def goOrientScan():
  gx  = readImage(gxfile)
  sem = readImage(semfile)
  sem = sub(1,sem)
  sem = pow(sem,0.4)
  sub(sem,min(sem),sem)
  div(sem,max(sem),sem)
  if not plotOnly:
    fs = LocalOrientScanner(1,sigmaPhi,sigmaTheta)
    fl,fp,ft = fs.scan(minPhi,maxPhi,minTheta,maxTheta,sem)
    print "fl min =",min(fl)," max =",max(fl)
    print "fp min =",min(fp)," max =",max(fp)
    print "ft min =",min(ft)," max =",max(ft)
    sub(fl,min(fl),fl)
    div(fl,max(fl),fl)
    writeImage(flefile,fl)
    writeImage(fpefile,fp)
    writeImage(ftefile,ft)
  else:
    fl = readImage(flefile)
    fp = readImage(fpefile)
    ft = readImage(ftefile)
  plot3(gx,k3=56,clab = "Amplitude",png="gx")
  plot3(gx,sem,k3=56,cmin=0.2,cmax=1,cmap=jetRamp(1.0),
        clab="Fault attribute",png="fa")
  plot3(gx,fl,k3=56,cmin=0.2,cmax=1,cmap=jetRamp(1.0),
        clab="Enhanced fault attribute",png="fae")
  plot3(gx,fp,k3=56,cmin=0,cmax=360,cmap=hueFill(1.0),
        clab="Fault strike (degrees)",cint=45,png="fpe")
  plot3(gx,convertDips(ft),k3=56,cmin=25,cmax=65,cmap=jetFill(1.0),
        clab="Fault dip (degrees)",png="fte")

def goThinE():
  print "goThin ..."
  gx = readImage(gxfile)
  fl = readImage(flefile)
  fp = readImage(fpefile)
  ft = readImage(ftefile)
  flt,fpt,ftt = FaultScanner.thin([fl,fp,ft])
  writeImage(fletfile,flt)
  writeImage(fpetfile,fpt)
  writeImage(ftetfile,ftt)
  plot3(gx,flt,k3=56,cmin=0.2,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flet")
  plot3(gx,fpt,k3=56,cmin=0,cmax=360,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=45,png="fpet")
  plot3(gx,convertDips(ftt),k3=56,cmin=15,cmax=55,cmap=jetFillExceptMin(1.0),
        clab="Fault dip (degrees)",png="ftet")

def goSkinF():
  print "goSkin ..."
  gx = readImage(gxfile)
  if not plotOnly:
    fl = readImage(flfile)
    fp = readImage(fpfile)
    ft = readImage(ftfile)
    fs = FaultSkinner()
    fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
    fs.setMaxDeltaStrike(10)
    fs.setMaxPlanarDistance(0.1)
    fs.setMinSkinSize(minSkinSize)
    cells = fs.findCells([fl,fp,ft])
    writeCells("cells",cells)
    skins = fs.findSkins(cells)
    for skin in skins:
      skin.smoothCellNormals(4)
    print "total number of cells =",len(cells)
    print "total number of skins =",len(skins)
    print "number of cells in skins =",FaultSkin.countCells(skins)
    fd = FaultDisplay()
    print "fault skins load finish..."
    removeAllSkinFiles(fskbaseF)
    writeSkins(fskbaseF,skins)
    plot3(gx,cells=cells,png="cellsF")
  else:
    skins = readSkins(fskbaseF)
  cells = readCells("cells")
  plot3(gx,cells=cells)
  plot3(gx,skins=skins,png="skinsF")

def goSkinE():
  print "goSkin ..."
  gx = readImage(gxfile)
  if not plotOnly:
    fl = readImage(flefile)
    fp = readImage(fpefile)
    ft = readImage(ftefile)
    fs = FaultSkinner()
    fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
    fs.setMaxDeltaStrike(10)
    fs.setMaxPlanarDistance(0.1)
    fs.setMinSkinSize(minSkinSize)
    cells = fs.findCells([fl,fp,ft])
    skins = fs.findSkins(cells)
    for skin in skins:
      skin.smoothCellNormals(4)
    print "total number of cells =",len(cells)
    print "total number of skins =",len(skins)
    print "number of cells in skins =",FaultSkin.countCells(skins)
    fd = FaultDisplay()
    print "fault skins load finish..."
    removeAllSkinFiles(fskbaseE)
    writeSkins(fskbaseE,skins)
    plot3(gx,cells=cells,png="cellsE")
  else:
    skins = readSkins(fskbaseE)
  plot3(gx,skins=skins,png="skinsE")
  plot3(gx,skins=[skins[1],skins[3],skins[4]],png="skinsCrossE")

def goTv():
  gx = readImage(gxfile)
  if not plotOnly:
    fl = readImage(flefile)
    fp = readImage(fpefile)
    ft = readImage(ftefile)
    skins = readSkins(fskbaseE)
    fcs = FaultSkin.getCells(skins)
    cells=[]
    for ic in range(0,len(fcs),4):
      cells.append(fcs[ic])
    plot3(gx,cells=cells,png="cellsTv")
    tv3 = TensorVoting3()
    tv3.setSigma(30)
    tv3.setVoteWindow(15,20,20)
    sm,cm,fp,ft = tv3.applyVote(n1,n2,n3,cells)
    sm = pow(sm,0.4)
    cm = pow(cm,0.8)
    sub(sm,min(sm),sm)
    div(sm,max(sm),sm)
    sub(cm,min(cm),cm)
    div(cm,max(cm),cm)
    writeImage(cmfile,cm)
    writeImage(smfile,sm)
    writeImage(fpvfile,fp)
    writeImage(ftvfile,ft)
  else: 
    sm = readImage(smfile)
    cm = readImage(cmfile)
    fp = readImage(fpvfile)
    ft = readImage(ftvfile)
  plot3(gx,sm,k3=56,cmin=0.20,cmax=1.0,cmap=jetRamp(1.0),clab="Surfaceness",png="sm")
  plot3(gx,cm,k3=56,cmin=0.20,cmax=1.0,cmap=jetRamp(1.0),
        clab="Intersection",png="cm")

def goSkinTv():
  gx = readImage(gxfile)
  if not plotOnly:
    fl = readImage(smfile)
    fc = readImage(cmfile)
    fp = readImage(fpvfile)
    ft = readImage(ftvfile)
    fs = FaultSkinner()
    fs.setGrowLikelihoods(0.4,upperLikelihood)
    fs.setMaxDeltaStrike(10)
    fs.setMaxPlanarDistance(0.1)
    fs.setMinSkinSize(minSkinSize)
    cells = fs.findCells([fl,fp,ft])
    skins = fs.findSkins(cells)
    fcs = FaultSkin.getCells(skins)
    cells = []
    for ic in range(0,len(fcs),4):
      cells.append(fcs[ic])
    fs = FaultScanner(sigmaPhi,sigmaTheta)
    sp = fs.makePhiSampling(minPhi,maxPhi)
    st = fs.makeThetaSampling(minTheta,maxTheta)
    fsx = FaultSkinnerX()
    fsx.resetCells(cells)
    fsx.setGrowLikelihoods(0.3,upperLikelihood)
    fsx.setMinSkinSize(minSkinSize)
    fsx.setMaxPlanarDistance(0.1)
    fsx.setGaussWeights(10,2,sp,st)
    skins = fsx.findSkins(n1,n2,n3,cells)
    removeAllSkinFiles(fsktv)
    writeSkins(fsktv,skins)
  else:
    skins = readSkins(fsktv)
  print len(skins)
  '''
  fd = FaultDisplay()
  cells = FaultSkin.getCells(skins)
  flt = fillfloat(-0.001,n1,n2,n3)
  fd.getFlImage(cells,flt)
  flt = fd.setSlices(56,25,436,flt)
  plot3(gx,flt,cmin=0.15,cmax=0.8,cmap=jetRamp(1.0),clab="Surfaceness",png="smt")
  '''
  plot3(gx,skins=skins,png="skinsTv")
  plot3(gx,skins=[skins[1],skins[3]],png="skinsCrossTv")




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
def slog(f):
  return mul(sgn(f),log(add(1.0,abs(f))))

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

def plot3(f,g=None,k3=92,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          xyz=None,cells=None,skins=None,smax=0.0,slices=None,
          links=False,curve=False,trace=False,htgs=None,png=None):
  n3 = len(f)
  n2 = len(f[0])
  n1 = len(f[0][0])
  s1,s2,s3=Sampling(n1),Sampling(n2),Sampling(n3)
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
    cmap = ColorMap(0.2,1.0,ColorMap.JET)
    xyz,uvw,rgb = FaultCell.getXyzUvwRgbForLikelihood(0.7,cmap,cells,False)
    #xyz,uvw,rgb = FaultCell.getXyzUvwRgbForLikelihood(1.0,cmap,cells,False)
    qg = QuadGroup(xyz,uvw,rgb)
    qg.setStates(ss)
    sf.world.addChild(qg)
  if htgs:
    for htg in htgs:
      sf.world.addChild(htg)
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
        cmap = ColorMap(0.2,1.0,ColorMap.JET)
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
  ipg.setSlices(98,5,k3)
  #ipg.setSlices(85,5,102)
  #ipg.setSlices(n1,0,n3) # use only for subset plots
  if cbar:
    sf.setSize(837,700)
  else:
    sf.setSize(700,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.45*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  ov.setWorldSphere(BoundingSphere(0.5*n1-10,0.5*n2,0.5*n3-6,radius))
  ov.setAzimuthAndElevation(-65.0,25.0)
  ov.setTranslate(Vector3(0.0241,0.0517,0.0103))
  ov.setScale(1.2)
  #ov.setScale(1.3) #use only for subset plots
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")


#############################################################################
run(main)
