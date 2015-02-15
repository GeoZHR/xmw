"""
Demonstrate 3D seismic image processing for faults and horizons
Author: Xinming Wu and Dave Hale, Colorado School of Mines
Version: 2014.02.02
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
flfile  = "fl" # fault likelihood
fpfile  = "fp" # fault strike (phi)
ftfile  = "ft" # fault dip (theta)
fltfile = "flt" # fault likelihood thinned
fptfile = "fpt" # fault strike thinned
fttfile = "ftt" # fault dip thinned
fs1file = "fs1" # fault slip (1st component)
fs2file = "fs2" # fault slip (2nd component)
fs3file = "fs3" # fault slip (3rd component)
fskbase = "fsk" # fault skin (basename only)
fskgood = "fsr" # fault skin (basename only)
fskslip = "fss" # fault skin (basename only)
frc1file = "rc1"
frc2file = "rc2"
frc3file = "rc3"
frs1file = "rs1"
frs2file = "rs2"
frs3file = "rs3"
fwpfile = "wp"
fcpfile = "cp"
fwcfile = "fwc"
fwsfile = "fws"
gwsfile = "gws"


# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
minPhi,maxPhi = 0,360
minTheta,maxTheta = 65,85
sigmaPhi,sigmaTheta = 4,20

# These parameters control the construction of fault skins.
# See the class FaultSkinner for more information.
lowerLikelihood = 0.2
upperLikelihood = 0.5
minSkinSize = 3000

# These parameters control the computation of fault dip slips.
# See the class FaultSlipper for more information.
minThrow = 0.01
maxThrow = 15.0

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
#pngDir = None
pngDir = "../../../png/uff/fake/"
plotOnly = False
# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  '''
  goFakeData()
  goSlopes()
  goScan()
  goThin()
  #goSmooth()
  goSkin()
  goReSkin()
  goSmooth()
  goSlip()
  '''
  #goUnfaultC()
  goUnfaultS()
  #goSlipTest()
def goSlipTest():
  smark=-99
  gx  = readImage(gxfile)
  p2  = readImage(p2file)
  p3  = readImage(p3file)
  gsx  = readImage(gsxfile)
  sks = readSkins(fskgood)
  sks = [sks[1]]
  fst = FaultSlipTest(sks[0])
  skt = fst.getSubSkin()
  fst.checkCellArrays(skt[0])
  skp = fst.applySlipper(skt,gsx,p2,p3)
  fsl = FaultSlipper(gsx,p2,p3)
  s1,s2,s3 = fsl.getDipSlips(skp,smark)
  plot3(gx,s1,cmin=-0.01,cmax=max(s1),cmap=jetFillExceptMin(1.0),
        clab="Fault throw (samples)",png="gxs1")
  plot3(gx,skins=skt,links=True,clab="subskin")
  plot3(gx,skins=skp,smax=10.0,clab="fault throw")
def goFakeData():
  #sequence = 'A' # 1 episode of faulting only
  #sequence = 'OA' # 1 episode of folding, followed by one episode of faulting
  #sequence = 'OOOOOAAAAA' # 5 episodes of folding, then 5 of faulting
  sequence = 'OAOAOAOAOA' # 5 interleaved episodes of folding and faulting
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
  '''
  plot3(gx,clab="Amplitude")
  plot3(gx,fl,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
        clab="Fault likelihood",png="fl")
  plot3(gx,fp,cmin=0,cmax=360,cmap=hueFill(1.0),
        clab="Fault strike (degrees)",cint=45,png="fp")
  plot3(gx,convertDips(ft),cmin=25,cmax=65,cmap=jetFill(1.0),
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
  '''
  plot3(gx,clab="Amplitude")
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  plot3(gx,fpt,cmin=0,cmax=360,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=45,png="fpt")
  plot3(gx,convertDips(ftt),cmin=25,cmax=65,cmap=jetFillExceptMin(1.0),
        clab="Fault dip (degrees)",png="ftt")
  '''

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
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  fs = FaultSkinner()
  fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
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
  plot3(gx,skins=skins,png="skins")
  for iskin,skin in enumerate(skins):
    plot3(gx,skins=[skin],links=True,png="skin"+str(iskin))

def goReSkin():
  print "goReSkin ..."
  flstop = 0.1
  fsigma = 8.0
  useOldCells = False
  gx = readImage(gxfile)
  fl = readImage(flfile)
  sk = readSkins(fskbase)
  fsx = FaultSkinnerX()
  fsx.setParameters(20,10,10)
  fsx.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  fsx.setMinSkinSize(minSkinSize)
  fsx.setMaxPlanarDistance(1.0)
  fsx.setSkinning(useOldCells)
  cells = FaultSkin.getCells(sk)
  fsx.resetCells(cells)
  skins = fsx.findSkinsXX(cells,fl)

  print "total number of cells =",len(cells)
  print "total number of skins =",len(skins)
  print "number of cells in skins =",FaultSkin.countCells(skins)
  removeAllSkinFiles(fskgood)
  writeSkins(fskgood,skins)
  plot3(gx,cells=cells,png="cells")
  plot3(gx,skins=skins,png="skins")
  for iskin,skin in enumerate(skins):
    plot3(gx,skins=[skin],links=True,png="skin"+str(iskin))


def goSmooth():
  print "goSmooth ..."
  flstop = 0.1
  fsigma = 8.0
  fl = readImage(flfile)
  gx = readImage(gxfile)
  skins = readSkins(fskgood)
  FaultSkin.setCells(skins,fl)
  flt = zerofloat(n1,n2,n3)
  fsx = FaultSkinnerX()
  fsx.getFl(skins,flt)

  p2,p3,ep = FaultScanner.slopes(8.0,1.0,1.0,5.0,gx)
  gsx = FaultScanner.smooth(flstop,fsigma,p2,p3,flt,gx)
  writeImage(p2file,p2)
  writeImage(p3file,p3)
  writeImage(epfile,ep)
  writeImage(gsxfile,gsx)
  plot3(gx)
  plot3(gsx)


def goSlip():
  print "goSlip ..."
  gx = readImage(gxfile)
  gsx = readImage(gsxfile)
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  #skins = readSkins(fskbase)
  skins = readSkins(fskgood)
  #skins = [skins[1]]
  '''
  cells = skins[0].getCellsLR()
  for fc in cells[25]:
    fc.fl = 1.0
  for iskin,skin in enumerate(skins):
    plot3(gx,skins=[skin],links=True,clab="skin"+str(iskin))
  '''

  fsl = FaultSlipper(gsx,p2,p3)
  fsl.setOffset(2.0) # the default is 2.0 samples
  fsl.setZeroSlope(False) # True only if we want to show the error
  fsl.computeDipSlips(skins,minThrow,maxThrow)
  print "  dip slips computed, now reskinning ..."
  print "  number of skins before =",len(skins),
  fsk = FaultSkinner() # as in goSkin
  fsk.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  fsk.setMinSkinSize(minSkinSize)
  fsk.setMinMaxThrow(minThrow,maxThrow)
  skins = fsk.reskin(skins)
  #skins = skins
  print ", after =",len(skins)
  #removeAllSkinFiles(fskbase)
  #writeSkins(fskbase,skins)
  removeAllSkinFiles(fskslip)
  writeSkins(fskslip,skins)
  smark = -999.999
  s1,s2,s3 = fsl.getDipSlips(skins,smark)
  writeImage(fs1file,s1)
  writeImage(fs2file,s2)
  writeImage(fs3file,s3)
  plot3(gsx)
  #plot3(gx,skins=skins,links=True,png="skinss1")
  plot3(gx,skins=skins,png="skinss1")
  plot3(gx,skins=skins,smax=10.0,png="skinss1")
  plot3(gx,s1,cmin=-0.01,cmax=10.0,cmap=jetFillExceptMin(1.0),
        clab="Fault throw (samples)",png="gxs1")
  s1,s2,s3 = fsl.interpolateDipSlips([s1,s2,s3],smark)
  plot3(gx,s1,cmin=0.0,cmax=10.0,cmap=jetFill(0.3),
        clab="Vertical shift (samples)",png="gxs1i")
  plot3(gx,s2,cmin=-2.0,cmax=2.0,cmap=jetFill(0.3),
        clab="Inline shift (samples)",png="gxs2i")
  plot3(gx,s3,cmin=-1.0,cmax=1.0,cmap=jetFill(0.3),
        clab="Crossline shift (samples)",png="gxs3i")
  gw = fsl.unfault([s1,s2,s3],gx)
  plot3(gx)
  plot3(gw,clab="Amplitude",png="gw")

def goUnfaultC():
  if not plotOnly:
    gx = readImage(gxfile)
    ep = readImage(epfile)
    fw = zerofloat(n1,n2,n3)
    cp = zerofloat(n1,n2,n3)
    skins = readSkins(fskslip)
    fsc = FaultSlipConstraints(skins)
    wp = pow(ep,2.0)
    ws = pow(ep,2.0)
    cs = fsc.controlPoints(ws,wp,cp)
    u1 = fillfloat(1.0,n1,n2,n3)
    u2 = fillfloat(0.0,n1,n2,n3)
    u3 = fillfloat(0.0,n1,n2,n3)
    lof = LocalOrientFilter(2.0,1.0,1.0)
    lof.applyForNormal(gx,u1,u2,u3)
    ps = array(u1,u2,u3,wp)
    flattener = FlattenerRTD(4.0,4.0)
    fsc.setNormals(ps)
    [r1,r2,r3] = flattener.computeShifts(True,cs,ws,ps)
    flattener.applyShifts([r1,r2,r3],gx,fw)
    writeImage(fwcfile,fw)
    writeImage(fwpfile,wp)
    writeImage(fcpfile,cp)
    writeImage(frc1file,r1)
    writeImage(frc2file,r2)
    writeImage(frc3file,r3)
  else :
    gx = readImage(gxfile)
    fw = readImage(fwcfile)
    wp = readImage(fwpfile)
    cp = readImage(fcpfile)
    r1 = readImage(frc1file)
    r2 = readImage(frc2file)
    r3 = readImage(frc3file)
  plot3(wp)
  plot3(gx,cp,cmin=0,cmax=6,cmap=jetFill(0.3),
        clab="Vertical shift (samples)",png="gxs1i")
  plot3(gx)
  plot3(fw)
  plot3(gx,r1,cmin=-5.0,cmax=8.0,cmap=jetFill(0.3),
        clab="Vertical shift (samples)",png="gxs1i")
  plot3(gx,r2,cmin=-2.0,cmax=2.0,cmap=jetFill(0.3),
        clab="Inline shift (samples)",png="gxs2i")
  plot3(gx,r3,cmin=-1.0,cmax=1.0,cmap=jetFill(0.3),
        clab="Crossline shift (samples)",png="gxs3i")

def goUnfaultS():
  if not plotOnly:
    gx = readImage(gxfile)
    fw = zerofloat(n1,n2,n3)
    gw = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    u1 = fillfloat(1.0,n1,n2,n3)
    u2 = fillfloat(0.0,n1,n2,n3)
    u3 = fillfloat(0.0,n1,n2,n3)
    lof = LocalOrientFilter(2.0,1.0,1.0)
    lof.applyForNormalPlanar(gx,u1,u2,u3,ep)
    ws = pow(ep,4.0)
    wp = pow(ep,4.0)
    skins = readSkins(fskslip)
    fsc = FaultSlipConstraints(skins)
    cs = fsc.screenPoints(wp)
    cn = fsc.screenPointsNearFaults(n2,n3)
    ps = array(u1,u2,u3,copy(wp))
    flattener = FlattenerRTS(6.0,6.0)
    flattener.setIters(20,20)
    flattener.setScreenNearFaults(cn)
    fsc.setNormals(ps)
    fl = mul(pow(cs[3][0],4.0),1.0)
    [r1,r2,r3] = flattener.findShifts(cs[0],cs[1],cs[2],fl,ps)
    flattener.applyShifts([r1,r2,r3],cs[0],cs[1],gx,gw)
    ps = array(u1,u2,u3,copy(wp))
    [t1,t2,t3] = flattener.unfaultShifts(cs[0],cs[1],cs[2],fl,ps)
    flattener.applyShifts([t1,t2,t3],None,None,gx,fw)
    writeImage(gwsfile,gw)
    writeImage(fwsfile,fw)
    writeImage(fwpfile,wp)
    writeImage(frs1file,r1)
    writeImage(frs2file,r2)
    writeImage(frs3file,r3)
  else :
    gx = readImage(gxfile)
    fw = readImage(fwsfile)
    gw = readImage(gwsfile)
    wp = readImage(fwpfile)
    cp = readImage(fcpfile)
    r1 = readImage(frs1file)
    r2 = readImage(frs2file)
    r3 = readImage(frs3file)
  plot3(gx)
  h1 = readImage(fs1file)
  plot3(fw,clab="unfaulted")
  plot3(gw,clab="unfaultAndUnfold")
  u1 = zerofloat(n1,n2,n3)
  u2 = zerofloat(n1,n2,n3)
  u3 = zerofloat(n1,n2,n3)
  lof = LocalOrientFilter(4.0,1.0,1.0)
  lof.applyForNormal(gx,u1,u2,u3)
  plot2(s1,s2,gw[47],png="flatten2d")
  u1,h1 = u1[47],h1[47]
  x1,x2,x3=getThrowArrs(h1)
  plot2(s1,s2,gx[47],x1=x1,x2=x2,x3=x3,
        clab="Vertical component of fault throws",cint=0.1,png="throws2d")
  plot2(s1,s2,gx[47],g=u1,gmin=0.7,gmax=1.0,
        clab="Vertical component of normal vectors",cint=0.1,png="normals2d")
  plot2(s1,s2,gx[47],g=h1,gmin=0.1,gmax=max(h1),
        clab="Vertical component of normal vectors",cint=1.0,png="throwCbar2d")

  '''
  plot3(gx,r1,cmin=-5.0,cmax=8.0,cmap=jetFill(0.3),
        clab="Vertical shift (samples)",png="gxs1i")
  plot3(gx,r2,cmin=-2.0,cmax=2.0,cmap=jetFill(0.3),
        clab="Inline shift (samples)",png="gxs2i")
  plot3(gx,r3,cmin=-1.0,cmax=1.0,cmap=jetFill(0.3),
        clab="Crossline shift (samples)",png="gxs3i")
  '''

def array(x1,x2,x3=None,x4=None):
  if x3 and x4:
    return jarray.array([x1,x2,x3,x4],Class.forName('[[[F'))
  elif x3:
    return jarray.array([x1,x2,x3],Class.forName('[[[F'))
  else:
    return jarray.array([x1,x2],Class.forName('[[[F'))
def getThrowArrs(h):
  n2 = len(h)
  n1 = len(h[0])
  smark = -999.999
  d1,d2=1.0,1.0
  x1,x2,x3=[],[],[]
  for i2 in range(1,n2-1,1):
    for i1 in range(1,n1-1,1):
      hi = h[i2][i1]
      if(hi>0):
        x1.extend([i1*d1])
        x2.extend([i2*d2])
        x3.extend([hi*d1])
  return x1,x2,x3

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

def plot3(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          xyz=None,cells=None,skins=None,smax=0.0,
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
  #ipg.setSlices(95,5,51)
  ipg.setSlices(95,5,90)
  #ipg.setSlices(95,5,n3-1)
  if cbar:
    sf.setSize(837,700)
  else:
    sf.setSize(700,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  ov.setWorldSphere(BoundingSphere(0.5*n1,0.5*n2,0.5*n3,radius))
  ov.setAzimuthAndElevation(-55.0,25.0)
  ov.setTranslate(Vector3(0.0241,0.0517,0.0103))
  ov.setScale(1.2)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(137,1,pngDir+png+"cbar.png")

def plot2(s1,s2,f,g=None,x1=None,x2=None,x3=None,gmin=None,gmax=None,
          clab=None,cint=None,title=None,png=None):
  n2 = len(f)
  n1 = len(f[0])
  panel = panel2()
  if clab:
    cb = panel.addColorBar(clab)
    cb.setInterval(cint)
  else:
    cb = panel.addColorBar()
    cb.setInterval(1.0)
    cb.setLabel("Amplitude")
    #cb.setLabel("Vector difference")
  if title:
    panel.setTitle(title)
  panel.setColorBarWidthMinimum(120)
  panel.setVInterval(20)
  panel.setHInterval(20)
  panel.setHLabel("Crossline (Samples)")
  panel.setVLabel("Time (Samples)")
  pv = panel.addPixels(s1,s2,f)
  #pv = panel.addPixels(f)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setColorModel(ColorMap.GRAY)
  pv.setClips(-3.0,3.0)
  if g:
    alpha = 0.3
    pv = panel.addPixels(s1,s2,g)
    #pv = panel.addPixels(g)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    if gmin==None: gmin = min(g)
    if gmax==None: gmax = max(g)
    pv.setClips(gmin,gmax)
    pv.setColorModel(ColorMap.getJet(alpha))
    if gmin==0:
      updateColorModel(pv,1.0)
  if x1:
    cmap = ColorMap(ColorMap.getJet(1.0))
    x3 = sub(x3,min(x3))
    x3 = div(x3,max(x3)-min(x3))
    for i1 in range(len(x1)):
      x1i,x2i,x3i=x1[i1],x2[i1],x3[i1]
      pv = panel.addPoints([x1i],[x2i])
      pv.setMarkSize(4.0)
      pv.setMarkColor(cmap.getColor(x3i))
      pv.setLineStyle(PointsView.Line.NONE)
      pv.setMarkStyle(PointsView.Mark.FILLED_SQUARE)
  frame2(panel,png)

def updateColorModel(pv,alpha):
    n = 256
    r = zerobyte(n)
    g = zerobyte(n)
    b = zerobyte(n)
    a = zerobyte(n)
    icm = pv.getColorModel()
    icm.getReds(r)
    icm.getGreens(g)
    icm.getBlues(b)
    for i in range(n):
      ai = int(255.0*alpha*i/n)
      if ai>127:
        ai -= 256
      a[i] = ai
    icm = IndexColorModel(8,n,r,g,b,a)
    pv.setColorModel(icm)

def panel2():
  #panel = PlotPanel(1,1,
  #  PlotPanel.Orientation.X1DOWN_X2RIGHT,
  #  PlotPanel.AxesPlacement.NONE)
  panel = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT)
  return panel

def frame2(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setSize(n2*6+120,n1*5);
  frame.setFontSizeForSlide(1.0,0.9,16.0/9.0)
  frame.setVisible(True)
  if png:
    frame.paintToPng(720,3.3,pngDir+png+".png")
  return frame
#############################################################################
run(main)
