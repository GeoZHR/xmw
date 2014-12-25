"""
Demonstrate 3D seismic image processing for faults and horizons
Author: Dave Hale, Colorado School of Mines
Version: 2014.06.17
"""
from uff import *
from schutils import *
setupForSubset("s2a")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

# Names and descriptions of files.
g0file  = "g0" # raw input image
gxfile  = "gx" # input image, after bilateral filtering
gsxfile = "gsx" # image after lsf with sharp faults
gwfile  = "gw" # image after unfaulting
epfile  = "ep" # eigenvalue-derived planarity
p2file  = "p2" # inline slopes
p3file  = "p3" # crossline slopes
flfile  = "fl" # fault likelihood
fifile  = "fic" # flattened image
gffile  = "gf" # flattened image
wsfile  = "ws" # weight image for flattening
fpfile  = "fp" # fault strike (phi)
ftfile  = "ft" # fault dip (theta)
fltfile = "flt" # fault likelihood thinned
fptfile = "fpt" # fault strike thinned
fttfile = "ftt" # fault dip thinned
fs1file = "fs1" # fault slip (1st component)
fs2file = "fs2" # fault slip (2nd component)
fs3file = "fs3" # fault slip (3rd component)
ft1file = "ft1" # fault slip interpolated (1st component)
ft2file = "ft2" # fault slip interpolated (2nd component)
ft3file = "ft3" # fault slip interpolated (3rd component)
fskbase = "fsk" # fault skin (basename only)
fslbase = "fsl" # fault skins after reskinning (basename only)
fsibase = "fsi" # fault skins after reskinning (basename only)
fsgbase = "fsg"
r1file = "r1"
r2file = "r2"
r3file = "r3"
hxfile = "hx"
ftcfile = "ftcm"
r1tfile = "r1c"
r2tfile = "r2t"
r3tfile = "r3t"
cpfile = "cp"
hxfile = "hx"
hxmfile = "hxm"
ftcfile = "ftcm"
ftcmfile = "ftcm"
r1tfile = "r1c"
uffile = "uf"
cpfile = "cp"


# These parameters control the scan over fault strikes and dips.
sigmaPhi,sigmaTheta = 8,40
minPhi,maxPhi = 0,360
minTheta,maxTheta = 65,85

# These parameters control the construction of fault skins.
lowerLikelihood = 0.01
upperLikelihood = 0.20
minSkinSize = 10000

# These parameters control the computation of fault dip slips.
minThrow =  0.0
maxThrow = 20.0

# Directory for saved pn images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
#pngDir = None
pngDir = "../../../png/sch/"

# We can avoid most computations entirely be setting plotOnly to True.
plotOnly = False

# Processing begins here. When experimenting with one part of this demo, we
# can comment out other parts that have already written results to files.
def main(args):
  #goDisplay()
  #goSlopes()
  #goScan()
  #goThin()
  #goSmooth()
  #goSkin()
  #goSlip()
  #goUnfault()
  #goUnfoldc()
  #goUnfaultc()
  #goUnfold()
  #goFlatten()
  #goDisplay()
  #goFS()
  goShow()

def goFS():
  print "goFaultSurfer ..."
  gx = readImage(gxfile)
  sk = readSkins(fskbase)
  cells = FaultSkin.getCells(sk)
  fs = FaultSurfer(n1,n2,n3,cells)
  sks = fs.applySurferM()
  writeSkins(fsgbase,sks)
  sks = readSkins(fsgbase)

  for i in range(len(sks)):
    skin=sks[i]
    cells=FaultSkin.getCells(skin)
    if(len(cells)>100000):
      plot3(gx,skins=[skin],clab=str(i))

def goShow():
  print "goFaultSurfer ..."
  gx = readImage(gxfile)

  sks = readSkins(fsgbase)

  '''
  sk = readSkins(fskbase)
  cells = FaultSkin.getCells(sk)
  fs = FaultSurfer(n1,n2,n3,cells)
  sks = fs.applySurferM()
  writeSkins(fskintp,sks)
  sks = readSkins(fskintp)

  plot3(gx,skins=sk,png="oldSkins")
  plot3(gx,skins=sks,png="newSkins")
  '''
  #for iskin,skin in enumerate(skss):
  for i in range(411):
    skin=sks[i]
    cells=FaultSkin.getCells(skin)
    if(len(cells)>40000):
      plot3(gx,skins=[skin],clab=str(i))
  ''' 

  for iskin,skin in enumerate(sk):
    plot3(gx,skins=[skin],links=True,png="oldSkin"+str(iskin))
  ''' 


def goDisplay():
  print "goDisplay ..."
  gx = readImage(gxfile)
  plot3(gx)
  #gx = slog(gx)
  #plot3(gx)

def goSlopes():
  print "goSlopes ..."
  gx = readImage(gxfile)
  if not plotOnly:
    sigma1,sigma2,sigma3,pmax = 16.0,2.0,2.0,2.0
    p2,p3,ep = FaultScanner.slopes(sigma1,sigma2,sigma3,pmax,gx)
    writeImage(p2file,p2)
    writeImage(p3file,p3)
    writeImage(epfile,ep)
  else:
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    ep = readImage(epfile)
  print "p2 min =",min(p2)," max =",max(p2)
  print "p3 min =",min(p3)," max =",max(p3)
  print "ep min =",min(ep)," max =",max(ep)
  plot3(gx,p2, cmin=-1,cmax=1,cmap=bwrNotch(1.0),
        clab="Inline slope (sample/sample)",png="p2")
  plot3(gx,p3, cmin=-1,cmax=1,cmap=bwrNotch(1.0),
        clab="Crossline slope (sample/sample)",png="p3")
  plot3(gx,sub(1,ep),cmin=0,cmax=1,cmap=jetRamp(1.0),
        clab="Planarity")

def goScan():
  print "goScan ..."
  def slog(f): # logarithmic gain to balance amplitudes
    return mul(sgn(f),log(add(1.0,abs(f))))
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  gx = readImage(gxfile)
  gx = slog(gx)
  if not plotOnly:
    gtx = FaultScanner.taper(50,0,0,gx);
    fsc = FaultScanner(sigmaPhi,sigmaTheta)
    fl,fp,ft = fsc.scan(minPhi,maxPhi,minTheta,maxTheta,p2,p3,gtx)
    writeImage(flfile,fl)
    writeImage(fpfile,fp)
    writeImage(ftfile,ft)
  else:
    fl = readImage(flfile)
    fp = readImage(fpfile)
    ft = readImage(ftfile)
  fp = convertStrikes(fp) # for display only
  ft = convertDips(ft) # for display only
  print "fl min =",min(fl)," max =",max(fl)
  print "fp min =",min(fp)," max =",max(fp)
  print "ft min =",min(ft)," max =",max(ft)
  plot3(gx,clab="Amplitude")
  plot3(gx,fl,cmin=0.01,cmax=1,cmap=jetRamp(1.0),
        clab="Fault likelihood",png="fl")
  plot3(gx,fp,cmin=0,cmax=360,cmap=hueFill(1.0),
        clab="Fault strike (degrees)",cint=45,png="fp")
  plot3(gx,ft,cmin=25,cmax=65,cmap=jetFill(1.0),
        clab="Fault dip (degrees)",png="ft")

def goThin():
  print "goThin ..."
  gx = readImage(gxfile)
  if not plotOnly:
    fl = readImage(flfile)
    fp = readImage(fpfile)
    ft = readImage(ftfile)
    flt,fpt,ftt = FaultScanner.thin([fl,fp,ft])
    writeImage(fltfile,flt)
    writeImage(fptfile,fpt)
    writeImage(fttfile,ftt)
  else:
    flt = readImage(fltfile)
    fpt = readImage(fptfile)
    ftt = readImage(fttfile)
  fpt = convertStrikes(fpt) # for display only
  ftt = convertDips(ftt) # for display only
  plot3(gx,clab="Amplitude")
  plot3(gx,flt,cmin=0.01,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  plot3(gx,fpt,cmin=-0.01,cmax=360,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=45,png="fpt")
  plot3(gx,ftt,cmin=25,cmax=65,cmap=jetFillExceptMin(1.0),
        clab="Fault dip (degrees)",png="ftt")

def goSmooth():
  print "goSmooth ..."
  gx = readImage(gxfile)
  if not plotOnly:
    flstop = 0.01
    fsigma = 8.0
    gx = readImage(gxfile)
    flt = readImage(fltfile)
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    gsx = FaultScanner.smooth(flstop,fsigma,p2,p3,flt,gx)
    writeImage(gsxfile,gsx)
  else:
    gsx = readImage(gsxfile)
  plot3(gx,clab="Amplitude",png="gx")
  plot3(gsx,clab="Amplitude",png="gsx")

def goSkin():
  print "goSkin ..."
  gx = readImage(gxfile)
  gsx = readImage(gsxfile)
  if not plotOnly:
    p2 = readImage(p2file)
    p3 = readImage(p3file)
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
  else:
    skins = readSkins(fskbase)
  plot3(gx,skins=skins,png="skins")
  #for iskin,skin in enumerate(skins):
  #  plot3(gx,skins=[skin],links=True,png="skin"+str(iskin))

def goSlip():
  print "goSlip ..."
  gx = readImage(gxfile)
  if not plotOnly:
    skins = readSkins(fskbase)
    gsx = readImage(gsxfile)
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    fsl = FaultSlipper(gsx,p2,p3)
    fsl.setOffset(2.0) # the default is 2.0 samples
    fsl.setZeroSlope(False) # True only to show the error
    fsl.computeDipSlips(skins,minThrow,maxThrow)
    print "  dip slips computed, now reskinning ..."
    print "  number of skins before =",len(skins),
    fsk = FaultSkinner() # as in goSkin
    fsk.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
    fsk.setMinSkinSize(minSkinSize)
    fsk.setMinMaxThrow(minThrow,maxThrow)
    skins = fsk.reskin(skins)
    print ", after =",len(skins)
    removeAllSkinFiles(fslbase)
    writeSkins(fslbase,skins)
    smark = -999.999
    s1,s2,s3 = fsl.getDipSlips(skins,smark)
    writeImage(fs1file,s1)
    writeImage(fs2file,s2)
    writeImage(fs3file,s3)
    t1,t2,t3 = fsl.interpolateDipSlips([s1,s2,s3],smark)
    writeImage(ft1file,t1)
    writeImage(ft2file,t2)
    writeImage(ft3file,t3)
  else:
    skins = readSkins(fslbase)
    s1 = readImage(fs1file)
    t1 = readImage(ft1file)
    t2 = readImage(ft2file)
    t3 = readImage(ft3file)
  plot3(gx,skins=skins,png="skinsfl")
  plot3(gx,skins=skins,smax=10.0,png="skinss1")
  plot3(gx,s1,cmin=0.0,cmax=10.0,cmap=jetFill(0.3),
        clab="Fault throw (samples)",png="gxs1")
  plot3(gx,t1,cmin=0.0,cmax=10.0,cmap=jetFill(0.3),
        clab="Vertical shift (samples)",png="gxt1")
  plot3(gx,t2,cmin=-2.0,cmax=2.0,cmap=jetFill(0.3),
        clab="Inline shift (samples)",png="gxt2")
  plot3(gx,t3,cmin=-2.0,cmax=2.0,cmap=jetFill(0.3),
        clab="Crossline shift (samples)",png="gxt3")

def goUnfault():
  print "goUnfault ..."
  gx = readImage(gxfile)
  if not plotOnly:
    ft1 = readImage(ft1file)
    ft2 = readImage(ft2file)
    ft3 = readImage(ft3file)
    gw = FaultSlipper.unfault([ft1,ft2,ft3],gx)
    writeImage(gwfile,gw)
  else:
    gw = readImage(gwfile)
  plot3(gx,clab="Amplitude")
  plot3(gw,clab="Amplitude",png="gw")
  slices = (370,159,34)
  plot3(gx,clab="Amplitude",slices=slices,png="gx159")
  plot3(gw,clab="Amplitude",slices=slices,png="gw159")
  slices = (370,208,34)
  plot3(gx,clab="Amplitude",slices=slices,png="gx208")
  plot3(gw,clab="Amplitude",slices=slices,png="gw208")
  slices = (370,288,34)
  plot3(gx,clab="Amplitude",slices=slices,png="gx288")
  plot3(gw,clab="Amplitude",slices=slices,png="gw288")
  """
  sf = SimpleFrame(AxesOrientation.XRIGHT_YOUT_ZDOWN)
  ipgw = sf.addImagePanels(s1,s2,s3,gw)
  ipgx = sf.addImagePanels(s1,s2,s3,gx)
  ipgw.setClips(-1.0,1.0)
  ipgx.setClips(-1.0,1.0)
  ov = sf.getOrbitView()
  ov.setScale(2.5)
  """

def goUnfaultc():
  if not plotOnly:
    ft = zerofloat(n1,n2,n3)
    gx = readImage(gxfile)
    cp  = zerofloat(n1,n2,n3)
    p2,p3,ep = FaultScanner.slopes(2.0,1.0,1.0,5.0,gx)
    skins = readSkins(fslbase)
    cfs = ConstraintsFromFaults(skins,ep)
    wp = pow(ep,2.0)
    cs = cfs.getWeightsAndConstraints(wp,cp)
    fm = cfs.getFaultMap()
    u1 = fillfloat(1.0,n1,n2,n3)
    u2 = fillfloat(0.0,n1,n2,n3)
    u3 = fillfloat(0.0,n1,n2,n3)
    p = array(u1,u2,u3,wp)
    flattener = FlattenerRTD(4.0,4.0)
    [r1,r2,r3] = flattener.computeShifts(True,fm,cs,p)
    flattener.applyShifts([r1,r2,r3],gx,ft)
    writeImage(r1tfile,r1)
    writeImage(r2tfile,r2)
    writeImage(r3tfile,r3)
    writeImage(ftcfile,ft)
    writeImage(cpfile,cp)
  else:
    r1 = readImage(r1tfile)
    r2 = readImage(r2tfile)
    r3 = readImage(r3tfile)
    ft = readImage(ftcfile)
    cp = readImage(cpfile)
    gx = readImage(gxfile)
  hmin,hmax,hmap = -3.0,3.0,ColorMap.GRAY

  plot3(cp,cmin=hmin,cmax=hmax,cmap=hmap,clab="ControlPointsM",png="cp")
  plot3(ft,cmin=hmin,cmax=hmax,cmap=hmap,clab="UnfaultC",png="ft")
  plot3(gx,r1,cmin=-5.0,cmax=8.0,cmap=jetFill(0.3),
        clab="Vertical shift (samples)",png="gxs1i")
  plot3(gx,r2,cmin=-2.0,cmax=2.0,cmap=jetFill(0.3),
        clab="Inline shift (samples)",png="gxs2i")
  plot3(gx,r3,cmin=-1.0,cmax=1.0,cmap=jetFill(0.3),
        clab="Crossline shift (samples)",png="gxs3i")

def goUnfold():
  gx = readImage(ftcfile)
  hx = zerofloat(n1, n2, n3)
  u1 = zerofloat(n1, n2, n3)
  u2 = zerofloat(n1, n2, n3)
  u3 = zerofloat(n1, n2, n3)
  ep = zerofloat(n1, n2, n3)
  lof = LocalOrientFilter(4.0,1.0)
  lof.applyForNormalPlanar(gx,u1,u2,u3,ep)
  pow(ep, 8.0, ep)
  p = array(u1, u2, u3, ep)
  flattener = FlattenerRT(6.0, 6.0)
  r = flattener.findShifts(p)
  flattener.applyShifts(r, gx, hx)
  writeImage(hxfile, hx)
  hmin, hmax, hmap = -1.0, 1.0, ColorMap.GRAY
  plot3(gx, cmin=hmin, cmax=hmax, cmap=hmap, clab="Amplitude", png="gx")
  plot3(hx, cmin=hmin, cmax=hmax, cmap=hmap, clab="Amplitude", png="hx")

def goFlatten():
  gx = readImage(ftcfile)
  sigma1,sigma2,sigma3,pmax = 4.0,2.0,2.0,5.0
  p2,p3,ep = FaultScanner.slopes(sigma1,sigma2,sigma3,pmax,gx)
  ws = pow(ep,6.0)
  fl = Flattener3()
  fl.setSmoothings(8.0,8.0);
  fl.setIterations(0.01,50);
  fm = fl.getMappingsFromSlopes(s1,s2,s3,p2,p3,ep)
  gf = fm.flatten(gx)
  writeImage(gffile,gf)
  plot3(gx,clab="Amplitude",png="gx")
  plot3(gf,clab="Amplitude",png="gf")

def goUnfoldc():
  uf = zerofloat(n1,n2,n3)
  cp = zerofloat(n1,n2,n3)
  gx = readImage(gxfile)
  gw = readImage(gwfile)
  u1 = zerofloat(n1,n2,n3)
  u2 = zerofloat(n1,n2,n3)
  u3 = zerofloat(n1,n2,n3)
  ep = zerofloat(n1,n2,n3)
  lof = LocalOrientFilter(4.0,2.0,2.0)
  lof.applyForNormalPlanar(gx,u1,u2,u3,ep)
  wp = copy(ep)
  skins = readSkins(fslbase)
  cfs = ConstraintsFromFaults(skins,wp)
  wp = pow(wp,4.0)
  cs = cfs.getWeightsAndConstraints(wp,cp)
  fm = cfs.getFaultMap()
  #u1 = fillfloat(1.0,n1,n2,n3)
  #u2 = fillfloat(0.0,n1,n2,n3)
  #u3 = fillfloat(0.0,n1,n2,n3)
  p = array(u1,u2,u3,wp)
  flattener = FlattenerRTD(4.0,4.0)
  r = flattener.computeShifts(fm,cs,p,cp)
  flattener.applyShifts(r,gx,uf)
  writeImage(r1file,r[0])
  writeImage(r2file,r[1])
  writeImage(r3file,r[2])
  writeImage(uffile,uf)
  hmin,hmax,hmap = -3.0,3.0,ColorMap.GRAY
  plot3(cp,cmin=hmin,cmax=hmax,cmap=hmap,clab="ControlPointsM",png="cp")
  plot3(uf,cmin=hmin,cmax=hmax,cmap=hmap,clab="Amplitude",png="uf")
  plot3(gx,r[0],cmin=0.0,cmax=10.0,cmap=jetFill(0.3),
        clab="Vertical shift (samples)",png="gxs1i")
  plot3(gx,r[1],cmin=-2.0,cmax=2.0,cmap=jetFill(0.3),
        clab="Inline shift (samples)",png="gxs2i")
  plot3(gx,r[2],cmin=-1.0,cmax=1.0,cmap=jetFill(0.3),
        clab="Crossline shift (samples)",png="gxs3i")

  

def goDisplay():
  gx = readImage(gxfile)
  cp = readImage(cpfile)
  gw = readImage(gwfile)
  ftc = readImage(ftcfile)
  r1 = readImage(r1tfile)
  '''
  r2 = readImage(r2tfile)
  r3 = readImage(r3tfile)
  '''
  ft1 = readImage(ft1file)
  #fs2 = readImage(fs1file)
  #fs3 = readImage(fs1file)
  hmin,hmax,hmap = -1.0,1.0,ColorMap.GRAY
  plot3(ftc,cmin=hmin,cmax=hmax,cmap=hmap,clab="UnfaultC",png="ftc")
  plot3(gw,cmin=hmin,cmax=hmax,cmap=hmap,clab="Unfault",png="gw")
  plot3(gx,cmin=hmin,cmax=hmax,cmap=hmap,clab="Amplitude",png="gx")
  plot3(cp,cmin=hmin,cmax=hmax,cmap=hmap,clab="ControlPoints",png="cp")
  plot3(gx,r1,cmin=-20.0,cmax=10.0,cmap=hueFill(0.3),
        clab="Vertical shift for unfaulting with constraints",png="gxs1i")
  '''
  plot3(gx,r2,cmin=-2.0,cmax=2.0,cmap=jetFill(0.3),
        clab="r2",png="gxs1i")
  plot3(gx,r3,cmin=-2.0,cmax=2.0,cmap=jetFill(0.3),
        clab="r3",png="gxs1i")
  '''
  plot3(gx,ft1,cmin=0.0,cmax=10.0,cmap=jetFill(0.3),
        clab="Vertical shift for unfaulting",png="gxs1i")



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

def convertStrikes(fp):
  return FaultScanner.convertStrikes(True,-90.0,fp)

def plot3(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          xyz=None,cells=None,skins=None,smax=0.0,slices=None,
          links=False,curve=False,trace=False,png=None):
  n1 = len(f[0][0])
  n2 = len(f[0])
  n3 = len(f)
  #sf = SimpleFrame(AxesOrientation.XRIGHT_YOUT_ZDOWN)
  sf = SimpleFrame(AxesOrientation.XRIGHT_YIN_ZDOWN)
  cbar = None
  if g==None:
    ipg = sf.addImagePanels(s1,s2,s3,f)
    if cmap!=None:
      ipg.setColorModel(cmap)
    if cmin!=None and cmax!=None:
      ipg.setClips(cmin,cmax)
    else:
      ipg.setClips(-1.0,1.0)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-1.0,1.0)
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
    xyz,uvw,rgb = FaultCell.getXyzUvwRgbForLikelihood(0.5,cmap,cells,True)
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
        xyz,uvw,rgb = skin.getCellXyzUvwRgbForThrow(size,cmap,True)
      else: # show fault likelihood
        cmap = ColorMap(0.0,1.0,ColorMap.JET)
        xyz,uvw,rgb = skin.getCellXyzUvwRgbForLikelihood(size,cmap,True)
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
  if slices:
    k1,k2,k3 = slices
  else:
    k1,k2,k3 = (370,105,34) # most plots use these
    #k1,k2,k3 = (370,150,0) # most plots use these
  ipg.setSlices(k1,k2,k3)
  if cbar:
    sf.setSize(985,700) # for sch data
    #sf.setSize(837,700) # for fake data
  else:
    sf.setSize(848,700) # for sch data
    #sf.setSize(700,700) # for fake data
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  ov.setEyeToScreenDistance(3018.87) # for consistency with brooks
  ov.setWorldSphere(BoundingSphere(0.5*n1,0.5*n2,0.5*n3,radius))
  #ov.setAzimuthAndElevation(25.0,20.0)
  ov.setAzimuthAndElevation(150.0,15.0)
  ov.setScale(1.5)
  #ov.setTranslate(Vector3(-0.182,-0.238,-0.012))
  ov.setTranslate(Vector3(-0.190,-0.168,-0.006))
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(137,1,pngDir+png+"cbar.png")

#############################################################################
run(main)
