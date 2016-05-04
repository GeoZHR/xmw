"""
Demonstrate 3D seismic image processing for faults and horizons
Author: Dave Hale, Colorado School of Mines
Version: 2014.07.17
"""

from utils2 import *
setupForSubset("oregan")
s1,s2 = getSamplings()
n1,n2 = s1.count,s2.count

# Names and descriptions of image files used below.
gxfile  = "gx" # input image (maybe after bilateral filtering)
gsxfile = "gsx" # image after lsf with fault likelihoods
epfile  = "ep" # eigenvalue-derived planarity
p2file  = "p2" # inline slopes
flfile  = "fl" # fault likelihood
ftfile  = "ft" # fault dip (theta)
fltfile = "flt" # fault likelihood thinned
fttfile = "ftt" # fault dip thinned

# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
minTheta,maxTheta = 80,89
sigmaTheta = 70

# These parameters control the construction of fault skins.
# See the class FaultSkinner for more information.
lowerLikelihood = 0.5
upperLikelihood = 0.7

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
pngDir = None
plotOnly = True
pngDir = "../../../png/oregan/"

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goSlopes()
  #goScan()
  goThin()

def goScan():
  print "goScan ..."
  gx = readImage(gxfile)
  if not plotOnly:
    gx = FaultScanner2.taper(10,0,gx)
    fs = FaultScanner2(sigmaTheta)
    sig1,sig2,smooth=8.0,36.0,48.0
    fl,ft = fs.scan(minTheta,maxTheta,sig1,sig2,smooth,gx)
    print "fl min =",min(fl)," max =",max(fl)
    print "ft min =",min(ft)," max =",max(ft)
    writeImage(flfile,fl)
    writeImage(ftfile,ft)
  else:
    fl = readImage(flfile)
    ft = readImage(ftfile)
  plot2(s1,s2,gx,g=fl,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
      label="Fault likelihood",png="fl")
  plot2(s1,s2,gx,g=abs(ft),cmin=minTheta,cmax=maxTheta,cmap=jetFill(1.0),
      label="Fault dip (degrees)",png="ft")

def goThin():
  print "goThin ..."
  gx = readImage(gxfile)
  if not plotOnly:
    fl = readImage(flfile)
    ft = readImage(ftfile)
    flt,ftt = FaultScanner2.thin([fl,ft])
    writeImage(fltfile,flt)
    writeImage(fttfile,ftt)
  else:
    flt = readImage(fltfile)
    ftt = readImage(fttfile)
  c1 = Sampling(n1,1,1200)
  for j2 in range(0,n2,1000):
    c2 = Sampling(1000,1,j2)
    gs  = copy(n1,1000,0,j2,gx)
    fls = copy(n1,1000,0,j2,flt)
    fts = copy(n1,1000,0,j2,ftt)
    plot2(c1,c2,gs,g=fls,cmin=0.9,cmax=1,cmap=jetFillExceptMin(1.0),
        png="section"+str(j2)+"~"+str(j2+1000))

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
  removeAllSkinFiles(fskbase)
  writeSkins(fskbase,skins)
  plot3(gx,cells=cells,png="cells")
  plot3(gx,skins=skins)
  for iskin,skin in enumerate(skins):
    plot3(gx,skins=[skin],links=True,)

def goReSkin():
  print "goReSkin ..."
  useOldCells = True
  gx = readImage(gxfile)
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
  plot3(gx,skins=skins,png="skinsNew")
  plot3(gx,skins=skins,links=True,png="skinsNewLinks")
  #plot3(gx,skins=[skins[2],skins[3]],png="skinsIntNew")
  '''
  for iskin,skin in enumerate(skins):
    plot3(gx,skins=[skin],links=True,png="skin"+str(iskin))
  '''

def goSmooth():
  print "goSmooth ..."
  flstop = 0.1
  fsigma = 8.0
  fl = readImage(flfile)
  gx = readImage(gxfile)
  skins = readSkins(fskgood)
  flt = zerofloat(n1,n2,n3)
  fsx = FaultSkinnerX()
  fsx.getFl(skins,flt)
  p2,p3,ep = FaultScanner.slopes(8.0,1.0,1.0,5.0,gx)
  gsx = FaultScanner.smooth(flstop,fsigma,p2,p3,flt,gx)
  writeImage(p2file,p2)
  writeImage(p3file,p3)
  writeImage(epfile,ep)
  writeImage(gsxfile,gsx)
  plot3(gx,flt,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
        clab="Fault likelihood",png="fli")
  plot3(gsx,png="gsx")

def goSlip():
  print "goSlip ..."
  gx = readImage(gxfile)
  gsx = readImage(gsxfile)
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  skins = readSkins(fskgood)
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
  print ", after =",len(skins)
  removeAllSkinFiles(fslbase)
  writeSkins(fslbase,skins)
  smark = -999.999
  s1,s2,s3 = fsl.getDipSlips(skins,smark)
  writeImage(fs1file,s1)
  writeImage(fs2file,s2)
  writeImage(fs3file,s3)
  plot3(gx,skins=skins,smax=10.0,slices=[85,5,60],png="skinss1")
  plot3(gx,s1,cmin=-10,cmax=10.0,cmap=jetFillExceptMin(1.0),
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

def goUnfaultS():
  if not plotOnly:
    gx = readImage(gxfile)
    fw = zerofloat(n1,n2,n3)
    lof = LocalOrientFilter(8.0,2.0,2.0)
    et = lof.applyForTensors(gx)
    et.setEigenvalues(0.001,1.0,1.0)

    wp = fillfloat(1.0,n1,n2,n3)
    skins = readSkins(fslbase)
    fsc = FaultSlipConstraints(skins)
    sp = fsc.screenPoints(wp)

    uf = UnfaultS(4.0,2.0)
    uf.setIters(100)
    uf.setTensors(et)
    mul(sp[3][0],10,sp[3][0])
    [t1,t2,t3] = uf.findShifts(sp,wp)
    uf.convertShifts(40,[t1,t2,t3])
    uf.applyShifts([t1,t2,t3],gx,fw)
    writeImage(fwsfile,fw)
    writeImage(sw1file,t1)
    writeImage(sw2file,t2)
    writeImage(sw3file,t3)
  else :
    gx = readImage(gxfile)
    fw = readImage(fwsfile)
  plot3(gx,png="gxuf")
  plot3(fw,png="fwuf")
  skins = readSkins(fslbase)
  mark = -999.99
  s1 = fillfloat(mark,n1,n2,n3)
  FaultSkin.getThrow(mark,skins,s1)
  plot3(gx,s1,cmin=-10,cmax=10.0,cmap=jetFillExceptMin(1.0),
        clab="Fault throw (samples)",png="gxs1")
  plot3(gx,t1,cmin=-6.0,cmax=6.0,cmap=jetFill(0.3),
        clab="Vertical shift (samples)",png="gxs1i")
  plot3(gx,t2,cmin=-2.0,cmax=2.0,cmap=jetFill(0.3),
        clab="Inline shift (samples)",png="gxs2i")
  plot3(gx,t3,cmin=-1.0,cmax=1.0,cmap=jetFill(0.3),
        clab="Crossline shift (samples)",png="gxs3i")

def goFlatten():
  fw = readImage(fwsfile)
  if not plotOnly:
    sig1,sig2,sig3,pmax=4.0,1.0,1.0,5.0
    p2,p3,ep = FaultScanner.slopes(sig1,sig2,sig3,pmax,fw)
    wp = pow(ep,6.0)
    fl = Flattener3()
    fl.setSmoothings(6.0,6.0)
    fl.setIterations(0.01,200)
    mp = fl.getMappingsFromSlopes(s1,s2,s3,p2,p3,wp)
    gu = mp.flatten(fw)
    x1 = mp.x1
    u1 = mp.u1
    writeImage(gufile,gu)
    writeImage(x1file,x1)
    writeImage(u1file,u1)
  else:
    gu = readImage(gufile)
    x1 = readImage(x1file)
    u1 = readImage(u1file)
  plot3(fw)
  plot3(gu,png="gu")
  plot3(fw,u1,cmin=10.0,cmax=n1,cmap=jetFill(1.0),
        clab="Relative geologic time (samples)",png="u1")

def goHorizonExtraction():
  gx = readImage(gxfile)
  u1 = readImage(u1file)
  x1 = readImage(x1file)
  w1 = readImage(sw1file)
  w2 = readImage(sw2file)
  w3 = readImage(sw3file)
  uf = UnfaultS(4.0,2.0)
  he = HorizonExtraction(s1,s2,s3,u1,x1)
  # extract a single horizon 
  u1i = 0.5*(s1.getFirst()+s1.getLast())
  hi = he.singleHorizon(u1i) # horizon extracted in unfaulted space
  uf.applyShiftsR([w1,w2,w3],hi,hi)# convert horizon back to original space
  stgs = he.applyForTgs(u1i,hi)
  plot3(gx,htgs=[stgs])
  # extract a set of horizons
  ft = s1.getFirst()+5
  dt = 10.0
  nt = (round((s1.getLast()-ft)/dt)-1)
  st = Sampling(nt,dt,ft)
  hs = he.multipleHorizons(st) # horizon extracted in unfaulted space
  for hk in hs:
    uf.applyShiftsR([w1,w2,w3],hk,hk)# convert horizon back to original space
  mtgs = he.applyForTgs(st,hs)
  plot3(gx,htgs=mtgs)

def goSlices():
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  gx  = readImage(gxfile)
  fw = readImage(fwsfile)
  gx  = gain(gx)
  fw  = gain(fw)
  flt = readImage(fltfile)
  sks = readSkins(fskgood)
  skl = readSkins(fslbase)
  fls = like(flt)
  fss = like(flt)
  FaultSkin.getLikelihood(sks,fls)
  FaultSkin.getThrow(skl,fss)
  fs = FaultSkinner()
  fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  fs.setMinSkinSize(minSkinSize)
  cls = fs.findCells([fl,fp,ft])
  plot3(gx,cells=cls,png="cells")
  plot3(gx,skins=sks,png="skins")
  plot3(gx,skins=skl,smax=max(fss)-10,png="throw")
  slt,sls,sfs = mul(100,flt),mul(100,fls),mul(100,fss)
  gxt,gxs,fws = sub(gx,slt), sub(gx,sls),sub(fw,sfs)
  gxw = sub(gx,sfs)
  d1 = 0.002
  mul(fss,d1*1000,fss)
  plot3f(gxt,a=flt,amin=0.01,amax=0.8,
        amap=jetFillExceptMin(1.0),alab="Fault likelihood",aint=0.1,png="flt")
  plot3f(gxs,a=fls,amin=0.01,amax=0.8,
        amap=jetFillExceptMin(1.0),alab="Fault likelihood",aint=0.1,png="fls")
  print max(fss)
  plot3f(gxw,a=fss,amin=0.01,amax=max(fss)-18,
        amap=jetFillExceptMin(1.0),alab="Vertical component of throw (ms)",
        aint=2.0,png="fss")
  plot3f(fws,a=fss,amin=0.01,amax=max(fss)-18,
        amap=jetFillExceptMin(1.0),alab="Vertical component of throw (ms)",
        aint=2.0,png="unfss")

def like(x):
  n2 = len(x)
  n1 = len(x[0])
  return zerofloat(n1,n2)



def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(80.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2)
  div(x,sqrt(g),y)
  return y

#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET
def jetFill(alpha):
  return ColorMap.setAlpha(ColorMap.JET,alpha)

def jetFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.JET,a)
def bwrNotch(alpha):
  a = zerofloat(256)
  for i in range(len(a)):
    if i<128:
      a[i] = alpha*(128.0-i)/128.0
    else:
      a[i] = alpha*(i-127.0)/128.0
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,a)


def bwrFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,a)

def jetRamp(alpha):
  return ColorMap.setAlpha(ColorMap.JET,rampfloat(0.0,alpha/256,256))

def bwrRamp(alpha):
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,rampfloat(0.0,alpha/256,256))

def grayRamp(alpha):
  return ColorMap.setAlpha(ColorMap.GRAY,rampfloat(0.0,alpha/256,256))


def plot2(s1,s2,f,g=None,cmin=None,cmax=None,cmap=None,label=None,png=None):
  n2 = len(f)
  n1 = len(f[0])
  f1,f2 = s1.getFirst(),s2.getFirst()
  d1,d2 = s1.getDelta(),s2.getDelta()
  panel = panel2Teapot()
  #panel.setHInterval(2.0)
  #panel.setVInterval(0.2)
  #panel.setHLabel("Inline (km)")
  #panel.setVLabel("Time (s)")
  #panel.setHInterval(100.0)
  #panel.setVInterval(100.0)
  #panel.setHLabel("Pixel")
  #panel.setVLabel("Pixel")
  '''
  if label:
    panel.addColorBar(label)
  else:
    panel.addColorBar()
  '''
  panel.setColorBarWidthMinimum(130)
  pv = panel.addPixels(s1,s2,f)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setColorModel(ColorMap.GRAY)
  if g:
    pv = panel.addPixels(s1,s2,g)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setColorModel(cmap)
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
  frame.setFontSize(12)
  frame.setSize(1180,950)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(400,3.2,pngDir+png+".png")
  return frame


#############################################################################
run(main)
