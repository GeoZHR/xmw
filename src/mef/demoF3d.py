"""
Demonstrate 3D seismic image processing for faults
Author: Xinming Wu, Colorado School of Mines
Version: 2016.01.28
"""

from utils import *
setupForSubset("f3d")
#setupForSubset("clyde")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

# Names and descriptions of image files used below.
gxfile  = "gx" # input image (maybe after bilateral filtering)
gsxfile = "gsx" # image after lsf with fault likelihoods
epfile  = "ep" # eigenvalue-derived planarity
eptfile  = "ept" # eigenvalue-derived planarity
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
flvfile  = "flv" # fault likelihood
fpvfile  = "fpv" # fault strike (phi)
ftvfile  = "ftv" # fault dip (theta)
fltvfile = "fltv" # fault likelihood thinned
fptvfile = "fptv" # fault strike thinned
fttvfile = "fttv" # fault dip thinned
fs1file = "fs1" # fault slip (1st component)
fs2file = "fs2" # fault slip (2nd component)
fs3file = "fs3" # fault slip (3rd component)
fskbase = "fsk" # fault skin (basename only)
fslbase = "fsl" # fault skin (basename only)
fskgood = "fsg" # fault skin (basename only)
fsktv = "fst" # fault skin (basename only)
fwsfile = "fws" # unfaulted image
sw1file = "sw1" # 1st component of unfaulting shifts
sw2file = "sw2" # 2nd component of unfaulting shifts
sw3file = "sw3" # 3rd component of unfaulting shifts
gufile = "gu" # flattened image
x1file = "x1" # horizon volume
u1file = "u1" # relateive geologic time volume
smfile = "sm"
cmfile = "cm"
semfile = "sem"


# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
minPhi,maxPhi = 0,360
minTheta,maxTheta = 75,85
sigmaPhi,sigmaTheta = 10,12

# These parameters control the construction of fault skins.
# See the class FaultSkinner for more information.
lowerLikelihood = 0.2
upperLikelihood = 0.5
minSkinSize = 2000

# These parameters control the computation of fault dip slips.
# See the class FaultSlipper for more information.
minThrow = 0.0
maxThrow = 25.0

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
pngDir = getPngDir()
pngDir = None
plotOnly = True

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  goPlanarity()
  #goSemblance()
  #goOrientScan()
  #goSlopes()
  #goScan()
  #goThin()
  #goTvThin()
  #goThinTv()
  #goSkin()
  #goCells()
  #goTv()
  #goSkinTv()
  #goTI()
  #goReSkin()
  #goSmooth()
  #goSlip()
  #goUnfaultS()
  #goFlatten()
  #goHorizonExtraction()
  #goComparison()
  '''
  gx = readImage(gxfile)
  fl = readImage(flfile)
  f2 = zerofloat(n1,n2,n3)
  f3 = zerofloat(n1,n2,n3)
  rgf = RecursiveGaussianFilterP(1.0)
  rgf.apply010(fl,f2)
  rgf.apply001(fl,f3)
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(n1):
        fli = fl[i3][i2][i1]
        if(fli<0.5):
          f3[i3][i2][i1]=-30
  plot3(gx, fbs=f3)
  '''

  '''
  gx = readImage(gxfile)
  plot3(gx)
  fx = copy(65,n2,n3,25,0,0,gx)
  writeImage(gxfile,fx)
  sem = readImage(semfile)
  gxs  = copy(300,300,300,50,400,0,gx)
  sems = copy(300,300,300,50,400,0,sem)
  writeImage("gxSub",gxs)
  writeImage("semSub",sems)
  '''
def goPlanarity():
  print "go planarity ..."
  gx = readImage(gxfile)
  lof = LocalOrientFilterP(3.0,1.0,1.0)
  ep = zerofloat(n1,n2,n3)
  u1 = zerofloat(n1,n2,n3)
  u2 = zerofloat(n1,n2,n3)
  u3 = zerofloat(n1,n2,n3)
  lof.applyForNormalPlanar(gx,u1,u2,u3,ep)
  plot3(gx,ep,cmin=0.1,cmax=1.0,cmap=jetRamp(1.0),clab="Semblance")
  plot3(ep)
  eps = zerofloat(n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      eps[i3][i2] = ep[i3][i2][56]
  writeImage("ep56",eps)

def goSemblance():
  print "go semblance ..."
  gx = readImage(gxfile)
  if not plotOnly:
    lof = LocalOrientFilterP(8.0,2.0,2.0)
    ets = lof.applyForTensors(gx)
    lsf = LocalSemblanceFilter(3,3)
    sem = lsf.semblance(LocalSemblanceFilter.Direction3.VW,ets,gx)
    writeImage(semfile,sem)
  else:
    sem = readImage(semfile)
  sem=sub(1,sem)
  plot3(gx,sem,cmin=0.1,cmax=1.0,cmap=jetRamp(1.0),clab="Semblance")
  plot3(sem)

def goOrientScan():
  gx = readImage(gxfile)
  sem = readImage(semfile)
  sem=sub(1,sem)
  if not plotOnly:
    fs = LocalOrientScanner(1,sigmaPhi,sigmaTheta)
    fl,fp,ft = fs.scan(minPhi,maxPhi,minTheta,maxTheta,sem)
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
  sub(fl,min(fl),fl)
  div(fl,max(fl),fl)
  plot3(gx,clab = "Amplitude",png="gx")
  plot3(gx,sem,cmin=0.1,cmax=1,cmap=jetRamp(1.0),
        clab="Fault attribute",png="fa")
  plot3(gx,fl,cmin=0.1,cmax=1,cmap=jetRamp(1.0),
        clab="Enhanced fault attribute",png="efa")
  plot3(gx,fp,cmin=0,cmax=360,cmap=hueFill(1.0),
        clab="Fault strike (degrees)",cint=45,png="fp")
  plot3(gx,convertDips(ft),cmin=25,cmax=65,cmap=jetFill(1.0),
        clab="Fault dip (degrees)",png="ft")

def goCells():
  gx = readImage(gxfile)
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  sub(fl,min(fl),fl)
  div(fl,max(fl),fl)
  fs = FaultSkinner()
  fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  fcs = fs.findCells([fl,fp,ft])
  cells=[]
  for ic in range(0,len(fcs),4):
    cells.append(fcs[ic])
  plot3(gx,cells=fcs,png="cells")
  plot3(gx,cells=cells,png="cellSub")

def goTv():
  gx = readImage(gxfile)
  if not plotOnly:
    fl = readImage(flfile)
    fp = readImage(fpfile)
    ft = readImage(ftfile)
    sub(fl,min(fl),fl)
    div(fl,max(fl),fl)
    fs = FaultSkinner()
    fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
    fcs = fs.findCells([fl,fp,ft])
    cells=[]
    for ic in range(0,len(fcs),4):
      cells.append(fcs[ic])
    tv3 = TensorVoting3()
    tv3.setSigma(20)
    tv3.setVoteWindow(10,20,20)
    sm,cm,fp,ft = tv3.applyVote(n1,n2,n3,cells)
    sm = pow(sm,0.4)
    sub(sm,min(sm),sm)
    div(sm,max(sm),sm)
    cm = pow(cm,0.4)
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
  plot3(gx,sm,cmin=0.15,cmax=0.8,cmap=jetRamp(1.0),clab="Surfaceness",png="sm")
  plot3(gx,cm,cmin=0.15,cmax=0.8,cmap=jetRamp(1.0),
        clab="Intersection",png="cm")

def goSkinTv():
  gx = readImage(gxfile)
  if not plotOnly:
    fl = readImage(smfile)
    fc = readImage(cmfile)
    fp = readImage(fpvfile)
    ft = readImage(ftvfile)
    fs = FaultScanner(sigmaPhi,sigmaTheta)
    sp = fs.makePhiSampling(minPhi,maxPhi)
    st = fs.makeThetaSampling(minTheta,maxTheta)
    fsx = FaultSkinnerX()
    fsx.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
    fsx.setMinSkinSize(minSkinSize)
    fsx.setMaxPlanarDistance(0.2)
    fcs = fsx.findCells([fl,fp,ft])
    cells = []
    for ic in range(0,len(fcs),8):
      cells.append(fcs[ic])
    fsx.resetCells(cells)
    fsx.setGaussWeights(sp,st)
    skins = fsx.findSkins(n1,n2,n3,cells)
    removeAllSkinFiles(fsktv)
    writeSkins(fsktv,skins)
  else:
    skins = readSkins(fsktv)
  print len(skins)
  fd = FaultDisplay()
  cells = FaultSkin.getCells(skins)
  flt = fillfloat(-0.001,n1,n2,n3)
  fd.getFlImage(cells,flt)
  flt = fd.setSlices(56,25,436,flt)
  plot3(gx,flt,cmin=0.15,cmax=0.8,cmap=jetRamp(1.0),clab="Surfaceness",png="smt")
  plot3(gx,skins=skins,png="skinsTv")

def goSlopes():
  print "goSlopes ..."
  gx = readImage(gxfile)
  sigma1,sigma2,sigma3,pmax = 8.0,4.0,4.0,5.0
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
  ep = sub(1,ep)
  plot3(gx,ep,cmin=min(ep),cmax=max(ep),cmap=jetRamp(1.0),
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
  plot3(gx,fl,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
        clab="Enhanced faults",png="fl")
  plot3(gx,fp,cmin=0,cmax=360,cmap=hueFill(1.0),
        clab="Fault strike (degrees)",cint=45,png="fp")
  plot3(gx,convertDips(ft),cmin=25,cmax=65,cmap=jetFill(1.0),
        clab="Fault dip (degrees)",png="ft")

def goThinTv():
  print "goThin ..."
  gx = readImage(gxfile)
  if not plotOnly:
    fl = readImage(smfile)
    fp = readImage(fpvfile)
    ft = readImage(ftvfile)
    flt,fpt,ftt = FaultScanner.thin([fl,fp,ft])
    writeImage(fltfile,flt)
    writeImage(fptfile,fpt)
    writeImage(fttfile,ftt)
  else:
    flt = readImage(fltvfile)
    fpt = readImage(fptvfile)
    ftt = readImage(fttvfile)
  plot3(gx,clab="Amplitude",png="gx")
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  plot3(gx,fpt,cmin=0,cmax=360,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=45,png="fpt")
  plot3(gx,convertDips(ftt),cmin=15,cmax=55,cmap=jetFillExceptMin(1.0),
        clab="Fault dip (degrees)",png="ftt")

def goThin():
  print "goThin ..."
  gx = readImage(gxfile)
  if not plotOnly:
    fl = readImage(flfile)
    sub(fl,min(fl),fl)
    div(fl,max(fl),fl)

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
  plot3(gx,clab="Amplitude",png="gx")
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  plot3(gx,fpt,cmin=0,cmax=360,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=45,png="fpt")
  plot3(gx,convertDips(ftt),cmin=15,cmax=55,cmap=jetFillExceptMin(1.0),
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
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  fs = FaultSkinner()
  sub(fl,min(fl),fl)
  div(fl,max(fl),fl)
  fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  fs.setMaxDeltaStrike(10)
  fs.setMaxPlanarDistance(0.2)
  fs.setMinSkinSize(minSkinSize)
  cells = fs.findCells([fl,fp,ft])
  skins = fs.findSkins(cells)
  print "total number of cells =",len(cells)
  print "total number of skins =",len(skins)
  print "number of cells in skins =",FaultSkin.countCells(skins)
  removeAllSkinFiles(fskbase)
  writeSkins(fskbase,skins)
  plot3(gx,cells=cells)
  plot3(gx,skins=skins)
  '''
  for iskin,skin in enumerate(skins):
    plot3(gx,skins=[skin],links=True,)
  '''

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
  #plot3(gx,skins=skins,png="skinsNew")
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

def plot3(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          xyz=None,cells=None,skins=None,smax=0.0,slices=None,
          links=False,curve=False,trace=False,htgs=None,fbs=None,png=None):
  n3 = len(f)
  n2 = len(f[0])
  n1 = len(f[0][0])
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
  if htgs:
    for htg in htgs:
      sf.world.addChild(htg)
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
        cmap = ColorMap(0.15,0.8,ColorMap.JET)
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
  ipg.setSlices(106,80,207)
  ipg.setSlices(56,25,436)
  #ipg.setSlices(115,25,167)
  if cbar:
    sf.setSize(987,720)
  else:
    sf.setSize(850,720)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  zscale = 0.3*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.4)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(0.0,0.20,0.12))
  ov.setAzimuthAndElevation(-56.0,35.0)
  #ov.setAzimuthAndElevation(-56.0,40.0)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")

#############################################################################
run(main)
