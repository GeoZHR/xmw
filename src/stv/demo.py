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
smtfile = "smt"


# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
minPhi,maxPhi = 0,360
minTheta,maxTheta = 65,85
sigmaPhi,sigmaTheta = 4,20

# These parameters control the construction of fault skins.
# See the class FaultSkinner for more information.
lowerLikelihood = 0.2
upperLikelihood = 0.5
minSkinSize = 2000

# These parameters control the computation of fault dip slips.
# See the class FaultSlipper for more information.
minThrow = -15.0
maxThrow =  15.0

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
pngDir = None
plotOnly = False

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goFakeData()
  #goSlopes()
  #goScan()
  #goThin()
  #goRescan()
  #goSkin()
  #goTest()
  #goKdTest()
  #goRescan()
  #goReSkin()
  '''
  goSmooth()
  goSlip()
  goUnfaultS()
  goFlatten()
  goHorizonExtraction()
  '''
  #goSubset()
  #gx = readImage(gxfile)
  #fl = readImage(flfile)
  #goSteer()
  #goRescanX()
  #phaseShift()
  #goHsurfer()
  #goTI()
  #goSemblance()
  #goBallVote()
  #goSemblanceThin()
  #goSemblanceTv()
  goVote()
  #voteScale()
def voteScale():
  tv3 = TensorVoting3()
  tv3.setVoteWindow(50,50,50)
  fc = FaultScanner(4,20)
  sp = fc.getPhiSampling(minPhi,maxPhi)
  st = fc.getThetaSampling(minTheta,maxTheta)
  scs = tv3.voteScale(sp,st,20)
  plot3(scs[0][0])

def goVote():
  gx = readImage(gxfile)
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  '''
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(n1):
        print ft[i3][i2][i1]
  '''
  fc = FaultScanner(4,20)
  sp = fc.getPhiSampling(minPhi,maxPhi)
  st = fc.getThetaSampling(minTheta,maxTheta)
  fs = FaultSkinner()
  fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  fs.setMaxDeltaStrike(10)
  fs.setMaxPlanarDistance(0.2)
  fs.setMinSkinSize(minSkinSize)
  fcs = fs.findCells([fl,fp,ft])
  cells = []
  for ic in range(0,len(fcs),5):
    cells.append(fcs[ic])
  tv3 = TensorVoting3()
  tv3.setVoteWindow(30,30,30)
  ss,cs,fp,ft = tv3.applyVoteFast(n1,n2,n3,15,sp,st,cells)
  plot3(gx,ss,cmin=0.0,cmax=1.0,cmap=jetRamp(1.0),
    clab="Surfaceness",png="sm")
  plot3(gx,fp,cmin=0,cmax=360,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=45,png="fpt")
  plot3(gx,convertDips(ft),cmin=15,cmax=55,cmap=jetFillExceptMin(1.0),
        clab="Fault dip (degrees)",png="ftt")

def goSemblance():
  print "go semblance ..."
  gx = readImage(gxfile)
  if not plotOnly:
    lof = LocalOrientFilterP(8.0,2.0,2.0)
    ets = lof.applyForTensors(gx)
    lsf = LocalSemblanceFilter(2,4)
    sem = lsf.semblance(LocalSemblanceFilter.Direction3.VW,ets,gx)
    sem=sub(1,sem)
    writeImage(semfile,sem)
  else:
    sem = readImage(semfile)
  plot3(gx,sem,cmin=min(sem),cmax=max(sem),cmap=jetRamp(1.0),clab="Semblance")


def goBallVote():
  gx = readImage(gxfile)
  fx = readImage(semfile)
  los = LocalOrientScanner(2,10)
  fl,fp,ft = los.scan(minPhi,maxPhi,minTheta,maxTheta,fx)
  writeImage(flvfile,fl)
  writeImage(fpvfile,fp)
  writeImage(ftvfile,ft)
  plot3(gx,fx,cmin=min(fx),cmax=max(fx),cmap=jetRamp(1.0),
    clab="Sem")
  plot3(gx,fl,cmin=min(fl),cmax=max(fl),cmap=jetRamp(1.0),
    clab="Fl")
  plot3(gx,fp,cmin=0,cmax=360,cmap=hueFill(1.0),
      clab="Fault strike (degrees)",cint=45,png="fp")
  ft = convertDips(ft)
  plot3(gx,ft,cmin=min(ft),cmax=max(ft),cmap=jetFill(1.0),
      clab="Fault dip (degrees)",png="ft")

def goSemblanceThin():
  print "goThin ..."
  gx = readImage(gxfile)
  fl = readImage(flvfile)
  fp = readImage(fpvfile)
  ft = readImage(ftvfile)
  flt,fpt,ftt = FaultScanner.thin([fl,fp,ft])
  writeImage(fltvfile,flt)
  writeImage(fptvfile,fpt)
  writeImage(fttvfile,ftt)
  plot3(gx,clab="Amplitude",png="gx")
  plot3(gx,fl,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
        clab="Fault likelihood",png="fl")
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  plot3(gx,fpt,cmin=0,cmax=360,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=45,png="fpt")
  plot3(gx,convertDips(ftt),cmin=15,cmax=55,cmap=jetFillExceptMin(1.0),
        clab="Fault dip (degrees)",png="ftt")

def goSemblanceTv():
  gx = readImage(gxfile)
  u1 = zerofloat(n1,n2,n3)
  u2 = zerofloat(n1,n2,n3)
  u3 = zerofloat(n1,n2,n3)
  sem = readImage(flfile)
  if not plotOnly:
    tv3 = TensorVoting3()
    fl = readImage(flfile)
    fp = readImage(fpfile)
    ft = readImage(ftfile)
    fs = FaultSkinner()
    fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
    fcs = fs.findCells([fl,fp,ft])
    plot3(gx,cells=fcs,png="cells")
    cells=[]
    for ic in range(0,len(fcs),5):
      cells.append(fcs[ic])
    tv3.setSigma(15)
    tv3.setWindow(30,20,20)
    sm,cm,u1,u2,u3 = tv3.applyVote(n1,n2,n3,cells)
    cells = tv3.findCells(0.2,sm,u1,u2,u3)
    plot3(gx,cells=cells,png="cells")
    fs = FaultSkinner()
    fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
    fs.setMaxDeltaStrike(10)
    fs.setMaxPlanarDistance(0.2)
    fs.setMinSkinSize(minSkinSize)
    skins = fs.findSkins(cells)
    removeAllSkinFiles(fsktv)
    writeSkins(fsktv,skins)
    writeImage(cmfile,cm)
    writeImage(smfile,sm)
  else: 
    sm = readImage(smfile)
    cm = readImage(cmfile)
    skins = readSkins(fsktv)
  plot3(gx,sem,cmin=0.0,cmax=1.0,cmap=jetRamp(1.0),clab="semblance",png="sm")
  plot3(gx,cm,cmin=0.2,cmax=1.0,skins=skins,cmap=jetRamp(1.0),png="skins")
  plot3(gx,sm,cmin=0.0,cmax=1.0,cmap=jetRamp(1.0),clab="Surfaceness",png="sm")
  #plot3(gx,cm,cmin=0.0,cmax=1.0,cmap=jetRamp(1.0),clab="Junction",png="cm")

def goTI():
  fx = randfloat(n1,n2,n3)
  gx = readImage(gxfile)
  sk = readSkins(fskbase)
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  fs = FaultSkinner()
  fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  fs.setMaxDeltaStrike(10)
  fs.setMaxPlanarDistance(0.5)
  fs.setMinSkinSize(minSkinSize)
  #fcs = fs.findCells([fl,fp,ft])
  fcs = FaultSkin.getCells(sk)
  plot3(gx,skins=sk,png="oldSkins")
  fls = zerofloat(n1,n2,n3)
  cells=[]
  for ic in range(0,len(fcs),5):
    cell = fcs[ic]
    cells.append(cell)
    ks = cell.getI()
    ms = cell.getIm()
    ps = cell.getIp()
    fls[ks[2]][ks[1]][ks[0]] = cell.getFl()
    fls[ms[2]][ms[1]][ms[0]] = cell.getFl()
    fls[ps[2]][ps[1]][ps[0]] = cell.getFl()
  plot3(gx,fls,cmin=min(fls),cmax=max(fls),cmap=jetRamp(1.0),
    clab="Fault likelihood",png="sm")
  ti = TensorInterp()
  ti.setParameters(80,5,0.4)
  sm,u1,u2,u3 = ti.apply(n1,n2,n3,cells)
  plot3(gx,fl,cmin=min(fl),cmax=max(fl),cmap=jetRamp(1.0),
    clab="Fault likelihood (resampled)",png="fl")
  plot3(gx,sm,cmin=min(sm),cmax=max(sm),cmap=jetRamp(1.0),
    clab="Salient map",png="sm")
  tv3 = TensorVoting3()
  cells = tv3.findCells(0.2,sm,u1,u2,u3)
  plot3(gx,cells=cells,png="cells")
  skins = fs.findSkins(cells)
  plot3(gx,skins=skins,png="skins")

def goFed():
  fx = readImage(gxfile)
  wp = fillfloat(1,n1,n2,n3)
  fs = fillfloat(1,n1,n2,n3)
  lof = LocalOrientFilter(4,4)
  ets = lof.applyForTensors(fx)
  ets.setEigenvalues(0.0001,1,1)
  fed = FastExplicitDiffusion()
  fed.setParameters(5,5,0.2)
  gx = fed.apply(ets,wp,fx)
  plot3(fx,clab="seismic")
  plot3(gx,cmin=min(gx),cmax=max(gx),clab="FED")
  lsf = LocalSmoothingFilter()
  lsf.apply(ets,20,fx,fs)
  plot3(fs,cmin=min(fs),cmax=max(fs),clab="LSF")

def goHsurfer():
  gx = readImage(gxfile)
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  '''
  plot3(gx,fp,cmin=0,cmax=360,cmap=hueFill(1.0),
      clab="Fault strike (degrees)",cint=45,png="fp")
  plot3(gx,convertDips(ft),cmin=15,cmax=55,cmap=jetFill(1.0),
      clab="Fault dip (degrees)",png="ft")
  '''
  flt = readImage(fltfile)
  fs = FaultSkinner()
  fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  fs.setMaxDeltaStrike(10)
  fs.setMaxPlanarDistance(0.2)
  fs.setMinSkinSize(minSkinSize)
  fcs = fs.findCells([fl,fp,ft])
  '''
  plot3(gx,flt,cmin=min(flt),cmax=max(flt),cmap=jetRamp(1.0),
    clab="flt",png="sm")
  '''
  sk = readSkins(fskbase)
  plot3(gx,skins=sk,png="skinOld")
  fcs = FaultSkin.getCells(sk)
  plot3(gx,cells=fcs,png="cells")
  fls = zerofloat(n1,n2,n3)
  for ic in range(0,len(fcs),1):
    cell = fcs[ic]
    ks = cell.getI()
    ms = cell.getIm()
    ps = cell.getIp()
    fls[ks[2]][ks[1]][ks[0]] = cell.getFl()
    fls[ms[2]][ms[1]][ms[0]] = cell.getFl()
    fls[ps[2]][ps[1]][ps[0]] = cell.getFl()
  plot3(gx,fls,cmin=min(fls),cmax=max(fls),cmap=jetRamp(1.0),
    clab="Fault likelihood",png="sm")

  cells=[]
  fl = zerofloat(n1,n2,n3)
  for ic in range(0,len(fcs),5):
    cell = fcs[ic]
    cells.append(cell)
    ks = cell.getI()
    ms = cell.getIm()
    ps = cell.getIp()
    fl[ks[2]][ks[1]][ks[0]] = cell.getFl()
    fl[ms[2]][ms[1]][ms[0]] = cell.getFl()
    fl[ps[2]][ps[1]][ps[0]] = cell.getFl()
  plot3(gx,cells=cells,png="cells")
  plot3(gx,fl,cmin=min(fl),cmax=max(fl),cmap=jetRamp(1.0),
    clab="Fault likelihood (resampled)",png="sm")
  tv3 = HermiteSurface()
  tv3.setSigma(20)
  tv3.setWindow(20,15,15)
  sm = tv3.applyVote(n1,n2,n3,cells)
  plot3(gx,sm,cmin=min(sm),cmax=max(sm),cmap=jetRamp(1.0),
    clab="Surfaceness",png="sm")

def goSteer():
  gx = readImage(gxfile)
  fl = readImage(fltfile)
  spyr = SteerablePyramid(0.5,1.0)
  x = copy(fl)
  y = copy(x)
  # Smooth locally-planar and plot.
  for i in range(0,10,1):
    pyr = spyr.makePyramid(y)
    attr = spyr.estimateAttributes(0,2.0,pyr)
    spyr.steerScale(0,0,50.0,0.5,attr,pyr)
    y = spyr.sumPyramid(1,pyr)
    plot3(gx,y,cmin=min(y),cmax=max(y),cmap=jetRamp(1.0),clab="fl"+str(i),png="sm")
  plot3(gx,fl,cmin=min(fl),cmax=max(fl),cmap=jetRamp(1.0),
    clab="linear",png="sm")
def phaseShift():
  gx = readImage(gxfile)
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  plot3(gx,clab="original")
  fs = FaultRescanner(sigmaPhi,sigmaTheta)
  for ph in range(10,180,10):
    px = fs.phaseShift(ph,gx)
    #se = fs.semblance(p2,p3,px)
    plot3(px,clab="phase="+str(ph))

def goRescanX():
  gx = readImage(gxfile)
  #sk = readSkins(fskbase)
  #fc = FaultSkin.getCells(sk)
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  fs = FaultSkinner()
  fc = fs.findCells([fl,fp,ft])
  cells=[]
  for ic in range(0,len(fc),1):
    cell = fc[ic]
    cells.append(cell)
  fs = FaultRescanner(sigmaPhi,sigmaTheta)
  sm,cm = fs.scan(n1,n2,n3,minPhi,maxPhi,minTheta,maxTheta,cells)
  plot3(gx,sm,cmin=min(sm),cmax=max(sm),cmap=jetRamp(1.0),
    clab="sm",png="sm")
  plot3(gx,cm,cmin=min(cm),cmax=max(cm),cmap=jetRamp(1.0),
    clab="cm",png="sm")



def goKdTest():
  gx = readImage(gxfile)
  sk = readSkins(fskbase)
  cl = FaultSkin.getCells(sk)
  g11 = fillfloat(0.0,n1,n2,n3)
  g12 = fillfloat(0.0,n1,n2,n3)
  g13 = fillfloat(0.0,n1,n2,n3)
  g22 = fillfloat(0.0,n1,n2,n3)
  g23 = fillfloat(0.0,n1,n2,n3)
  g33 = fillfloat(0.0,n1,n2,n3)
  tv = TensorVoting3()
  et = tv.mark(g11,g12,g13,g22,g23,g33,cl)
  lsf = LocalSmoothingFilter()
  lsf.applySmoothS(g11,g11)
  lsf.applySmoothS(g12,g12)
  lsf.applySmoothS(g13,g13)
  lsf.applySmoothS(g22,g22)
  lsf.applySmoothS(g23,g23)
  lsf.applySmoothS(g33,g33)
  lsf.apply(et,440,g11,g11)
  lsf.apply(et,440,g12,g12)
  lsf.apply(et,440,g13,g13)
  lsf.apply(et,440,g22,g22)
  lsf.apply(et,440,g23,g23)
  lsf.apply(et,440,g33,g33)
  sm,cm,jm = tv.solveEigenproblems(g11,g12,g13,g22,g23,g33)
  plot3(gx,sm,cmin=min(sm),cmax=max(sm),cmap=jetRamp(1.0),
    clab="sm",png="sm")
  plot3(gx,cm,cmin=min(cm),cmax=max(cm),cmap=jetRamp(1.0),
    clab="sm",png="sm")
  plot3(gx,sm,cmin=min(sm),cmax=max(sm),cells=cl,cmap=jetRamp(1.0),
    clab="sm",png="sm")


def goRescan():
  gx = readImage(gxfile)
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  plot3(gx)
  hbt = HilbertTransformFilter()
  fx = zerofloat(n1,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      hbt.apply(n1,gx[i3][i2],fx[i3][i2])
  plot3(fx)
  fs = FaultRescanner(sigmaPhi,sigmaTheta)
  ss = fs.scan(n1,n2,n3,minPhi,maxPhi,minTheta,maxTheta,p2,p3,fx,gx)
  '''
  ss = sub(1,ss)
  plot3(gx,ss,cmin=min(ss),cmax=max(ss),cmap=jetRamp(1.0),
    clab="sm",png="sm")
  '''
  plot3(gx,ss[0],cmin=min(ss[0]),cmax=max(ss[0]),cmap=jetRamp(1.0),
    clab="fl",png="sm")
  plot3(gx,ss[1],cmin=min(ss[1]),cmax=max(ss[1]),cmap=jetRamp(1.0),
    clab="sm",png="sm")
  plot3(gx,ss[2],cmin=min(ss[2]),cmax=max(ss[2]),cmap=jetRamp(1.0),
    clab="cm",png="sm")
  plot3(gx,ss[3],cmin=min(ss[3]),cmax=max(ss[3]),cmap=jetRamp(1.0),
    clab="jm",png="sm")
def coherence(sigma, et, fx):
  lsf = LocalSemblanceFilter(sigma, 2*sigma)
  return lsf.semblance(LocalSemblanceFilter.Direction3.VW, et, fx)
def goTest():
  gx = readImage(gxfile)
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  '''
  plot3(gx,fp,cmin=0,cmax=360,cmap=hueFill(1.0),
      clab="Fault strike (degrees)",cint=45,png="fp")
  plot3(gx,convertDips(ft),cmin=15,cmax=55,cmap=jetFill(1.0),
      clab="Fault dip (degrees)",png="ft")
  '''
  flt = readImage(fltfile)
  fs = FaultSkinner()
  fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  fs.setMaxDeltaStrike(10)
  fs.setMaxPlanarDistance(0.2)
  fs.setMinSkinSize(minSkinSize)
  fcs = fs.findCells([fl,fp,ft])
  '''
  plot3(gx,flt,cmin=min(flt),cmax=max(flt),cmap=jetRamp(1.0),
    clab="flt",png="sm")
  '''
  sk = readSkins(fskbase)
  plot3(gx,skins=sk,png="skinOld")
  #fcs = FaultSkin.getCells(sk)
  plot3(gx,cells=fcs,png="cells")
  fls = zerofloat(n1,n2,n3)
  for ic in range(0,len(fcs),5):
    cell = fcs[ic]
    ks = cell.getI()
    ms = cell.getIm()
    ps = cell.getIp()
    fls[ks[2]][ks[1]][ks[0]] = cell.getFl()
    fls[ms[2]][ms[1]][ms[0]] = cell.getFl()
    fls[ps[2]][ps[1]][ps[0]] = cell.getFl()
  plot3(gx,fls,cmin=min(fls),cmax=max(fls),cmap=jetRamp(1.0),
    clab="Fault likelihood",png="sm")

  tv3 = TensorVoting3()
  tv3.setSigma(20)
  tv3.setWindow(20,15,15)
  #cells = tv3.getFaultCells(n1,n2,n3,fcs)
  cells = tv3.randCells(len(fcs)/5,988,n1,n2,n3,fcs)
  plot3(gx,cells=cells,png="cells")
  sm,cm,u1,u2,u3 = tv3.applyVote(n1,n2,n3,cells)
  cells = tv3.findCells(0.1,sm,u1,u2,u3)
  plot3(gx,cells=cells,png="cells")
  skins = fs.findSkins(cells)
  for skin in skins:
    skin.smoothCellNormals(4)
  plot3(gx,cm,cmin=0.2,cmax=1.0,skins=skins,cmap=jetRamp(1.0),
        clab="Faults",png="skins")
  plot3(gx,sm,cmin=0.0,cmax=1.0,cmap=jetRamp(1.0),
    clab="Surfaceness",png="sm")
  plot3(gx,cm,cmin=0.0,cmax=1.0,cmap=jetRamp(1.0),
    clab="Junction",png="sm")
  '''
  plot3(gx,fm,cmin=0.0,cmax=1.0,cmap=jetRamp(1.0),
    clab="fm",png="sm")
  '''


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
    fl = readImage(flfile)
    fp = readImage(fpfile)
    ft = readImage(ftfile)
  plot3(gx,fl,cmin=min(fl),cmax=max(fl),cmap=jetRamp(1.0),
      clab="Fault likelihood",png="fl")
  plot3(gx,fp,cmin=0,cmax=360,cmap=hueFill(1.0),
      clab="Fault strike (degrees)",cint=45,png="fp")
  ft = convertDips(ft)
  plot3(gx,ft,cmin=min(ft),cmax=max(ft),cmap=jetFill(1.0),
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
  plot3(gx,clab="Amplitude",png="gx")
  plot3(gx,fl,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
        clab="Fault likelihood",png="fl")
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
  fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  fs.setMaxDeltaStrike(15)
  fs.setMaxPlanarDistance(0.5)
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

def plot3(f,g=None,isosurface=False,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          xyz=None,cells=None,skins=None,smax=0.0,slices=None,
          links=False,curve=False,trace=False,htgs=None,png=None):
  n3=len(f)
  n2=len(f[0])
  n1=len(f[0][0])
  s1,s2,s3=Sampling(n1),Sampling(n2),Sampling(n3)
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
  if isosurface:
    mc = MarchingCubes(s1,s2,s3,g)
    mc.setNormals(True)
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
    tg.setStates(states)
    sf.world.addChild(tg)
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
  ipg.setSlices(85,5,56)
  ipg.setSlices(90,5,83)
  if cbar:
    sf.setSize(837,700)
  else:
    sf.setSize(700,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.45*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  ov.setWorldSphere(BoundingSphere(0.5*n1-10,0.5*n2,0.5*n3-6,radius))
  ov.setAzimuthAndElevation(-70.0,25.0)
  ov.setTranslate(Vector3(0.0241,0.0517,0.0103))
  # for subset plots
  #ov.setWorldSphere(BoundingSphere(0.5*n1,0.5*n2,0.5*n3,radius))
  #ov.setAzimuthAndElevation(-40.0,25.0)
  #ov.setTranslate(Vector3(0.0241,-0.0400,0.0103))
  ov.setScale(1.2)
  #ov.setScale(1.3) #use only for subset plots
  
def plot3X(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          xyz=None,cells=None,skins=None,smax=0.0,slices=None,
          links=False,curve=False,trace=False,png=None):
  n1 = len(f[0][0])
  n2 = len(f[0])
  n3 = len(f)
  s1,s2,s3=Sampling(n1),Sampling(n2),Sampling(n3)
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
      ipg.setClips(-2.0,2.0)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-2.0,2.0)
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
    cmap = ColorMap(0.0,0.8,ColorMap.JET)
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
        cmap = ColorMap(0.0,0.8,ColorMap.JET)
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
    #k1,k2,k3 = (370,105,34) # most plots use these
    #k1,k2,k3 = (370,0,0) # most plots use these
    #k1,k2,k3 = (200,208,228) # most plots use these
    #k1,k2,k3 = (150,194,212) # most plots use these
    k1,k2,k3 = (135,194,216) # most plots use these
    k1,k2,k3 = (144,195,216) # most plots use these
  ipg.setSlices(k1,k2,k3)
  if cbar:
    sf.setSize(985,700) # for sch data
    #sf.setSize(837,700) # for fake data
  else:
    sf.setSize(848,700) # for sch data
    #sf.setSize(1048,900) # for sch data
    #sf.setSize(700,700) # for fake data
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  ov.setEyeToScreenDistance(3018.87) # for consistency with brooks
  ov.setWorldSphere(BoundingSphere(0.5*n1+20,0.5*n2,0.5*n3,radius))
  #ov.setAzimuthAndElevation(25.0,20.0)
  ov.setAzimuthAndElevation(-35.0,45.0)
  #ov.setAzimuthAndElevation(150.0,15.0)
  #ov.setAzimuthAndElevation(160.0,65.0)
  ov.setScale(1.2)
  #ov.setTranslate(Vector3(-0.182,-0.238,-0.012))
  #ov.setTranslate(Vector3(-0.190,-0.168,-0.006))
  ov.setTranslate(Vector3(-0.190,0.200,0.356))
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(137,1,pngDir+png+"cbar.png")


#############################################################################
run(main)
