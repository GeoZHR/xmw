"""
Demonstrate 3D seismic image processing for faults and horizons
Author: Dave Hale, Colorado School of Mines
Version: 2014.07.17
"""

from utils import *
#setupForSubset("bpSub1")
setupForSubset("jake")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count
# Names and descriptions of image files used below.
gxfile  = "gx" # input image (maybe after bilateral filtering)
gwfile  = "gw" # input image (maybe after bilateral filtering)
hxfile  = "horizon"
gsxfile = "gsx" # image after lsf with fault likelihoods
epfile  = "ep" # eigenvalue-derived planarity
wpfile  = "wp" # eigenvalue-derived planarity
p2file  = "p2" # inline slopes
p3file  = "p3" # crossline slopes
p2kfile = "p2k" # inline slopes (known)
p3kfile = "p3k" # crossline slopes (known)
flfile  = "fl" # fault likelihood
fpfile  = "fp" # fault strike (phi)
ftfile  = "ft" # fault dip (theta)
flvfile  = "flv" # fault likelihood
fpvfile  = "fpv" # fault strike (phi)
ftvfile  = "ftv" # fault dip (theta)
fltfile = "flt" # fault likelihood thinned
fptfile = "fpt" # fault strike thinned
fttfile = "ftt" # fault dip thinned
fs1file = "fs1" # fault slip (1st component)
fs2file = "fs2" # fault slip (2nd component)
fs3file = "fs3" # fault slip (3rd component)
fskbase = "fsk" # fault skin (basename only)
fslbase = "fsl" # fault skin (basename only)
fskgood = "fsg" # fault skin (basename only)
fwsfile = "fws" # unfaulted image
sw1file = "sw1" # 1st component of unfaulting shifts
sw2file = "sw2" # 2nd component of unfaulting shifts
sw3file = "sw3" # 3rd component of unfaulting shifts
gufile = "gu" # flattened image
gtfile = "gt" # flattened image
ghfile = "gh" # flattened image
grfile = "gr" # flattened image
x1file = "x1" # horizon volume
u1file = "u1" # first component of normal
u2file = "u2" # second component of normal
u3file = "u3" # third component of normal
smfile = "sm"
cmfile = "cm"


# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
minPhi,maxPhi = 0,360
minTheta,maxTheta = 65,85
sigmaPhi,sigmaTheta = 8,30

# These parameters control the construction of fault skins.
# See the class FaultSkinner for more information.
lowerLikelihood = 0.1
upperLikelihood = 0.3
minSkinSize = 10000

# These parameters control the computation of fault dip slips.
# See the class FaultSlipper for more information.
minThrow = 0.0
maxThrow = 85.0

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
pngDir = None
#pngDir = "../../../png/beg/hongliu/"
#pngDir = "../../../png/beg/bp/sub1/"
plotOnly = True

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goSlopes()
  #goScan()
  #goSkin()
  #goThin()
  #goTv()
  #goSkinTv()
  #goReSkin()
  #goSmooth()
  #goSlip()
  #goUnfaultS()
  #goFlattenWeights()
  #goHorizonExtraction1()
  #goHorizonExtraction2()
  #goHorizonExtraction3()
  #goFlattenC()
  #goDisplay()
  #gu = readImage(gufile)
  #plot3(gu)
  #sk = readSkins(fskgood)
  #plot3(gx,skins=sk)

def goDisplay():
  gx  = readImage(gxfile)
  hz1 = readHorizon("hz1")
  hz2 = readHorizon("hz2")
  hz3 = readHorizon("hz3")
  plot3(gx,horizon=hz1)
  plot3(gx,horizon=hz2)
  plot3(gx,horizon=hz3)

def goSlopes():
  print "goSlopes ..."
  gx = readImage(gxfile)
  sigma1,sigma2,sigma3,pmax = 8.0,2.0,2.0,5.0
  p2,p3,ep = FaultScanner.slopes(sigma1,sigma2,sigma3,pmax,gx)
  writeImage(p2file,p2)
  writeImage(p3file,p3)
  writeImage(epfile,ep)
  '''
  print "p2  min =",min(p2)," max =",max(p2)
  print "p3  min =",min(p3)," max =",max(p3)
  print "ep min =",min(ep)," max =",max(ep)
  plot3(gx,p2, cmin=-1,cmax=1,cmap=bwrNotch(1.0),
        clab="Inline slope (sample/sample)",png="p2")
  plot3(gx,p3, cmin=-1,cmax=1,cmap=bwrNotch(1.0),
        clab="Crossline slope (sample/sample)",png="p3")
  plot3(gx,sub(1,ep),cmin=0,cmax=1,cmap=jetRamp(1.0),
        clab="Planarity")
  '''

def goFlattenWeights():
  hp = Helper()
  gx = readImage(gxfile)
  ep = readImage(epfile)
  sk = readSkins(fskgood)
  wp = pow(ep,10)
  hp.setWeights(sk,wp)
  writeImage(wpfile,wp)
  '''
  plot3(gx,wp,cmin=0,cmax=1,cmap=jetRamp(1.0),
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
    fl = readImage(flfile)
    fp = readImage(fpfile)
    ft = readImage(ftfile)
    '''
  plot3(gx,fl,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
      clab="Fault likelihood",png="fl")
  plot3(gx,fp,cmin=0,cmax=360,cmap=hueFill(1.0),
      clab="Fault strike (degrees)",cint=45,png="fp")
  plot3(gx,ft,cmin=65,cmax=85,cmap=jetFill(1.0),
      clab="Fault dip (degrees)",png="ft")
    '''

def goThin():
  print "goThin ..."
  gx = readImage(gxfile)
  flt = readImage(fltfile)
  ftt = readImage(fttfile)
  '''
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
  '''
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  plot3(gx,ftt,cmin=70,cmax=85,cmap=jetFillExceptMin(1.0),
        clab="Fault dip (degrees)",png="ftt")
  '''
  plot3(gx,fpt,cmin=0,cmax=360,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=45,png="fpt")
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
    #plot3(gx,cells=cells,png="cells")
  else:
    skins = readSkins(fskgood)
  plot3(gx)
  plot3(gx,skins=skins)
  '''
  for iskin,skin in enumerate(skins):
    plot3(gx,skins=[skin],links=True,)
  '''
def goTv():
  gx = readImage(gxfile)
  if not plotOnly:
    tv3 = TensorVoting3()
    sk = readSkins(fskbase)
    fcs = FaultSkin.getCells(sk)
    cells=[]
    for ic in range(0,len(fcs),5):
      cells.append(fcs[ic])
    tv3.setSigma(10)
    tv3.setVoteWindow(20,20,20)
    fsc = FaultScanner(4,20)
    sp = fsc.getPhiSampling(minPhi,maxPhi)
    st = fsc.getThetaSampling(minTheta,maxTheta)
    fl,cm,fp,ft = tv3.applyVote(n1,n2,n3,cells)
    writeImage(flvfile,fl)
    writeImage(fpvfile,fp)
    writeImage(ftvfile,ft)
  else: 
    fl = readImage(flvfile)
    fp = readImage(fpvfile)
    ft = readImage(ftvfile)

def goSkinTv():
  print "goSkinTv ..."
  gx = readImage(gxfile)
  if not plotOnly:
    fl = readImage(flvfile)
    fp = readImage(fpvfile)
    ft = readImage(ftvfile)
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
    removeAllSkinFiles(fskgood)
    writeSkins(fskgood,skins)
  else:
    skins = readSkins(fskgood)
  plot3(gx,skins=skins)

def goReSkinX():
  print "goReSkin ..."
  gx = readImage(gxfile)
  if not plotOnly:
    fl = readImage(flfile)
    sk = readSkins(fskbase)
    fsx = FaultSkinnerX()
    fsx.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
    fsx.setMinSkinSize(minSkinSize)
    skins = fsx.reskin(sk,fl)
    removeAllSkinFiles(fskgood)
    writeSkins(fskgood,skins)
  skins = readSkins(fskgood)
  '''
  for skin in skins:
    skin.smoothCellNormals(4)
  '''
  plot3(gx,skins=skins,png="skinsNew")
  #plot3(gx,skins=skins,links=True,png="skinsNewLinks")

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
  '''
  plot3(gx,flt,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
        clab="Fault likelihood",png="fli")
  plot3(gsx,png="gsx")
  '''

def goSlip():
  print "goSlip ..."
  gx = readImage(gxfile)
  if not plotOnly:
    gsx = readImage(gsxfile)
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    skins = readSkins(fskgood)
    fsl = FaultSlipper(gsx,p2,p3)
    fsl.setOffset(3.0) # the default is 2.0 samples
    fsl.setZeroSlope(False) # True only if we want to show the error
    fsl.computeDipSlips(skins,minThrow,maxThrow)
    print "  dip slips computed, now reskinning ..."
    print "  number of skins before =",len(skins),
    fsk = FaultSkinner() # as in goSkin
    fsk.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
    fsk.setMinSkinSize(minSkinSize)
    fsk.setMinMaxThrow(minThrow,maxThrow)
    #skins = fsk.reskin(skins)
    print ", after =",len(skins)
    removeAllSkinFiles(fslbase)
    writeSkins(fslbase,skins)
    '''
    smark = -999.999
    s1,s2,s3 = fsl.getDipSlips(skins,smark)
    s1,s2,s3 = fsl.interpolateDipSlips([s1,s2,s3],smark)
    gw = fsl.unfault([s1,s2,s3],gx)
    writeImage(gwfile,gw)
    writeImage(fs1file,s1)
    writeImage(fs2file,s2)
    writeImage(fs3file,s3)
    '''
  else:
    gw = readImage(gwfile)
    #s1 = readImage(fs1file)
    #skins = readSkins(fslbase)
  '''
  plot3(gx,skins=skins,smax=10.0,png="skinss1")
  plot3(gx,s1,cmin=-10,cmax=10.0,cmap=jetFillExceptMin(1.0),
        clab="Fault throw (samples)",png="gxs1")
  plot3(gx,s1,cmin=0.0,cmax=10.0,cmap=jetFill(0.3),
        clab="Vertical shift (samples)",png="gxs1i")
  plot3(gx,s2,cmin=-2.0,cmax=2.0,cmap=jetFill(0.3),
        clab="Inline shift (samples)",png="gxs2i")
  plot3(gx,s3,cmin=-1.0,cmax=1.0,cmap=jetFill(0.3),
        clab="Crossline shift (samples)",png="gxs3i")
  plot3(gx)
  plot3(gw,png="gw")
  '''

def goUnfaultS():
  gx = readImage(gxfile)
  if not plotOnly:
    fw = zerofloat(n1,n2,n3)
    lof = LocalOrientFilter(8.0,4.0,4.0)
    et = lof.applyForTensors(gx)
    et.setEigenvalues(0.001,1.0,1.0)

    wp = fillfloat(1.0,n1,n2,n3)
    skins = readSkins(fslbase)
    fsc = FaultSlipConstraints(skins)
    sp = fsc.screenPoints(wp)
    mul(sp[3][0],10,sp[3][0])

    uf = UnfaultS(8.0,8.0)
    uf.setIters(100)
    uf.setTensors(et)
    [t1,t2,t3] = uf.findShifts(sp,wp)
    #uf.convertShifts(40,[t1,t2,t3])
    uf.applyShifts([t1,t2,t3],gx,fw)
    writeImage(fwsfile,fw)
    writeImage(sw1file,t1)
    writeImage(sw2file,t2)
    writeImage(sw3file,t3)
  else :
    fw = readImage(fwsfile)
    gw = readImage(gwfile)
  plot3(gx,png="gxuf")
  plot3(fw,png="fwuf")
  plot3(gw,png="fwuf")
  '''
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
  '''

def goHorizonExtraction1():
  k11 = [ 142,126,128, 75,48, 48, 28, 26, 66, 36, 43,121,138, 133,
          134, 143, 91, 60, 73, 67,120,117,131,157, 159,137, 77,
          182,101, 176, 198, 196,102, 98, 84, 191, 66]
  k12 = [1403,985,858,388,72,  8,  8,  8,161,161,161,794,895,1029,
         1223,1475,332,332,428,428,428,428,723,969,1045,507,507,
          644,661,1047,1194,1446,408,399,370,1299,383]
  k13 = [  21, 34, 29, 29,29,139,228,651,195,254,882, 60, 60,  59,
           60,  60,275,378,950,634,530,351,142,142, 142,585,698,
          644,816, 215, 317, 308,480,471,143, 292,559]
  gx = readImage(gxfile)
  if not plotOnly:
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    wp = readImage(wpfile)
    lmt = n1-1
    se = SurfaceExtractorC()
    se.setWeights(0.0)
    se.setSmoothings(6.0,6.0)
    se.setCG(0.01,100)
    surf = se.surfaceInitialization(n2,n3,lmt,k11,k12,k13)
    se.surfaceUpdateFromSlopes(wp,p2,p3,k11,k12,k13,surf)
    writeImage("hz1",surf) 
  else:
    surf = readHorizon("hz1")
  '''
  plot3(gx)
  plot3(gx,horizon=surf)
  '''

def goHorizonExtraction2():
  k11 = [ 99, 74, 70, 76,106, 81,133,152, 156, 159,130, 87,
          83, 82,135, 98,160,120,121, 147, 176, 183,149]

  k12 = [  8,  8, 71, 71, 68, 68,582,777,1003,1382,266,253,
         236,260,303,322,483,393,285,1163,1177,1478,894]

  k13 = [132,208,382,901,133,263, 24, 24,  39,  39,178,264,
         567,929,212,342,312,415,  6,   6, 139, 188, 14]
  gx = readImage(gxfile)
  if not plotOnly:
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    wp = readImage(wpfile)
    lmt = n1-1
    se = SurfaceExtractorC()
    se.setWeights(0.0)
    se.setSmoothings(6.0,6.0)
    se.setCG(0.01,100)
    surf = se.surfaceInitialization(n2,n3,lmt,k11,k12,k13)
    se.surfaceUpdateFromSlopes(wp,p2,p3,k11,k12,k13,surf)
    writeImage("hz2",surf) 
  else:
    surf = readHorizon("hz2")
  '''
  plot3(gx)
  plot3(gx,horizon=surf)
  '''

def goHorizonExtraction3():
  k11 = [23, 28, 10, 29, 48, 75, 125, 133, 45, 15, 40, 42, 75,
         67, 80, 53,122, 66,108,184,143,137,200, 160, 41, 49,110,106,110,115]
  k12 = [61, 61, 61, 91,313,573,1041,1493,202,209,222,390,420,
        420,409,429,563,563,563,748,765,841,879,1024,357,475,498,466,498,502]
  k13 = [53,146,227,901,102,102,  97, 124,178,355,105,631,201,
        160,440,866,585,751,348,617,326,208,735, 260,470,676,535,393,488,504]
  gx = readImage(gxfile)
  if not plotOnly:
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    wp = readImage(wpfile)
    lmt = n1-1
    se = SurfaceExtractorC()
    se.setWeights(0.0)
    se.setSmoothings(6.0,6.0)
    se.setCG(0.01,100)
    surf = se.surfaceInitialization(n2,n3,lmt,k11,k12,k13)
    se.surfaceUpdateFromSlopes(wp,p2,p3,k11,k12,k13,surf)
    writeImage("hz3",surf) 
  else:
    surf = readHorizon("hz3")
  plot3(gx)
  plot3(gx,horizon=surf)

def goFlattenC():
  print "Flatten with control points..."
  gx = readImage(gxfile)
  if not plotOnly:
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    wp = readImage(wpfile)
    hz1 = readHorizon("hz1")
    hz2 = readHorizon("hz2")
    hz3 = readHorizon("hz3")
    sc = SetupConstraints()
    kk1 = sc.constraintsFromSurface(n1-2,hz1)
    kk2 = sc.constraintsFromSurface(n1-2,hz2)
    kk3 = sc.constraintsFromSurface(n1-2,hz3)
    k1 = [kk1[0],kk2[0],kk3[0]]
    k2 = [kk1[1],kk2[1],kk3[1]]
    k3 = [kk1[2],kk2[2],kk3[2]]
    k4 = [kk1[3],kk2[3],kk3[3]]
    fl = Flattener3C()
    fl.setIterations(0.01,200)
    fl.setSmoothings(6.0,6.0)
    fl.setWeight1(0.05)  
    fl.setScale(0.001)
    fm = fl.getMappingsFromSlopes(s1,s2,s3,p2,p3,wp,k4,k1,k2,k3)
    gu = fm.flatten(gx) # flattened image
    gt = fm.u1 # rgt volume
    gh = fm.x1 # horizon volume
    gr = fl.resampleRgt(s1,gt)
    writeImage(gufile,gu)
    writeImage(gtfile,gt)
    writeImage(ghfile,gh)
    writeImage(grfile,gr)
  else:
    gu = readImage(gufile)
    gt = readImage(gtfile)
    gh = readImage(ghfile)
    gr = readImage(grfile)
    '''
  plot3(gx,png="seismic")
  plot3(gu,png="flattened")
  plot3(gx,gt,cmin=min(gt)+20,cmax=max(gt)-10,cmap=jetRamp(1.0),
        clab="Relative geologic time",png="rgt")
  plot3(gx,p2, cmin=-1,cmax=1,cmap=jetRamp(1.0),
      clab="Inline slope (sample/sample)",png="p2")
  plot3(gx,p3, cmin=-1,cmax=1,cmap=jetRamp(1.0),
      clab="Crossline slope (sample/sample)",png="p3")
  #fl = Flattener3C()
  #gr = fl.resampleRgt(s1,gt)
  f1 = min(gr)+50
  d1 = 5*s1.getDelta()
  n1 = round((max(gr)-f1)/d1)-1
  st = Sampling(n1,d1,f1)
  hfr = HorizonExtraction(s1,s2,s3,None,gr)
  k2 = 320
  k3s = [11,61,91,151,181,251,291]
  for k3 in k3s:
    #hls  = hfr.horizonCurves(st,k2,k3)
    plot3X(gx,k2=k2,k3=k3,png="seismic"+str(k3)+"new")
    #plot3X(gx,k2=k2,k3=k3,curve=True,hs=hls,png="horizonLines"+str(k3)+"new")
    '''


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
          horizon=None,xyz=None,cells=None,skins=None,smax=0.0,slices=None,
          links=False,curve=False,trace=False,png=None):
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
      ipg.setClips(-0.5,0.5)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-0.5,0.5)
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
  if horizon:
    sd = SurfaceDisplay()
    ts = sd.horizonWithAmplitude([-0.5,0.5],horizon,f)
    tg = TriangleGroup(True,ts[0],ts[1])
    #tg = TriangleGroup(True,s3,s2,horizon)
    #tg.setColor(Color.CYAN)
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
  #ipg.setSlices(85,5,43)
  #ipg.setSlices(85,5,102)
  #ipg.setSlices(n1,0,n3) # use only for subset plots
  ipg.setSlices(n1,376,350)
  if cbar:
    sf.setSize(1037,700)
  else:
    sf.setSize(900,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  zscale = 0.5*max(n2*d2,n3*d3)/(n1*d1)
  #zscale = 1.5*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.5)
  #ov.setScale(2.5)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(0.0,-0.00,-0.05))
  ov.setAzimuthAndElevation(45.0,35.0)
  #ov.setAzimuthAndElevation(-55.0,35.0)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")


#############################################################################
run(main)
