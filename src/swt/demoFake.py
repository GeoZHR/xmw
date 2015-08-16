"""
Demonstrate 3D seismic image processing for faults and horizons
Author: Xinming Wu, Colorado School of Mines
Version: 2015.07.16
"""

from fakeutils import *
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

# Names and descriptions of image files used below.
gxfile  = "gx" # input image (maybe after bilateral filtering)
vxfile  = "vx" # input image (maybe after bilateral filtering)
dxfile  = "dx" # input image (maybe after bilateral filtering)
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
fslbase = "fsl" # fault skin (basename only)
fskgood = "fsg" # fault skin (basename only)
fwsfile = "fws" # unfaulted image
sw1file = "sw1" # 1st component of unfaulting shifts
sw2file = "sw2" # 2nd component of unfaulting shifts
sw3file = "sw3" # 3rd component of unfaulting shifts
gufile = "gu" # flattened image
x1file = "x1" # horizon volume
u1file = "u1" # relateive geologic time volume
ulfile = "ul"
uncfile = "unc"
fgfile = "fg"
rgtfile = "rgt"
vqxfile = "vqx"
dqxfile = "dqx"
vqwfile = "vqw"
dqwfile = "dqw"
vqufile = "vqu"
dqufile = "dqu"

logType = "vel"
logType = "den"

# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
minPhi,maxPhi = 0,360
minTheta,maxTheta = 65,85
sigmaPhi,sigmaTheta = 8,20

# These parameters control the construction of fault skins.
# See the class FaultSkinner for more information.
lowerLikelihood = 0.2
upperLikelihood = 0.5
minSkinSize = 2000

# These parameters control the computation of fault dip slips.
# See the class FaultSlipper for more information.
minThrow = 0.0
maxThrow = 20.0

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
pngDir = None
pngDir = "../../../png/swt/fake/"
plotOnly = True

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goFakeData()
  #goSlopes()
  #goScan()
  #goThin()
  #goSkin()
  #goReSkin()
  #goSmooth()
  #goSlip()
  #goUnfaultS()
  goUnfaultX()
  #goUncScan()
  #goFlatten()
  #goInterp()
  #goHorizonExtraction()
  #goTest()
  #test()
def test():
  fg = readImage(fgfile)
  ri = RgtInterp3()
  fk1 = zerofloat(n2,n3)
  fk2 = zerofloat(n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      fk1[i3][i2] = fg[i3][i2][110]
      fk2[i3][i2] = fg[i3][i2][111]
  writeImage("fake110",fk1)
  writeImage("fake111",fk2)
  plot3(fg)
def goTest():
  uncs = readUncs(uncfile)
  sc = SetupConstraints()
  uncs = add(uncs,1.0)
  cs = sc.uncConstraints(uncs)
  rgt = readImage(rgtfile)
  fw  = readImage(fwsfile)
  fl3 = Flattener3Unc()
  fus,x1s = fl3.rgtSampling(s1,cs,rgt)
  for fu in fus:
    print fu
  g = fl3.flatten(s1,x1s,fw)
  n1s = len(g[0][0])
  dn = (n1s-n1)/2
  gs = copy(n1,n2,n3,0,0,0,g)
  fu1 = round(min(fus))
  lu1 = round(max(fus))
  nu1 = len(fus)
  su1 = Sampling(nu1,1.0,fu1)
  gt = fl3.rgtFromHorizonVolume(n1,su1,x1s)
  mp = fl3.getMappingsFromRgt(s1,gt)
  fg = mp.flatten(fw)
  plot3(gs)
  plot3(fw)
  plot3(fg)
  plot3(fw,rgt,cmin=10.0,cmax=n1,cmap=jetFill(1.0),
        clab="Relative geologic time (samples)",png="rgt")
  plot3(fw,gt,cmin=10.0,cmax=n1,cmap=jetFill(1.0),
        clab="Relative geologic time (samples)",png="rgt")
  '''
  gx = readImage(gxfile)
  gw = readImage(fwsfile)
  gt = readImage(rgtfile)
  fg = readImage(fgfile)
  sw1 = readImage(sw1file)
  sw2 = readImage(sw2file)
  sw3 = readImage(sw3file)
  skins = readSkins(fskgood)
  logs = getLogs(logType)
  cp = ConvertPoints()
  ps = cp.setUnfaultCoord(logs,skins,sw1,sw2,sw3)
  ps = cp.setFlattenedCoord(x1s,ps)
  fw,w1,w2,w3=cp.getSamplesW(ps)
  fx,x1,x2,x3=cp.getSamplesX(ps)
  fu,u1,u2,u3=cp.getSamplesU(ps)
  samplesW = fw,w1,w2,w3
  samplesX = fx,x1,x2,x3
  samplesU = fu,u1,u2,u3
  if not plotOnly:
    ri = RgtInterp3(ps)
    ri.setScales(0.001,1.0)
    ri.setRgt(gt)
    uf  = UnfaultS(4.0,2.0)
    fqw = ri.gridXX(n1,s2,s3,x1s,g)
    fqx = zerofloat(n1,n2,n3)
    uf.applyShiftsX([sw1,sw2,sw3],fqw,fqx)
    if logType=="vel":
      writeImage(vqxfile,fqx)
      writeImage(vqwfile,fqw)
    if logType=="den":
      writeImage(dqxfile,fqx)
      writeImage(dqwfile,fqw)
  else:
    if logType=="vel":
      fqw = writeImage(vqwfile)
      fqx = writeImage(vqxfile)
    if logType=="den":
      fqw = writeImage(dqwfile)
      fqx = writeImage(dqxfile)
  if logType=="vel":
    fqx = div(fqx,1000)
    fqw = div(fqw,1000)
    vmin,vmax,vmap= 3.0,5.5,ColorMap.JET
    clab = "Velocity (km/s)"
    pngx = "vqx"
    pngw = "vqw"
  if logType=="den":
    vmin,vmax,vmap= 2.2,3.1,ColorMap.JET
    clab = "Density (g/cc)"
    pngx = "dqx"
    pngw = "dqw"

  plot3(gx,samples=samplesX)
  plot3(gw,samples=samplesW)
  plot3(gs,samples=samplesU)
  plot3(fqx,cmin=vmin,cmax=vmax,cmap=vmap,clab=clab,png=pngx)
  plot3(fqw,cmin=vmin,cmax=vmax,cmap=vmap,clab=clab,png=pngw)
  '''

def goFakeData():
  sequence = 'OA' # 1 episode of folding, followed by one episode of faulting
  #sequence = 'OOOOOAAAAA' # 5 episodes of folding, then 5 of faulting
  nplanar = 3 # number of planar faults
  conjugate = False # if True, two large planar faults will intersect
  conical = False # if True, may want to set nplanar to 0 (or not!)
  impedance = False # if True, data = impedance model
  wavelet = True # if False, no wavelet will be used
  lateralViriation = True
  noise = 0.5 # (rms noise)/(rms signal) ratio
  if not plotOnly:
    gx,p2,p3,vx,dx = FakeData.seismicAndSlopes3d2015A(sequence,
    nplanar,conjugate,conical,impedance,wavelet,lateralViriation,noise)
    writeImage(gxfile,gx)
    writeImage(vxfile,vx)
    writeImage(dxfile,dx)
    writeImage(p2kfile,p2)
    writeImage(p3kfile,p3)
  else:
    gx = readImage(gxfile)
    p2 = readImage(p2kfile)
    p3 = readImage(p3kfile)
    vx = readImage(vxfile)
    dx = readImage(dxfile)
  print "gx min =",min(gx)," max =",max(gx)
  print "p2 min =",min(p2)," max =",max(p2)
  print "p3 min =",min(p3)," max =",max(p3)
  dmin = 2.1
  dmax = 3.2
  vmin = 2.8
  vmax = 5.8
  dmap = ColorMap.JET
  vmap = ColorMap.JET
  print "dmin=",dmin, "dmax=",dmax
  gmin,gmax,gmap = -2.0,1.5,ColorMap.GRAY
  if impedance:
    gmin,gmax,gmap = 0.0,1.4,ColorMap.JET
  vx = div(vx,1000)
  plot3(gx,cmin=gmin,cmax=gmax,cmap=gmap,clab="Amplitude",png="gx")
  plot3(vx,cmin=vmin,cmax=vmax,cmap=vmap,clab="Velocity (km/s)",png="vx")
  plot3(dx,cmin=dmin,cmax=dmax,cmap=dmap,clab="Density (g/cc)",png="dx")

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
  plot3(gx,fl,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
      clab="Fault likelihood",png="fl")
  plot3(gx,fp,cmin=0,cmax=360,cmap=hueFill(1.0),
      clab="Fault strike (degrees)",cint=45,png="fp")
  plot3(gx,convertDips(ft),cmin=15,cmax=55,cmap=jetFill(1.0),
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
  #plot3(gx,skins=skins,links=True,png="skinsNew")
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
  fsl.setOffset(1.0) # the default is 2.0 samples
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
  plot3(gx,skins=skins,smax=15.0,slices=[85,5,60],png="skinss1")
  plot3(gx,s1,cmin=0,cmax=15.0,cmap=jetFillExceptMin(1.0),
        clab="Fault throw (samples)",png="gxs1")
  '''
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
  '''

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

    uf = UnfaultS(8.0,4.0)
    uf.setIters(100)
    uf.setTensors(et)
    mul(sp[3][0],10,sp[3][0])
    [t1,t2,t3] = uf.findShifts(sp,wp)
    #[t1,t2,t3] = uf.convertShifts(40,[t1,t2,t3])
    uf.applyShifts([t1,t2,t3],gx,fw)
    writeImage(fwsfile,fw)
    writeImage(sw1file,t1)
    writeImage(sw2file,t2)
    writeImage(sw3file,t3)
  else :
    t1 = readImage(sw1file)
    t2 = readImage(sw2file)
    t3 = readImage(sw3file)
    gx = readImage(gxfile)
    fw = readImage(fwsfile)
  plot3(gx,png="gxuf")
  plot3(fw,png="fwuf")
  skins = readSkins(fslbase)
  mark = -999.99
  s1 = fillfloat(mark,n1,n2,n3)
  FaultSkin.getThrow(mark,skins,s1)
  plot3(gx,s1,cmin=0.0,cmax=15.0,cmap=jetFillExceptMin(1.0),
        clab="Fault throw (samples)",png="gxs1")
  plot3(gx,t1,cmin=-10.0,cmax=10.0,cmap=jetFill(0.3),
        clab="Vertical shift (samples)",png="gxs1i")
  plot3(gx,t2,cmin=-3.0,cmax=3.0,cmap=jetFill(0.3),
        clab="Inline shift (samples)",png="gxs2i")
  plot3(gx,t3,cmin=-1.0,cmax=1.0,cmap=jetFill(0.3),
        clab="Crossline shift (samples)",png="gxs3i")

def goUnfaultX():
  gx = readImage(gxfile)
  skins = readSkins(fslbase)
  mark = -999.99
  ufx = UnfaultX()
  t1,t2,t3 = ufx.getDipSlips(n1,n2,n3,skins,mark)
  x1,x2,x3 = ufx.interpolateDipSlips([t1,t2,t3],mark)
  fw = ufx.unfault([x1,x2,x3],gx)
  plot3(gx)
  plot3(fw)
  plot3(gx,x1,cmin=-10.0,cmax=10.0,cmap=jetFill(0.3),
        clab="Vertical shift (samples)")#,png="gxs1i")
  plot3(gx,x2,cmin=-3.0,cmax=3.0,cmap=jetFill(0.3),
        clab="Inline shift (samples)")#,png="gxs2i")
  plot3(gx,x3,cmin=-1.0,cmax=1.0,cmap=jetFill(0.3),
        clab="Crossline shift (samples)")#,png="gxs3i")

def goUncScan():
  sig1s,sig2s=1.0,2.0
  fw = readImage(fwsfile)
  if not plotOnly:
    fws = smoothF(fw)
    ip = InsPhase()
    cs = like(fws)
    ip.applyForCosine(fws,cs)
    unc = UncSurfer()
    unc.setSampling(2,2)
    unc.setForLof(sig1s,sig2s)
    ul=unc.likelihood(cs)
    uli = unc.interp(n1,n2,n3,ul)
    writeImage(ulfile,uli)
  ul = readImage(ulfile)
  unc = UncSurfer()
  ult = like(ul)
  unc.thin(0.6,ul,ult)
  sfs = unc.surfer2(n2,n3,0.2,3000,ult)
  sfs = unc.extractUncs(sfs,fw)
  removeAllUncFiles(uncfile)
  writeUncs(uncfile,sfs)
  plot3(fw,uncs=sfs,png="uncs")
  plot3(fw,ul,cmin=0.3,cmax=1.0,cmap=jetRamp(1.0),
        clab="Unconformity likelihood",png="ul")
def goFlatten():
  gw = readImage(fwsfile)
  if not plotOnly:
    sigma1,sigma2,sigma3,pmax = 2.0,1.0,1.0,5.0
    ip = InsPhase()
    cs = like(gw)
    ip.applyForCosine(gw,cs)
    p2,p3,ep = FaultScanner.slopes(sigma1,sigma2,sigma3,pmax,cs)
    zm = ZeroMask(0.2,1.0,1.0,1.0,gw)
    zero = 0.00;
    tiny = 0.01;
    zm.setValue(zero,p2)#set inline slopes for samples above water bottom
    zm.setValue(zero,p3)#set crossline slopes for samples above water bottom
    zm.setValue(tiny,ep)#set planarities for samples above water bottom
    wp = pow(ep,6)
    uncs = readUncs(uncfile)
    sc = SetupConstraints()
    uncs = add(uncs,2.0)
    cs = sc.constraintsFromSurfaces(uncs)
    sfs = copy(uncs)
    for i3 in range(0,35):
      for i2 in range(115,n2):
        sfs[0][i3][i2] = -100
    sfs = sc.uncConstraints(sfs)
    rs = zerofloat(n1,n2,n3)
    fl3 = Flattener3Unc()
    sig1,sig2=10.0,10.0
    fl3.setSmoothings(sig1,sig2)
    fl3.setIterations(0.01,300);
    mp = fl3.getMappingsFromSlopes(s1,s2,s3,p2,p3,wp,None,sfs,rs)
    gt = mp.u1
    su1 = Sampling(n1,1,1)
    #gu  = mp.flatten(gw)
    gu = fl3.flatten(s1,su1,gt,gw)
    writeImage(fgfile,gu)
    writeImage(rgtfile,gt)
  gx  = readImage(gxfile)
  gt = readImage(rgtfile)
  plot3(gx)
  plot3(gw)
  plot3(gu,png="gu")
  plot3(gw,gt,cmin=10.0,cmax=n1+20,cmap=ColorMap.JET,
        clab="Relative geologic time (samples)",png="rgt")
  plot3(gw,gt,cmin=10.0,cmax=n1+20,cmap=jetFill(0.9),
        clab="Relative geologic time (samples)",png="gwt")

def goInterp():
  #p=getImageWithLogPoints(logType)
  #samplesW=getLogPointsW(logType)
  #samplesX=getLogPointsX(logType)
  su1 = Sampling(n1,1,1)
  gx = readImage(gxfile)
  gw = readImage(fwsfile)
  gt = readImage(rgtfile)
  fl = Flattener3Unc()
  gu = fl.flatten(s1,su1,gt,gw)
  #gu = readImage(gufile)
  sw1 = readImage(sw1file)
  sw2 = readImage(sw2file)
  sw3 = readImage(sw3file)
  skins = readSkins(fskgood)
  logs = getLogs(logType)
  cp = ConvertPoints()
  ps = cp.setUnfaultCoord(logs,skins,sw1,sw2,sw3)
  ps = cp.setFlattenedCoord(s1,s2,s3,gt,ps)
  fw,w1,w2,w3=cp.getSamplesW(ps)
  fx,x1,x2,x3=cp.getSamplesX(ps)
  fu,u1,u2,u3=cp.getSamplesU(ps)
  samplesW = fw,w1,w2,w3
  samplesX = fx,x1,x2,x3
  samplesU = fu,u1,u2,u3
  if not plotOnly:
    ri = RgtInterp3(ps)
    ri.setScales(0.001,1.0)
    ri.setRgt(gt)
    uf  = UnfaultS(4.0,2.0)
    fqu,fqw = ri.gridX(s1,s2,s3,gw)
    fqu = copy(n1,n2,n3,fqu)
    fqx = zerofloat(n1,n2,n3)
    uf.applyShiftsX([sw1,sw2,sw3],fqw,fqx)
    if logType=="vel":
      writeImage(vqxfile,fqx)
      writeImage(vqwfile,fqw)
      writeImage(vqufile,fqu)
      fqx = div(fqx,1000)
      fqw = div(fqw,1000)
      fqu = div(fqu,1000)
    if logType=="den":
      writeImage(dqxfile,fqx)
      writeImage(dqwfile,fqw)
      writeImage(dqufile,fqu)
  else:
    if logType=="vel":
      fqu = readImage(vqufile)
      fqw = readImage(vqwfile)
      fqx = readImage(vqxfile)
      fqx = div(fqx,1000)
      fqw = div(fqw,1000)
    if logType=="den":
      fqu = readImage(dqufile)
      fqw = readImage(dqwfile)
      fqx = readImage(dqxfile)
  if logType=="vel":
    fx = readImage(vxfile)
    vmin,vmax,vmap= 3.0,5.5,jetFill(0.5)
    clab = "Velocity (km/s)"
    pngx = "vqx"
    pngw = "vqw"
    pngu = "vqu"
    pngf = "vgx"
    pngxw = "gxvw"
    pngww = "gwvw"
    pnguw = "guvw"
    pngi = "vxi"
  if logType=="den":
    fx = readImage(dxfile)
    vmin,vmax,vmap= 2.1,3.2,jetFill(0.5)
    clab = "Density (g/cc)"
    pngx = "dqx"
    pngw = "dqw"
    pngu = "dqu"
    pngf = "dgx"
    pngxw = "gxdw"
    pngww = "gwdw"
    pnguw = "gudw"
    pngi = "dxi"
  plot3(gx,samples=samplesX,png=pngxw)
  plot3(gw,samples=samplesW,png=pngww)
  plot3(gu,samples=samplesU,png=pnguw)
  plot3(gu,fqu,cmin=vmin,cmax=vmax,cmap=vmap,clab=clab,samples=samplesU,png=pngu)
  plot3(gw,fqw,cmin=vmin,cmax=vmax,cmap=vmap,clab=clab,samples=samplesW,png=pngw)
  plot3(gx,fqx,cmin=vmin,cmax=vmax,cmap=vmap,clab=clab,samples=samplesX,png=pngx)
  plot3(gx,fqx,cmin=vmin,cmax=vmax,cmap=vmap,clab=clab,png=pngi)
  plot3(gx,fx, cmin=vmin,cmax=vmax,cmap=vmap,clab=clab,png=pngf)
#def getImageWithLogPoints(logType):
def getLogPointsW(logType):
  k1u = [109,109,109,109,109, 21, 21, 21, 21, 21, 21, 21]
  k2u = [128,100, 81, 59, 21,142, 70, 20,  8, 96, 50,140]
  k3u = [ 86, 79, 89, 94, 68,145,130,135,  8, 20, 30, 10]
  np = len(k1u)
  r1 = readImage(sw1file)
  r2 = readImage(sw2file)
  r3 = readImage(sw3file)
  f = zerofloat(n1,n2,n3)
  if logType == "vel":
    f = readImage(vxfile)
  if logType == "den":
    f = readImage(dxfile)
  si = SincInterpolator()
  si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT)
  fx,x1,x2,x3=[],[],[],[]
  #p = zerofloat(n1,n1,n2)
  for ip in range(np):
    i1 = round(k1u[ip])
    i2 = round(k2u[ip])
    i3 = round(k3u[ip])
    k2 = round(k2u[ip]+r2[i3][i2][i1])
    k3 = round(k3u[ip]+r3[i3][i2][i1])
    for k1 in range(n1):
      p1 = k1-r1[k3][k2][k1]
      p2 = k2-r2[k3][k2][k1]
      p3 = k3-r3[k3][k2][k1]
      fp = f[k3][k2][k1]
      x1.append(p1)
      x2.append(p2)
      x3.append(p3)
      fx.append(fp)
  return fx,x1,x2,x3

def getLogs(logType):
  k1u = [105,105,105,105,105, 21, 21, 21, 21, 21, 21, 21]
  k2u = [127,100, 81, 56, 20,142, 70, 20,  8, 96, 50,140]
  k3u = [ 90, 78, 90, 93, 65,145,130,135,  8, 20, 30, 10]
  np = len(k1u)
  r1 = readImage(sw1file)
  r2 = readImage(sw2file)
  r3 = readImage(sw3file)
  f = zerofloat(n1,n2,n3)
  if logType == "vel":
    f = readImage(vxfile)
  if logType == "den":
    f = readImage(dxfile)
  fx = zerofloat(n1,np,4)
  for ip in range(np):
    i1 = round(k1u[ip])
    i2 = round(k2u[ip])
    i3 = round(k3u[ip])
    k2 = round(k2u[ip]+r2[i3][i2][i1])
    k3 = round(k3u[ip]+r3[i3][i2][i1])
    for k1 in range(n1):
      fx[1][ip][k1]=k1
      fx[2][ip][k1]=k2
      fx[3][ip][k1]=k3
      fx[0][ip][k1]=f[k3][k2][k1]
  return fx

def getLogPointsX(logType):
  k1u = [109,109,109,109,109, 21, 21, 21, 21, 21, 21, 21]
  k2u = [128,100, 81, 59, 21,142, 70, 20,  8, 96, 50,140]
  k3u = [ 86, 79, 89, 93, 68,145,130,135,  8, 20, 30, 10]
  np = len(k1u)
  r1 = readImage(sw1file)
  r2 = readImage(sw2file)
  r3 = readImage(sw3file)
  f = zerofloat(n1,n2,n3)
  if logType == "vel":
    f = readImage(vxfile)
  if logType == "den":
    f = readImage(dxfile)
  si = SincInterpolator()
  si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT)
  fx,x1,x2,x3=[],[],[],[]
  #p = zerofloat(n1,n1,n2)
  for ip in range(np):
    i1 = round(k1u[ip])
    i2 = round(k2u[ip])
    i3 = round(k3u[ip])
    k2 = round(k2u[ip]+r2[i3][i2][i1])
    k3 = round(k3u[ip]+r3[i3][i2][i1])
    for k1 in range(n1):
      x1.append(k1)
      x2.append(k2)
      x3.append(k3)
      fx.append(f[k3][k2][k1])
  return fx,x1,x2,x3

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

def gain2(x,sigma):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(sigma)
  ref.apply1(g,g)
  y = like(x)
  div(x,sqrt(g),y)
  return y

def smoothF(x):
  fsigma = 4.0
  flstop = 0.9
  flt = fillfloat(0.0,n1,n2,n3)
  sigma1,sigma2,sigma3,pmax = 8.0,1.0,1.0,1.0
  p2,p3,ep = FaultScanner.slopes(sigma1,sigma2,sigma3,pmax,x)
  return FaultScanner.smooth(flstop,fsigma,p2,p3,flt,x)

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
          links=False,curve=False,trace=False,htgs=None,
          uncs=None,samples=None,png=None):
  n3 = len(f)
  n2 = len(f[0])
  n1 = len(f[0][0])
  s1,s2,s3=Sampling(n1),Sampling(n2),Sampling(n3)
  d1,d2,d3=s1.getDelta(),s2.getDelta(),s3.getDelta()
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
  ipg.setSlices(109,138,59)
  #ipg.setSlices(92,140,59)
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
    us = UncSurfer()
    ul=readImage(ulfile)
    #ul = div(exp(ul),exp(1.0))
    for unc in uncs:
      [xyz,rgb]=us.buildTrigs(n1,s3,s2,-0.1,unc,ul)
      #[xyz,rgb]=us.buildTrigs(n1,s3,s2,0.01,unc,ul)
      tg  = TriangleGroup(True,xyz,rgb)
      sg.addChild(tg)
    sf.world.addChild(sg)
  if samples:
    fx,x1,x2,x3 = samples
    if logType=="vel":
      fx = div(fx,1000)
      vmin,vmax,vmap= 3.0,5.5,ColorMap.JET
    if logType=="den":
      vmin,vmax,vmap= 2.1,3.2,ColorMap.JET
    pg = makePointGroup(fx,x1,x2,x3,vmin,vmax,None)
    sf.world.addChild(pg)
  if cbar:
    sf.setSize(887,700)
  else:
    sf.setSize(750,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.48*sqrt(n1*n1+n2*n2+n3*n3)
  zscale = 0.80*max(n2*d2,n3*d3)/(n1*d1)
  ov = sf.getOrbitView()
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setWorldSphere(BoundingSphere(0.5*n1,0.4*n2,0.4*n3,radius))
  ov.setAzimuthAndElevation(120.0,25.0)
  ov.setTranslate(Vector3(0.02,0.16,-0.27))
  ov.setScale(1.25)
  # for subset plots
  #ov.setWorldSphere(BoundingSphere(0.5*n1,0.5*n2,0.5*n3,radius))
  #ov.setAzimuthAndElevation(-40.0,25.0)
  #ov.setTranslate(Vector3(0.0241,-0.0400,0.0103))
  #ov.setScale(1.3) #use only for subset plots
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")


#############################################################################
run(main)
