"""
Demo of surface reconstruction from fault cells/oriented points
Author: Xinming Wu and Dave Hale, Colorado School of Mines
Version: 2015.02.09
"""

from fakeutils import *
#setupForSubset("fake")
#setupForSubset("tp")
setupForSubset("pnz")
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
gffile = "gf"
gsfile = "gs"
ghfile = "gh"
grfile = "gr"
dwfile = "dw"
grifile = "gri"

# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
minPhi,maxPhi = 0,145
minTheta,maxTheta = 65,85
sigmaPhi,sigmaTheta = 8,30

# These parameters control the construction of fault skins.
# See the class FaultSkinner for more information.
lowerLikelihood = 0.5
upperLikelihood = 0.7
minSkinSize = 3000

# These parameters control the computation of fault dip slips.
# See the class FaultSlipper for more information.
minThrow = 0.01
maxThrow = 15.0

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
pngDir = None
plotOnly = False
pngDir = "../../../png/flc/tp/"
pngDir = "../../../png/flc/fake/"

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goFakeData()
  #goFlatten2d()
  #goFlatten3d()
  #goSlopes()
  #goScan()
  #goThin()
  #goSkin()
  goReferImage()
  #goRefine3d()
  '''
  gx = readImage(gxfile)
  plot3p(gx)
  gx = readImage("pnz00")
  n1,n2,n3=300,450,450
  j1,j2,j3=120,0 ,400
  gxs = copy(n1,n2,n3,j1,j2,j3,gx)
  gxs = div(gxs,50000)
  writeImage("gx",gxs)
  #plot3(gx)
  '''

def goFlatten2d():
  gx = readImage(gxfile)
  gx = gain(gx)
  g3 = gx[50]
  p2 = zerofloat(n1,n2)
  ep = zerofloat(n1,n2)
  lsf = LocalSlopeFinder(4.0,2.0,5.0)
  lsf.findSlopes(g3,p2,ep)
  ep = pow(ep,8)
  fl = Flattener2()
  fl.setWeight1(0.02)
  fl.setIterations(0.001,400)
  fl.setSmoothings(6.0,6.0)
  fm = fl.getMappingsFromSlopes(s1,s2,p2,ep)
  gf = fm.flatten(g3)
  flr = FlattenerR()
  gs = flr.getReferImageX(gf)
  smin,smax = -5.0,5.0
  r1min,r1max = -0.5,0.5
  r2min,r2max = -0.5,0.5
  dwk = DynamicWarpingK(8,smin,smax,s1,s2)
  dwk.setStrainLimits(r1min,r1max,r2min,r2max)
  dwk.setSmoothness(4,2)
  ss = dwk.findShifts(s1,gs,s1,gf)
  gh = dwk.applyShifts(s1,gf,ss)
  plot2(s1,s2,g3,clab="Amplitude",cmin=-2,cmax=2,png="gx")
  plot2(s1,s2,gs,clab="Amplitude",cmin=-2,cmax=2,png="gs")
  plot2(s1,s2,gf,clab="Amplitude",cmin=-2,cmax=2,png="gf")
  plot2(s1,s2,gh,clab="Amplitude",cmin=-2,cmax=2,png="gh")
  plot2(s1,s2,ss,cmap=ColorMap.JET,clab="Shifts (samples)",cmin=-3,cmax=3,png="ss")

def goFlatten3d():
  gx = readImage(gxfile)
  p2 = zerofloat(n1,n2,n3)
  p3 = zerofloat(n1,n2,n3)
  ep = zerofloat(n1,n2,n3)
  lsf = LocalSlopeFinder(4.0,2.0,2.0,5.0)
  lsf.findSlopes(gx,p2,p3,ep)
  ep = pow(ep,8)
  fl = Flattener3()
  fl.setWeight1(0.02)
  fl.setIterations(0.001,200)
  fl.setSmoothings(6.0,6.0)
  fm = fl.getMappingsFromSlopes(s1,s2,s3,p2,p3,ep)
  gf = fm.flatten(gx)
  writeImage(gffile,gf)
  plot3p(gx)
  plot3p(gf)
def goReferImage():
  gf = readImage(gffile)
  flr = FlattenerR()
  dr2,dr3 = 10,10
  fr2,fr3 = 10,10
  nr2 = (n2-fr2-10)/dr2
  nr3 = (n3-fr3-10)/dr3
  sr2 = Sampling(nr2,dr2,fr2)
  sr3 = Sampling(nr3,dr3,fr3)
  gr = flr.resample(s2,s3,sr2,sr3,gf)
  g2 = flr.imageToTraces(gr)
  sx2 = Sampling(len(g2)) 
  wlw = WellLogWarpingD()
  wlw.setMaxShift(10)
  wlw.setPowError(0.05)
  g2s = wlw.findShifts(g2)
  g2d = wlw.applyShiftsX(g2,g2s)
  '''
  smax = 10
  r1min,r1max = -0.2,0.2
  r2min,r2max = -0.2,0.2
  r3min,r3max = -0.2,0.2
  g3d = flr.tracesToImage(nr2,nr3,g2d)
  g3f = flr.flattenImage(smax,r1min,r1max,r2min,r2max,r3min,r3max,g3d)
  g2f = flr.imageToTraces(g3f)
  '''
  g2f = flr.flattenTraces(10,-0.2,0.2,g2d)
  plot2(s1,sx2,g2,clab="Amplitude",cmin=-2,cmax=2,png="g2")
  plot2(s1,sx2,g2d,clab="Amplitude",cmin=-2,cmax=2,png="g2d")
  plot2(s1,sx2,g2f,clab="Amplitude",cmin=-2,cmax=2,png="g2f")
  g3f = flr.tracesToImage(nr2,nr3,g2f)
  g3r = flr.sincInterp(sr2,sr3,s2,s3,g3f)
  plot3p(g3r,clab="Amplitude",png="g3r")
  writeImage(grifile,g3r)
def goRefine3d():
  gf = readImage(gffile)
  if not plotOnly:
    flr = FlattenerR()
    gr = readImage(grifile)
    smin,smax = -10.0,10.0
    dwk = DynamicWarpingK(4,smin,smax,s1,s2,s3)
    dwk.setStrainLimits(-0.2,0.2,-0.2,0.2,-0.2,0.2)
    dwk.setSmoothness(4,2,2)
    dw = dwk.findShifts(s1,gr,s1,gf)
    gh = dwk.applyShifts(s1,gf,dw)
    writeImage(dwfile,dw)
    writeImage(ghfile,gh)
    writeImage(grfile,gr)
  else:
    gh = readImage(ghfile)
    gr = readImage(grfile)
    dw = readImage(dwfile)
  plot3p(gf,clab="Amplitude",png="gf1")
  plot3p(gh,clab="Amplitude",png="gh1")
  plot3p(gr,clab="Amplitude",png="gr1")
  plot3p(gf,dw,cmin=-8,cmax=6,clab="Shifts (samples)",cmap=jetFill(1.0),png="dw1")

def goRefine3dV():
  gf = readImage(gffile)
  if not plotOnly:
    sk = readSkins(fskbase)
    flr = FlattenerR()
    #gr = flr.getReferImage(gf)
    k2,k3=98,175
    #gr = flr.getReferImageX(k2,k3,gf)
    gr = readImage(grifile)
    smin,smax = -10.0,10.0
    r1mins = fillfloat(-0.2,n1,n2,n3)
    r1maxs = fillfloat( 0.2,n1,n2,n3)
    r2mins = fillfloat(-0.2,n1,n2,n3)
    r2maxs = fillfloat( 0.2,n1,n2,n3)
    r3mins = fillfloat(-0.2,n1,n2,n3)
    r3maxs = fillfloat( 0.2,n1,n2,n3)
    FaultSkin.setValuesOnFaults(-10.0,sk,r1mins)
    FaultSkin.setValuesOnFaults( 10.0,sk,r1maxs)
    FaultSkin.setValuesOnFaults(-10.0,sk,r2mins)
    FaultSkin.setValuesOnFaults( 10.0,sk,r2maxs)
    FaultSkin.setValuesOnFaults(-10.0,sk,r3mins)
    FaultSkin.setValuesOnFaults( 10.0,sk,r3maxs)
    dwc = DynamicWarpingC(8,smin,smax,s1,s2,s3)
    dwc.setStrains(r1mins,r2mins,r3mins,r1maxs,r2maxs,r3maxs)
    dwc.setSmoothness(4,2,2)
    dw = dwc.findShifts(s1,gr,s1,gf)
    gh = dwc.applyShifts(s1,gf,dw)
    #dwk = DynamicWarpingK(8,smin,smax,s1,s2,s3)
    #dwk.setStrainLimits(-0.2,0.2,-0.2,0.2,-0.2,0.2)
    #dwk.setSmoothness(4,2,2)
    #dw = dwk.findShifts(s1,gr,s1,gf)
    #gh = dwk.applyShifts(s1,gf,dw)
    #writeImage(dwfile,dw)
    #writeImage(ghfile,gh)
    #writeImage(grfile,gr)
  else:
    gh = readImage(ghfile)
    gr = readImage(grfile)
    dw = readImage(dwfile)
  plot3(gf,clab="Amplitude",png="gf1")
  plot3(gh,clab="Amplitude",png="gh1")
  plot3(gr,clab="Amplitude",png="gr1")
  plot3(gf,dw,cmin=-8,cmax=6,clab="Shifts (samples)",cmap=jetFill(1.0),png="dw1")
'''
def goFlatten3d():
  gx = readImage(gxfile)
  p2 = zerofloat(n1,n2,n3)
  p3 = zerofloat(n1,n2,n3)
  ep = zerofloat(n1,n2,n3)
  lsf = LocalSlopeFinder(4.0,2.0,2.0,5.0)
  lsf.findSlopes(gx,p2,p3,ep)
  ep = pow(ep,8)
  fl = Flattener3()
  fl.setWeight1(0.02)
  fl.setIterations(0.001,400)
  fl.setSmoothings(6.0,6.0)
  fm = fl.getMappingsFromSlopes(s1,s2,s3,p2,p3,ep)
  gf = fm.flatten(gx)
  flr = FlattenerR()
  gs = flr.getReferImageM(gf)
  smin,smax = -5.0,5.0
  r1min,r1max = -0.5,0.5
  r2min,r2max = -0.5,0.5
  r3min,r3max = -0.5,0.5
  dwr = DynamicWarpingK(8,smin,smax,s1,s2,s3)
  dwr.setStrainLimits(r1min,r1max,r2min,r2max,r3min,r3max)
  dwr.setSmoothness(2,5,5)
  ss = dwr.findShifts(s1,gs,s1,gf)
  gh = dwr.applyShifts(s1,gf,ss)
  #plot3(gx)
  #plot3(gs)
  plot3(gf)
  plot3(gh)
'''
def goSPS():
  print "screened Poisson surface method ..."
  gx = readImage(gxfile)
  p2 = zerofloat(n1,n2,n3)
  p3 = zerofloat(n1,n2,n3)
  ep = zerofloat(n1,n2,n3)
  lsf = LocalSlopeFinder(4.0,1.0,1.0,5.0)
  lsf.findSlopes(gx,p2,p3,ep)
  ep = pow(ep,8)
  fl = Flattener3()
  fl.setWeight1(0.02)
  fl.setIterations(0.001,400)
  fl.setSmoothings(6.0,6.0)
  fm = fl.getMappingsFromSlopes(s1,s2,s3,p2,p3,ep)
  ut = fm.u1
  gf = fm.flatten(gx)
  fx = fm.unflatten(gf)
  plot3(gx,png="gx")
  plot3(gf,png="gf")
  plot3(fx,png="fx")
  plot3(gx,ut,cmin=min(ut),cmax=max(ut),cmap=jetFill(1.0),
    clab="ScreenPoissonSurface",png="pss")

def goFakeData():
  #sequence = 'A' # 1 episode of faulting only
  sequence = 'OA' # 1 episode of folding, followed by one episode of faulting
  #sequence = 'OOOOOAAAAA' # 5 episodes of folding, then 5 of faulting
  #sequence = 'OAOAOAOAOA' # 5 interleaved episodes of folding and faulting
  nplanar = 0 # number of planar faults
  conjugate = False # if True, two large planar faults will intersect
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
  gx = readImage(gffile)
  sigma1,sigma2,sigma3,pmax = 16.0,1.0,1.0,5.0
  p2,p3,ep = FaultScanner.slopes(sigma1,sigma2,sigma3,pmax,gx)
  zm = ZeroMask(0.1,1,1,1,gx)
  zero,tiny=0.0,0.01
  zm.setValue(zero,p2)
  zm.setValue(zero,p3)
  zm.setValue(tiny,ep)
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
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  gx = readImage(gffile)
  gx = FaultScanner.taper(10,0,0,gx)
  fs = FaultScanner(sigmaPhi,sigmaTheta)
  fl,fp,ft = fs.scan(minPhi,maxPhi,minTheta,maxTheta,p2,p3,gx)
  zm = ZeroMask(0.1,1,1,1,gx)
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
  plot3(gx,clab="Amplitude")
  plot3(gx,fl,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
        clab="Fault likelihood",png="fl")
  plot3(gx,fp,cmin=0,cmax=145,cmap=jetFill(1.0),
        clab="Fault strike (degrees)",cint=45,png="fp")
  plot3(gx,convertDips(ft),cmin=25,cmax=65,cmap=jetFill(1.0),
        clab="Fault dip (degrees)",png="ft")

def goThin():
  print "goThin ..."
  gx = readImage(gffile)
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  flt,fpt,ftt = FaultScanner.thin([fl,fp,ft])
  writeImage(fltfile,flt)
  writeImage(fptfile,fpt)
  writeImage(fttfile,ftt)
  plot3(gx,clab="Amplitude")
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  plot3(gx,fpt,cmin=0,cmax=145,cmap=jetFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=45,png="fpt")
  plot3(gx,convertDips(ftt),cmin=25,cmax=65,cmap=jetFillExceptMin(1.0),
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
  gx = readImage(gffile)
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  fs = FaultSkinner()
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
  plot3(gx,cells=cells,png="cells")
  plot3(gx,skins=skins,png="skins")
  for iskin,skin in enumerate(skins):
    plot3(gx,skins=[skin],links=True,png="skin"+str(iskin))

def like(x):
  n2 = len(x)
  n1 = len(x[0])
  return zerofloat(n1,n2)

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

def plot1s(s1,ss,ys,rs=None,vmin=None,vmax=None,color=Color.RED,
  hlabel="Seismic traces",vlabel="time (ms)",png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sf = 10.0
  yf = sf
  sp.setVLimits(0,n1)
  if vmin and vmax:
    sp.setVLimits(vmin,vmax)
  sp.setHLimits(5.0,115)
  for il,y in enumerate(ys):
    ya = sum(y)/len(y)
    y = sub(y,ya)
    y = add(y,yf)
    pv = sp.addPoints(ss[il],y)
    pv.setLineColor(color)
    yf = yf+sf
  rf = sf
  if rs:
    for il,r in enumerate(rs):
      ra = sum(r)/len(r)
      r = sub(r,ra)
      r = add(r,rf)
      pv = sp.addPoints(s1,r)
      pv.setLineColor(Color.BLACK)
      rf = rf+sf
  sp.setSize(600,800)
  sp.setHLabel(hlabel)
  sp.setVLabel(vlabel)
  if png and pngDir:
    sp.paintToPng(300,7.0,pngDir+png+".png")

def plot2(s1,s2,x,cmap=ColorMap.GRAY,clab=None,cmin=0,cmax=0,
         title=None,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  if title:
    sp.setTitle(title)
  sp.addColorBar(clab)
  sp.setSize(680,600)
  sp.plotPanel.setColorBarWidthMinimum(80)
  pv = sp.addPixels(s1,s2,x)
  pv.setColorModel(cmap)
  #pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if pngDir and png:
    sp.paintToPng(300,3.333,pngDir+png+".png")


def plot3(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
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

def plot3p(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
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
    sf.setSize(937,600)
  else:
    sf.setSize(800,600)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  ov.setWorldSphere(BoundingSphere(0.5*n1,0.5*n2,0.5*n3,radius))
  ov.setAzimuthAndElevation(55.0,25.0)
  ov.setTranslate(Vector3(0.7,0.33,0.7))
  ov.setScale(1.1)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(137,1,pngDir+png+"cbar.png")

#############################################################################
run(main)
