"""
Demonstrate simultaneous multiple-well ties
Author: Xinming Wu, Colorado School of Mines
Version: 2015.06.08
"""

from utils import *
setupForSubset("subt1")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count
d1,d2,d3 = s1.delta,s2.delta,s3.delta

# Names and descriptions of image files used below.
gxfile = "gx" # input seismic image 
gtfile = "gt" # RGT volume
gtcfile = "gtc" # corrected RGT volume
gsfile = "gt" # flatten shifts
ghfile = "gh" # horizon volume
grfile = "gr" # horizon volume
dwfile = "dw" # horizon volume
gffile = "gf" # flattened image 
gufile = "gu" # flattened image 
txfile = "tx"
tufile = "tu"
tsfile = "ts"
wshifts = "wellShifts"
wshiftd = "wellShiftsDimensions"
wsample = "wellSampling"
p2file = "p2"
p3file = "p3"
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
fpifile  = "fpi" # nearest interp
fqifile  = "fqi" # blended interp

# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
minPhi,maxPhi = 0,90
minTheta,maxTheta = 70,89
sigmaPhi,sigmaTheta = 10,40

# These parameters control the construction of fault skins.
# See the class FaultSkinner for more information.
lowerLikelihood = 0.5
upperLikelihood = 0.7
minSkinSize = 20000

# These parameters control the computation of fault dip slips.
# See the class FaultSlipper for more information.
minThrow = 0.01
maxThrow = 15.0


# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
pngDir = None
pngDir = "../../../png/swt/print/"
plotOnly = False

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goDisplay()
  #goSynSeis()
  #goSynsFlatten()
  #goSeisFlatten()
  #goSynsSeisTie()
  goTimeUpdate()
  #goSlopes()
  #goScan()
  #goThin()
  #goSkin()
  #goImageFlatten()
  #goRefine3dV()
  #goTest()
def goTest():
  gx = readImage(gxfile)
  lgs = getLogs()
  nl = len(lgs)
  fx = zerofloat(n1,nl)
  for il, lg in enumerate(lgs):
    model = SynSeis.getModel(lg)
    i2 = s2.indexOfNearest(model.x2)
    i3 = s3.indexOfNearest(model.x3)
    fx[il] = gx[i3][i2]
  swt = SeismicWellTie()
  mds = swt.updateTimeDepthS(s1,s2,s3,lgs,gx)
  ndfx = zerodouble(3,nl)
  ndfu = zerodouble(3,nl)
  wxs = swt.getSyns(True,lgs,ndfx)
  wus = swt.getSyns(True,mds,ndfu)
  swx = []
  swu = []
  for il in range(nl):
    sxi = Sampling((int)(ndfx[il][0]),ndfx[il][1],ndfx[il][2])
    sui = Sampling((int)(ndfu[il][0]),ndfu[il][1],ndfu[il][2])
    swx.append(sxi)
    swu.append(sui)
  plot1s(s1,swx,wxs,rs=fx,vmin=0.1,vmax=1.15, 
         vlabel="Time (s)",png="synsSeisBS")
  plot1s(s1,swu,wus,rs=fx,vmin=0.1,vmax=1.15,
         vlabel="Time (s)",png="synsSeisAS")
  svw,sdw=swt.getSamples(s1,lgs)
  svs,sds=swt.getSamples(s1,mds)
  vlabel2 = "RGT"
  vlabel1 = "Time (s)"
  hlabel = "Synthetic seismic traces"
  plot3(gx,sps=svw,curve="vel",wmin=2.4,wmax=5.0,png="seisVelBS")
  plot3(gx,sps=svs,curve="vel",wmin=2.4,wmax=5.0,png="seisVelAS")
  goRgtInterp(svw,"wS")
  goRgtInterp(svs,"sS")

def goRgtInterp(sps,fname):
  gx = readImage(gxfile)
  gt = readImage(gtcfile)
  swt = SeismicWellTie()
  spc = swt.convertPoints(sps)
  ri = RgtInterp3(spc[0],spc[1],spc[2],spc[3])
  ri.setRgt(gt)
  ri.setScales(0.001,1.0)
  fti,fpi,fqi = ri.grid(s1,s2,s3)
  writeImage(fpifile+fname,fpi)
  #writeImage(fqifile+fname,fqi)
  plot3(gx,fpi,cmin=2.4,cmax=5.0,cmap=jetRamp(1.0),
        clab="Velocity (km/s)",png=fpifile+fname)
  plot3(gx,fpi,sps=sps,cmin=2.4,cmax=5.0,wmin=2.4,wmax=5.0,cmap=jetRamp(1.0),
        clab="Velocity (km/s)",png=fpifile+fname+"Wells")

def goDisplay():
  gx = readImage(gxfile)
  lgs = getLogs()
  swt = SeismicWellTie()
  sps = swt.getSamples(s1,lgs)
  plot3(gx,sps=sps[1],wmin=2.2,wmax=2.8,clab="Density (g/cc)",png="seisDen")
  plot3(gx,sps=sps[0],wmin=2.3,wmax=5.4,clab="Velocity (km/s)",png="seisVel")

def goTimeUpdate():
  gx = readImage(gxfile)
  lgs = getLogs()
  nl = len(lgs)
  swt = SeismicWellTie()
  ndfs = zerodouble(3,nl)
  ndfw = zerodouble(3,nl)
  fx,fs,fu=goSeisFlatten()
  swsu,wsu,mds= swt.updateTimeDepthM(5,s1,lgs,fx,fs,fu)
  wwx = swt.getSyns(True,lgs,ndfw)
  wsx = swt.getSyns(True,mds,ndfs)
  swsx = []
  swwx = []
  for il in range(nl):
    swi = Sampling((int)(ndfw[il][0]),ndfw[il][1],ndfw[il][2])
    ssi = Sampling((int)(ndfs[il][0]),ndfs[il][1],ndfs[il][2])
    swwx.append(swi)
    swsx.append(ssi)
  svw,sdw=swt.getSamples(s1,lgs)
  svs,sds=swt.getSamples(s1,mds)
  hlabel = "Log index"
  vlabel1 = "Time (s)"
  vlabel2 = "Relative geologic time (s)"
  plot1s(s1,swwx,wwx,rs=fx,vmin=0.1,vmax=1.15,vlabel=vlabel1,png="synsSeisMB")
  plot1s(s1,swsx,wsx,rs=fx,vmin=0.1,vmax=1.15,vlabel=vlabel1,png="synsSeisMA")
  plot1s(s1,swsu,wsu,rs=fu,vmin=0.1,vmax=1.15,vlabel=vlabel2,png="synsSeisF")
  plot3(gx,sps=svw,curve="vel",wmin=2.3,wmax=5.4,png="seisVelB")
  plot3(gx,sps=svs,curve="vel",wmin=2.3,wmax=5.4,png="seisVelA")
  goRgtInterp(svw,"w")
  goRgtInterp(svs,"s")


def goSynSeis():
  simple = True
  logs = getLogs()
  nl = len(logs)
  ndft = zerofloat(3,nl)
  ndfz = zerofloat(3,nl)
  sw = SeismicWellTie()
  sa = sw.computeSyns(simple,logs,ndft,ndfz)
  ss = []
  for il in range(nl):
    sa[il] = normalize(sa[il])
    sl = Sampling((int)(ndft[il][0]),ndft[il][1],ndft[il][2])
    ss.append(sl)
  gx = readImage(gxfile)
  fx = zerofloat(n1,nl)
  for il, log in enumerate(logs):
    model = SynSeis.getModel(log)
    i2 = s2.indexOfNearest(model.x2)
    i3 = s3.indexOfNearest(model.x3)
    fx[il] = gx[i3][i2]
  plot1s(s1,ss,sa,hlabel="Seismic traces",vlabel="time (s)")
  plot1s(s1,ss,sa,rs=fx,hlabel="Seismic traces",vlabel="time (s)")

def goSynsFlatten():
  smax = 40
  epow = 0.0625
  simple = True
  logs = getLogs()
  nl = len(logs)
  ndfx = zerodouble(3,nl)
  ndfu = zerodouble(3,nl)
  swt = SeismicWellTie()
  wx,ws,wu = swt.synsFlatten(simple,smax,epow,logs,ndfx,ndfu)
  ssx = []
  ssu = []
  for il in range(nl):
    sxi = Sampling((int)(ndfx[il][0]),ndfx[il][1],ndfx[il][2])
    sui = Sampling((int)(ndfu[il][0]),ndfu[il][1],ndfu[il][2])
    ssx.append(sxi)
    ssu.append(sui)
  hlabel = "Log index"
  vlabel1 = "Time (s)"
  vlabel2 = "Relative geologic time (s)"
  png1,png2="syns","synsF"
  plot1s(s1,ssx,wx,vmin=0.1,vmax=1.15,hlabel=hlabel,vlabel=vlabel1,png=png1)
  plot1s(s1,ssu,wu,vmin=0.1,vmax=1.15,hlabel=hlabel,vlabel=vlabel2,png=png2)
  return wu,ssu

def goSeisFlatten():
  logs = getLogs()
  nl = len(logs)
  gx = readImage(gxfile)
  fx = zerofloat(n1,nl)
  for il, log in enumerate(logs):
    model = SynSeis.getModel(log)
    i2 = s2.indexOfNearest(model.x2)
    i3 = s3.indexOfNearest(model.x3)
    fx[il] = gx[i3][i2]
  smax,epow = 40,0.125
  swt = SeismicWellTie()
  fs,fu = swt.seisFlatten(smax,epow,fx)
  sst = []
  for il in range(nl):
    sst.append(s1)
  writeImage(tufile,fu)
  writeImage(txfile,fx)
  hlabel = "Log index"
  vlabel1 = "Time (s)"
  vlabel2 = "Relative geologic time (s)"
  plot1s(s1,sst,fx,vmin=0.1,vmax=1.15,color=Color.BLACK,
        hlabel=hlabel,vlabel=vlabel1,png="seis")
  plot1s(s1,sst,fu,vmin=0.1,vmax=1.15,color=Color.BLACK,
        hlabel=hlabel,vlabel=vlabel2,png="seisF")
  return fx,fs,fu

def goSynsSeisTie():
  gx,gs,gu=goSeisFlatten() 
  wu,sw=goSynsFlatten()
  dmin = 20*d1
  smin,smax = -0.10,0.30
  rmin,rmax = -0.05,0.05
  swt = SeismicWellTie()
  ndf = zerofloat(3,len(wu))
  wgu = swt.synsToSeis(sw,s1,smin,smax,rmin,rmax,dmin,wu,gu,ndf)
  swg = []
  for il in range(len(wu)):
    si = Sampling((int)(ndf[il][0]),ndf[il][1],ndf[il][2])
    swg.append(si)
  vlabel = "Relative geologic time (s)"
  plot1s(s1,sw,wu,rs=gu,vmin=0.1,vmax=1.15,vlabel=vlabel,png="synSeisFB")
  plot1s(s1,swg,wgu,rs=gu,vmin=0.1,vmax=1.15,vlabel=vlabel,png="synSeisFA")

def goImageFlatten():
  gx = readImage(gxfile)
  if not plotOnly:
    p2 = zerofloat(n1,n2,n3)
    p3 = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    lsf = LocalSlopeFinder(4.0,2.0,2.0,5.0)
    lsf.findSlopes(gx,p2,p3,ep)
    zm = ZeroMask(0.1,1,1,1,gx)
    zero,tiny=0.0,0.01
    zm.setValue(zero,p2)
    zm.setValue(zero,p3)
    zm.setValue(tiny,ep)
    ep = pow(ep,8)
    p2 = mul(d1/d2,p2)
    p3 = mul(d1/d3,p3)
    fl = Flattener3()
    fl.setWeight1(0.04)
    fl.setIterations(0.001,200)
    fl.setSmoothings(6.0,6.0)
    fm = fl.getMappingsFromSlopes(s1,s2,s3,p2,p3,ep)
    u1 = fm.u1
    gs = fm.getShiftsR();
    gs = div(gs,d1)
    writeImage(gtfile,u1)
    writeImage(gsfile,gs)
    gf = fm.flatten(gx)
    writeImage(gffile,gf)
  else:
    gf = readImage(gffile)
  plot3(gx,png="seisImage")
  plot3(gf,png="seisImageFlat")
  plot3(gf,gs,cmin=min(gs),cmax=max(gs),clab="Shifts (s)",
        cmap=jetFill(0.3),png="dw1")

def goRefine3dV():
  gf = readImage(gffile)
  if not plotOnly:
    sk = readSkins(fskbase)
    flr = FlattenerR()
    k2,k3=57,37
    gr = flr.getReferImageX(k2,k3,gf)
    smin,smax = -15,15
    s1 = Sampling(n1)
    s2 = Sampling(n2)
    s3 = Sampling(n3)
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
    writeImage(dwfile,dw)
    writeImage(ghfile,gh)
  else:
    gh = readImage(ghfile)
    gr = readImage(grfile)
    dw = readImage(dwfile)
  plot3(gf,clab="Amplitude",png="gf1")
  plot3(gh,clab="Amplitude",png="gh1")
  plot3(gr,clab="Amplitude",png="gr1")
  plot3(gf,dw,cmin=min(dw),cmax=max(dw),clab="Shifts (samples)",
        cmap=jetFill(1.0),png="dw1")
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
  print "p2  min =",min(p2)," max =",max(p2)
  print "p3  min =",min(p3)," max =",max(p3)
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
  plot3F(gx,cells=cells,png="cells")
  plot3F(gx,skins=skins,png="skins")
  for iskin,skin in enumerate(skins):
    plot3(gx,skins=[skin],links=True,png="skin"+str(iskin))

def normalize(f):
  sigma = 100
  n = len(f)
  g = zerofloat(n)
  ref = RecursiveExponentialFilter(sigma)
  g = mul(f,f)
  ref.apply(g,g)
  g = div(f,sqrt(g))
  return g

def like(x):
  n2 = len(x)
  n1 = len(x[0])
  return zerofloat(n1,n2)

def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(10.0)
  ref.apply1(g,g)
  y = like(x)
  div(x,sqrt(g),y)
  return y

def plot1(s1,ys,hlabel="Seismic traces",vlabel="depth (km)",png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  for y in ys:
    pv = sp.addPoints(s1,y)
    pv.setLineColor(Color.BLACK)
  #sp.setVLimits(0.1,1.1)
  sp.setSize(800,800)
  sp.setHLabel(hlabel)
  sp.setVLabel(vlabel)
  if png and pngDir:
    sp.paintToPng(300,7.0,pngDir+png+".png")

def plot1X(s1,y,vmin=None,vmax=None, 
    color=Color.BLACK,hlabel="Stacked trace",vlabel="Time (s)",png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPoints(s1,y)
  pv.setLineColor(color)
  sp.setVLimits(vmin,vmax)
  sp.setHLimits(-3.5,3.5)
  sp.setSize(190,800)
  sp.setHLabel(hlabel)
  sp.setVLabel(vlabel)
  if png and pngDir:
    sp.paintToPng(300,7.0,pngDir+png+".png")

def plot1s(s1,ss,ys,rs=None,vmin=None,vmax=None,color=Color.RED,
  hlabel="Log index",vlabel="Time (s)",png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sf = 1.0
  yf = sf
  sp.setVLimits(0.1,1.0)
  if vmin and vmax:
    sp.setVLimits(vmin,vmax)
  sp.setHLimits(0.5,11.5)
  sp.setHInterval(2)
  for il,y in enumerate(ys):
    ya = sum(y)/len(y)
    y = sub(y,ya)
    y = div(y,10)
    y = add(y,yf)
    pv = sp.addPoints(ss[il],y)
    pv.setLineColor(color)
    pv.setLineWidth(1.5)
    yf = yf+sf
  rf = sf
  if rs:
    for il,r in enumerate(rs):
      ra = sum(r)/len(r)
      r = sub(r,ra)
      r = div(r,10)
      r = add(r,rf)
      pv = sp.addPoints(s1,r)
      pv.setLineColor(Color.BLACK)
      pv.setLineWidth(1.5)
      rf = rf+sf
  sp.setSize(600,650)
  sp.setHLabel(hlabel)
  sp.setVLabel(vlabel)
  sp.setFontSize(20)
  sp.setVInterval(0.2)
  if png and pngDir:
    sp.paintToPng(300,7.0,pngDir+png+".png")


def plot2(w,sz,sl,wmin=0.0,wmax=0.0,vlabel="Time (s)",cbar=None,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(500,900)
  sp.setVLabel(vlabel)
  sp.setHLabel("Log index")
  sp.addColorBar(cbar)
  sp.plotPanel.setColorBarWidthMinimum(90)
  pv = sp.addPixels(sz,sl,w)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ColorMap.GRAY)
  pv.setClips(wmin,wmax)
  if png and pngDir:
    sp.paintToPng(300,7.0,pngDir+png+".png")

def plot3(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          slices=None,surf=None,hs=None,logs=None,sps=None,curve=None,
          wmin=0,wmax=0,png=None):
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  f1 = f1+d1*50
  s1s = Sampling(n1-50,d1,f1)
  f = copy(n1-50,n2,n3,50,0,0,f)
  if g!=None:
    g = copy(n1-50,n2,n3,50,0,0,g)

  sf = SimpleFrame(AxesOrientation.XRIGHT_YOUT_ZDOWN)
  cbar = None
  if g==None:
    ipg = sf.addImagePanels(s1s,s2,s3,f)
    if cmap!=None:
      ipg.setColorModel(cmap)
    if wmin!=0 and wmax!=0:
      ipg.setClips(wmin,wmax)
    if cmin!=None and cmax!=None:
      ipg.setClips(cmin,cmax)
    else:
      #ipg.setClips(-2.0,2.0)
      ipg.setClips(-2.0,1.5) # use for subset plots
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1s,s2,s3,f,g)
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
    #cbar.setWidthMinimum(120) # for slides
    cbar.setWidthMinimum(80)
  if logs:
    wg = wellGroup(logs,curve,wmin,wmax)
    sf.world.addChild(wg)
  if sps:
    #samples = sps[0],sps[1],sps[2],sps[3]
    wg = makeLogPoints(sps,wmin,wmax,cbar)
    sf.world.addChild(wg)
  if hs:
    x1 = readImage(ghfile)
    u1 = readImage(gtfile)
    hfr = HorizonFromRgt(s1s,s2,s3,x1,u1)
    for hi in hs:
      [xyz,rgb] = hfr.singleHorizon(hi)
      tg = TriangleGroup(True,xyz,rgb)
      sf.world.addChild(tg)
  if surf:
    tgs = Triangle()
    xyz = tgs.trianglesForSurface(surf,0,n1-1)
    tg  = TriangleGroup(True,xyz)
    sf.world.addChild(tg)
  ipg.setSlices(924,224,68)
  #ipg.setSlices(n1,0,n3) # use only for subset plots
  if cbar:
    sf.setSize(837,700)
  else:
    #sf.setSize(700,700) # for slides
    sf.setSize(740,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  ov = sf.getOrbitView()
  zscale = 0.8*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.1)
  ov.setAzimuthAndElevation(235,25)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(0.0,0.05,0.08))
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      #cbar.setFont(Font("Arial", Font.PLAIN, 36)) #for slides
      cbar.setFont(Font("Arial", Font.PLAIN, 24)) #for print
      cbar.setInterval(0.5)
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")

def plot3F(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
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

def wellGroup(logs,curve,cmin=0,cmax=0,cbar=None):
  print "number of logs =",len(logs)
  #s1 = Sampling(2762,0.002,0.000)
  #s2 = Sampling(357,0.025,0.000)
  #s3 = Sampling(161,0.025,0.000)
  fl,x1l,x2l,x3l = [],[],[],[]
  for log in logs:
    samples = log.getSamples(curve,s1,s2,s3)
    f,x1,x2,x3 = samples
    fl.append(f)
    x1l.append(x1)
    x2l.append(x2)
    x3l.append(x3)
  samples = fl,x1l,x2l,x3l
  lg = makeLogPoints(samples,cmin,cmax,cbar)
  return lg

def makeLogPoints(samples,cmin,cmax,cbar):
  lg = Group()
  fl,x1l,x2l,x3l = samples
  for i,f in enumerate(fl):
    f = fl[i]
    x1 = x1l[i]
    x2 = x2l[i]
    x3 = x3l[i]
    pg = makePointGroup(f,x1,x2,x3,cmin,cmax,cbar)
    lg.addChild(pg)
  return lg

def makePoint(f,x1,x2,x3,cmin,cmax,cbar):
  xyz = zerofloat(3)
  xyz[0],xyz[1],xyz[2]=x3,x2,x1
  rgb = None
  if cmin<cmax:
    cmap = ColorMap(cmin,cmax,ColorMap.JET)
    if cbar:
      cmap.addListener(cbar)
    rgb = cmap.getRgbFloats([f])
  pg = PointGroup(xyz,rgb)
  ps = PointState()
  ps.setSize(4)
  ps.setSmooth(False)
  ss = StateSet()
  ss.add(ps)
  pg.setStates(ss)
  return pg

def makePointGroup(f,x1,x2,x3,cmin,cmax,cbar):
  n = len(x1)
  xyz = zerofloat(3*n)
  copy(n,0,1,x3,0,3,xyz)
  copy(n,0,1,x2,1,3,xyz)
  copy(n,0,1,x1,2,3,xyz)
  rgb = None
  if cmin<cmax:
    cmap = ColorMap(cmin,cmax,ColorMap.JET)
    if cbar:
      cmap.addListener(cbar)
    rgb = cmap.getRgbFloats(f)
  pg = PointGroup(xyz,rgb)
  ps = PointState()
  ps.setSize(4)
  ps.setSmooth(False)
  ss = StateSet()
  ss.add(ps)
  pg.setStates(ss)
  return pg

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

#############################################################################
run(main)
