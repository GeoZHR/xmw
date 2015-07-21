"""
Demonstrate simultaneous multiple-well ties
Author: Xinming Wu, Colorado School of Mines
Version: 2015.06.08
"""

from utils import *
setupForSubset("subt")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count
d1,d2,d3 = s1.delta,s2.delta,s3.delta

# Names and descriptions of image files used below.
gxfile = "gx" # input seismic image 
gtfile = "gt" # RGT volume
ghfile = "gh" # horizon volume
gufile = "gu" # flattened image 
txfile = "tx"
tufile = "tu"
tsfile = "ts"
wshifts = "wellShifts"
wshiftd = "wellShiftsDimensions"
wsample = "wellSampling"

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
pngDir = None
pngDir = "../../../png/swt/"
plotOnly = True

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goDisplay()
  #goSynSeis()
  #goSynSeisFlat()
  #goSynSeisFlatten()
  #goFlattenD()
  goTie()

def goTie():
  vnull = -999.25
  #gu=goFlattenD()
  gu = readImage2(n1,11,tufile)
  sfm,sfs,fu=goSynSeisFlatten()
  esum,e,u=findShifts(sfm,sfs,fu,gu)
  sw = SeismicWellTie()
  ndfs = zerofloat(3,len(fu))
  fus = sw.applyShifts(vnull,sfs,sfm,u,fu,ndfs)
  sss = []
  for il in range(len(fu)):
    ssi = Sampling((int)(ndfs[il][0]),ndfs[il][1],ndfs[il][2])
    sss.append(ssi)
  plot1s(s1,sfs,fu,rs=gu,vmin=0.1,vmax=1.2,vlabel="Relative geologic time",
         png="synsFlattenXR")
  plot1s(s1,sss,fus,rs=gu,vmin=0.1,vmax=1.2,vlabel="Relative geologic time",
         png="synsFlattenXS")

def findShifts(sfm,sfs,fu,gu):
  vnull = -999.25
  umin,umax = -0.00,0.30
  #rmin,rmax = -0.02,0.02
  rmin,rmax = -0.15,0.15
  dmin = 0.1
  dw = DynamicWarpingX(umin,umax,sfm)
  dw.setStrainLimits(rmin,rmax)
  dw.setSmoothness(dmin)
  ss = dw.samplingS
  e = zerofloat(ss.count,sfm.count)
  e = add(e,dw.computeErrorsX(sfs,fu,s1,gu))
  u = dw.findShifts(e)
  esum = 0.0
  si = SincInterpolator()
  for jf in range(sfm.count):
    esum += si.interpolate(ss,e[jf],u[jf])
  plot = True
  if plot:
    sp = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
    pv = sp.addPixels(sfm,ss,pow(transpose(e),0.5))
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv = sp.addPoints(sfm,u)
    pv.setLineColor(Color.WHITE)
    sp.paintToPng(300,7.0,pngDir+"alignError.png")
  return esum,e,u

def invertShifts(ss,u):
  """ 
  Given a uniformly sampled u(s) such that s+u(s) increases monotonically, 
  computes a uniformly sampled r(t) = u(t-r(t)) such that t-r(t) increases 
  monotonically. Uses the following 3-step algorithm:
  (1) computes t(s) = s+u(s)
  (2) computes s(t) = by inverse linear interpolation of t(s)
  (3) computes r(t) = t-s(t)
  Returns the sampling of time st and the sequence of sampled r(t).
  """
  ns,ds,fs,ls = ss.count,ss.delta,ss.first,ss.last
  dt = ds # make sampling intervals equal
  ft = fs+u[0] # ft = time t of first sample
  lt = ls+u[ns-1] # lt = time of last sample
  ft = int(ft/dt)*dt # force ft to be a multiple of interval dt
  nt = 1+int((lt-ft)/dt) # number of t samples
  st = Sampling(nt,dt,ft) # sampling of t
  t = add(rampfloat(fs,ds,ns),u) # t(s) = s+u(s)
  s = zerofloat(nt)
  ii = InverseInterpolator(ss,st)
  ii.invert(t,s) # both t(s) and s(t) increase monotonically
  r = sub(rampfloat(ft,dt,nt),s) # r(t) = t-s(t)
  return st,r


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
  rs = getRealSeis()
  plot1s(s1,ss,sa,hlabel="Seismic traces",vlabel="time (ms)")
  plot1s(s1,ss,sa,rs=rs,hlabel="Seismic traces",vlabel="time (ms)")

def goSynSeisFlatten():
  simple = True  
  logs = getLogs()
  nl = len(logs)
  ndf = zerofloat(3)
  sw = SeismicWellTie()
  gx = sw.computeSyns(simple,logs,ndf)
  sz = Sampling((int)(ndf[0]),ndf[1],ndf[2])
  sl = Sampling(nl,1,0)
  n1 = len(gx)
  n2 = len(gx[0])
  ga = zerofloat(n1,n2,1)
  ga[0] = gx

  #maxShift = 40
  #errorPow = 0.05
  weight = 1.0
  maxShift = 40
  errorPow = 0.05

  wlw = WellLogWarping()
  wlw.setMaxShift(maxShift)
  wlw.setPowError([errorPow])
  ss = wlw.findShifts([weight],ga)
  gu = wlw.applyShiftsW(ga[0],ss)

  cbar   = "Amplitude"
  vlabel = "Relative geologic time"
  wmin,wmax = -2.0,2.0
  #plot2(ga[0],sz,sl,wmin=wmin,wmax=wmax,cbar=cbar,png="syns")
  #plot2(gu,sz,sl,wmin=wmin,wmax=wmax,vlabel=vlabel,cbar=cbar,png="synsFlatten")

  ndfx = zerofloat(3,nl)
  ndfu = zerofloat(3,nl)
  fx = sw.getValidSamples(sz,ga[0],ndfx)
  fu = sw.getValidSamples(sz,gu,ndfu)
  mk = 0
  nu = 0
  ssx = []
  ssu = []
  for il in range(nl):
    sxi = Sampling((int)(ndfx[il][0]),ndfx[il][1],ndfx[il][2])
    sui = Sampling((int)(ndfu[il][0]),ndfu[il][1],ndfu[il][2])
    ssx.append(sxi)
    ssu.append(sui)
    if ndfx[il][0]>nu:
      nu=ndfx[il][0]
      mk=il
  smu = Sampling((int)(ndfu[mk][0]),ndfu[mk][1],ndfu[mk][2])
  plot1s(s1,ssx,fx,vmin=0.1,vmax=1.0,vlabel="Time (s)",png="synsX")
  plot1s(s1,ssu,fu,vmin=0.1,vmax=1.0,vlabel="Relative geologic time",
         png="synsFlattenX")
  return smu,ssu,fu

def goFlattenD():
  gx = readImage(gxfile)
  logs = getLogs()
  gxs = zerofloat(n1,len(logs),1)
  for il, log in enumerate(logs):
    model = SynSeis.getModel(log)
    i2 = s2.indexOfNearest(model.x2)
    i3 = s3.indexOfNearest(model.x3)
    gxs[0][il] = gx[i3][i2]
  weight = 1.0
  maxShift = 40
  errorPow = 0.05

  wlw = WellLogWarping()
  wlw.setMaxShift(maxShift)
  wlw.setPowError([errorPow])
  s = wlw.findShifts([weight],gxs)
  gus = wlw.applyShiftsX(gxs[0],s)
  sst = []
  for i2 in range(len(gus)):
    sst.append(s1)
    for i1 in range(len(gus[0])):
      if(gus[i2][i1]<-10):
        gus[i2][i1] = 0.0
  writeImage(tsfile,s)
  writeImage(tufile,gus)
  writeImage(txfile,gxs[0])
  plot1s(s1,sst,gxs[0],vmin=0.15,vmax=1.5,color=Color.BLACK,
        vlabel="depth (km)",png="originalSeisTraces")
  plot1s(s1,sst,gus,vmin=0.15,vmax=1.5,color=Color.BLACK,
        vlabel="Relative geologic time",png="seisTracesFlattenD")
  return gus

def goSynSeisFlat():
  simple = True
  logs = getLogs()
  nl = len(logs)
  ndft = zerofloat(3,nl)
  ndfz = zerofloat(3,nl)
  ndfw = zerofloat(3,nl)
  sw = SeismicWellTie()
  td = sw.getTimeDepth(logs)
  sa = sw.computeSyns(simple,logs,ndft,ndfz)
  sst = []
  ssz = []
  ssw = []
  for il in range(nl):
    sa[il] = normalize(sa[il])
    sti = Sampling((int)(ndft[il][0]),ndft[il][1],ndft[il][2])
    szi = Sampling((int)(ndfz[il][0]),ndfz[il][1],ndfz[il][2])
    sst.append(sti)
    ssz.append(szi)
  zs = readImage1(3,wsample)
  sn = readImage1(2,wshiftd)
  ss = readImage2((int)(sn[1]),(int)(sn[0]),wshifts)
  sz = Sampling((int)(zs[0]),zs[1],zs[2])
  gw = sw.applySynsFlat(sz,ssz,sst,ndfw,ss,td,sa)
  for il in range(nl):
    swi = Sampling((int)(ndfw[il][0]),ndfw[il][1],ndfw[il][2])
    ssw.append(swi)
  plot1s(s1,sst,sa,hlabel="Seismic traces",vlabel="time (ms)",png="syns")
  #plot1s(s1,ssw,gw,vmin=0.25,vmax=1.85,vlabel="Relative geologic time")
  plot1s(s1,ssw,gw,vmin=0.04,vmax=1.82,vlabel="Relative geologic time",png="synsFlat")

def getRealSeis():
  logs = getLogs()
  gx = readImage(gxfile)
  gxs = zerofloat(n1,len(logs))
  for il, log in enumerate(logs):
    model = SynSeis.getModel(log)
    i2 = s2.indexOfNearest(model.x2)
    i3 = s3.indexOfNearest(model.x3)
    gxs[il] = gx[i3][i2]
  return gxs

def getRealSeisFlat():
  logs = getLogs()
  gu = readImage(gufile)
  gus = zerofloat(n1,len(logs))
  for il, log in enumerate(logs):
    model = SynSeis.getModel(log)
    i2 = s2.indexOfNearest(model.x2)
    i3 = s3.indexOfNearest(model.x3)
    gus[il] = gu[i3][i2]
  return gus

def normalize(f):
  sigma = 100
  n = len(f)
  g = zerofloat(n)
  ref = RecursiveExponentialFilter(sigma)
  g = mul(f,f)
  ref.apply(g,g)
  g = div(f,sqrt(g))
  return g

def goDisplay():
  gx = readImage(gxfile)
  gu = readImage(gufile)
  logs = getLogs()
  yxs,yus=[],[]
  yf = 10
  logs = logs
  gxs = zerofloat(n1,12,1)
  for il, log in enumerate(logs):
    model = SynSeis.getModel(log)
    i2 = s2.indexOfNearest(model.x2)
    i3 = s3.indexOfNearest(model.x3)
    gxs[0][il] = gx[i3][i2]
    yx = add(gx[i3][i2],yf)
    yu = add(gu[i3][i2],yf)
    yf = yf+10
    yxs.append(yx)
    yus.append(yu)
  wlw = WellLogWarping()
  wlw.setMaxShift(50)
  wlw.setPowError([2.0])
  s = wlw.findShifts([1.0],gxs)
  gus = wlw.applyShifts(gxs[0],s)
  wus = []
  uf = 10
  for i in range(12):
    for i1 in range(n1):
      if(gus[i][i1]<-10.0):
        gus[i][i1] = 0.0
    wus.append(add(gus[i],uf))
    uf = uf+10
  plot1(s1,yxs,png="originalSeisTraces")
  plot1(s1,yus,vlabel="Relative geologic time",png="flattenedSeisTraces")
  plot1(s1,wus,vlabel="Relative geologic time (dw)",png="flattenedSeisTraces")


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
  hlabel="Seismic traces",vlabel="time (ms)",png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sf = 10.0
  yf = sf
  sp.setVLimits(0.1,1.0)
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


#############################################################################
run(main)
