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
plotOnly = False

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goDisplay()
  #goSynSeis()
  #goSynsFlatten()
  #goSeisFlatten()
  goSynsSeisTie()

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

def goSynsFlatten():
  smax = 40
  epow = 0.125
  simple = True
  logs = getLogs()
  nl = len(logs)
  ndfx = zerofloat(3,nl)
  ndfu = zerofloat(3,nl)
  swt = SeismicWellTie()
  wx,wu = swt.synsFlatten(simple,smax,epow,logs,ndfx,ndfu)
  ssx = []
  ssu = []
  for il in range(nl):
    sxi = Sampling((int)(ndfx[il][0]),ndfx[il][1],ndfx[il][2])
    sui = Sampling((int)(ndfu[il][0]),ndfu[il][1],ndfu[il][2])
    ssx.append(sxi)
    ssu.append(sui)
  plot1s(s1,ssx,wx,vmin=0.1,vmax=1.0,vlabel="Time (s)",png="synsX")
  plot1s(s1,ssu,wu,vmin=0.1,vmax=1.0,vlabel="Relative geologic time",
         png="synsFlattenX")
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
  fu = swt.seisFlatten(smax,epow,fx)
  sst = []
  for il in range(nl):
    sst.append(s1)
  writeImage(tufile,fu)
  writeImage(txfile,fx)
  plot1s(s1,sst,fx,vmin=0.15,vmax=1.5,color=Color.BLACK,
        vlabel="depth (km)",png="originalSeisTraces")
  plot1s(s1,sst,fu,vmin=0.15,vmax=1.5,color=Color.BLACK,
        vlabel="Relative geologic time",png="seisTracesFlattenD")
  return fu

def goSynsSeisTie():
  gu=goSeisFlatten() 
  wu,sw=goSynsFlatten()
  dmin = 0.1
  smin,smax = -0.10,0.30
  rmin,rmax = -0.15,0.15
  swt = SeismicWellTie()
  ndf = zerofloat(3,len(wu))
  wgu = swt.synsToSeis(sw,s1,smin,smax,rmin,rmax,dmin,wu,gu,ndf)
  swg = []
  for il in range(len(wu)):
    si = Sampling((int)(ndf[il][0]),ndf[il][1],ndf[il][2])
    swg.append(si)
  plot1s(s1,sw,wu,rs=gu,vmin=0.1,vmax=1.2,vlabel="Relative geologic time",
         png="synsFlattenXR")
  plot1s(s1,swg,wgu,rs=gu,vmin=0.1,vmax=1.2,vlabel="Relative geologic time",
         png="synsFlattenXS")

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
