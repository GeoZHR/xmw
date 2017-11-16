"""
Demo of generating avo models
Author: Xinming Wu, University of Texas at Austin
Version: 2017.05.05
"""

from utils import *
setupForSubset("bh")
#setupForSubset("tp")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

# Names and descriptions of image files used below.
gxfile = "gx"
gdfile = "rou"
gsfile = "gxs"
dsfile = "gds"
gcfile = "gxc"
p2file = "p2"
p3file = "p3"
epfile = "ep"

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
pngDir = None
pngDir = "../../../png/bh/"
plotOnly = False
minPhi,maxPhi = 0,360
minTheta,maxTheta = 65,89
sigmaPhi,sigmaTheta = 10,35

# These parameters control the construction of fault skins.
# See the class FaultSkinner for more information.
lowerLikelihood = 0.02
upperLikelihood = 0.6
minSkinSize = 4000

# These parameters control the computation of fault dip slips.
# See the class FaultSlipper for more information.
minThrow = 0.0
maxThrow = 30.0

def main(args):
  #goTopBottomHorizons()
  goSlopes()
  goScan()
def goDen():
  gx = readImage3D(n1,n2,n3,gxfile)
  gd = readImage3D(n1,n2,n3,gdfile)
  for i3 in range(n3):
    gt = copy(gd[i3])
    for i2 in range(n2):
      gd[i3][i2] = gt[n2-i2-1]
  ds = copy(550,n2,n3,780,0,0,gd)
  gs = copy(550,n2,n3,780,0,0,gx)
  writeImage("gds",ds)
  writeImage("gxs",gs)
  plot3(gs)
  plot3(gs,ds,cmin=2.2,cmax=2.7,cmap=jetRamp(0.4))
def goTopBottomHorizons():
  gs = readImage3D(n1,n2,n3,"gxs")
  ds = readImage3D(n1,n2,n3,"gds")
  hp = Helper()
  hs = hp.getTopBottomHorizons(ds)
  gc = hp.getSection(hs,gs)
  writeImage(gcfile,gc)
  plot3(gs)
  plot3(gs,ds,cmin=2.2,cmax=2.7,cmap=jetRamp(1.0))
  plot3(gc,ds,cmin=2.2,cmax=2.7,cmap=jetRamp(1.0))

def goSlopes():
  print "goSlopes ..."
  gx = readImage(gcfile)
  sigma1,sigma2,sigma3,pmax = 16.0,2.0,2.0,5.0
  p2,p3,ep = FaultScanner.slopes(sigma1,sigma2,sigma3,pmax,gx)
  writeImage(p2file,p2)
  writeImage(p3file,p3)
  writeImage(epfile,ep)
  zm = ZeroMask(0.10,1,1,1,gx)
  zero,tiny=0.0,0.01
  zm.setValue(zero,p2)
  zm.setValue(zero,p3)
  zm.setValue(tiny,ep)
  print "p2  min =",min(p2)," max =",max(p2)
  print "p3  min =",min(p3)," max =",max(p3)
  print "ep min =",min(ep)," max =",max(ep)
  '''
  plot3(gx,p2, cmin=-1,cmax=1,cmap=bwrNotch(1.0),
        clab="Inline slope (sample/sample)",png="p2")
  plot3(gx,p3, cmin=-1,cmax=1,cmap=bwrNotch(1.0),
        clab="Crossline slope (sample/sample)",png="p3")
  plot3(gx,sub(1,ep),cmin=0,cmax=1,cmap=jetRamp(1.0),
        clab="Planarity")
  '''

def goScan():
  print "goScan ..."
  if not plotOnly:
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    gx = readImage(gcfile)
    zm = ZeroMask(0.10,1,1,1,gx)
    zero,tiny=0.0,0.01
    zm.setValue(tiny,gx)
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
  plot3(gx,ft,cmin=65,cmax=89,cmap=jetFill(1.0),
      clab="Fault dip (degrees)",png="ft")
  '''

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
  plot3(gx,flt,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",png="flt")
  plot3(gx,ftt,cmin=65,cmax=85,cmap=jetFillExceptMin(1.0),
        clab="Fault dip (degrees)",png="ftt")
  plot3(gx,fpt,cmin=0,cmax=360,cmap=hueFillExceptMin(1.0),
        clab="Fault strike (degrees)",cint=45,png="fpt")

def goSlopes():
  gx = readImage3D(n1,n2,n3,gsfile)
  p2 = zerofloat(n1,n2,n3)
  p3 = zerofloat(n1,n2,n3)
  ep = zerofloat(n1,n2,n3)
  sigma1,sigma2=4.0,2.0
  lsf = LocalSlopeFinder(sigma1,sigma2,sigma2,5) 
  lsf.findSlopes(gx,p2,p3,ep);
  writeImage(p2file,p2)
  writeImage(p3file,p3)
  writeImage(epfile,ep)

def goModelSmooth():
  n1,n2=501,501
  rh = readImage2D(n1,n2,"Rho")
  vp = readImage2D(n1,n2,"Vp")
  vs = readImage2D(n1,n2,"Vs")
  vp = div(304800.0,vp)
  vs = div(304800.0,vs)
  sx = readImage2DL(n1,n2,sxFile)
  lof = LocalOrientFilter(3.0,1.0)
  tensors = lof.applyForTensors(sx)
  tensors.setEigenvalues(0.8,1.0)
  lsf = LocalSmoothingFilter()
  rh1 = zerofloat(n1,n2)
  rh2 = zerofloat(n1,n2)
  vp1 = zerofloat(n1,n2)
  vp2 = zerofloat(n1,n2)
  vs1 = zerofloat(n1,n2)
  vs2 = zerofloat(n1,n2)
  lsf.applySmoothS(rh,rh1)
  lsf.applySmoothS(vp,vp1)
  lsf.applySmoothS(vs,vs1)
  lsf.apply(tensors,600,rh1,rh2)
  lsf.apply(tensors,600,vp1,vp2)
  lsf.apply(tensors,600,vs1,vs2)
  s1,s2=Sampling(n1),Sampling(n2)
  writeImageL("rhs600",rh2)
  writeImageL("vps600",vp2)
  writeImageL("vss600",vs2)

  plotLogs(s1,s2,rh,wh=700,wv=600,cmin=2.0,cmax=2.7,
    hint=2,vlab="Depth (Samples)",cbar="Density",png="Rh")
  plotLogs(s1,s2,vp,wh=700,wv=600,cmin=2500,cmax=5000,
    hint=2,vlab="Depth (Samples)",cbar="Vp",png="Vp")
  plotLogs(s1,s2,vs,wh=700,wv=600,cmin=900,cmax=2300,
    hint=2,vlab="Depth (Samples)",cbar="Vs",png="Vs")
  plotLogs(s1,s2,rh2,wh=700,wv=600,cmin=2.0,cmax=2.7,
    hint=2,vlab="Depth (Samples)",cbar="Density",png="Rh2")
  plotLogs(s1,s2,vp2,wh=700,wv=600,cmin=2500,cmax=5000,
    hint=2,vlab="Depth (Samples)",cbar="Vp",png="Vp2")
  plotLogs(s1,s2,vs2,wh=700,wv=600,cmin=900,cmax=2300,
    hint=2,vlab="Depth (Samples)",cbar="Vs",png="Vs2")

def goLogCorrelation():
  lmin,lmax=-350,0
  rv2 = readTxtLogs(logFile2)
  rv3 = readTxtLogs(logFile3)
  rv5 = readTxtLogs(logFile5)
  rv6 = readTxtLogs(logFile6)
  print len(rv2[0])
  print len(rv3[0])
  print len(rv5[0])
  print len(rv6[0])
  writeImage(rvFile2,rv2) #np2 = 4450
  writeImage(rvFile3,rv3) #np3 = 3760
  writeImage(rvFile5,rv5) #np5 = 4367
  writeImage(rvFile6,rv6) #np6 = 3800
  hp = Helper()
  las = hp.sortLogs(rv2,rv5,rv6,rv3)
  las = copy(2500,4,3,0,0,0,las)
  nz = len(las[0][0])
  logs = las
  logs = zerofloat(nz,3,3)
  logs[0][0] = las[0][0]
  logs[0][1] = las[0][2]
  logs[0][2] = las[0][3]
  logs[1][0] = las[1][0]
  logs[1][1] = las[1][2]
  logs[1][2] = las[1][3]
  logs[2][0] = las[2][0]
  logs[2][1] = las[2][2]
  logs[2][2] = las[2][3]
  nz = len(logs[0][0])
  nw = len(logs[0])
  sz = Sampling(nz,1,1)
  sw = Sampling(nw,1,1)
  ww = WellFlattener(lmin,lmax)
  ww.setErrorExponent(0.125)
  ww.setStrainMax(1.0)
  nl = lmax-lmin+1
  flogs = ww.flatten(logs)
  sl = Sampling(nl)
  s1 = Sampling(nz)
  clab1 = "Density (g/cc)"
  clab2 = "Vp"
  clab3 = "Vs"
  plotLogs(sz,sw,logs[0],cmin=2.,cmax=2.7,cbar=clab1,png="den")
  plotLogs(sz,sw,logs[1],cmin=60,cmax=115,cbar=clab2,png="vp")
  plotLogs(sz,sw,logs[2],cmin=140,cmax=250,cbar=clab3,png="vs")
  vlab = "Relative geologic time"
  plotLogs(sz,sw,flogs[0],cmin=2.,cmax=2.7,cbar=clab1,png="denf")
  plotLogs(sz,sw,flogs[1],cmin=60,cmax=115,cbar=clab2,png="vpf")
  plotLogs(sz,sw,flogs[2],cmin=140,cmax=250,cbar=clab3,png="vsf")

def goModels():
  rv2 = readTxtLogs(logFile2)
  rv3 = readTxtLogs(logFile3)
  rv5 = readTxtLogs(logFile5)
  rv6 = readTxtLogs(logFile6)
  print len(rv2[0])
  print len(rv3[0])
  print len(rv5[0])
  print len(rv6[0])
  writeImage(rvFile2,rv2) #np2 = 4450
  writeImage(rvFile3,rv3) #np3 = 3760
  writeImage(rvFile5,rv5) #np5 = 4367
  writeImage(rvFile6,rv6) #np6 = 3800
  hp = Helper()
  rv21 = copy(rv2)
  rv22 = copy(rv2)
  rv23 = copy(rv2)
  mul(rv21[0],1.005,rv21[0])
  mul(rv22[0],1.000,rv22[0])
  mul(rv23[0],0.99,rv23[0])

  mul(rv21[1],1.01,rv21[1])
  mul(rv22[1],0.98,rv22[1])
  mul(rv23[1],0.94,rv23[1])

  mul(rv21[2],1.02,rv21[2])
  mul(rv22[2],0.96,rv22[2])
  mul(rv23[2],0.93,rv23[2])
  la = hp.sortLogs(rv2,rv21,rv22,rv23)
  fx = zerofloat(4450,4)
  nz = 4450
  sz = Sampling(len(la[0][0]))
  sw = Sampling(4)
  plotLogs(sz,sw,la[0],wh=400,wv=800,cmin=2.0,cmax=2.7,
    hint=2,vlab="Depth (Samples)",cbar="Density")
  plotLogs(sz,sw,la[1],wh=400,wv=800,cmin=60,cmax=115,
    hint=2,vlab="Depth (Samples)",cbar="Vp")
  plotLogs(sz,sw,la[2],wh=400,wv=800,cmin=140,cmax=250,
    hint=2,vlab="Depth (Samples)",cbar="Vs")
  lar = hp.resampleLogs(501,7,la)
  rvs = FakeData.densityAndVelocity2d(0.0,lar)
  lof = LocalOrientFilter(3.0,1.0)
  tensors = lof.applyForTensors(rvs[0])
  tensors.setEigenvalues(0.2,1.0)
  lsf = LocalSmoothingFilter()
  #lsf.applySmoothS(rvs[0],rvs[0])
  lsf.apply(tensors,6,rvs[0],rvs[0])
  lsf.apply(tensors,6,rvs[1],rvs[1])
  lsf.apply(tensors,6,rvs[2],rvs[2])
  s1,s2=Sampling(501),Sampling(501)
  plotLogs(s1,s2,rvs[0],wh=700,wv=600,cmin=2.0,cmax=2.7,
    hint=2,vlab="Depth (Samples)",cbar="Density",png="Rho")
  plotLogs(s1,s2,rvs[1],wh=700,wv=600,cmin=70,cmax=115,
    hint=2,vlab="Depth (Samples)",cbar="Vp",png="Vp")
  plotLogs(s1,s2,rvs[2],wh=700,wv=600,cmin=140,cmax=250,
    hint=2,vlab="Depth (Samples)",cbar="Vs",png="Vs")
  writeImage("Rho",rvs[0])
  writeImage("Vp",rvs[1])
  writeImage("Vs",rvs[2])
def goInterpolation():
  n1,n2,nt=501,501,15
  s1,s2=Sampling(n1),Sampling(n2)
  sx = readImage2D(n1,n2,sxFile)
  ra = readImage3D(n1,4,nt,raFile)
  rb = readImage3D(n1,4,nt,rbFile)
  rp = readImage3D(n1,4,nt,rpFile)
  x2s = [60,180,320,420]
  if not plotOnly:
    lof = LocalOrientFilter(3.0,1.0)
    tensors = lof.applyForTensors(sx)
    tensors.setEigenvalues(0.001,1.0)
    sa = zerofloat(n1,n2,nt)
    sb = zerofloat(n1,n2,nt)
    sp = zerofloat(n1,n2,nt)
    fnull = -1.0
    for it in range(nt):
      fa,x1,x2=getPoints(ra[it],x2s)
      fb,x1,x2=getPoints(rb[it],x2s)
      fp,x1,x2=getPoints(rp[it],x2s)
      pa = goSimpleGridder(s1,s2,fa,x1,x2,fnull)
      pb = goSimpleGridder(s1,s2,fb,x1,x2,fnull)
      pp = goSimpleGridder(s1,s2,fp,x1,x2,fnull)
      qa = goBlendGridder(tensors,pa,fa,x1,x2,fnull)
      qb = goBlendGridder(tensors,pb,fb,x1,x2,fnull)
      qp = goBlendGridder(tensors,pp,fp,x1,x2,fnull)
      sa[it],sb[it],sp[it] = qa,qb,qp
    writeImageL(saFile,sa)
    writeImageL(sbFile,sb)
    writeImageL(spFile,sp)
  else:
    sa = readImage3D(n1,n2,nt,saFile)
    sb = readImage3D(n1,n2,nt,sbFile)
    sp = readImage3D(n1,n2,nt,spFile)
  kt = 13
  fa,x1,x2=getPoints(ra[kt],x2s)
  fb,x1,x2=getPoints(rb[kt],x2s)
  fp,x1,x2=getPoints(rp[kt],x2s)
  plot2(fa,x1,x2,sx,s1,s2,gmin=0.5,gmax=1,label="Aas",png="ra")
  plot2(fb,x1,x2,sx,s1,s2,gmin=-0.8,gmax=-0.1,label="Abs",png="rb")
  plot2(fp,x1,x2,sx,s1,s2,gmin=-1,gmax=1,label="Aps",png="rp")
  plot2(fa,x1,x2,sx,s1,s2,g=sa[kt],gmin=0.5,gmax=1,label="Aas",png="rai")
  plot2(fb,x1,x2,sx,s1,s2,g=sb[kt],gmin=-0.8,gmax=-0.1,label="Abs",png="rbi")
  plot2(fp,x1,x2,sx,s1,s2,g=sp[kt],gmin=-1,gmax=1,label="Aps",png="rpi")

def getPoints(rx,x2s):
  n2 = len(rx)
  n1 = len(rx[0])
  fx = zerofloat(n1*n2)
  x1 = zerofloat(n1*n2)
  x2 = zerofloat(n1*n2)
  for i2 in range(n2):
    for i1 in range(n1):
      x1[n1*i2+i1] = i1
      x2[n1*i2+i1] = x2s[i2]
      fx[n1*i2+i1] = rx[i2][i1]
  return fx,x1,x2

def goBlendGridder(tensors,px,fx,x1,x2,fnull):
  smooth = 100.0
  bg = BlendedGridder2(fx,x1,x2)
  bg.setSmoothness(smooth)
  bg.setTimeMax(FLT_MAX)
  bg.setTensors(tensors)
  dx = bg.gridNearest(fnull,px)
  qx = copy(px)
  bg.gridBlended(dx,px,qx)
  return qx

def goSimpleGridder(s1,s2,fx,x1,x2,fnull):
  sg = SimpleGridder2(fx,x1,x2)
  sg.setNullValue(fnull)
  return sg.grid(s1,s2)

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
  ref = RecursiveExponentialFilter(100.0)
  ref.apply1(g,g)
  y = like(x)
  div(x,sqrt(g),y)
  return y

#############################################################################
# graphics
#############################################################################
# plotting
backgroundColor = Color.WHITE
cjet = ColorMap.JET
alpha = fillfloat(1.0,256); alpha[0] = 0.0
ajet = ColorMap.setAlpha(cjet,alpha)
def plotLogs(sz,sw,fx,wh=500,wv=800,cmin=None,cmax=None,
    hint=100,vlab="Depth (samples)",cbar=None,png=None):
  fz = sz.first
  dz = sz.delta
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(wh,wv)
  #sp.setHInterval(hint)
  sp.setVLabel(vlab)
  sp.setHLabel("Lateral position (samples)")
  sp.addColorBar(cbar)
  sp.plotPanel.setColorBarWidthMinimum(85)
  #sp.setHLimits(sw.first-sw.delta/2,sw.last)
  #sp.setVLimits(sz.first,sz.last)
  pv = sp.addPixels(sz,sw,fx)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cjet)
  n1 = len(fx[0])
  sp.setFontSize(20)
  if cmin and cmax:
    pv.setClips(cmin,cmax)
  if png and pngDir:
    sp.paintToPng(1080,3.33333,pngDir+png+".png")


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

def addTensorsInImage(ip,et,esize):
  tp = TensorsPanel(s1,s2,s3,et)
  tp.setEllipsoidSize(esize)
  ip.getFrame().addChild(tp)
  return tp

def addColorBar(frame,clab=None,cint=None):
  cbar = ColorBar(clab)
  if cint:
    cbar.setInterval(cint)
  cbar.setFont(Font("Arial",Font.PLAIN,24)) # size by experimenting
  cbar.setWidthMinimum
  cbar.setBackground(Color.WHITE)
  frame.add(cbar,BorderLayout.EAST)
  return cbar

def convertDips(ft):
  return FaultScanner.convertDips(0.2,ft) # 5:1 vertical exaggeration

def makePointGroup(f,x1,x2,x3,cmin,cmax,cbar):
  n = len(x1)
  xyz = zerofloat(3*n)
  copy(n,0,1,x3,0,3,xyz)
  copy(n,0,1,x2,1,3,xyz)
  copy(n,0,1,x1,2,3,xyz)
  rgb = None
  if cmin<cmax:
    cmap = ColorMap(cmin,cmax,ColorMap.getJet(0.3))
    if cbar:
      cmap.addListener(cbar)
    rgb = cmap.getRgbFloats(f)
  pg = PointGroup(xyz,rgb)
  ps = PointState()
  ps.setSize(8)
  ps.setSmooth(False)
  ss = StateSet()
  ss.add(ps)
  pg.setStates(ss)
  return pg

def plot1s(s1,ys,rs=None,vmin=None,vmax=None,color=Color.RED,
  hlabel="Seismic traces",vlabel="time (ms)",png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sf = 1.0
  yf = sf
  sp.setVLimits(0,n1)
  if vmin and vmax:
    sp.setVLimits(vmin,vmax)
  sp.setHLimits(0,len(ys))
  for il,y in enumerate(ys):
    ya = sum(y)/len(y)
    y = sub(y,ya)
    y = add(y,yf)
    pv = sp.addPoints(s1,y)
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
  sp.setSize(600,500)
  sp.setHLabel(hlabel)
  sp.setVLabel(vlabel)
  if png and pngDir:
    sp.paintToPng(300,7.0,pngDir+png+".png")

def plot2(f,x1,x2,s,s1,s2,g=None,gmin=None,gmax=None,
                label=None,png=None,et=None):
  n2 = len(s)
  n1 = len(s[0])
  panel = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.LEFT_TOP)
  panel.setVLimits(0,500)
  panel.setHInterval(100)
  panel.setVInterval(100)
  panel.setHLabel("Lateral position (sample)")
  panel.setVLabel("Time (sample)")
  #panel.setHInterval(100.0)
  #panel.setVInterval(100.0)
  #panel.setHLabel("Pixel")
  #panel.setVLabel("Pixel")
  if label:
    panel.addColorBar(label)
  else:
    panel.addColorBar()
  panel.setColorBarWidthMinimum(90)
  pv = panel.addPixels(s1,s2,s)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setColorModel(ColorMap.GRAY)
  pv.setClips(-0.1,0.1)
  if g:
    alpha = 1.0
  else:
    g = zerofloat(s1.count,s2.count)
    alpha = 0.0
  pv = panel.addPixels(s1,s2,g)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setColorModel(ColorMap.getJet(alpha))
  if label and label[0]=="T":
    pv.setClips(0.0,1000.0)
    pv.setClips(min(g),max(g))
  else:
    if gmin:
      pv.setClips(gmin,gmax)
    else:
      pv.setClips(2.0,2.8)
  cmap = ColorMap(gmin,gmax,ColorMap.JET)
  fs,x1s,x2s = makePointSets(cmap,f,x1,x2)
  for i in range(len(fs)):
    #color = cmap.getColor((fs[i][0]-min(fs))/max(fs))
    color = cmap.getColor((fs[i][0]))
    #color = Color(color.red,color.green,color.blue)
    pv = panel.addPoints(x1s[i],x2s[i])
    pv.setLineStyle(PointsView.Line.NONE)
    pv.setMarkStyle(PointsView.Mark.FILLED_SQUARE)
    pv.setMarkSize(4)
    pv.setMarkColor(color)
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  #frame.setFontSizeForPrint(8,240)
  #frame.setSize(1240,774)
  #frame.setFontSizeForSlide(1.0,0.8)
  frame.setSize(1290,777)
  frame.setSize(900,777)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(400,3.2,pngDir+png+".png")
  return frame
def makePointSets(cmap,f,x1,x2):
  sets = {}
  for i in range(len(f)):
    if f[i] in sets:
      points = sets[f[i]]
      points[0].append(f[i])
      points[1].append(x1[i])
      points[2].append(x2[i])
    else:
      points = [[f[i]],[x1[i]],[x2[i]]] # lists of f, x1, x2
      sets[f[i]] = points
  ns = len(sets)
  fs = zerofloat(1,ns)
  x1s = zerofloat(1,ns)
  x2s = zerofloat(1,ns)
  il = 0
  for points in sets:
    fl = sets[points][0]
    x1l = sets[points][1]
    x2l = sets[points][2]
    nl = len(fl)
    fs[il] = zerofloat(nl)
    x1s[il] = zerofloat(nl)
    x2s[il] = zerofloat(nl)
    copy(fl,fs[il])
    copy(x1l,x1s[il])
    copy(x2l,x2s[il])
    il += 1
  return fs,x1s,x2s

def plot3(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          hs=None,skins=None,smax=0.0,slices=None, htgs=None,
          et=None,size=20,samples=None,png=None):
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
      ipg.setClips(-1.,1.) # use for subset plots
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-1.,1.)
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
    cbar.setWidthMinimum(137)
  if htgs:
    for htg in htgs:
      sf.world.addChild(htg)
  if hs:
    for hi in hs:
      tg = TriangleGroup(True,s3,s2,hi)
      tg.setColor(Color.CYAN)
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
    sf.world.addChild(sg)
  ipg.setSlices(106,138,59)
  #ipg.setSlices(92,140,59)
  if et:
    addTensorsInImage(ipg.getImagePanel(Axis.X),et,size)
    addTensorsInImage(ipg.getImagePanel(Axis.Y),et,size)
    addTensorsInImage(ipg.getImagePanel(Axis.Z),et,size)
  if samples:
    fx,x1,x2,x3 = samples
    vmin,vmax,vmap= 6000,16000,ColorMap.JET
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
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")

def plot3c(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          skins=None,smax=0.0,slices=None, htgs=None,hz=None,mk=-1,
          et=None,w=True,samples=None,png=None):
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
    cbar.setWidthMinimum(85)
  if htgs:
    for htg in htgs:
      sf.world.addChild(htg)
  if et:
    tv = TensorView()
    if w:
      hs = tv.applyForSegmentsW(2,et,hz)
    else:
      hs = tv.applyForSegmentsV(2,et,hz,skins)
    cp = ColorMap(0,1,ColorMap.JET)
    vi = fillfloat(0.9,6)
    cb = cp.getRgbFloats(vi)
    for hi in hs:
      lg = LineGroup(hi,cb)
      ss = StateSet()
      lg.setStates(ss)
      ls = LineState()
      ls.setWidth(8)
      ls.setSmooth(False)
      ss.add(ls)
      sf.world.addChild(lg)
  '''
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
    sf.world.addChild(sg)
  '''
  ipg.setSlices(n1,138,59)
  #ipg.setSlices(92,140,59)
  if hz:
    sd = SurfaceDisplay()
    ts = sd.horizonWithAmplitude(mk,[cmin,cmax],hz,f)
    tg = TriangleGroup(True,ts[0],ts[1])
    sf.world.addChild(tg)
  if samples:
    fx,x1,x2,x3 = samples
    vmin,vmax,vmap= 6000,16000,ColorMap.JET
    pg = makePointGroup(fx,x1,x2,x3,vmin,vmax,None)
    sf.world.addChild(pg)
  if cbar:
    sf.setSize(852,700)
  else:
    sf.setSize(750,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.48*sqrt(n1*n1+n2*n2+n3*n3)
  zscale = 0.80*max(n2*d2,n3*d3)/(n1*d1)
  ov = sf.getOrbitView()
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setWorldSphere(BoundingSphere(0.5*n1,0.4*n2,0.4*n3,radius))
  ov.setAzimuthAndElevation(115.0,40.0)
  ov.setTranslate(Vector3(0.02,0.16,-0.3))
  ov.setScale(1.5)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")

#############################################################################
run(main)
