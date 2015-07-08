"""
Interpolates scattered data, such as data from well logs.
"""
from tputils import *

#setupForSubset("subz_51_4_1400")
#setupForSubset("subz_401_4_400")
setupForSubset("subz_401_4_600")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.getCount(),s2.getCount(),s3.getCount()
d1,d2,d3 = s1.getDelta(),s2.getDelta(),s3.getDelta()
f1,f2,f3 = s1.getFirst(),s2.getFirst(),s3.getFirst()
#s1,s2,s3 = Sampling(n1),Sampling(n2),Sampling(n3)
#d1,d2,d3 = s1.getDelta(),s2.getDelta(),s3.getDelta()
fmin,fmax=-5,5 # seismic values
method = "b" # blended
logSet = "d" # deep logs only
logType="v"; logLabel="Velocity (km/s)"; vmin,vmax = 2.4,5.6;k1,k2,k3 = 366,15,95
#logType="d"; logLabel="Density (g/cc)";  vmin,vmax = 2.0,2.8;k1,k2,k3 = 366,15,77
#logType = "p"; logLabel = "Porosity"; vmin,vmax = 0.0,0.4
#logType = "g"; logLabel = "Gamma ray (API units)"; vmin,vmax = 0.0,200.0
smin,smax = -5.5,5.5
smooth = 50 # half-width of smoothing filter for logs
smooth = 0 # no smoothing

sfile = "tpsz" # seismic image
efile = "tpet" # eigen-tensors (structure tensors)
esfile = "tpets" # eigen-tensors scaled by semblances
s1file = "tps1" # semblance w,uv
s2file = "tps2" # semblance vw,u
s3file = "tps3" # semblance uvw,
gfile = "tpg"+logType # simple gridding with null for unknown samples
pfile = "tpp"+logType+method # values of nearest known samples
qfile = "tpq"+logType+method # output of blended gridder
tfile = "tpt"+logType+method # times to nearest known samples

flfile  = "fl" # fault likelihood
fpfile  = "fp" # fault strike (phi)
ftfile  = "ft" # fault dip (theta)
gtfile  = "gt" # RGT volume
ghfile  = "gh" # horizon volume
gufile  = "gu" # flattened image 
p2file  = "p2" # inline slopes
p3file  = "p3" # crossline slopes
epfile  = "ep" # eigenvalue-derived planarity
wpfile  = "wp" # weight image for flattening

fqfile = "fq"+logType # output of blended gridder
fpfile = "fp"+logType # output of blended gridder
ftfile = "ft"+logType # output of blended gridder

fw1file = "fw1"+logType # one well log data
fq1file = "fq1"+logType # output of blended gridder
fp1file = "fp1"+logType # output of blended gridder
ft1file = "ft1"+logType # output of blended gridder

horizon = "CrowMountainCRMT"
"""
horizons = [
  "CrowMountainCRMT",
  "TensleepASand",
  "TensleepBbaseC1Dolo"]
"""

pngDir = "../../../png/rgi/tp3d/"
pngDir = None

plotOnly = True

def main(args):
  #goSlopes()
  #goFaultLikelihoods()
  #goControlSurfs()
  #goRgt()
  goRgtInterp()
  goInterpO()
  #goOneWell()
  #goRgtInterpOne()

def goSlopes():
  print "goSlopes ..."
  if not plotOnly:
    # set half-width of smoother for computing structure tensors
    sig1 = 8.0 #half-width in vertical direction
    sig2 = 2.0 #half-width in one literal direction
    sig3 = 2.0 #half-width in another literal direction
    pmax = 5.0 #maximum slope returned by this slope finder
    gx = readImage(sfile)
    p2 = zerofloat(n1,n2,n3)
    p3 = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    lsf = LocalSlopeFinder(sig1,sig2,sig3,pmax)
    lsf.findSlopes(gx,p2,p3,ep);

    zm = ZeroMask(0.20,1,1,1,gx)
    zero,tiny=0.0,0.01
    zm.setValue(zero,p2)
    zm.setValue(zero,p3)
    zm.setValue(tiny,ep)

    writeImage(p2file,p2)
    writeImage(p3file,p3)
    writeImage(epfile,ep)
    print "p2  min =",min(p2)," max =",max(p2)
    print "p3  min =",min(p3)," max =",max(p3)
    print "ep  min =",min(ep)," max =",max(ep)
  else:
    gx = readImage(sfile)
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    ep = readImage(epfile)
    print min(gx)
    print max(gx)
  plot3(gx)
  plot3(gx,p2, cmin=-1,cmax=1,cmap=jetRamp(1.0),
      clab="Inline slope (sample/sample)",png="p2")
  plot3(gx,p3, cmin=-1,cmax=1,cmap=jetRamp(1.0),
      clab="Crossline slope (sample/sample)",png="p3")
  plot3(gx,pow(ep,4.0),cmin=0,cmax=1,cmap=jetRamp(1.0),
      clab="Planarity")

def goFaultLikelihoods():
  if not plotOnly:
    minPhi,maxPhi = -20,20
    minTheta,maxTheta = 80,85
    sigmaPhi,sigmaTheta = 8,30
    gx = readImage(sfile)
    sig1,sig2,sig3,pmax = 16.0,1.0,1.0,5.0
    p2,p3,ep = FaultScanner.slopes(sig1,sig2,sig3,pmax,gx)
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
    gx = readImage(sfile)
    fl = readImage(flfile)
    fp = readImage(fpfile)
    ft = readImage(ftfile)
  plot3(gx,fl,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
      clab="Fault likelihood",png="fl")
  plot3(gx,fp,cmin=-20,cmax=20,cmap=jetFill(1.0),
      clab="Fault strike (degrees)",cint=5,png="fp")
  plot3(gx,ft,cmin=80,cmax=85,cmap=jetFill(1.0),
      clab="Fault dip (degrees)",png="ft")
def goControlSurfs():
  k01 = [306,309,338,317,313]
  k02 = [ 86,157,331, 17, 56]
  k03 = [ 96, 96, 55, 55, 55]

  k11 = [113,117,116,115]
  k12 = [122,149, 59,280]
  k13 = [ 96, 96, 55, 50]

  k21 = [ 56, 67, 61, 91]
  k22 = [104,137, 68,332]
  k23 = [ 96, 96, 58, 52]

  k31 = [ 27, 34, 33, 31, 68]
  k32 = [ 99,156,155, 62,309]
  k33 = [ 96, 96, 85, 62, 23]

  k41 = [346,341,374,356,352]
  k42 = [202,231,232, 71, 52]
  k43 = [ 67, 68,115,113, 93]

  k51 = [ 90,101,110, 93]
  k52 = [102,135,309, 57]
  k53 = [ 96, 95, 53, 64]

  k61 = [185,189,194,208]
  k62 = [121,148, 56,330]
  k63 = [ 96, 96, 46, 46]

  k71 = [142,147,146,156,175]
  k72 = [122,162, 59,141,346]
  k73 = [ 96, 96, 55, 55, 65]

  k81 = [229,230,237,260]
  k82 = [144,157, 48,335]
  k83 = [ 96, 96, 42, 41]

  k91 = [280,276,304,276]
  k92 = [117,200,326, 51]
  k93 = [101,101, 67, 54]


  gx = readImage(sfile)
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  ep = readImage(epfile)
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  fs = FaultSkinner()
  cells = fs.findCells([fl,fp,ft])
  fci = FaultCellsToImage(cells)
  fli = fci.getLikelihoodThick(0.6,fl)
  wp = sub(1.0,fli)
  wp = pow(wp,6.0)
  ep = pow(ep,4.0)
  wp = mul(wp,ep)
  sc = SetupConstraints()
  kk0 = sc.extend(k01,k02,k03,n2,n3,p2,p3,wp,gx)
  kk1 = sc.extend(k11,k12,k13,n2,n3,p2,p3,wp,gx)
  kk2 = sc.extend(k21,k22,k23,n2,n3,p2,p3,wp,gx)
  kk3 = sc.extend(k31,k32,k33,n2,n3,p2,p3,wp,gx)
  kk4 = sc.extend(k41,k42,k43,n2,n3,p2,p3,wp,gx)
  kk5 = sc.extend(k51,k52,k53,n2,n3,p2,p3,wp,gx)
  kk6 = sc.extend(k61,k62,k63,n2,n3,p2,p3,wp,gx)
  kk7 = sc.extend(k71,k72,k73,n2,n3,p2,p3,wp,gx)
  kk8 = sc.extend(k81,k82,k83,n2,n3,p2,p3,wp,gx)
  kk9 = sc.extend(k91,k92,k93,n2,n3,p2,p3,wp,gx)
  writeImage("kk0",kk0)
  writeImage("kk1",kk1)
  writeImage("kk2",kk2)
  writeImage("kk3",kk3)
  writeImage("kk4",kk4)
  writeImage("kk5",kk5)
  writeImage("kk6",kk6)
  writeImage("kk7",kk7)
  writeImage("kk8",kk8)
  writeImage("kk9",kk9)
  nks = zerofloat(10,1)
  nks[0][0]=len(kk0[0])
  nks[0][1]=len(kk1[0])
  nks[0][2]=len(kk2[0])
  nks[0][3]=len(kk3[0])
  nks[0][4]=len(kk4[0])
  nks[0][5]=len(kk5[0])
  nks[0][6]=len(kk6[0])
  nks[0][7]=len(kk7[0])
  nks[0][8]=len(kk8[0])
  nks[0][9]=len(kk9[0])
  writeImage("nks",nks)
def goRgt():
  print "goRgt ..."
  nks = readImage2D(10,1,"nks")
  np0=int(nks[0][0])
  np1=int(nks[0][1])
  np2=int(nks[0][2])
  np3=int(nks[0][3])
  np4=int(nks[0][4])
  np5=int(nks[0][5])
  np6=int(nks[0][6])
  np7=int(nks[0][7])
  np8=int(nks[0][8])
  np9=int(nks[0][9])
  if not plotOnly:
    kk0 = readImage2D(np0,4,"kk0")
    kk1 = readImage2D(np1,4,"kk1")
    kk2 = readImage2D(np2,4,"kk2")
    kk3 = readImage2D(np3,4,"kk3")
    kk4 = readImage2D(np4,4,"kk4")
    kk5 = readImage2D(np5,4,"kk5")
    kk6 = readImage2D(np6,4,"kk6")
    kk7 = readImage2D(np7,4,"kk7")
    kk8 = readImage2D(np8,4,"kk8")
    kk9 = readImage2D(np9,4,"kk9")
    k1 = [kk0[0],kk1[0],kk2[0],kk3[0],kk4[0],kk5[0],kk6[0],kk7[0],kk8[0],kk9[0]]
    k2 = [kk0[1],kk1[1],kk2[1],kk3[1],kk4[1],kk5[1],kk6[1],kk7[1],kk8[1],kk9[1]]
    k3 = [kk0[2],kk1[2],kk2[2],kk3[2],kk4[2],kk5[2],kk6[2],kk7[2],kk8[2],kk9[2]]
    k4 = [kk0[3],kk1[3],kk2[3],kk3[3],kk4[3],kk5[3],kk6[3],kk7[3],kk8[3],kk9[3]]

    gx = readImage(sfile)
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    ep = readImage(epfile)
    fl = readImage(flfile)
    fp = readImage(fpfile)
    ft = readImage(ftfile)
    fs = FaultSkinner()
    cells = fs.findCells([fl,fp,ft])
    fci = FaultCellsToImage(cells)
    fli = fci.getLikelihoodThick(0.6,fl)
    wp = sub(1.0,fli)
    wp = pow(wp,6.0)
    ep = pow(ep,12.0)
    wp = mul(wp,ep)
    p2 = mul(d1/d2,p2)
    p3 = mul(d1/d3,p3)
    fl = Flattener3C()
    fl.setWeight1(0.08)
    fl.setScale(0.001)
    fl.setIterations(0.01,200)
    fl.setSmoothings(12.0,12.0,6.0)
    fm = fl.getMappingsFromSlopes(s1,s2,s3,p2,p3,wp,k4,k1,k2,k3)
    gu = fm.flatten(gx) # flattened image
    gt = fm.u1 # rgt volume
    gh = fm.x1 # horizon volume
    writeImage(gufile,gu)
    writeImage(gtfile,gt)
    writeImage(ghfile,gh)
    writeImage(wpfile,wp)
  else:
    gx = readImage(sfile)
    gu = readImage(gufile)
    gt = readImage(gtfile)
    gh = readImage(ghfile)
    wp = readImage(wpfile)
  plot3(gx,png="seismic")
  plot3(gu,png="flattenedC")
  plot3(gx,gt,cmin=min(gt),cmax=max(gt),cmap=jetRamp(1.0),
        clab="Relative geologic time",png="rgtC")
  plot3(gx,wp,cmin=0.0,cmax=1.0,cmap=jetRamp(1.0),
        clab="Weights",png="weights")

def goRgtInterp():
  print "goRgtInterp ..."
  s1,s2,s3 = Sampling(n1),Sampling(n2),Sampling(n3)
  if not plotOnly:
    p  = readImage(gfile)
    gx = readImage(sfile)
    gt = readImage(gtfile)
    ri = RgtInterp3(0.0,p)
    ri.setScales(0.001,1.0)
    ri.setRgt(gt)
    samples=ri.getPoints(0.0,p)
    ft,fp,fq = ri.grid(s1,s2,s3)
    writeImage(fpfile,fp)
    writeImage(fqfile,fq)
    writeImage(ftfile,ft)
  else:
    gx = readImage(sfile)
    fp = readImage(fpfile)
    fq = readImage(fqfile)
    ft = readImage(ftfile)
  '''
  display(gx,ft,0.0,20.0,logType)
  display(gx,fp,vmin,vmax,logType+"Wells")
  display(gx,fq,vmin,vmax,logType+"Wells")
  '''
  samples=readLogSamples(logSet,logType,smooth)
  plot3(gx,samples=samples,horizon=horizon,clab="Amplitude")
  plot3(gx,fq,samples=samples,horizon=horizon,
        cmin=vmin,cmax=vmax,cmap=jetFill(0.3),clab=logLabel)

  '''
  plot3(gx,fq,horizon=horizon,
        cmin=vmin,cmax=vmax,cmap=jetFill(0.3),clab=logLabel)
  '''

def goRgtInterpOne():
  print "goRgtInterpOne ..."
  s1,s2,s3 = Sampling(n1),Sampling(n2),Sampling(n3)
  if not plotOnly:
    gx = readImage(sfile)
    gt = readImage(gtfile)
    p  = readImage(fw1file)
    ri = RgtInterp3(0.0,p)
    ri.setScales(0.001,1.0)
    ri.setRgt(gt)
    samples=ri.getPoints(0.0,p)
    ft,fp,fq = ri.grid(s1,s2,s3)
    writeImage(fp1file,fp)
    writeImage(fq1file,fq)
    writeImage(ft1file,ft)
  else:
    gx = readImage(sfile)
    fp = readImage(fp1file)
    fq = readImage(fq1file)
    ft = readImage(ft1file)
    p  = readImage(fw1file)
  ri = RgtInterp3()
  fx,x1,x2,x3=ri.getPoints(0.0,p)
  x1=mul(x1,d1)
  x2=mul(x2,d2)
  x3=mul(x3,d3)
  x1=add(x1,f1)
  x2=add(x2,f2)
  x3=add(x3,f3)
  samples=fx,x1,x2,x3
  plot3(gx,ft,cmin=0,cmax=50,cmap=jetFill(0.3),
        clab="Time")
  plot3(gx,fp,samples=[samples],horizon=horizon,
        cmin=vmin,cmax=vmax,cmap=jetFill(0.3),clab=logLabel)

def goInterpO():
  global k1,k2,k3
  k1,k2,k3 = 366,15,96
  k1,k2,k3 = 366,15,95
  #gridBlendedP()
  #gridBlendedQ()
  s = readImage(sfile); print "s min =",min(s)," max =",max(s)
  #display1(s,True,cmin=vmin,cmax=vmax,png=logType+"Wells")
  #display1(s,False)
  #display1(s,False,["CrowMountainCRMT","TensleepASand"])
  #display1(s,True,["CrowMountainCRMT","TensleepASand"])
  #p = readImage(pfile); print "p min =",min(p)," max =",max(p)
  q = readImage(qfile); print "q min =",min(q)," max =",max(q)
  #t = readImage(tfile); print "t min =",min(t)," max =",max(t)
  #display(s,p,vmin,vmax,logType)
  #display(s,q,vmin,vmax,logType,png=logType+"Old")
  #display(s,t,0.0,100.0,logType)
  #display(s,q,vmin,vmax,logType,["CrowMountainCRMT"])
  plot3(s,q,horizon=horizon,
        cmin=vmin,cmax=vmax,cmap=jetFill(0.3),clab=logLabel)

  #display(s,q,vmin,vmax,logType,["TensleepASand"])
def goOneWell():
  p = readImage(gfile)
  s = readImage(sfile)
  q=zerofloat(n1,n2,n3)
  q[77][190]=p[77][190]
  display(s,q,vmin,vmax,logType)
  display(s,p,vmin,vmax,logType)
  writeImage(fw1file,q)


def goFigures():
  global k1,k2,k3
  k1,k2,k3 = 228,170,96 # intersect low-velocity layers
  s = readImage(sfile); print "s min =",min(s)," max =",max(s)
  #p = readImage(pfile); print "p min =",min(p)," max =",max(p)
  #q = readImage(qfile); print "q min =",min(q)," max =",max(q)
  #display3(s,None,0.0,0.0,"tpsz")
  #display3(s,p,vmin,vmax,"tppvb")
  #display3(s,q,vmin,vmax,"tpqvb")
  p = readImage("ig6/tppvb"); print "p min =",min(p)," max =",max(p)
  display3(s,p,vmin,vmax,"tppvb")
  p = readImage("tppvo09"); print "p min =",min(p)," max =",max(p)
  display3(s,p,vmin,vmax,"tppvb")

def gridBlendedP():
  e = getEigenTensors()
  bi = BlendedGridder3(e)
  p = readImage(gfile)
  t = bi.gridNearest(0.0,p)
  writeImage(pfile,p)
  writeImage(tfile,t)

def gridBlendedQ():
  e = getEigenTensors()
  bg = BlendedGridder3(e)
  bg.setSmoothness(1.0)
  p = readImage(pfile)
  t = readImage(tfile)
  t = clip(0.0,50.0,t)
  q = copy(p)
  bg.gridBlended(t,p,q)
  writeImage(qfile,q)

def getEigenTensors():
  e = readTensors(esfile)
  return e

def goImpedance():
  global k1,k2,k3,logLabel
  k1,k2,k3 = 366,15,96
  logType = None
  logLabel = "Impedance (g/cc x km/s)"
  d = readImage("tpqdb"); print "d min =",min(d)," max =",max(d)
  v = readImage("tpqvb"); print "v min =",min(v)," max =",max(v)
  s = readImage("tpsz"); print "s min =",min(s)," max =",max(s)
  i = mul(d,v)
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply2XX(mul(0.5,log(i)),i)
  #imin,imax = min(i),max(i)
  imin,imax = -0.05,0.05
  display2S(s,i,imin,imax,logType)

def display2S(s,g,cmin,cmax,logType,horizons=[]):
  world = World()
  ipg = addImageToWorld(world,s)
  ipg.setClips(smin,smax)
  ipg.setSlices(k1,k2,k3)
  ipg = addImageToWorld(world,g)
  ipg.setClips(cmin,cmax)
  ipg.setSlices(k1,k2,k3)
  if logType:
    addLogsToWorld(world,logSet,logType,cmin,cmax,smooth=smooth)
  for horizon in horizons:
    addHorizonToWorld(world,horizon)
  frame = makeFrame(world)

def display(s,g,cmin,cmax,logType,horizons=[],png=None):
  world = World()
  ipg = addImage2ToWorld(world,s,g)
  ipg.setClips1(smin,smax)
  ipg.setClips2(cmin,cmax)
  ipg.setSlices(k1,k2,k3)
  '''
  if logType:
    addLogsToWorld(world,logSet,logType,cmin,cmax,smooth=smooth)
  '''
  for horizon in horizons:
    addHorizonToWorld(world,horizon)
  frame = makeFrame(world)
  cbar = addColorBar(frame,logLabel)
  cbar.setWidthMinimum(120)
  ipg.addColorMap2Listener(cbar)
  if png and pngDir:
    frame.paintToFile(pngDir+png+".png")

def display1(s,wells=True,horizons=[],cmin=0,cmax=0,png=None):
  world = World()
  ipg = addImageToWorld(world,s)
  ipg.setClips(smin,smax)
  ipg.setSlices(k1,k2,k3)
  frame = makeFrame(world)
  if wells:
    cbar = addColorBar(frame,logLabel)
    addLogsToWorld(world,logSet,logType,cmin,cmax,cbar,smooth=smooth)
  for horizon in horizons:
    addHorizonToWorld(world,horizon)
  if png and pngDir:
    frame.paintToFile(pngDir+png+".png")

def addColorBar(frame,label):
  cbar = ColorBar(logLabel)
  cbar.setFont(cbar.getFont().deriveFont(36.0))
  frame.add(cbar,BorderLayout.EAST)
  cbar.setWidthMinimum
  #frame.viewCanvas.setBackground(frame.getBackground())
  return cbar

def display2(s,g=None,cmin=0,cmax=0):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.addColorBar()
  sp.getPlotPanel().setColorBarWidthMinimum(80)
  pv = sp.addPixels(s)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if g!=None:
    pv = sp.addPixels(g)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setColorModel(ColorMap.getJet(0.3))
    if cmin!=cmax:
      pv.setClips(cmin,cmax)

def display3(s,g=None,cmin=0,cmax=0,png=None):
  pp = PlotPanelPixels3(
    PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT,
    PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM,
    s1,s2,s3,s)
  pp.setSlices(k1,k2,k3)
  pp.setLabel1("Depth (km)")
  pp.setLabel2("Crossline (km)")
  pp.setLabel3("Inline (km)")
  pp.setClips(smin,smax)
  if g:
    pp.setLineColor(Color.BLACK)
    cb = pp.addColorBar(logLabel)
    cb.setInterval(1.0)
  else:
    pp.setLineColor(Color.YELLOW)
    cb = pp.addColorBar("Amplitude")
    cb.setInterval(5.0)
  pp.setColorBarWidthMinimum(100)
  pp.setInterval1(0.5)
  pp.setInterval2(2.0)
  pp.setInterval3(2.0)
  pp.mosaic.setHeightElastic(0,100)
  pp.mosaic.setHeightElastic(1,200)
  if g:
    pv12 = PixelsView(s1,s2,slice12(k3,g))
    pv12.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv13 = PixelsView(s1,s3,slice13(k2,g))
    pv13.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv23 = PixelsView(s2,s3,slice23(k1,g))
    pv23.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
    for pv in [pv12,pv13,pv23]:
      pv.setColorModel(ColorMap.getJet(0.5))
      if cmin!=cmax:
        pv.setClips(cmin,cmax)
    pp.pixelsView12.tile.addTiledView(pv12)
    pp.pixelsView13.tile.addTiledView(pv13)
    pp.pixelsView23.tile.addTiledView(pv23)
  pf = PlotFrame(pp)
  pf.setFontSizeForSlide(1.0,1.0)
  pf.setSize(996,815)
  pf.setVisible(True)
  if png and pngDir:
    pf.paintToPng(300,6,pngDir+png+".png")

#############################################################################

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

def plot3(f,g=None,samples=None,horizon=None, 
          cmin=None,cmax=None,cmap=None,clab=None,cint=None,png=None):
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  world = World()
  sf = SimpleFrame(world)
  cbar = None
  if g==None:
    ipg = sf.addImagePanels(s1,s2,s3,f)
    if cmap!=None:
      ipg.setColorModel(cmap)
    if cmin!=None and cmax!=None:
      ipg.setClips(cmin,cmax)
    else:
      ipg.setClips(fmin,fmax)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setColorModel1(ColorMap.getGray())
    ipg.setClips1(fmin,fmax)
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
  ipg.setSlices(k1,k2,k3)
  if horizon:
    surf = readHorizon(horizon)
    ijk = surf.getIABC();
    xyz = surf.getX321();
    if g==None:
      for k in range(2,len(xyz),3):
        xyz[k] = xyz[k]+2*d1
      tg = TriangleGroup(ijk,xyz)
      tg.setColor(Color.GRAY)
    else:
      cp = ColorMap(vmin,vmax,ColorMap.JET)
      xyz,uvw,rgb = surf.xyzUvwRgb(s1,s2,s3,cp,xyz,f,g)
      tg = TriangleGroup(ijk,xyz,uvw,rgb)
    sf.world.addChild(tg)

  if samples:
    fl,x1l,x2l,x3l = samples
    for i,f in enumerate(fl):
      f = fl[i]
      x1 = x1l[i]
      x2 = x2l[i]
      x3 = x3l[i]
      pg = makePointGroup(f,x1,x2,x3,vmin,vmax,cbar)
      sf.world.addChild(pg)
  if cbar:
    sf.setSize(837,700)
    sf.setSize(1250,900)
  else:
    sf.setSize(700,700)
    sf.setSize(1250,900)

  view = sf.getOrbitView()
  zscale = 0.75*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.3)
  #view.setAzimuth(75.0)
  #view.setAzimuth(-75.0)
  view.setAzimuth(-65.0)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  sf.viewCanvas.setBackground(sf.getBackground())
  sf.setSize(1250,900)

  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")

#############################################################################
run(main)
