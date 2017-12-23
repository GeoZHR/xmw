"""
Demonstrate velocity estimation
Author: Xinming Wu, Colorado School of Mines
Version: 2016.07.20
"""
from utils2 import *
setupForSubset("fakeNew")
s1,s2,s3= getSamplings()
n1,n2= s1.count,s2.count

# Names and descriptions of image files used below.
gxfile  = "image" # migrated image 
gffile  = "gf"
smfile  = "sm" # migrated image 
pkfile  = "pk" # picked velocity 
epfile  = "ep" # eigenvalue-derived planarity
p2file  = "p2" # inline slopes
flfile  = "fl" # fault likelihood
ftfile  = "ft" # fault dip (theta)
fltfile = "flt" # fault likelihood thinned
fttfile = "ftt" # fault dip thinned

# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
minTheta,maxTheta = 65,89
sigmaTheta = 20

# These parameters control the construction of fault skins.
# See the class FaultSkinner for more information.
lowerLikelihood = 0.1
upperLikelihood = 0.50
minSize = 10

minThrow = -35*8
maxThrow =  0

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
plotOnly = False
pngDir = None
pngDir = "../../../png/tjxd/fake/"

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goFault()
  #goScan()
  #goFaultThrow()
  #goInterp()
  #goMigVel()
  #goHorizonsB()
  goRgtInterp()
  #goRgtInterpSub()
  #goFlatten()
  #goLogVelocity()
def goHorizonsB():
  sigma1,sigma2=8.0,2.0
  pmax = 5
  m2 = 1570
  gx = readImage2D(n1,m2,"fxs")
  p2 = zerofloat(n1,m2)
  el = zerofloat(n1,m2)
  lsf = LocalSlopeFinder(sigma1,sigma2,pmax)
  lsf.findSlopes(gx,p2,el) # estimate slopes and linearity
  wp = fillfloat(1,n1,m2)
  for i1 in range(135):
    wp[370][i1] = 0.02
    wp[371][i1] = 0.05
    wp[369][i1] = 0.05
    wp[820][i1] = 0.02
    wp[821][i1] = 0.05
    wp[819][i1] = 0.05
  he2 = MultigridHorizon2()
  he2.setSmoothing(12)
  he2.setWeight(0.001)
  he2.setExternalIterations(20)
  # pick a horizon with 2 control points
  k11,k12=[42,9,1,27,47],[90,188,200,700,1180]  
  k21,k22=[50,21,47,64],[125,320,659,1180]  
  k31,k32=[75,30,59,76],[60,275,642,1250]  
  k41,k42=[ 96,48,66,97],[60,245,600,1280]  
  k51,k52=[78,94,116,114],[190,555,1090,1560]  
  k61,k62=[130,130,122],[80,485,1000]  
  k71,k72=[140,140,140],[10,620,900]  
  k1s = [k11,k21,k31,k41,k51,k61,k71]
  k2s = [k12,k22,k32,k42,k52,k62,k72]
  cv1s,cv2s,cv3s = [],[],[]
  p1s,p2s = [],[]
  np = len(k1s)
  for ip in range(np):
    k1 = k1s[ip]
    k2 = k2s[ip]
    cv1i = he2.curveInitialization(m2,n1-1,k1,k2)
    cv2i = copy(cv1i)
    cv1i = he2.curveUpdateFromSlopes(wp,p2,k1,k2,cv1i);
    cv1s.append(cv1i)
    for ik in range(len(k1)):
      p1s.append(k1[ik])
      p2s.append(k2[ik])
  hp = Helper()
  np = len(cv1s[0])
  hs = []
  for iv in range(len(cv1s)):
    hs.append(cv1s[iv][np-1])
  for i2 in range(900):
    hs[len(hs)-1][i2] = 140
  tx = hp.rgtInterpolate(n1,m2,hs)
  c2 = Sampling(m2)
  fl = readImage2D(n1,n2,flfile)
  fx = readImage2D(n1,n2,gxfile)
  txs = hp.unfaultBlocker(2,fl,tx)
  txs = hp.unfaultBlocker(1,fl,txs)
  writeImage("tx",txs)
  print len(txs)
  print len(txs[0])
  plots(s1,c2,gx,hs=cv1s,k2=p2s,k1=p1s,label="Amplitude",png="hs")
  plot2x(gx,s1,c2,g=tx,gmin=0,gmax=180,label="Relative geologic time",png="rgt")
  plot2x(txs,s1,s2,g=txs,gmin=0,gmax=180,label="Relative geologic time",png="rgt")

def goLogVelocity():
  vt = readImage2D(n1,n2,"vt")
  vm = readImage2D(n1,n2,"vm")
  vt = mul(vt,0.001)
  vm = mul(vm,0.001)
  vs = zerofloat(n1)
  rgf = RecursiveGaussianFilterP(2)
  rgf.apply0(vt[50],vs)
  vd = sub(vs,vm[50])
  plot1s(s1,vt[50],vm[50],png="logVel")
  plot1s(s1,vs,vm[50],png="logVels")
  plot1s(s1,vd,vd,vmin=-0.5,vmax=2.0,hlabel="Velocity update (km/s)",png="logVelUpdate")

def goFlatten():
  gx = readImage2D(n1,n2,gxfile)
  gx = gain(gx)
  fx = zerofloat(n1,n2)
  c2 = 50
  dw = DynamicWarping(-50,0)
  dw.setStrainMax(0.25)
  for k2 in range(50,n2,1):
    sk = dw.findShifts(gx[c2],gx[k2])
    fx[k2] = dw.applyShifts(sk,gx[k2])
  plot2x(gx,s1,s2,label="Amplitude",png="seis")
  plot2x(fx,s1,s2,label="Amplitude",png="seis")
def goHorizons():
  sigma1,sigma2=8.0,2.0
  pmax = 5
  gx = readImage2D(n1,n2,gxfile)
  gx = gain(gx)
  p2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lsf = LocalSlopeFinder(sigma1,sigma2,pmax)
  lsf.findSlopes(gx,p2,el) # estimate slopes and linearity
  fl = readImage2D(n1,n2,flfile)
  wp = sub(1,fl)
  he2 = MultigridHorizon2()
  he2.setSmoothing(12)
  he2.setWeight(0.001)
  he2.setExternalIterations(20)
  # pick a horizon with 2 control points
  k11,k12=[44,9,0,24,45],[60,188,340,380,580]  
  k21,k22=[50,21,44,61],[125,320,338,580]  
  k31,k32=[75,30,60,75],[60,275,322,580]  
  k41,k42=[ 96,48,93,96],[60,245,275,620]  
  k51,k52=[78,116],[190,620]  
  k61,k62=[137],[620]  
  k1s = [k11,k21,k31,k41,k51,k61]
  k2s = [k12,k22,k32,k42,k52,k62]
  cv1s,cv2s,cv3s = [],[],[]
  p1s,p2s = [],[]
  np = len(k1s)
  for ip in range(np):
    k1 = k1s[ip]
    k2 = k2s[ip]
    cv1i = he2.curveInitialization(n2,n1-1,k1,k2)
    cv2i = copy(cv1i)
    cv1i = he2.curveUpdateFromSlopes(18,15,60,gx,wp,p2,k1,k2,cv1i);
    cv1s.append(cv1i)
    for ik in range(len(k1)):
      p1s.append(k1[ik])
      p2s.append(k2[ik])
  hp = Helper()
  np = len(cv1s[0])
  hs = []
  for iv in range(len(cv1s)):
    hs.append(cv1s[iv][np-1])
  tx = hp.rgtInterpolate(n1,n2,hs)
  writeImage("tx",tx)
  plots(s1,s2,gx,hs=cv1s,k2=p2s,k1=p1s,label="Amplitude",png="hs")
  plot2x(gx,s1,s2,g=tx,gmin=0,gmax=180,label="Relative geologic time",png="rgt")
def goMigVel():
  gx = readImage2D(n1,n2,gxfile)
  gx = gain(gx)
  plot2x(gx,s1,s2,label="Amplitude",png="seis")
  #plot2x(vx,c1,c2,g=vx,gmin=2500,gmax=6000,label="Velocity",png="migVel")
def goFault():
  p11 = [  0, 10, 35, 49,110,125,n1-1]
  p12 = [370,354,300,250,160,125,20]
  p21 = [  0, 57, 76,115,139,n1-1]
  p22 = [470,341,266,190,131,50]
  gx = readImage2D(n1,n2,gxfile)
  gx = gain(gx)
  hpr = Helper()
  fl = zerofloat(n1,n2)
  hpr.fault(1,n1,n2,[p21],[p22],fl)
  hpr.fault(2,n1,n2,[p11],[p12],fl)
  fxs = hpr.faultBlocker(1,fl,gx)
  fxs = hpr.faultBlocker(2,fl,fxs)
  m2 = len(fxs)
  c2 = Sampling(m2)
  writeImage(flfile,fl)
  plot2(s1,s2,gx,g=fl,cmin=0.20,cmax=1,cmap=jetRamp(1.0),
      label="Fault likelihood",png="fl")
  writeImage("fxs",fxs)
  print m2
  plot2(s1,s2,gx)
  plot2(s1,c2,fxs)

def goInterp():
  vx = readImage2D(n1,n2,"vel")
  gx = readImage2D(n1,n2,gxfile)
  gx = gain(gx)
  x1,x2,fx = getWells(50,750)
  et = getTensors()
  bg = BlendedGridder2(et,fx,x1,x2)
  bg.setSmoothness(50)
  vi = bg.grid(s1,s2)
  plot2x(gx,s1,s2,f=fx,x1=x1,x2=x2,gmin=2500,gmax=6000,label="Velocity",png="seisAndLogs")
  plot2x(gx,s1,s2,g=vi,f=fx,x1=x1,x2=x2,gmin=2500,gmax=6000,label="Velocity",png="interpVel")
  plot2x(gx,s1,s2,g=vx,gmin=2500,gmax=6000,label="Velocity",png="trueVel")
  au = zerofloat(n1,n2)
  av = zerofloat(n1,n2)
  et.getEigenvalues(au,av)
  au = clip(0.003,1,au)
  et.setEigenvalues(au,av)
  e1 = Sampling(19,10,1)
  e2 = Sampling(21,40,1)
  plotTensors(gx,s1,s2,d=et,dscale=3.5,ne=20,e1=e1,e2=e2,cmin=0,cmax=0,png="tensors")

def goRgtInterp():
  tx = readImage2D(n1,n2,"tx")
  vt = readImage2D(n1,n2,"vt")
  vm = readImage2D(n1,n2,"vm")
  vt = mul(vt,0.001)
  vm = mul(vm,0.001)
  vs = zerofloat(n1)
  rgf = RecursiveGaussianFilterP(8)
  rgf.apply0(vt[50],vs)
  vd = sub(vs,vm[50])

  gx = readImage2D(n1,n2,gxfile)
  x1,x2,fx = getWell(50)
  gx = gain(gx)
  hp = Helper()
  vi = hp.interpWithRgt(50,vd,tx)
  lof = LocalOrientFilter(4,2)
  ets = lof.applyForTensors(gx)
  ets.setEigenvalues(0.02,1.0)
  plot2x(gx,s1,s2,f=vd,x1=x1,x2=x2,gmin=-0.5,gmax=2.0,label="Velocity (km/s)",png="seisAndLogs")
  plot2x(gx,s1,s2,g=vi,f=vd,x1=x1,x2=x2,gmin=-0.5,gmax=2.0,label="Velocity update (km/s)",png="interpVelS8")
  plot2x(gx,s1,s2,g=vt,gmin=2.5,gmax=6.5,label="Velocity (km/s)",png="trueVel")
  vm = readImage2D(n1,n2,"vm")
  vm = mul(0.001,vm)
  vu = add(vi,vm)
  writeImage("vu8",vu)
  plot2x(gx,s1,s2,g=vu,gmin=2.5,gmax=6.5,label="Updated migration velocity (km/s)",png="updatedVelS8")
  plot2x(gx,s1,s2,g=vm,gmin=2.5,gmax=6.5,label="Initial migration velocity (km/s)",png="initialVel")
def goRgtInterpSub():
  vt = readImage2D(n1,n2,"vt")
  vm = readImage2D(n1,n2,"vm")
  vt = mul(vt,0.001)
  vm = mul(vm,0.001)
  vs = zerofloat(n1)
  rgf = RecursiveGaussianFilterP(2)
  rgf.apply0(vt[130],vs)
  vd = sub(vs,vm[130])

  tx = readImage2D(n1,n2,"tx")
  fl = readImage2D(n1,n2,flfile)
  gx = readImage2D(n1,n2,gxfile)
  b1,e1=30,110
  x1,x2,fx = getWellSub(130,b1,e1)
  gx = gain(gx)
  hp = Helper()
  vi = hp.interpWithRgt(130,vd,tx)
  wpx = sub(1,fl)
  vd = copy(70,20,vd)
  vss = hp.maskout(130,130,b1,e1,vi,tx)
  print max(vss)
  plot2x(gx,s1,s2,f=vd,x1=x1,x2=x2,gmin=-0.5,gmax=2.0,label="Velocity (km/s)",png="seisAndLogsSub")
  plot2x(gx,s1,s2,g=vss,f=vd,x1=x1,x2=x2,gmin=-0.5,gmax=2.0,cmap=jetFillExceptMin(1.0),
          label="Velocity update (km/s)",png="interpVelSub2")
  plot2x(gx,s1,s2,g=vt,gmin=2.5,gmax=6.5,label="Velocity (km/s)",png="trueVel")
  vm = readImage2D(n1,n2,"vm")
  vm = mul(0.001,vm)
  vu = add(vi,vm)
  vus = hp.maskout(130,400,b1,e1,vu,tx)
  plot2x(gx,s1,s2,g=vus,gmin=2.5,gmax=6.5,cmap=jetFillExceptMin(1.0),label="Updated migration velocity (km/s)",png="updatedVelSub2")
  plot2x(gx,s1,s2,g=vm,gmin=2.5,gmax=6.5,label="Initial migration velocity (km/s)",png="initialVel")

def getWell(k2):
  vx = readImage2D(n1,n2,"vt")
  x1 = zerofloat(n1)
  x2 = zerofloat(n1)
  fx = zerofloat(n1)
  k = 0
  for i1 in range(n1):
    x1[k] = i1
    x2[k] = k2
    fx[k] = vx[k2][i1]
    k = k+1
  return x1,x2,fx

def getWellSub(k2,b1,e1):
  vx = readImage2D(n1,n2,"vt")
  x1 = zerofloat(e1)
  x2 = zerofloat(e1)
  fx = zerofloat(e1)
  k = 0
  for i1 in range(b1,e1,1):
    x1[k] = i1
    x2[k] = k2
    fx[k] = vx[k2][i1]
    k = k+1
  return x1,x2,fx

def getWells(k1,k2):
  vx = readImage2D(n1,n2,"vt")
  x1 = zerofloat(n1*2)
  x2 = zerofloat(n1*2)
  fx = zerofloat(n1*2)
  k = 0
  for i1 in range(n1):
    x1[k] = i1
    x2[k] = k1
    fx[k] = vx[k1][i1]
    k = k+1
    x1[k] = i1
    x2[k] = k2
    fx[k] = vx[k2][i1]
    k = k+1
  return x1,x2,fx
  

def getTensors():
  gx = readImage2D(n1,n2,gxfile)
  fl = readImage2D(n1,n2,flfile)
  gx = gain(gx)
  fl = sub(1,fl)
  av = clip(0.00001,1.0,fl)
  au = fillfloat(0.0001,n1,n2)
  lof = LocalOrientFilter(3,2)
  ets = lof.applyForTensors(gx)
  ets.setEigenvalues(au,av)
  return ets

def goScan():
  print "goScan ..."
  gx = readImage2D(n1,n2,gxfile)
  p11 = [  0, 10*8, 42*8, 45*8,125*8]
  p12 = [370,354,285,270,125]
  p21 = [  2*8, 57*8, 79*8,115*8,129*8]
  p22 = [463,346,260,192,152]
  gx = readImage2D(n1,n2,gxfile)
  gx = gain(gx)
  hpr = Helper()
  fxs = hpr.fault(n1*8,n2,[p11,p21],[p12,p22])
  gxi = zerofloat(n1*8,n2)
  s1i = Sampling(n1*8,0.125,0)
  si = SincInterpolator()
  s1p = Sampling(n1*8,1,0)
  for i2 in range(n2):
    si.interpolate(s1,gx[i2],s1i,gxi[i2])
  writeImage("gxi",gxi)
  if not plotOnly:
    fs = FaultScanner2(sigmaTheta)
    fl,ft = fs.scan(minTheta,maxTheta,fxs)
    print "fl min =",min(fl)," max =",max(fl)
    print "ft min =",min(ft)," max =",max(ft)
    writeImage(flfile,fl)
    writeImage(ftfile,ft)
  else:
    fl = readImage2D(n1,n2,flfile)
    ft = readImage2D(n1,n2,ftfile)
  fl = sub(fl,min(fl))
  fl = div(fl,max(fl))
  plot2(s1p,s2,gxi,g=fl,cmin=0.10,cmax=1,cmap=jetRamp(1.0),
      label="Fault likelihood",png="fl")
  '''
  plot2(s1,s2,gx,g=abs(ft),cmin=minTheta,cmax=maxTheta,cmap=jetFill(1.0),
      label="Fault dip (degrees)",png="ft")
  '''
def goThin():
  print "goThin ..."
  gx = readImage2D(n1,n2,gxfile)
  gx = gain(gx)
  if not plotOnly:
    fl = readImage2D(n1,n2,flfile)
    ft = readImage2D(n1,n2,ftfile)
    fs = FaultScanner2(sigmaTheta)
    flt,ftt = fs.thin([fl,ft])
    writeImage(fltfile,flt)
    writeImage(fttfile,ftt)
  else:
    flt = readImage2D(n1,n2,fltfile)
    ftt = readImage2D(n1,n2,fttfile)
  plot2(s1,s2,gx)
  plot2(s1,s2,gx,g=flt,cmin=0.2,cmax=1,cmap=jetFillExceptMin(1.0))
  plot2(s1,s2,gx,g=abs(ftt),cmin=minTheta,cmax=maxTheta,cmap=jetFill(1.0),
      label="Fault dip (degrees)",png="ft")

def goSmooth():
  gx = readImage2D(n1,n2,gxfile)
  gx = gain(gx)
  fl = readImage2D(n1,n2,flfile)
  fl = sub(fl,min(fl))
  fl = div(fl,max(fl))
  ft = readImage2D(n1,n2,ftfile)
  fc = FaultCurver()
  fc.setMinCurveSize(minSize)
  fc.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  ps = fc.findPoints([fl,ft])
  cc = fc.findCurves(ps)
  ftt = zerofloat(n1,n2)
  flt = zerofloat(n1,n2)
  FaultCurve.getFlsImage(cc,flt)
  wp = sub(1,flt);
  wp = pow(wp,10)
  gs = zerofloat(n1,n2)
  lof = LocalOrientFilterP(8,2)
  ets = lof.applyForTensors(gx)
  ets.setEigenvalues(0.001,1.0)
  lsf = LocalSmoothingFilter()
  lsf.apply(ets,400,wp,gx,gs)
  plot2(s1,s2,gx)
  plot2(s1,s2,gs)
  plot2(s1,s2,gx,g=flt,cmin=0.1,cmax=1,cmap=jetFillExceptMin(1.0))
  return gs
def goFaultCurve():
  gx = readImage2D(n1,n2,gxfile)
  gx = gain(gx)
  fl = readImage2D(n1,n2,flfile)
  ft = readImage2D(n1,n2,ftfile)
  fc = FaultCurver()
  fc.setMinCurveSize(minSize)
  fc.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  ps = fc.findPoints([fl,ft])
  print len(ps)
  cs = fc.findCurves(ps)
  ftt = zerofloat(n1,n2)
  flt = fillfloat(-1000,n1,n2)
  FaultCurve.getFtImage(cs,ftt)
  FaultCurve.getFlImage(cs,flt)
  writeImage(fltfile,flt)
  plot2(s1,s2,gx,g=flt,cmin=0.1,cmax=1,cmap=jetFillExceptMin(1.0))

def goFaultThrow():
  m1 = n1*8
  gx = readImage2D(m1,n2,"gxi")
  fl = readImage2D(m1,n2,flfile)
  ft = readImage2D(m1,n2,ftfile)
  fc = FaultCurver()
  fc.setMinCurveSize(minSize)
  fc.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  ps = fc.findPoints([fl,ft])
  cc = fc.findCurves(ps)
  cs = [cc[1],cc[0]]
  ftt = zerofloat(m1,n2)
  flt = zerofloat(m1,n2)
  FaultCurve.getFlsImage(cc,flt)
  s1i = Sampling(m1)
  plot2(s1i,s2,gx,g=flt,cmin=0.1,cmax=1,cmap=jetFillExceptMin(1.0))
  lsp = LocalSlopeFinder(64,2,5)
  el = zerofloat(m1,n2)
  p2 = zerofloat(m1,n2)
  lsp.findSlopes(gx,p2,el)
  fcr = FaultCorrelater(gx,p2)
  fcr.setZeroSlope(False) # True only if we want to show the error
  fcr.setOffset(2)
  fcr.computeThrow([cs[0]],minThrow,maxThrow)
  '''
  fcr.setOffset(4)
  fcr.computeThrow([cs[1]],minThrow,maxThrow)
  '''
  fst = zerofloat(m1,n2)
  FaultCurve.getFsImage(cs,fst)
  fs1 = fillfloat(-1000,m1,n2)
  fs2 = fillfloat(-1000,m1,n2)
  FaultCurve.getFsImage(cs,fs1,fs2)
  print min(fst)
  print max(fst)
  smark = -999.999
  gw = copy(gx)
  plot2(s1i,s2,gx,label="Amplitude",png="gx")
  for ci in cs:
    p1,p2 = fcr.getDipSlips(m1,n2,[cs[0]],smark)
    p1,p2 = fcr.interpolateDipSlips([p1,p2],smark)
    gw = fcr.unfault([p1,p2],gw)
  plot2(s1i,s2,gw,label="Amplitude",png="gw")
  return cs

def goTensors():
  gx = readImage2D(n1,n2,gxfile)
  fl = readImage2D(n1,n2,fltfile)
  lof = LocalOrientFilterP(4,2)
  et = lof.applyForTensors(gx)
  eu = zerofloat(n1,n2)
  ev = zerofloat(n1,n2)
  et.getEigenvalues(eu,ev)
  eu=clip(0.008,max(eu),eu)
  ev=clip(0.008,max(ev),ev)
  et.setEigenvalues(eu,ev)
  et.invertStructure(1.0,1.0)
  et.getEigenvalues(eu,ev)
  print min(eu)
  print max(eu)
  #plotTensors(gx,s1,s2,d=et,dscale=20,mk=mk,cmin=-2,cmax=2,png="tensors")
  plot2(s1,s2,gx,cmin=-2,cmax=2,label="Amplitdue",png="seis")
  plotTensors(gx,s1,s2,d=et,dscale=1,ne=25,cmin=-2,cmax=2,png="tensors")

def goUnfault():
  cc = goFaultThrow()
  cc = [cc[0]]
  fl = zerofloat(n1,n2)
  fw = zerofloat(n1,n2)
  wp = fillfloat(1.0,n1,n2)
  gx = readImage2D(n1,n2,gxfile)
  gx = gain(gx)
  lof = LocalOrientFilterP(8,2)
  et = lof.applyForTensors(gx)
  et.setEigenvalues(0.001,1.0)
  slp = FaultSlipConstraints2(cc)
  sp = slp.screenPointsX(wp)
  uf = UnfaultS2(10,10)
  uf.setIters(100)
  uf.setTensors(et)
  np =  len(sp[0][0])
  scale = (n1*n2/np)
  mul(sp[3][0],scale,sp[3][0])
  [t1,t2] = uf.findShifts(sp,wp)
  #[t1,t2] = uf.convertShifts(40,[t1,t2])
  uf.applyShifts([t1,t2],gx,fw)
  plot2(s1,s2,gx)
  plot2(s1,s2,fw)
  plot2(s1,s2,wp)


def goVelocity():
  cc = goFaultThrow()
  fl = zerofloat(n1,n2)
  wp = fillfloat(1.0,n1,n2)
  slp = FaultSlipConstraints2(cc)
  sp = slp.screenPointsX(wp)
  #sp = None
  FaultCurve.getFlsImage(cc,fl)
  gx = readImage2D(n1,n2,gxfile)
  sm = readImage2D(n1,n2,smfile)
  pk = readImage2D(n1,n2,pkfile)
  lof = LocalOrientFilterP(4,2)
  et = lof.applyForTensors(gx)
  eu = zerofloat(n1,n2)
  ev = zerofloat(n1,n2)
  et.getEigenvalues(eu,ev)
  eu=clip(0.005,max(eu),eu)
  ev=clip(0.005,max(ev),ev)
  et.setEigenvalues(eu,ev)
  et.invertStructure(1.0,1.0)
  sm = clip(0.0,1.0,sm)
  sm = pow(sm,2)
  sm = mul(sm,wp)
  ve = VelocityEstimator(10,10)
  ve.setIterations(0.00001,400)
  ve.setTensors(et)
  ve.setSmoothness(200.0)
  ps = pow(pk,2)
  vi = ve.applyForVelocity(sp,sm,wp,ps)
  vi = sqrt(vi)
  vp = ve.predictVelocity(vi)
  plot2(s1,s2,wp)
  clab = "Velocity (km/s)"
  clab1 = "Picked migration velocity (km/s)"
  clab2 = "Predicted migration velocity (km/s)"
  clab3 = "Interval velocity (km/s)"
  plot2(s1,s2,gx,g=pk,cmin=1.4,cmax=2.3,cmap=jetFill(1.0),label=clab1,png="pk")
  plot2(s1,s2,gx,g=vp,cmin=1.4,cmax=2.3,cmap=jetFill(1.0),label=clab2,png="vp")
  plot2(s1,s2,gx,g=vi,cmin=1.6,cmax=2.7,cmap=jetFill(1.0),label=clab3,png="vi")
  plot2(s1,s2,gx,g=vi,cmin=1.6,cmax=2.7,cmap=jetFill(0.6),label=clab3,png="seisvi")

def goShapping():
  dix = readImage2D(n1,n2,"dix")
  shp = readImage2D(n1,n2,"shp")
  shpp = readImage2D(n1,n2,"shpp")
  gx = readImage2D(n1,n2,gxfile)
  clab1 = "Picked migration velocity (km/s)"
  clab2 = "Predicted migration velocity (km/s)"
  clab3 = "Interval velocity (km/s)"
  plot2(s1,s2,gx,g=dix,cmin=1.6,cmax=2.7,cmap=jetFill(1.0),label=clab3,png="dix")
  plot2(s1,s2,gx,g=shp,cmin=1.6,cmax=2.7,cmap=jetFill(1.0),label=clab3,png="shp")
  plot2(s1,s2,gx,g=shp,cmin=1.6,cmax=2.7,cmap=jetFill(0.6),label=clab3,png="seishp")
  plot2(s1,s2,gx,g=shpp,cmin=1.4,cmax=2.3,cmap=jetFill(1.0),label=clab2,png="shpp")
def like(x):
  n2 = len(x)
  n1 = len(x[0])
  return zerofloat(n1,n2)



def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(80.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2)
  div(x,sqrt(g),y)
  return y

#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET
backgroundColor = Color(0xfd,0xfe,0xff) # easy to make transparent
def jetFill(alpha):
  return ColorMap.setAlpha(ColorMap.JET,alpha)

def jetFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.JET,a)
def bwrNotch(alpha):
  a = zerofloat(256)
  for i in range(len(a)):
    if i<128:
      a[i] = alpha*(128.0-i)/128.0
    else:
      a[i] = alpha*(i-127.0)/128.0
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,a)


def bwrFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,a)

def jetRamp(alpha):
  return ColorMap.setAlpha(ColorMap.JET,rampfloat(0.0,alpha/256,256))

def bwrRamp(alpha):
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,rampfloat(0.0,alpha/256,256))

def grayRamp(alpha):
  return ColorMap.setAlpha(ColorMap.GRAY,rampfloat(0.0,alpha/256,256))

def plot1s(s1,vt,vm,vmin=2.5,vmax=6.5,hlabel="Velocity (km/s)", vlabel="Depth (sample)",png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setVLimits(0,187)
  sp.setHLimits(vmin,vmax)
  sp.setHInterval(1)
  sp.setVInterval(50)
  pv1 = sp.addPoints(s1,vt)
  pv1.setLineColor(Color.RED)
  pv1.setLineWidth(3.0)
  pv2 = sp.addPoints(s1,vm)
  pv2.setLineColor(Color.BLACK)
  pv2.setLineWidth(3.0)
  sp.setSize(300,500)
  sp.setHLabel(hlabel)
  sp.setVLabel(vlabel)
  sp.setFontSize(18) #for print
  if png and pngDir:
    sp.paintToPng(720,3.3,pngDir+png+".png")


def plots(s1,s2,fx,hs=None,k2=None,k1=None,v1=None,v2=None,
         w1=1000,w2=500,cmap=ColorMap.GRAY,cmin=0,cmax=0,label=None,png=None):
  np = len(hs[0])
  m2 = len(fx)
  for ip in range(np-1,np):
    sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
    pv = sp.addPixels(s1,s2,fx)
    sp.setHLimits(0,m2-1)
    sp.setVLimits(0,n1-1)
    sp.setHLabel("Lateral position (sample)")
    sp.setVLabel("Depth (sample)")
    if label:
      sp.addColorBar(label)
    else:
      sp.addColorBar()
    pv.setColorModel(cmap)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setClips(-2,2)
    for ik in range(len(hs)):
      hi = hs[ik][ip]
      x = rampfloat(0,1,m2)
      pvh = sp.addPoints(hi,x)
      pvh.setLineColor(Color.YELLOW)
      pvh.setLineWidth(3.0)
      sp.add(pvh)
    '''
    for ik in range(len(k2)):
      pv = PointsView([k1[ik]],[k2[ik]])
      pv.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT)
      pv.setLineStyle(PointsView.Line.NONE)
      pv.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE)
      pv.setMarkColor(Color.RED)
      pv.setMarkSize(8)
      pv.setLineWidth(4)
      sp.add(pv)
    '''
    sp.setSize(900,500)
    sp.setFontSize(18)
    sp.plotPanel.setColorBarWidthMinimum(80)
    if pngDir and png:
      sp.paintToPng(720,2.2222,pngDir+png+str(ip)+".png")

def plotTensors(g,s1,s2,d=None,dscale=1,ne=20,mk=None,e1=None,e2=None,cmin=0,cmax=0,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setBackground(backgroundColor)
  sp.setHLabel("Lateral position (sample)")
  sp.setVLabel("Depth (sample)")
  sp.setHInterval(100)
  sp.setVInterval(50)
  sp.setFontSize(24)
  #sp.setFontSizeForPrint(8,240)
  #sp.setFontSizeForSlide(1.0,0.9)
  sp.setSize(793,500)
  pv = sp.addPixels(s1,s2,g)
  pv.setColorModel(ColorMap.GRAY)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setClips(-1,1)
  if d:
    tv = TensorsView(s1,s2,d)
    tv.setOrientation(TensorsView.Orientation.X1DOWN_X2RIGHT)
    tv.setLineColor(Color.YELLOW)
    tv.setLineWidth(2)
    if(e1 and e2):
      tv.setEllipsesDisplayed(e1,e2)
    else:
      tv.setEllipsesDisplayed(ne)
    tv.setScale(dscale)
    tile = sp.plotPanel.getTile(0,0)
    tile.addTiledView(tv)
  sp.setFontSize(18)
  if pngDir and png:
    sp.paintToPng(720,3.3,pngDir+png+".png")
    #sp.paintToPng(720,3.3,pngDir+png+".png")


def plot2x(s,s1,s2,g=None,f=None,x1=None,x2=None,gmin=None,gmax=None,cmap=jetFill(1.0),label=None,png=None,et=None):
  n2 = len(s)
  n1 = len(s[0])
  s1,s2=Sampling(n1),Sampling(n2)
  panel = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.LEFT_TOP)
  panel.setVLimits(0,n1-1)
  panel.setHInterval(100)
  panel.setVInterval(50)
  panel.setHLabel("Lateral position (sample)")
  panel.setVLabel("Depth (sample)")
  #panel.setHInterval(100.0)
  #panel.setVInterval(100.0)
  #panel.setHLabel("Pixel")
  #panel.setVLabel("Pixel")
  if label:
    panel.addColorBar(label)
  else:
    panel.addColorBar()
  panel.setColorBarWidthMinimum(80)
  pv = panel.addPixels(s1,s2,s)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setColorModel(ColorMap.GRAY)
  pv.setClips(-1.0,1.0)
  if g:
    alpha = 1.0
  if g:
    pv = panel.addPixels(s1,s2,g)
    pv.setInterpolation(PixelsView.Interpolation.LINEAR)
    pv.setColorModel(cmap)
    if gmin:
      pv.setClips(gmin,gmax)
  if x1 and x2:
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
  #frame.setFontSizeForPrint(8,240)
  #frame.setSize(1240,774)
  #frame.setFontSizeForSlide(1.0,0.8)
  frame.setSize(900,500)
  frame.setVisible(True)
  frame.setFontSize(18)
  if png and pngDir:
    frame.paintToPng(720,3.3,pngDir+png+".png")
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

def plot2(s1,s2,f,g=None,cmin=None,cmax=None,cmap1=ColorMap.GRAY,cmap=None,label=None,png=None):
  n2 = len(f)
  n1 = len(f[0])
  f1,f2 = s1.getFirst(),s2.getFirst()
  d1,d2 = s1.getDelta(),s2.getDelta()
  panel = panel2Teapot()
  panel.setHInterval(50.0)
  panel.setVInterval(50.0)
  panel.setHLabel("Lateral position (km)")
  panel.setVLabel("Time (s)")
  #panel.setHInterval(100.0)
  #panel.setVInterval(100.0)
  #panel.setHLabel("Pixel")
  #panel.setVLabel("Pixel")
  if label:
    panel.addColorBar(label)
  else:
    panel.addColorBar()
  panel.setColorBarWidthMinimum(80)
  pv = panel.addPixels(s1,s2,f)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setColorModel(ColorMap.GRAY)
  pv.setClips(-2,2)
  if g:
    pv = panel.addPixels(s1,s2,g)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setColorModel(cmap)
    if label:
      panel.addColorBar(label)
    else:
      panel.addColorBar()
  if cmin and cmax:
    pv.setClips(cmin,cmax)
  frame2Teapot(panel,png)
def panel2Teapot():
  panel = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT)#,PlotPanel.AxesPlacement.NONE)
  return panel
def frame2Teapot(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  #frame.setFontSizeForPrint(8,240)
  #frame.setSize(1240,774)
  #frame.setFontSizeForSlide(1.0,0.9)
  frame.setFontSize(12)
  frame.setSize(660,800)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(400,3.2,pngDir+png+".png")
  return frame


#############################################################################
run(main)
