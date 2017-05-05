"""
Demo of generating avo models
Author: Xinming Wu, University of Texas at Austin
Version: 2017.05.05
"""

from utils import *
setupForSubset("fake")
#setupForSubset("tp")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

# Names and descriptions of image files used below.
logFile2 = "BZ35-2-2_logs_E3"
logFile3 = "BZ35-2-3_logs_E3"
logFile5 = "BZ35-2-5_logs_E3"
logFile6 = "BZ35-2-6_logs_E3"
rvFile2 = "rv2"
rvFile3 = "rv3"
rvFile5 = "rv5"
rvFile6 = "rv6"

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
maxThrow = 20.0

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
pngDir = None
pngDir = "../../../png/avo/fake/"
plotOnly = True
# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  goModels()
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
  la = hp.sortLogs(rv2,rv5,rv6,rv3)
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

def goImpedanceHorizon():
  pc = readImage(pcfile)
  hz = readImage2D(n2,n3,hzfile)
  hz = add(hz,-1)
  ps = readImage(pxfile+str(0.9)+str(0.0)+str(1.0))
  pk = readImage(pxfile+str(0.9)+str(0.0)+str(0.001))
  pmin = 6000
  pmax = 18000
  cmap = ColorMap.JET
  plot3c(pc,cmin=pmin,cmax=pmax,cmap=cmap,hz=hz,mk=1,png="pchz")
  plot3c(ps,cmin=pmin,cmax=pmax,cmap=cmap,hz=hz,mk=1,png="pshz")
  plot3c(pk,cmin=pmin,cmax=pmax,cmap=cmap,hz=hz,mk=1,png="pkhz")
def goOrientations():
  sk = readSkins(fslbase)
  gc = readImage(gcfile)
  gn = readImage(gncfile)
  hz = readImage2D(n2,n3,hzfile)
  sk = readSkins(fslbase)
  lof = LocalOrientFilter(2,3)
  et = lof.applyForTensors(gn)
  et.invertStructure(0.0,1.0,1.0)
  '''
  plot3(gn,cmin=-2,cmax=1.5,clab="Amplitude",png="gn")
  plot3(gn,skins=sk,smax=15.0,clab="Fault throw (samples)",png="gnskin")
  plot3(gn,cmin=0,cmax=15,cmap=ColorMap.JET,
        clab="Fault throw (samples)",png="throwBar")
  plot3(gn,cmin=-2,cmax=1.5,clab="Amplitude",et=et,size=15,png="gnet")
  plot3c(gn,cmin=-2,cmax=1.5,clab="Amplitude",hz=hz,mk=-1,png="gnhz")
  '''
  plot3c(gn,cmin=-2,cmax=1.5,hz=hz,mk=-1,et=et,png="gnhzw")
  plot3c(gn,cmin=-2,cmax=1.5,hz=hz,mk=-1,et=et,skins=sk,w=False,png="gnhzv")
def goHorizon():
  hz = readImage2D(n2,n3,hzfile)
  gc = readImage(gnfile)
  gn = readImage(gnfile)
  d1 = .05
  srd = SurfaceRefinerDp()
  srd.setStrainMax(1.0,1.0)
  srd.setErrorSmoothing(3)
  gmin,gmax,gmap = -2.0,1.5,ColorMap.GRAY
  plot3(gc,hz=hz,cmin=gmin,cmax=gmax,cmap=gmap,clab="Amplitude",png="gnhz")
  b1,e1=-40,20
  em = srd.getErrorMatrix(b1,e1,d1,gc,hz)
  sfr = zerofloat(n2,n3)
  srd.findSurface(em,sfr)
  sfr = mul(sfr,d1)
  sfr = add(sfr,b1*d1)
  add(hz,sfr,hz)

  b1,e1=-40,0
  em = srd.getErrorMatrix(b1,e1,d1,gc,hz)
  sfr = zerofloat(n2,n3)
  srd.findSurface(em,sfr)
  sfr = mul(sfr,d1)
  sfr = add(sfr,b1*d1)
  add(hz,sfr,hz)

  b1,e1=-40,0
  em = srd.getErrorMatrix(b1,e1,d1,gc,hz)
  sfr = zerofloat(n2,n3)
  srd.findSurface(em,sfr)
  sfr = mul(sfr,d1)
  sfr = add(sfr,b1*d1)
  add(hz,sfr,hz)

  b1,e1=-40,0
  em = srd.getErrorMatrix(b1,e1,d1,gc,hz)
  sfr = zerofloat(n2,n3)
  srd.findSurface(em,sfr)
  sfr = mul(sfr,d1)
  sfr = add(sfr,b1*d1)
  add(hz,sfr,hz)

  writeImage(hzfile,hz)
  plot3(gc,hz=hz,cmin=gmin,cmax=gmax,cmap=gmap,clab="Amplitude",png="gnhz")
  plot3(gn,hz=hz,cmin=gmin,cmax=gmax,cmap=gmap,clab="Amplitude",png="gnhz")

def goInitial():
  k1,k2,k3,fp = getImpLogs()
  ai3 = AcousticImpedanceInv3(6.0,6.0)
  pt = zerofloat(n1,n2,n3)
  ai3.setInitial(pt,k1,k2,k3,fp)
  print min(pt)
  print max(pt)
  samples=fp,k1,k2,k3
  plot3(pt,cmin=6000,cmax=18000,cmap=ColorMap.JET,clab="Impedance",
      samples=samples,png="pInitial")

def goFaults():
  fs = readSkins(fskgood) 
  fh = FaultHelper()
  fl = zerofloat(n1,n2,n3)
  fh.getFlOnFaults(fs,fl)
  gx = readImage(gnfile)
  plot3(gx,fl,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0), 
          clab="Fault likelihoods", png="seisFl")

  
def goImpedance3(smooth,wps,lambda2):
  gx = readImage(gcfile)
  pc = readImage(pcfile)
  rx = readImage(rnfile)
  k1,k2,k3,fp = getImpLogs()
  fc = []
  for fpi in fp:
    fc.append(0.5*log(fpi))
  if not plotOnly:
    fh = FaultHelper()
    fk = readSkins(fskgood)
    wp = fillfloat(1,n1,n2,n3)
    fh.setValuesOnFaults(wps,fk,wp)
    lof = LocalOrientFilter(1.0,1.0)
    et3 = lof.applyForTensors(gx)
    et3.setEigenvalues(0.000001,lambda2,1.0)
    ai3 = AcousticImpedanceInv3(6.0,6.0)
    ai3.setIterations(0.001,1000)
    ai3.setTensors(et3)
    ai3.setSmoothness(smooth)
    pt = zerofloat(n1,n2,n3)
    ai3.setInitial(pt,k1,k2,k3,fc)
    px = ai3.applyForImpedance(copy(pt),rx,wp,k1,k2,k3,fc)
    px = mul(px,2.0)
    px = exp(px)
    writeImage(pxfile+str(smooth)+str(wps)+str(lambda2),px)
  else:
    px = readImage(pxfile+str(smooth)+str(wps)+str(lambda2))
  samples = fp,k1,k2,k3
  pmin = 6000
  pmax = 18000
  plot3(gx,clab="Amplitude",png="seismic")
  plot3(rx,cmin=-0.15,cmax=0.15,clab="Reflectivity",png="refx")
  plot3(pc,cmin=pmin,cmax=pmax,cmap=ColorMap.JET,clab="Impedance",png="pTrue")
  plot3(px,cmin=pmin,cmax=pmax,cmap=ColorMap.JET,clab="Impedance",
        samples=samples,png="px"+str(smooth)+str(wps)+str(lambda2))

def goImpFlatten(smooth,wps,lambda2):
  pc = readImage(pcfile)
  px = readImage(pxfile+str(smooth)+str(wps)+str(lambda2))
  t1 = readImage(sw1file)
  t2 = readImage(sw2file)
  t3 = readImage(sw3file)
  gt = readImage(rgtfile)
  uf = UnfaultS(8.0,4.0)
  pxw = zerofloat(n1,n2,n3)
  pcw = zerofloat(n1,n2,n3)
  uf.applyShifts([t1,t2,t3],px,pxw)
  uf.applyShifts([t1,t2,t3],pc,pcw)
  fl3 = Flattener3Unc()
  pxu = fl3.flatten(s1,s1,gt,pxw)
  pcu = fl3.flatten(s1,s1,gt,pcw)
  pmin = 6000
  pmax = 18000
  plot3(px,cmin=pmin,cmax=pmax,cmap=ColorMap.JET,clab="Impedance")
  plot3(pcu,cmin=pmin,cmax=pmax,cmap=ColorMap.JET,clab="Impedance",
         png="pcu")
  plot3(pxu,cmin=pmin,cmax=pmax,cmap=ColorMap.JET,clab="Impedance",
         png="pxu"+str(smooth)+str(wps)+str(lambda2))


def goImpedance2():
  gx = readImage(gxfile)
  px = readImage(pxfile)
  rx = readImage(rnfile)
  #rx = readImage(rxfile)
  #smooth = 0.5 # for clean data
  rx = readImage(rnfile)
  smooth = 0.9 # for noisy data
  r3,p3,g3 = rx[50],px[50],gx[50]
  k1,k2,fp = setWells(p3)
  ep = fillfloat(1.0,n1,n2)
  u1 = fillfloat(1.0,n1,n2)
  u2 = fillfloat(1.0,n1,n2)
  lof = LocalOrientFilter(8.0,2.0)
  et2 = lof.applyForTensors(g3)
  lof.applyForNormalLinear(r3,u1,u2,ep)
  wp = pow(ep,6.0)
  et2.setEigenvalues(0.000001,1.0)
  ai2 = AcousticImpedanceInv2(6.0,6.0)
  ai2.setIterations(0.001,800)
  ai2.setTensors(et2)
  ai2.setSmoothness(smooth)
  pt = zerofloat(n1,n2)
  ai2.setInitial(u1,pt,k1,k2,fp)
  pi2 = ai2.applyForImpedance(r3,wp,k1,k2,fp)
  plot2(s1,s2,r3,clab="Reflectivity",cmin=-0.05,cmax=0.05,png="rx")
  plot2(s1,s2,p3,cmap=ColorMap.JET,clab="True impedance",
        cmin=min(p3),cmax=max(p3),png="p3")
  plot2(s1,s2,pt,cmap=ColorMap.JET,clab="Initial impedance",
        cmin=min(p3),cmax=max(p3),interp=False,png="pt")
  plot2(s1,s2,pi2,cmap=ColorMap.JET,clab="Recovered impedance",
        cmin=min(p3),cmax=max(p3),png="pi")

def getImpLogs():
  k1u = [ 95,105,105,105,105, 21, 21]
  k2u = [138,127, 81, 56, 20, 70,124]
  k3u = [100, 90, 90, 93, 65, 30, 10]
  np = len(k1u)
  r2 = readImage(sw2file)
  r3 = readImage(sw3file)
  pc = readImage(pcfile)
  x1,x2,x3,fx=[],[],[],[]
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
      fx.append(pc[k3][k2][k1])
  return x1,x2,x3,fx

def getLogs(px):
  k2u = [81,56,20,8,96,50]
  k3u = [90,93,65,8,20,30]
  np = len(k2u)
  k1,k2,k3,fp = [],[],[],[]
  for ip in range(np):
    i2 = round(k2u[ip])
    i3 = round(k3u[ip])
    for i1 in range(n1):
      k1.append(i1)
      k2.append(i2)
      k3.append(i3)
      fp.append(px[i3][i2][i1])
  return k1,k2,k3,fp

def setWells3D(p3):
  k1 = zerofloat(n1*2)
  k2 = zerofloat(n1*2)
  k3 = zerofloat(n1*2)
  fp = zerofloat(n1*2)
  k = 0
  for i1 in range(n1):
    k1[k] = i1
    k2[k] = 15
    fp[k] = p3[15][i1]
    k = k+1
  for i1 in range(n1):
    k1[k] = i1
    k2[k] = 61
    fp[k] = p3[61][i1]
    k = k+1
  return k1,k2,fp

def setWells(p3):
  k1 = zerofloat(n1*2)
  k2 = zerofloat(n1*2)
  fp = zerofloat(n1*2)
  k = 0
  for i1 in range(n1):
    k1[k] = i1
    k2[k] = 15
    fp[k] = p3[15][i1]
    k = k+1
  for i1 in range(n1):
    k1[k] = i1
    k2[k] = 61
    fp[k] = p3[61][i1]
    k = k+1
  return k1,k2,fp

def goFakeData():
  sequence = 'OA' # 1 episode of folding, followed by one episode of faulting
  nplanar = 3 # number of planar faults
  conjugate = False # if True, two large planar faults will intersect
  conical = False # if True, may want to set nplanar to 0 (or not!)
  impedance = False # if True, data = impedance model
  wavelet = True # if False, no wavelet will be used
  lateralViriation = False
  noise = 0.6 # (rms noise)/(rms signal) ratio
  if not plotOnly:
    gc,gn,rn,pc,hs = FakeData.seismicAndSlopes3d2015A(sequence,
    nplanar,conjugate,conical,impedance,wavelet,lateralViriation,noise)
    #hz = hs[0]
    #writeImage(hzfile,hz)
    writeImage(gcfile,gc)
    writeImage(gnfile,gn)
    writeImage(rnfile,rn)
    writeImage(pcfile,pc)
  else:
    gc = readImage(gcfile)
    gn = readImage(gnfile)
    rn = readImage(rnfile)
    pc = readImage(pcfile)
    #hz = readImage2D(n2,n3,hzfile)
  dmin = 2.1
  dmax = 3.2
  rmin = -0.15
  rmax = 0.15
  pmin = 6000
  pmax = 18000

  pmap = ColorMap.JET
  gmin,gmax,gmap = -2.0,1.5,ColorMap.GRAY
  if impedance:
    gmin,gmax,gmap = 0.0,1.4,ColorMap.JET
  plot3(gc,cmin=gmin,cmax=gmax,cmap=gmap,clab="Amplitude",png="gc")
  plot3(gn,cmin=gmin,cmax=gmax,cmap=gmap,clab="Amplitude",png="gn")
  plot3(rn,cmin=rmin,cmax=rmax,cmap=gmap,clab="Reflectivity",png="rn")
  plot3(pc,cmin=pmin,cmax=pmax,cmap=pmap,clab="Impedance",png="pc")


def goSlip():
  print "goSlip ..."
  gsx = readImage(gcfile)
  p2,p3,ep = FaultScanner.slopes(4.0,2.0,2.0,5.0,gsx)
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
  plot3(gsx,skins=skins,smax=15.0,slices=[85,5,60],png="skinss1")



def goUnfaultS():
  if not plotOnly:
    gx = readImage(gcfile)
    fw = zerofloat(n1,n2,n3)
    lof = LocalOrientFilter(4.0,2.0,2.0)
    et = lof.applyForTensors(gx)
    et.setEigenvalues(0.001,1.0,1.0)

    wp = fillfloat(1.0,n1,n2,n3)
    skins = readSkins(fslbase)
    fsc = FaultSlipConstraints(skins)
    sp = fsc.screenPoints(wp)

    uf = UnfaultS(8.0,4.0)
    uf.setIters(200)
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
    gx = readImage(gncfile)
    fw = zerofloat(n1,n2,n3)
    uf = UnfaultS(8.0,4.0)
    uf.applyShifts([t1,t2,t3],gx,fw)
  plot3(gx)
  plot3(fw,clab="Amplitude",png="fw")

def goFlatten():
  gw = readImage(fwsfile)
  if not plotOnly:
    sigma1,sigma2,sigma3,pmax = 2.0,1.0,1.0,5.0
    p2,p3,ep = FaultScanner.slopes(sigma1,sigma2,sigma3,pmax,gw)
    zm = ZeroMask(0.2,1.0,1.0,1.0,gw)
    zero = 0.00;
    tiny = 0.01;
    zm.setValue(zero,p2)#set inline slopes for samples above water bottom
    zm.setValue(zero,p3)#set crossline slopes for samples above water bottom
    zm.setValue(tiny,ep)#set planarities for samples above water bottom
    wp = pow(ep,10)
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
    fl3.setIterations(0.01,400);
    mp = fl3.getMappingsFromSlopes(s1,s2,s3,p2,p3,wp,None,sfs,rs)
    gt = mp.u1
    su1 = Sampling(n1,1,1)
    #gu  = mp.flatten(gw)
    gu = fl3.flatten(s1,su1,gt,gw)
    writeImage(fgfile,gu)
    writeImage(rgtfile,gt)
  gx = readImage(gcfile)
  gu = readImage(fgfile)
  gt = readImage(rgtfile)
  plot3(gx)
  plot3(gw)
  plot3(gu,png="gu")



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
  sp.plotPanel.setColorBarWidthMinimum(72)
  sp.setHLimits(sw.first-sw.delta/2,sw.last)
  sp.setVLimits(sz.first,sz.last)
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

def plot2(s1,s2,x,cmap=ColorMap.GRAY,clab=None,cmin=0,cmax=0,
         title=None,interp=True,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  if title:
    sp.setTitle(title)
  sp.addColorBar(clab)
  sp.setSize(680,600)
  sp.plotPanel.setColorBarWidthMinimum(100)
  pv = sp.addPixels(s1,s2,x)
  pv.setColorModel(cmap)
  if interp:
    pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  else:
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if pngDir and png:
    sp.paintToPng(300,3.333,pngDir+png+".png")


def plot3(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          skins=None,smax=0.0,slices=None, htgs=None,
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
    cbar.setWidthMinimum(137)
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
