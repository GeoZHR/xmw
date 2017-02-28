#############################################################################
"""
Demo of dynamic warping for automatic picking
Author: Xinming Wu, University of Texas at Austin
Version: 2016.06.01
"""


from utils2 import * 
#setupForSubset("nwc1")
#setupForSubset("nwc2")
#setupForSubset("lulia")
#setupForSubset("bahamas")
#setupForSubset("fd2")
#setupForSubset("curt")
#setupForSubset("curt2")
#setupForSubset("tp2")
#setupForSubset("dgb")
#s1,s2,s3 = getSamplings()
#n1,n2,n3 = s1.count,s2.count,s3.count
#f1,f2,f3 = s1.getFirst(),s2.getFirst(),s3.getFirst()
#d1,d2,d3 = s1.getDelta(),s2.getDelta(),s3.getDelta()


#############################################################################
fxfile = "gx" # 
fxfile = "gx185" # 
fxfile = "gx569" # 
fxfile = "fxx582" # 
fxfile = "gx286" # 
fxfile = "gx666" # 
fxfile = "fx172" # 
fxfile = "gxsb" # 
fxfile = "gxsb2" # 
fxfile = "tp73s" # 
fxsfile = "fxs" # 
gxfile = "gx" # 
flfile = "fl" # 
ftfile = "ft" # 
fltfile = "flt" # 
fttfile = "ftt" # 
ls1file = "ls1" # 
ls2file = "ls2" # 
hsfile = "hs"
elfile = "el"
p2file = "p2" # seismic slopes
p3file = "p3" # seismic slopes
#pngDir = getPngDir()
plotOnly = False
pngDir =  "../../../png/hdw/"
minTheta,maxTheta = 75,85
sigmaTheta = 40

# These parameters control the construction of fault skins.
# See the class FaultSkinner for more information.
lowerLikelihood = 0.25
upperLikelihood = 0.55
minSize = 50

def main(args):
  #goTpd()
  #goCurt1()
  #goCurt2()
  #goTpdPaint()
  #goTpdDerek()
  #goTpdFlatten()
  #goCurt1Flatten()
  #goLulia()
  #goLuliaFlatten()
  goCurt2Flatten()
  #goCurt3Flatten()
  #goFd()
  #goFdFlatten()
def goFaultScan(gx):
  print "goScan ..."
  n2 = len(gx)
  n1 = len(gx[0])
  s1 = Sampling(n1)
  s2 = Sampling(n2)
  if not plotOnly:
    gx = FaultScanner2.taper(10,0,gx)
    fs = FaultScanner2(sigmaTheta)
    sig1,sig2,smooth=16.0,2.0,4.0
    fl,ft = fs.scan(minTheta,maxTheta,sig1,sig2,smooth,gx)
    print "fl min =",min(fl)," max =",max(fl)
    print "ft min =",min(ft)," max =",max(ft)
    writeImage(flfile,fl)
    writeImage(ftfile,ft)
  else:
    fl = readImage(flfile)
    ft = readImage(ftfile)

def goCurt():
  setupForSubset("curt")
  fxfile = "gxsb" # fx = copy(230,2700,280,1300,fx)
  s1,s2,s3 = getSamplings()
  n1,n2,n3 = s1.count,s2.count,s3.count
  dl = 60
  nl = 2*dl+1
  gx = readImage(fxfile)
  fx = zerofloat(n1,n2)
  for i2 in range(n2):
    fx[i2] = gx[n2-i2-1]
  fs = copy(210,n2,0,0,fx)
  writeImage("gxsb1",fs)
  plot2(s1,s2,fx,cmin=min(fx)/2,cmax=max(fx)/2)

def goLulia():
  setupForSubset("lulia")
  fxfile = "gx286" 
  s1,s2,s3 = getSamplings()
  n1,n2,n3 = s1.count,s2.count,s3.count
  fx = readImage(fxfile)
  fx = copy(360,n2,0,0,fx)
  n1 = 360 
  s1 = Sampling(n1)
  dl = 150
  nl = 2*dl+1
  fx = gain(fx)
  dp = DynamicPicking(-dl,dl)
  if not plotOnly:
    lsf = LocalSlopeFinder(8,2,5) 
    p = zerofloat(n1,n2)
    el= zerofloat(n1,n2)
    lsf.findSlopes(fx,p,el);
    dp.setGate(-1,1)
    dp.setWeights(1.0,0.2);
    gx = zerofloat(n1,n2)
    us = zerofloat(n2,n1)
    ks = rampint(0,1,n1)
    ds = zerofloat(nl,n1,n2)
    us = dp.pick(25,dl,0.5,ks,p,fx,ds)
    writeImage(hsfile,us)
    writeImage(elfile,el)
  else:
    us = readImageX(n2,n1,hsfile)
    el = readImage(elfile)
  '''
  dp.smoothHorizons(8,el,us)
  dp.refine(fx,us)
  dp.smoothHorizons(8,el,us)
  dp.refine(fx,us)
  dp.smoothHorizons(4,el,us)
  '''
  gx=dp.flattenWithHorizons(us,fx)
  cmin = min(fx)*0.6
  cmax = max(fx)*0.6
  plot2(s1,s2,fx,cmin=cmin,cmax=cmax)
  plot2(s1,s2,gx,cmin=cmin,cmax=cmax)
  plot2(s1,s2,fx,u=us,cmin=cmin,cmax=cmax,color=Color.RED)
  #plot2(s1,s2,fx,u=uss,cmin=cmin,cmax=cmax,color=Color.RED)
  #plot2(s1,s2,fx,u=ut,cmin=min(fx)/2,cmax=max(fx)/2,color=Color.YELLOW)
  #plot3(ds,k2=20,cmin=min(ds),cmax=max(ds),cmap=ColorMap.JET)

def goCurt1():
  setupForSubset("curt1")
  fxfile = "gxsb1" # fx = copy(210,3000,280,1300,fx)
  s1,s2,s3 = getSamplings()
  n1,n2,n3 = s1.count,s2.count,s3.count
  dl = 55
  nl = 2*dl+1
  fx = readImage(fxfile)
  fx = gain(fx)
  dp = DynamicPicking(-dl,dl)
  if not plotOnly:
    lsf = LocalSlopeFinder(8,2,5) 
    p = zerofloat(n1,n2)
    el= zerofloat(n1,n2)
    lsf.findSlopes(fx,p,el);
    dp.setGate(-1,1)
    dp.setWeights(1.0,0.2);
    gx = zerofloat(n1,n2)
    us = zerofloat(n2,n1)
    ks = rampint(0,1,n1)
    ds = zerofloat(nl,n1,n2)
    us = dp.pick(25,dl,0.5,ks,p,fx,ds)
    writeImage(hsfile,us)
    writeImage(elfile,el)
  else:
    us = readImageX(n2,n1,hsfile)
    el = readImage(elfile)
  dp.smoothHorizons(8,el,us)
  dp.refine(fx,us)
  dp.smoothHorizons(8,el,us)
  dp.refine(fx,us)
  dp.smoothHorizons(4,el,us)
  gx=dp.flattenWithHorizons(us,fx)
  cmin = min(fx)*0.6
  cmax = max(fx)*0.6
  plot2(s1,s2,fx,cmin=cmin,cmax=cmax)
  plot2(s1,s2,gx,cmin=cmin,cmax=cmax)
  plot2(s1,s2,fx,u=us,cmin=cmin,cmax=cmax,color=Color.RED)
  #plot2(s1,s2,fx,u=uss,cmin=cmin,cmax=cmax,color=Color.RED)
  #plot2(s1,s2,fx,u=ut,cmin=min(fx)/2,cmax=max(fx)/2,color=Color.YELLOW)
  plot3(ds,k2=20,cmin=min(ds),cmax=max(ds),cmap=ColorMap.JET)


def goCurt2():
  setupForSubset("curt2")
  fxfile = "gxsb2" # fx = copy(230,2700,280,1300,fx)
  s1,s2,s3 = getSamplings()
  n1,n2,n3 = s1.count,s2.count,s3.count
  dl = 30
  nl = 2*dl+1
  fx = readImage(fxfile) #fs = copy(230,3200,280,0,fx)
  fx = gain(fx)
  lsf = LocalSlopeFinder(4,2,5) 
  p2= zerofloat(n1,n2)
  el= zerofloat(n1,n2)
  lsf.findSlopes(fx,p2,el);
  dp = DynamicPicking(-dl,dl)
  ux = zerofloat(n1,n2)
  ut = zerofloat(n2,n1)
  gx = dp.flatten(p2,fx,ux)
  cmin = min(fx)*0.6
  cmax = max(fx)*0.6
  for i2 in range(n2): 
    for i1 in range(n1): 
      ut[i1][i2] = ux[i2][i1]
  us = copy(ut)
  dp.smoothHorizons(8,el,us)
  gs=dp.flattenWithHorizons(us,fx)
  plot2(s1,s2,fx,cmin=cmin,cmax=cmax)
  plot2(s1,s2,gx,cmin=cmin,cmax=cmax)
  plot2(s1,s2,gs,cmin=cmin,cmax=cmax)
  plot2(s1,s2,fx,u=ut,cmin=cmin,cmax=cmax,color=Color.RED)
  plot2(s1,s2,fx,u=us,cmin=cmin,cmax=cmax,color=Color.RED)
  '''
  if not plotOnly:
    lsf = LocalSlopeFinder(4,2,5) 
    p = zerofloat(n1,n2)
    el= zerofloat(n1,n2)
    lsf.findSlopes(fx,p,el);
    dp.setGate(-1,1)
    dp.setWeights(1.0,0.2);
    gx = zerofloat(n1,n2)
    us = zerofloat(n2,n1)
    ks = rampint(0,1,n1)
    ds = zerofloat(nl,n1,n2)
    us = dp.pick(30,dl,0.2,ks,p,fx,ds)
    #writeImage(hsfile,us)
    #writeImage(elfile,el)
  else:
    us = readImageX(n2,n1,hsfile)
    el = readImage(elfile)
  dp.smoothHorizons(8,el,us)
  dp.refine(fx,us)
  dp.smoothHorizons(8,el,us)
  dp.refine(fx,us)
  dp.smoothHorizons(4,el,us)
  gx=dp.flattenWithHorizons(us,fx)
  cmin = min(fx)*0.6
  cmax = max(fx)*0.6
  plot2(s1,s2,fx,cmin=cmin,cmax=cmax)
  plot2(s1,s2,gx,cmin=cmin,cmax=cmax)
  plot2(s1,s2,fx,u=us,cmin=cmin,cmax=cmax,color=Color.RED)
  #plot2(s1,s2,fx,u=uss,cmin=cmin,cmax=cmax,color=Color.RED)
  #plot2(s1,s2,fx,u=ut,cmin=min(fx)/2,cmax=max(fx)/2,color=Color.YELLOW)
  #plot3(ds,k2=20,cmin=min(ds),cmax=max(ds),cmap=ColorMap.JET)
  '''
  dp = DynamicPicking(-dl,dl)
def goFdFlatten():
  setupForSubset("fd2")
  fxfile = "fx172" # fx = copy(230,2700,280,1300,fx)
  s1,s2,s3 = getSamplings()
  n1,n2,n3 = s1.count,s2.count,s3.count
  fx = readImage(fxfile)
  fx = copy(260,800,20,0,fx)
  n1 = 260
  n2 = 800
  s1 = Sampling(n1)
  s2 = Sampling(n2)
  dl = 50
  fx = gain(fx)
  lsf = LocalSlopeFinder(4,1,5) 
  p2= zerofloat(n1,n2)
  el= zerofloat(n1,n2)
  lsf.findSlopes(fx,p2,el);
  df = DynamicFlattener2(-dl,dl)
  df.setStrainMax(1.0)
  df.setWindow(1,10)
  #df.setErrorExponent(4)
  df.setShiftSmoothing(1)
  fs = copy(fx)
  rgf = RecursiveGaussianFilterP(2)
  #rgf.apply00(fx,fs)
  ux = zerofloat(n1,n2)
  gx = df.flatten(el,p2,fs,ux)
  ut = zerofloat(n2,n1)
  for i2 in range(n2): 
    for i1 in range(n1): 
      ut[i1][i2] = ux[i2][i1]
  dp = DynamicPicking(-dl,dl)
  gx = dp.flattenWithHorizons(ut,fx)
  hw,vw=n2*2,round(n1*2.5)
  plot2(s1,s2,fx,hw=hw,vw=vw,cmin=min(fx)/2,cmax=max(fx)/2)
  plot2(s1,s2,gx,hw=hw,vw=vw,cmin=min(fx)/2,cmax=max(fx)/2)
  plot2(s1,s2,fx,u=ut,hw=hw,vw=vw,cmin=min(fx)/2,cmax=max(fx)/2)

def goLuliaFlatten():
  setupForSubset("lulia")
  fxfile = "gx286" # fx = copy(230,2700,280,1300,fx)
  s1,s2,s3 = getSamplings()
  n1,n2,n3 = s1.count,s2.count,s3.count
  fx = readImage(fxfile)
  fx = copy(360,n2,0,0,fx)
  n1 = 360
  s1 = Sampling(n1)
  dl = 150
  fx = gain(fx)
  lsf = LocalSlopeFinder(8,2,5) 
  p2= zerofloat(n1,n2)
  el= zerofloat(n1,n2)
  lsf.findSlopes(fx,p2,el);
  df = DynamicFlattener2(-dl,dl)
  df.setStrainMax(0.2)
  df.setWindow(200)
  df.setShiftSmoothing(1)
  rgf = RecursiveGaussianFilterP(1)
  #rgf.apply00(fx,fx)
  ux = zerofloat(n1,n2)
  gx = df.flatten(el,fx,ux)
  ut = zerofloat(n2,n1)
  for i2 in range(n2): 
    for i1 in range(n1): 
      ut[i1][i2] = ux[i2][i1]
  hw,vw=n2,round(n1*2.5)
  plot2(s1,s2,fx,cmin=min(fx)/2,cmax=max(fx)/2)
  plot2(s1,s2,gx,cmin=min(fx)/2,cmax=max(fx)/2)
  plot2(s1,s2,fx,u=ut,cmin=min(fx)/2,cmax=max(fx)/2)

def goCurt1Flatten():
  setupForSubset("curt1")
  fxfile = "gxsb1" # fx = copy(230,2700,280,1300,fx)
  s1,s2,s3 = getSamplings()
  n1,n2,n3 = s1.count,s2.count,s3.count
  dl = 55
  fx = readImage(fxfile)
  '''
  fx = zerofloat(n1,n2)
  for i2 in range(n2):
    fx[i2] = gx[n2-i2-1]
  '''
  fx = gain(fx)
  lsf = LocalSlopeFinder(8,2,5) 
  p2= zerofloat(n1,n2)
  el= zerofloat(n1,n2)
  lsf.findSlopes(fx,p2,el);
  df = DynamicFlattener(-dl,dl)
  df.setStrainMax(0.5)
  df.setWindow(2,200)
  df.setShiftSmoothing(1)
  df.setErrorSmoothing(3)
  fs = zerofloat(n1,n2)
  rgf = RecursiveGaussianFilterP(1)
  rgf.apply00(fx,fs)
  ux = zerofloat(n1,n2)
  gx = df.flatten(el,fs,ux)
  ut = zerofloat(n2,n1)
  for i2 in range(n2): 
    for i1 in range(n1): 
      ut[i1][i2] = ux[i2][i1]
  hw,vw=n2,round(n1*2.5)
  vlabel = "Relative geologic time (samples)"
  plot2(s1,s2,fx,cmin=min(fx)/2,cmax=max(fx)/2,png="curt1Seis")
  plot2(s1,s2,gx,cmin=min(gx)/2,cmax=max(gx)/2,vlabel=vlabel,png="curt1Flatten")
  us = []
  dp = DynamicPicking(-dl,dl)
  el = pow(el,0)
  dp.smoothHorizons(1,el,ut)
  for i1 in range(1,n1,10):
    us.append(ut[i1])
  plot2(s1,s2,fx,u=us,cmin=min(fx)/2,cmax=max(fx)/2,png="curt1Horizons")

def goCurt2Flatten():
  setupForSubset("curt2")
  fxfile = "gxsb2" # fx = copy(230,2700,280,1300,fx)
  s1,s2,s3 = getSamplings()
  n1,n2,n3 = s1.count,s2.count,s3.count
  dl = 40
  fx = readImage(fxfile)
  fx = gain(fx)
  #goFaultScan(fx)
  #goFaultThin(fx)
  #goFaultCurve(fx)
  lsf = LocalSlopeFinder(8,2,5) 
  p2= zerofloat(n1,n2)
  el= zerofloat(n1,n2)
  lsf.findSlopes(fx,p2,el);
  df = DynamicFlattener(-dl,dl)
  df.setStrainMax(0.2)
  df.setWindow(2,200)
  df.setGate(10)
  df.setErrorSmoothing(3)
  df.setShiftSmoothing(1)
  ux = zerofloat(n1,n2)
  rgf = RecursiveGaussianFilterP(2)
  fs = zerofloat(n1,n2)
  rgf.apply00(fx,fs)
  fl = readImage(flfile)
  gx = df.flatten(fl,fs,ux)
  ut = zerofloat(n2,n1)
  for i2 in range(n2): 
    for i1 in range(n1): 
      ut[i1][i2] = ux[i2][i1]
  gh = df.flattenWithHorizons(ut,fx)
  dp = DynamicPicking(-dl,dl)
  el = pow(el,10)
  dp.smoothHorizons(2,el,ut)
  hw,vw=n2,round(n1*2.5)
  plot2(s1,s2,fx,cmin=min(fx)/2,cmax=max(fx)/2,png="fxCurt2")
  plot2(s1,s2,gx,cmin=min(fx)/2,cmax=max(fx)/2,png="gxCurt2")
  plot2(s1,s2,gh,cmin=min(fx)/2,cmax=max(fx)/2,png="ghCurt2")
  us = []
  for i1 in range(1,n1,15):
    us.append(ut[i1])
  plot2(s1,s2,fx,u=us,cmin=min(fx)*0.6,cmax=max(fx)*0.6,png="usCurt2")

def goCurt3Flatten():
  setupForSubset("curt3")
  fxfile = "fx400" # fx = copy(230,2700,280,1300,fx)
  s1,s2,s3 = getSamplings()
  n1,n2,n3 = s1.count,s2.count,s3.count
  dl = 40
  fx = readImage(fxfile)
  fx = gain(fx)
  #goFaultScan(fx)
  #goFaultThin(fx)
  #goFaultCurve(fx)
  lsf = LocalSlopeFinder(8,2,5) 
  p2= zerofloat(n1,n2)
  el= zerofloat(n1,n2)
  lsf.findSlopes(fx,p2,el);
  df = DynamicFlattener2(-dl,dl)
  df.setStrainMax(0.2)
  df.setWindow(2,200)
  df.setErrorSmoothing(3)
  df.setShiftSmoothing(1)
  ux = zerofloat(n1,n2)
  fs = zerofloat(n1,n2)
  rgf = RecursiveGaussianFilterP(2)
  rgf.apply00(fx,fs)
  gx = df.flatten(el,p2,fs,ux)
  ut = zerofloat(n2,n1)
  for i2 in range(n2): 
    for i1 in range(n1): 
      ut[i1][i2] = ux[i2][i1]
  #df.refine(fx,ut)
  gh = df.flattenWithHorizons(ut,fx)
  dp = DynamicPicking(-dl,dl)
  rgf = RecursiveGaussianFilterP(2)
  u2 = zerofloat(n2,n1)
  rgf.applyX1(ut,u2)
  hw,vw=n2,round(n1*2.5)
  plot2(s1,s2,fx,cmin=min(fx)/2,cmax=max(fx)/2)
  plot2(s1,s2,gx,cmin=min(fx)/2,cmax=max(fx)/2)
  plot2(s1,s2,gh,cmin=min(fx)/2,cmax=max(fx)/2)
  plot2(s1,s2,fx,u=ut,cmin=min(fx)/2,cmax=max(fx)/2)

def goTpdPaint():
  setupForSubset("tpdPaint")
  ppfile = "paint" # 
  fxfile = "tp73s" # 
  s1,s2,s3 = getSamplings()
  n1,n2,n3 = s1.count,s2.count,s3.count
  dl = 60
  fx = readImageM(fxfile)
  pp = readImageM(ppfile)
  fh = FlattenHelper(s1)
  gx = fh.flattenByRgt(pp,fx)
  ut = fh.rgtToHorizonVolume(pp)
  m1 = 240
  fx = copy(m1,n2,0,0,fx)
  gx = copy(m1,n2,0,0,gx)
  c1 = Sampling(m1)
  hw,vw=n2,round(n1*2.5)
  plot2(c1,s2,fx,vint=20,hint=50,hw=hw,vw=vw,cmin=min(fx)/2,cmax=max(fx)/2,png="tpfx")
  plot2(c1,s2,gx,vint=20,hint=50,hw=hw,vw=vw,cmin=min(fx)/2,cmax=max(fx)/2,png="tpgxPaint")
  us = []
  for i1 in range(2,m1,10):
    ui = zerofloat(n2)
    for i2 in range(n2):
      ui[i2] = ut[i2][i1]
    us.append(ui)
  plot2(c1,s2,fx,u=us,vint=20,hint=50,hw=hw,vw=vw,cmin=min(fx)*0.6,cmax=max(fx)*0.6,png="tphsPaint")
  #plot2(s1,s2,fx,u=uss,vint=20,hint=50,hw=hw,vw=vw,cmin=min(fx)/2,cmax=max(fx)/2)
  #plot2(s1,s2,fx,u=ut,cmin=min(fx)/2,cmax=max(fx)/2,color=Color.YELLOW)
  #plot3(ds,k2=20,cmin=min(ds),cmax=max(ds),cmap=ColorMap.JET)
def goTpdDerek():
  setupForSubset("tpdPaint")
  ppfile = "paint" # 
  fxfile = "tp73s" # 
  fx = readImageM(fxfile)
  s1,s2,s3 = getSamplings()
  n1,n2,n3 = s1.count,s2.count,s3.count
  sigma1,sigma2=8.0,2.0 # good for Teapot Dome image tp73
  pmax = 10.0
  lsf = LocalSlopeFinder(sigma1,sigma2,pmax)
  p2 = zerofloat(n1,n2)
  wp = zerofloat(n1,n2)
  lsf.findSlopes(fx,p2,wp)
  wp = pow(wp,8)
  fl = Flattener2()
  fl.setWeight1(0.000)
  fl.setIterations(0.01,1000)
  fl.setSmoothings(4.0,8.0)
  fm = fl.getMappingsFromSlopes(s1,s2,p2,wp)
  gx = fm.flatten(fx)
  m1 = 240
  fx = copy(m1,n2,0,0,fx)
  gx = copy(m1,n2,0,0,gx)
  c1 = Sampling(m1)
  hw,vw=n2,round(n1*2.5)
  plot2(c1,s2,fx,vint=20,hint=50,hw=hw,vw=vw,cmin=min(fx)/2,cmax=max(fx)/2,png="tpfx")
  plot2(c1,s2,gx,vint=20,hint=50,hw=hw,vw=vw,cmin=min(fx)/2,cmax=max(fx)/2,png="tpgxDerek")
  us = []
  for i1 in range(2,m1,10):
    ui = zerofloat(n2)
    for i2 in range(n2):
      ui[i2] = fm.x1[i2][i1]
    us.append(ui)
  plot2(c1,s2,fx,u=us,vint=20,hint=50,hw=hw,vw=vw,cmin=min(fx)*0.6,cmax=max(fx)*0.6,png="tphsDerek")
  #plot2(s1,s2,fx,u=uss,vint=20,hint=50,hw=hw,vw=vw,cmin=min(fx)/2,cmax=max(fx)/2)
  #plot2(s1,s2,fx,u=ut,cmin=min(fx)/2,cmax=max(fx)/2,color=Color.YELLOW)
  #plot3(ds,k2=20,cmin=min(ds),cmax=max(ds),cmap=ColorMap.JET)


def goTpdFlatten():
  setupForSubset("tp2")
  fxfile = "tp73s" # 
  s1,s2,s3 = getSamplings()
  n1,n2,n3 = s1.count,s2.count,s3.count
  dl = 60
  fx = readImage(fxfile)
  #fx = gain(fx)
  lsf = LocalSlopeFinder(8,2,5) 
  p2= zerofloat(n1,n2)
  el= zerofloat(n1,n2)
  lsf.findSlopes(fx,p2,el);
  writeImageM("dip",p2)
  df = DynamicFlattener(-dl,dl)
  df.setStrainMax(0.25)
  df.setWindow(1,200)
  df.setGate(5)
  df.setErrorSmoothing(3)
  df.setShiftSmoothing(1)
  rgf = RecursiveGaussianFilterP(1)
  #rgf.apply00(fx,fx)
  ux = zerofloat(n1,n2)
  fl = readImage(flfile)
  gx = df.flatten(fl,fx,ux)
  ut = zerofloat(n2,n1)
  for i2 in range(n2): 
    for i1 in range(n1): 
      ut[i1][i2] = ux[i2][i1]
  dp = DynamicPicking(-dl,dl)
  el = pow(el,0)
  dp.smoothHorizons(1,el,ut)
  m1 = 240
  fx = copy(m1,n2,0,0,fx)
  gx = copy(m1,n2,0,0,gx)
  c1 = Sampling(m1)
  hw,vw=n2,round(n1*2.5)
  plot2(c1,s2,fx,vint=20,hint=50,hw=hw,vw=vw,cmin=min(fx)/2,cmax=max(fx)/2,png="tpfx")
  plot2(c1,s2,gx,vint=20,hint=50,hw=hw,vw=vw,cmin=min(fx)/2,cmax=max(fx)/2,png="tpgx")
  us = []
  for i1 in range(2,m1,10):
    us.append(ut[i1])
  plot2(c1,s2,fx,u=us,vint=20,hint=50,hw=hw,vw=vw,cmin=min(fx)*0.6,cmax=max(fx)*0.6,png="tphs")
  #plot2(s1,s2,fx,u=uss,vint=20,hint=50,hw=hw,vw=vw,cmin=min(fx)/2,cmax=max(fx)/2)
  #plot2(s1,s2,fx,u=ut,cmin=min(fx)/2,cmax=max(fx)/2,color=Color.YELLOW)
  #plot3(ds,k2=20,cmin=min(ds),cmax=max(ds),cmap=ColorMap.JET)

def goFd():
  setupForSubset("fd2")
  fxfile = "fx172" # fx = copy(230,2700,280,1300,fx)
  s1,s2,s3 = getSamplings()
  n1,n2,n3 = s1.count,s2.count,s3.count
  fx = readImage(fxfile)
  fx = copy(160,400,240,0,fx)
  n1 = 160
  n2 = 400
  s1 = Sampling(n1)
  s2 = Sampling(n2)
  dl = 30
  nl = 2*dl+1
  fx = gain(fx)
  dp = DynamicPicking(-dl,dl)
  if not plotOnly:
    lsf = LocalSlopeFinder(8,2,5) 
    p = zerofloat(n1,n2)
    el= zerofloat(n1,n2)
    lsf.findSlopes(fx,p,el);
    dp.setGate(-1,1)
    dp.setWeights(1.0,0.1);
    gx = zerofloat(n1,n2)
    us = zerofloat(n2,n1)
    ut = zerofloat(n2,n1)
    ks = rampint(0,1,n1)
    ds = zerofloat(nl,n1,n2)
    fs = zerofloat(n1,n2)
    rgf = RecursiveGaussianFilterP(2)
    rgf.apply00(fx,fs)
    us = dp.pick(30,dl,1.0,ks,p,fs,ds)
    writeImage(hsfile,us)
    writeImage(elfile,el)
  else:
    us = readImageX(n2,n1,hsfile)
    el = readImage(elfile)
    el = pow(el,0)
  hw,vw=n2,round(n1*2.5)
  '''
  dp.smoothHorizons(2,el,us)
  dp.refine(fx,us)
  dp.smoothHorizons(2,el,us)
  '''
  gx=dp.flattenWithHorizons(us,fx)
  plot2(s1,s2,fx,vint=20,hint=50,hw=hw,vw=vw,cmin=min(fx)/2,cmax=max(fx)/2)
  plot2(s1,s2,gx,vint=20,hint=50,hw=hw,vw=vw,cmin=min(fx)/2,cmax=max(fx)/2)
  plot2(s1,s2,fx,u=us,vint=20,hint=50,hw=hw,vw=vw,cmin=min(fx)/2,cmax=max(fx)/2)
  #plot2(s1,s2,fx,u=uss,vint=20,hint=50,hw=hw,vw=vw,cmin=min(fx)/2,cmax=max(fx)/2)
  #plot2(s1,s2,fx,u=ut,cmin=min(fx)/2,cmax=max(fx)/2,color=Color.YELLOW)
  #plot3(ds,k2=20,cmin=min(ds),cmax=max(ds),cmap=ColorMap.JET)

def goTpd():
  setupForSubset("tp2")
  fxfile = "tp73s" # 
  s1,s2,s3 = getSamplings()
  n1,n2,n3 = s1.count,s2.count,s3.count
  dl = 50
  nl = 2*dl+1
  fx = readImage(fxfile)
  fx = gain(fx)
  dp = DynamicPicking(-dl,dl)
  if not plotOnly:
    lsf = LocalSlopeFinder(8,2,5) 
    p = zerofloat(n1,n2)
    el= zerofloat(n1,n2)
    lsf.findSlopes(fx,p,el);
    dp.setGate(-1,1)
    dp.setWeights(1.0,0.1);
    gx = zerofloat(n1,n2)
    us = zerofloat(n2,n1)
    ut = zerofloat(n2,n1)
    ks = rampint(0,1,n1)
    ds = zerofloat(nl,n1,n2)
    us = dp.pick(30,dl,0.1,ks,p,fx,ds)
    writeImage(hsfile,us)
    writeImage(elfile,el)
  else:
    us = readImageX(n2,n1,hsfile)
    el = readImage(elfile)
    el = pow(el,0)
  hw,vw=n2,round(n1*2.5)
  '''
  dp.smoothHorizons(2,el,us)
  dp.refine(fx,us)
  dp.smoothHorizons(2,el,us)
  '''
  gx=dp.flattenWithHorizons(us,fx)
  plot2(s1,s2,fx,vint=20,hint=50,hw=hw,vw=vw,cmin=min(fx)/2,cmax=max(fx)/2)
  plot2(s1,s2,gx,vint=20,hint=50,hw=hw,vw=vw,cmin=min(fx)/2,cmax=max(fx)/2)
  plot2(s1,s2,fx,u=us,vint=20,hint=50,hw=hw,vw=vw,cmin=min(fx)/2,cmax=max(fx)/2)
  #plot2(s1,s2,fx,u=uss,vint=20,hint=50,hw=hw,vw=vw,cmin=min(fx)/2,cmax=max(fx)/2)
  #plot2(s1,s2,fx,u=ut,cmin=min(fx)/2,cmax=max(fx)/2,color=Color.YELLOW)
  #plot3(ds,k2=20,cmin=min(ds),cmax=max(ds),cmap=ColorMap.JET)
def goDipPick():
  dl = 50
  fx = readImage(fxfile)
  fx = gain(fx)
  fs = zerofloat(n1,n2)
  dp = DynamicPicking(-dl,dl)
  lsf = LocalSlopeFinder(8,4,5) 
  p = zerofloat(n1,n2)
  el= zerofloat(n1,n2)
  lsf.findSlopes(fx,p,el);
  dp.setGate(-1,1)
  dp.setWeights(0.0,1.0);
  gx = zerofloat(n1,n2)
  us = zerofloat(n2,n1)
  ut = zerofloat(n2,n1)
  ef = zerofloat(dl*2+1,n1,n2)
  eb = zerofloat(dl*2+1,n1,n2)
  k1 = 92
  k1 = 299
  k1 = 209
  k1 = 262
  k1 = 130
  k1 = 150
  #es = dp.pickThrough(k1,15,dl,el,p,fx,us)
  for k1 in range(0,n1,4):
    dp.pickDip(k1,dl,el,p,us,ef)
    dp.trackDip(k1,p,ut)
  #es = add(ef,eb)
  #es = mul(es,0.5)
  plot2(s1,s2,fx,cmin=min(fx)/2,cmax=max(fx)/2,)
  plot2(s1,s2,fx,u=us,cmin=min(fx)/2,cmax=max(fx)/2,color=Color.RED)
  plot2(s1,s2,fx,u=ut,cmin=min(fx)/2,cmax=max(fx)/2,color=Color.YELLOW)
  plot3(ef,k2=k1,cmin=-1000,cmax=max(ef)/2,cmap=ColorMap.JET)
  #plot3(eb,k2=k1,cmin=-100,cmax=max(eb)/2,cmap=ColorMap.JET)
  #plot3(es,k2=k1,cmin=-100,cmax=max(es)/2,cmap=ColorMap.JET)
def goDw():
  dl = 50
  fx = readImage(fxfile)
  fx = gain(fx)
  fs = zerofloat(n1,n2)
  lof = LocalOrientFilter(2,1)
  ets = lof.applyForTensors(fx)
  lsf = LocalSmoothingFilter()
  ets.setEigenvalues(0.05,1.0)
  dp = DynamicPicking(-dl,dl)
  lsfp = LocalSlopeFinder(4,2,5) 
  p = zerofloat(n1,n2)
  el= zerofloat(n1,n2)
  lsfp.findSlopes(fx,p,el);
  rgf = RecursiveGaussianFilterP(2)
  #fs = gain(fs)
  lsf.apply(ets,10,fx,fs)
  dp.setGate(-1,1)
  dp.setWeights(1.0,0.2);
  gx = zerofloat(n1,n2)
  us = zerofloat(n2,n1)
  ut = zerofloat(n2,n1)
  for k1 in range(35,100,5):
    es = dp.pickThrough(k1,35,dl,el,p,fs,us)
    es = dp.trackDip(k1,p,ut)
  plot2(s1,s2,fs,cmin=min(fs)/2,cmax=max(fs)/2)
  plot2(s1,s2,fx,u=us,cmin=min(fx)/2,cmax=max(fx)/2,color=Color.RED)
  plot2(s1,s2,fx,u=ut,cmin=min(fx)/2,cmax=max(fx)/2,color=Color.RED)
def goPick():
  fx = readImage(fxfile)
  df = DynamicFlattening(1,-2,2)
  dl = 100
  es = df.computeErrors(10,dl,fx)
  nl = len(es[0][0])
  k1 = 201
  e1 = zerofloat(nl,n2)
  for i2 in range(n2):
    e1[i2] = es[i2][k1]
  opp = OptimalPathPicker(10,1)
  ft = opp.applyTransform(e1)
  m2,m1 = len(e1),len(e1[0])
  wht = opp.applyForWeight(ft)
  tms1 = zerofloat(m2,m1)
  tms2 = zerofloat(m2,m1)
  pik1 = opp.forwardPick(dl,wht,tms1)
  #pik2 = opp.backwardPick(round(pik1[m2-1]),wht,tms2)
  sl = Sampling(m1)
  plot2(sl,s2,e1,u=pik1,cmin=min(e1),cmax=max(e1),color=Color.YELLOW)
  plot3(es,cmin=min(es),cmax=max(es),cmap=ColorMap.JET)

def goFlatten():
  fx = readImage(fxfile)
  #fx = gain(fx)
  lsf = LocalSlopeFinder(8,2,10) 
  p = zerofloat(n1,n2)
  el= zerofloat(n1,n2)
  lsf.findSlopes(fx,p,el);
  df = DynamicFlattening(1,-1,1)
  df.setWeights(1.0,0.0,0.0)
  dl = 50
  e = df.computeErrors(20,dl,fx)
  nl = len(e[0][0])
  eu = zerofloat(nl,n1,n2)
  #es = df.smooth1(0.1,e)
  es = copy(e)
  #us = df.pickX(dl,p,el,es,eu)
  us = df.pick(0.1,dl,p,el,es,eu)
  e1 = zerofloat(nl,n2)
  c1 = zerofloat(nl,n2)
  k1 = 425
  for i2 in range(n2):
    e1[i2] = e[i2][k1]
  plot2(s1,s2,fx,cmin=-1,cmax=1)
  print min(e)
  print max(e)
  sl = Sampling(nl)
 # plot2(sl,s2,e1,u=us[k1],cmin=min(e1),cmax=max(e1),color=Color.YELLOW)
  plot2(s1,s2,fx,u=add(sub(us[k1],dl),k1),cmin=min(fx)/2,cmax=max(fx)/2,color=Color.YELLOW)
  plot2(s1,s2,fx,u=add(sub(us[k1+2],dl),k1+2),cmin=min(fx)/2,cmax=max(fx)/2,color=Color.YELLOW)
  plot3(e,cmin=min(e),cmax=max(e),cmap=ColorMap.JET)
  plot3(es,cmin=min(es),cmax=max(es),cmap=ColorMap.JET)
def timeMark():
  fx = readImage(fxfile)
  #fx = gain(fx)
  lof = LocalOrientFilter(4,2)
  ets = lof.applyForTensors(fx)
  au = fillfloat(0.0001,n1,n2)
  av = mul(-1,fx)
  av = sub(av,min(av))
  av = sub(av,min(av))
  av = div(av,max(av))
  av = fillfloat(1.000,n1,n2)
  ets.setEigenvalues(au,av)
  f = [1]; k1 = [384,325,382]; k2 = [500,700,900]
  f = [1]; k1 = [162]; k2 = [600]
  t = fillfloat(1,n1,n2)
  p = zerofloat(n1,n2)
  for i in range(len(k1)):
    p[k2[i]][k1[i]] = 1
    t[k2[i]][k1[i]] = 0
  bd = BlendedGridder2(ets,f,k1,k2)
  bd.gridNearest(t,p)
  dpp = DynamicProgrammingPicker(5,-2,2)
  u = zerofloat(n2)
  dpp.track(t,u)
  v = fillfloat(0.0001,n1,n2)
  for i2 in range(n2):
    v[i2][round(u[i2])] = 1
  ets.setEigenvalues(v,v)
  for i in range(len(k1)):
    p[k2[i]][k1[i]] = 1
    t[k2[i]][k1[i]] = 0
  bd = BlendedGridder2(ets,f,k1,k2)
  bd.gridNearest(t,p)
  dpp.track(t,u)
  plot(fx,g=t,cmap=jetFillExceptMin(1.0),cmin=0,cmax=1000)
  plot2(s1,s2,fx,u=u,vint=20,hint=50,cmin=min(fx),cmax=max(fx),color=Color.YELLOW)
def goHorizonPik():
  fx = readImage(fxfile)
  fx = gain(fx)
  k = 10
  hp = HorizonPicker(k,-3,3)
  fp = hp.findPeaks(fx)
  plot2(s1,s2,fx)
  plot2(s1,s2,fp)
  '''
  dpp = DynamicProgrammingPicker(k,-3,3)
  dpp.setWeights(1.0,0.01,0.0)
  lsf = LocalSlopeFinder(0,0,10) 
  p = zerofloat(n1,n2)
  el= zerofloat(n1,n2)
  lsf.findSlopes(fx,p,el);
  k1 = [40]; k2 = [0] #fd
  #k1 = [48]; k2 = [n2-1] #fd
  k1 = [138]; k2 = [0] #fd
  k1 = [319]; k2 = [0]
  k1 = [330]; k2 = [0]
  #k1 = [124,112]; k2 = [650,750]
  #k1 = [382,380,328,358]; k2 = [0,900,700,580]
  #k1 = [382]; k2 = [n2-1]
  #k1 = [225]; k2 = [n2-1]
  fx = add(fx,-min(fx))
  fx = div(fx,max(fx))
  el = pow(el,8)
  el = sub(el,min(el))
  el = div(el,max(el))
  u = dpp.pick(k1,k2,p,fx)
  #u = dpp.pickForward(k1,k2,p,el,fx)
  ut  = dpp.trackForward(k1[0],k2[0],p)
  #u = dpp.pickBackward(k1,k2,p,el,fx)
  #ut  = dpp.trackBackward(k1[0],k2[0],p)
  #plot2(s1,s2,el,cmin=0.0,cmax=1)
  #d = dpp.accumulateForward(p,fx)
  #m2 = len(d)
  #m1 = len(d[0])
  #c1 = Sampling(m1)
  #c2 = Sampling(m2)
  #print min(d)
  #print max(d)
  #plot2(c1,c2,d,vint=200,hint=500,cmin=min(d)/2,cmax=max(d)/2)
  plot2(s1,s2,fx,u=u,vint=20,hint=50,cmin=min(fx),cmax=max(fx),color=Color.YELLOW)
  plot2(s1,s2,fx,u=ut,vint=20,hint=50,cmin=min(fx),cmax=max(fx))
  ei = dpp.mapResample(fx)
  eli = dpp.mapResample(el)
  pi = dpp.slopeResample(p)
  m2 = len(pi)
  m1 = len(pi[0])
  emax = sum(abs(ei))
  di = fillfloat(emax,m1,m2)
  #dpp.setControlPoints(k1,k2,di)
  ei = add(ei,-min(ei))
  dpp.accumulateForward(k1[0]*k,k2[0],pi,eli,ei,di)
  dpp.accumulateBackward(k1[0]*k,k2[0],pi,eli,ei,di)
  for i2 in range(m2):
    for i1 in range(m1):
      dii = di[i2][i1]
      if(dii>=emax):
        di[i2][i1] = 10
  c1 = Sampling(m1)
  c2 = Sampling(m2)
  plot2(c1,c2,di,vint=20,hint=50,cmin=0,cmax=max(di))
  plot2(c1,c2,di,u=mul(u,k),vint=20,hint=50,color=Color.YELLOW,cmin=0,cmax=max(di))
  '''
def goCorrelation():
  '''
  fx = readImage1L(fxfile)
  gx = readImage1L(gxfile)
  ls = readImage2L(lsfile)
  '''
  fxs = readImage2L(n1,n2,fxsfile)
  fx = fxs[1]
  gx = fxs[6]
  tp = LocalCorrelationFilter.Type.SYMMETRIC
  wd = LocalCorrelationFilter.Window.GAUSSIAN
  lcf = LocalCorrelationFilter(tp,wd,20)
  lcf.setInputs(fx,gx)
  nl = 81
  lc = zerofloat(n1,nl)
  sl = Sampling(nl)
  for il in range(nl):
    lcf.correlate(il-40,lc[il])
  opp = OptimalPathPicker(3,4)
  #lc = pow(lc,2)
  lc = sub(lc,min(lc))
  lc = div(lc,max(lc))
  wht = opp.applyForWeight(lc)
  tms1 = zerofloat(n1,nl)
  tms2 = zerofloat(n1,nl)
  tmsd = zerofloat(n1,nl)
  pik1 = opp.forwardPick(40,wht,tms1)
  #pik2 = opp.backwardPick(round(pik1[n1-1]),wht,tms2)
  pik2 = opp.backwardPick(40,wht,tms2)
  tmss = add(tms1,tms2)
  tmss = div(1,tmss)
  tmss = sub(tmss,min(tmss))
  tmss = div(tmss,max(tmss))
  print min(tmss)
  plot2(s1,sl,lc,vint=200,hint=10,cmin=min(lc),cmax=max(lc))
  plot2(s1,sl,tmss,vint=200,hint=10,cmin=max(tmss)-0.1,cmax=max(tmss))
  plot2(s1,sl,lc,u=pik1,vint=200,hint=10,cmin=0.2,cmax=1.0)
  plot2(s1,s2,ls,vint=0,hint=.2,cmin=min(ls),cmax=max(ls))

def goAbSemblancePik():
  sem = readImage(semfile)
  opp = OptimalPathPicker(3,2)
  wht = opp.applyForWeight(sem)
  tms1 = zerofloat(n1,n2)
  tms2 = zerofloat(n1,n2)
  pik1 = opp.forwardPick(100,wht,tms1)
  pik2 = opp.backwardPick(round(pik1[n1-1]),wht,tms2)
  pik3 = opp.applyForPath(100,6,sem)

  pik1 = add(mul(pik1,d2),s2.first)
  pik2 = add(mul(pik2,d2),s2.first)
  pik3 = add(mul(pik3,d2),s2.first)
  plot2(s1,s2,tms1,u=pik1,vint=1,hint=200,cmin=min(tms1),cmax=max(tms1))
  plot2(s1,s2,tms2,u=pik2,vint=1,hint=200,cmin=min(tms2),cmax=max(tms2))
  plot2(s1,s2,sem,u=pik1,vint=1,hint=200,cmin=0.2,cmax=1.0)
  plot2(s1,s2,sem,u=pik2,vint=1,hint=200,cmin=0.2,cmax=1.0)
  plot2(s1,s2,sem,u=pik3,vint=1,hint=200,cmin=0.2,cmax=1.0)
def goEnv3dPik():
  fx = readImageL("env")
  fx = pow(fx,0.5)
  fx = sub(fx,min(fx))
  fx = div(fx,max(fx))
  hp = Helper()
  gx = hp.transpose13(fx) 
  opp = OptimalPathPicker(3,1)
  sf1,sf2,sfs = opp.applyForSurface(12,12,gx)
  sf1,sf1s = opp.applyForSurfaceInline(12,12,gx)
  plot3(gx,cmin=0.0,cmax=0.5,cmap=ColorMap.JET)
  #plot3(gx,cmin=0.0,cmax=0.5,cmap=ColorMap.JET,surf=u)
  #plot3(gx,cmin=0.0,cmax=0.5,cmap=ColorMap.JET,surf=su)
  plot3(gx,cmin=0.0,cmax=0.5,cmap=ColorMap.JET,surf=sf1)
  plot3(gx,cmin=0.0,cmax=0.5,cmap=ColorMap.JET,surf=sf2)
  plot3(gx,cmin=0.0,cmax=0.5,cmap=ColorMap.JET,surf=sf1s)
  plot3(gx,cmin=0.0,cmax=0.5,cmap=ColorMap.JET,surf=sfs)

def goABsemblance():
  strainMax = 0.35
  sem = readImage(semfile)
  pik = readImage1D(pikfile)
  print min(pik)
  print max(pik)
  et = sub(1,sem)
  e = dtran(et)
  dw = DynamicProgramming()
  dw.setStrainMax(strainMax)
  #dw.setShiftSmoothing(4)
  u = zerofloat(n1)
  dw.findPath(e,u)
  u = smooth(8,u)
  u = add(mul(u,d2),s2.first)
  fc1 = 95*d1
  nc1 = n1-95
  c1 = Sampling(nc1,d1,fc1)
  uc = copy(nc1,95,u)
  pikc = copy(nc1,95,pik)
  semc = copy(nc1,n2,95,0,sem)
  plot2(c1,s2,semc,uc,vint=1,hint=200,cmin=0.2,cmax=1.0)
  plot2(c1,s2,semc,pikc,vint=1,hint=200,cmin=0.2,cmax=1.0)

def goChannel():
  stms = [1,1,0.9]
  cp1 = [182]
  cs1 = [202]
  cp2 = [79,388]
  cs2 = [81, 17]
  cp3 = [43, 143,302]
  cs3 = [167,158,116]
  cps = [cp1,cp2,cp3]
  css = [cs1,cs2,cs3]
  gx = readImage(chfile)
  e = dtran(gx)
  dw = DynamicProgramming()
  nc = len(cps)
  us = zerofloat(n1,nc)
  for k in range(nc):
    cpk = cps[k]
    csk = css[k]
    dw.setStrainMax(stms[k])
    dw.setControlPoints(n1,n2,1,cpk,csk)
    dw.findPath(e,us[k])
    us[k] = smooth(4,us[k])
    #d = dw.accumulateForward(e)
  c2 = Sampling(cp2[1],1,0)
  c3 = Sampling(cp3[2],1,0)
  ss = [s1,c2,c3]
  plot2(s1,gx,cps=cps,css=css,vint=100,hint=100,cmap=ColorMap.GRAY)
  plot2(s1,gx,us=us,ss=ss,cps=cps,css=css,vint=100,hint=100,cmap=ColorMap.GRAY)
  #plot2(s1,dtran(d),u,vint=200,hint=200,cmap=ColorMap.GRAY)
def goSlices():
  #fx = readImage2L(n1,n2,"sliceNew")
  hx = readImage2L(n1,n2,"sliceNewx")
  gx = readImage2L(n1,n2,"sliceOld")
  #fx = gain(fx)
  hx = gain(hx)
  gx = gain(gx)
  #plot2(s1,s2,fx,vint=200,hint=200,cmin=-2,cmax=2,cmap=ColorMap.GRAY,title="new")
  plot2(s1,s2,hx,vint=200,hint=200,cmin=-2,cmax=2,cmap=ColorMap.GRAY,title="new")
  plot2(s1,s2,gx,vint=200,hint=200,cmin=-2,cmax=2,cmap=ColorMap.GRAY,title="old")
  hxs = copy(150,150,550,550,hx)
  gxs = copy(150,150,550,550,gx)
  c1 = Sampling(150,1,550)
  c2 = Sampling(150,1,550)
  plot2(c1,c2,hxs,vint=10,hint=10,cmin=-2,cmax=2,cmap=ColorMap.GRAY,title="newSub")
  plot2(c1,c2,gxs,vint=10,hint=10,cmin=-2,cmax=2,cmap=ColorMap.GRAY,title="oldSub")

def goSemblance3d():
  fx = readImageL("semb")
  if not plotOnly:
    fx = sub(fx,min(fx))
    fx = div(fx,max(fx))
    fx = pow(fx,0.5)
    fx = sub(fx,min(fx))
    fx = div(fx,max(fx))
    hp = Helper()
    gx = hp.transpose(fx) 
    print min(gx)
    print max(gx)
    dw = DynamicProgramming()
    dw.setStrainMax(0.2,0.2)
    u = zerofloat(n1,n3)
    dw.setErrorSmoothing(2)
    dw.findSurface(sub(1,gx),u)
    eu = dw.getError(gx,u)
    su = dw.smooth(8,8,eu,u)
    u = smooth2(4,8,u)
    writeImageL("su",su)
    writeImageL("u",u)
    os = sub(n2-1,su)
    os = mul(os,0.01)
    os = add(os,1.4)
    writeImageL("offsetPik",os)
  else:
    su = readImage2L(n1,n3,"su")
    u  = readImage2L(n1,n3,"u")
    fx = sub(fx,min(fx))
    fx = div(fx,max(fx))
    fx = pow(fx,0.5)
    fx = sub(fx,min(fx))
    fx = div(fx,max(fx))
    hp = Helper()
    gx = hp.transpose(fx) 
    #os = sub(n2-1,su)
    os = su
    os = mul(os,0.01)
    os = add(os,1.4)
    writeImageL("offsetPik",os)
    print len(os)
    print len(os[0])
  plot3(gx,cmin=0.0,cmax=0.5,cmap=ColorMap.JET,surf=u)
  plot3(gx,cmin=0.0,cmax=0.5,cmap=ColorMap.JET,surf=su)

def goSurface():
  #gx = readImage3D(gxfile)
  gx = readImage3D("gxs")
  for i3 in range (50,60,1):
    for i2 in range (n2):
      for i1 in range (45,n1,1):
        gx[i3][i2][i1] = gx[i3][i2][i1]*1.5
  for i3 in range (100,105,1):
    for i2 in range (n2):
      for i1 in range (40,n1,1):
        gx[i3][i2][i1] = gx[i3][i2][i1]*2.0
  #gs = copy(55,254,137,187,357,223,gx)
  #writeImage("gxs",gs)
  dw = DynamicProgramming()
  dw.setStrainMax(0.8,0.8)
  u = zerofloat(n2,n3)
  dw.setErrorSmoothing(2)
  gx = mul(-1,gx)
  fx = copy(gx)
  dw.findSurface(fx,u)
  u = smooth2(2,u)
  plot3(gx)
  fx = pow(fx,0.4)
  plot3(fx,fx,cmin=0.0,cmax=0.2)
  plot3(gx,surf=u,png="saltSl")
def goSlopes():
  gx = readImage3D(gxfile)
  sigma1,sigma2=8.0,4.0
  lsf = LocalSlopeFinder(sigma1,sigma2,5) 
  p2 = zerofloat(n1,n2,n3)
  p3 = zerofloat(n1,n2,n3)
  ep = zerofloat(n1,n2,n3)
  lsf.findSlopes(gx,p2,p3,ep);
  writeImage(p2file,p2)
  writeImage(p3file,p3)
  plot3(gx,p2,cmin=-1,cmax=1,cmap=jetFill(1.0))
  plot3(gx,p3,cmin=-1,cmax=1,cmap=jetFill(1.0))

def goHorizon():
  k3 = [84,  81,227,387,400,204,189]
  k2 = [212,157,409,269,377,333,137]
  k1 = [142,164,154,202,193,142,167]
  gx = readImage3D(gxfile)
  p2 = readImage3D(p2file)
  p3 = readImage3D(p3file)
  p2 = abs(p2)
  p3 = abs(p3)
  dw = DynamicProgrammingS()
  dw.setControlPoints(n1,n2,n3,2.0,k1,k2,k3,gx)
  u = zerofloat(n2,n3)
  dw.setStrainMax(p2,p3)
  dw.setErrorSmoothing(2)
  dw.findSurface(copy(gx),u)
  u = smooth2(2,u)
  plot3(gx)
  plot3(gx,surf=u,png="saltSl")

def gain(x):
  n2 = len(x)
  n1 = len(x[0])
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(5.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2)
  div(x,sqrt(g),y)
  return y

def smooth(sig,u):
  v = copy(u)
  rgf = RecursiveGaussianFilterP(sig)
  rgf.apply0(u,v)
  return v

def smooth2(sig1,sig2,u):
  v = copy(u)
  rgf1 = RecursiveGaussianFilterP(sig1)
  rgf2 = RecursiveGaussianFilterP(sig2)
  rgf1.apply0X(u,v)
  rgf2.applyX0(v,v)
  return v


def normalize(e):
  emin = min(e)
  emax = max(e)
  return mul(sub(e,emin),1.0/(emax-emin))

def etran(e):
  #return transpose(pow(e,0.25))
  return transpose(e)

def dtran(d):
  return transpose(d)

def makeSequences():
  n = 500
  fpeak = 0.125
  shift = 2.0/fpeak
  #w = Warp1Function.constant(shift,n)
  w = WarpFunction1.sinusoid(shift,n)
  #f = makeCosine(fpeak,n)
  f = makeRandomEvents(n,seed=seed); 
  g = w.warp(f)
  f = addRickerWavelet(fpeak,f)
  g = addRickerWavelet(fpeak,g)
  f = addNoise(nrms,fpeak,f,seed=10*seed+1)
  g = addNoise(nrms,fpeak,g,seed=10*seed+2)
  s = zerofloat(n)
  for i in range(n):
    s[i] = w.ux(i)
  return f,g,s

def makeCosine(freq,n):
  return cos(mul(2.0*PI*freq,rampfloat(0.0,1.0,n)))

def makeRandomEvents(n,seed=0):
  if seed!=0:
    r = Random(seed)
  else:
    r = Random()
  return pow(mul(2.0,sub(randfloat(r,n),0.5)),15.0)

def addRickerWavelet(fpeak,f):
  n = len(f)
  ih = int(3.0/fpeak)
  nh = 1+2*ih
  h = zerofloat(nh)
  for jh in range(nh):
    h[jh] = ricker(fpeak,jh-ih)
  g = zerofloat(n)
  Conv.conv(nh,-ih,h,n,0,f,n,0,g)
  return g

def ricker(fpeak,time):
  x = PI*fpeak*time
  return (1.0-2.0*x*x)*exp(-x*x)

def addNoise(nrms,fpeak,f,seed=0):
  n = len(f)
  if seed!=0:
    r = Random(seed)
  else:
    r = Random()
  nrms *= max(abs(f))
  g = mul(2.0,sub(randfloat(r,n),0.5))
  g = addRickerWavelet(fpeak,g)
  #rgf = RecursiveGaussianFilter(3.0)
  #rgf.apply1(g,g)
  frms = sqrt(sum(mul(f,f))/n)
  grms = sqrt(sum(mul(g,g))/n)
  g = mul(g,nrms*frms/grms)
  return add(f,g)

#############################################################################
# plotting
def jetFill(alpha):
  return ColorMap.setAlpha(ColorMap.JET,alpha)
def jetFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.JET,a)
def jetRamp(alpha):
  return ColorMap.setAlpha(ColorMap.JET,rampfloat(0.0,alpha/256,256))
def bwrRamp(alpha):
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,rampfloat(0.0,alpha/256,256))
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


backgroundColor = Color.WHITE

def plot2f(s1,s2,f,g=None,cmin=None,cmax=None,cmap=None,label=None,png=None):
  n2 = len(f)
  n1 = len(f[0])
  f1,f2 = s1.getFirst(),s2.getFirst()
  d1,d2 = s1.getDelta(),s2.getDelta()
  panel = panel2Teapot()
  panel.setHInterval(1.0)
  panel.setVInterval(1.0)
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
  frame.setFontSize(24)
  frame.setSize(450+80,700)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(400,3.2,pngDir+png+".png")
  return frame


def plot(f,g=None,t=None,cmap=None,cmin=None,cmax=None,cint=None,
        label=None,neareast=False,png=None): 
  orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT;
  panel = PlotPanel(1,1,orientation)#,PlotPanel.AxesPlacement.NONE)
  panel.setVInterval(50)
  panel.setHInterval(200)
  #panel.setHLabel("Inline (traces)")
  #panel.setVLabel("Time (samples)")
  pxv = panel.addPixels(0,0,s1,s2,f);
  pxv.setColorModel(ColorMap.GRAY)
  #pxv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if g:
    pxv.setClips(-1,1)
  else:
    if cmin and cmax:
      pxv.setClips(cmin,cmax)
  if g:
    pv = panel.addPixels(s1,s2,g)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setColorModel(cmap)
    if cmin and cmax:
      pv.setClips(cmin,cmax)
  moc = panel.getMosaic();
  frame = PlotFrame(panel);
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  #frame.setTitle("normal vectors")
  frame.setVisible(True);
  frame.setSize(1400,700)
  frame.setFontSize(24)
  if pngDir and png:
    frame.paintToPng(720,3.333,pngDir+png+".png")

def plot2(s1,s2,c,u=None,vint=30,hint=200,hw=None,vw=None,
          cmin=0.0,cmax=0.0,vlabel="Time (samples)", cmap=ColorMap.GRAY,color=Color.RED,
          title=None,perc=None,png=None):
  n2 = s2.getCount()
  n1 = s1.getCount()
  panel = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT)
          #PlotPanel.AxesPlacement.NONE)
  panel.setHLabel("Inline (traces)")
  panel.setVLabel(vlabel)
  panel.setHLimits(0,s2.first,s2.last)
  panel.setVLimits(0,s1.first,s1.last)
  panel.setVInterval(0,vint)
  panel.setHInterval(0,hint)
  if title:
    panel.addTitle(title)
  cv = panel.addPixels(0,0,s1,s2,c)
  cv.setInterpolation(PixelsView.Interpolation.NEAREST)
  cv.setColorModel(cmap)
  if perc:
    cv.setPercentiles(100-perc,perc)
  elif cmin<cmax:
    cv.setClips(cmin,cmax)
  if u:
    nu = len(u)
    x2 = rampfloat(0,1,n2)
    cp = ColorMap(0,nu,ColorMap.PRISM)
    k = 1
    for iu in range(nu):
      color = Color.MAGENTA
      if(k%4==0): color=Color.MAGENTA
      if(k%4==1): color=Color.GREEN
      if(k%4==2): color=Color.RED
      if(k%4==3): color=Color.BLUE
      uv = panel.addPoints(0,0,u[iu],x2)
      #uv.setLineColor(cp.getColor(iu))
      uv.setLineColor(color)
      uv.setLineWidth(3.0)
      k = k+1
  #panel.addColorBar()
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(backgroundColor)
  #frame.setFontSizeForPrint(8,240)
  #frame.setSize(470,1000)
  frame.setFontSize(16)
  if hw and vw:
    frame.setSize(hw,vw)
  else:
    frame.setSize(round(n2*0.4),round(n1*2))
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(720,3.33333,pngDir+png+".png")

def plotfg(f,g,png=None):
  n = len(f)
  panel = PlotPanel(2,1,PlotPanel.Orientation.X1RIGHT_X2UP)
  panel.mosaic.setHeightElastic(0,25)
  panel.mosaic.setHeightElastic(1,25)
  panel.setHLimits(0,0,n-1)
  fv = panel.addPoints(0,0,f)
  gv = panel.addPoints(1,0,g)
  fv.setLineWidth(2)
  gv.setLineWidth(2)
  panel.setVLimits(0,-1.3,1.3)
  panel.setVLimits(1,-1.3,1.3)
  #panel.setHLabel("sample index i")
  #panel.setVLabel(0,"f")
  #panel.setVLabel(1,"g")
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(backgroundColor)
  frame.setFontSizeForPrint(8,240)
  frame.setSize(1000,650)
  frame.setVisible(True)
  if png and pngDir:
    png += "n"+str(int(10*nrms))
    frame.paintToPng(720,3.33333,pngDir+"/"+png+".png")

def plotc(c,s=None,u=None,cmin=0.0,cmax=0.0,perc=None,png=None):
  n,nlag = len(c[0]),len(c)
  s1 = Sampling(n,1.0,0.0)
  slag = Sampling(nlag,1.0,-(nlag-1)/2)
  panel = PlotPanel(1,1,PlotPanel.Orientation.X1RIGHT_X2UP)
  panel.setHLimits(0,0,n-1)
  panel.setVLimits(0,slag.first,slag.last)
  cv = panel.addPixels(0,0,s1,slag,c)
  cv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if perc:
    cv.setPercentiles(100-perc,perc)
  elif cmin<cmax:
    cv.setClips(cmin,cmax)
  if s:
    sv = panel.addPoints(0,0,s)
    sv.setLineColor(Color.WHITE)
    sv.setLineStyle(PointsView.Line.DOT)
    sv.setLineWidth(3)
  if u:
    uv = panel.addPoints(0,0,u)
    uv.setLineColor(Color.WHITE)
    uv.setLineWidth(3)
  panel.addColorBar()
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(backgroundColor)
  frame.setFontSizeForPrint(8,240)
  frame.setSize(1000,470)
  frame.setVisible(True)
  if png and pngDir:
    png += "n"+str(int(10*nrms))
    png += "s"+str(int(10*strainMax))
    frame.paintToPng(720,3.33333,pngDir+"/"+png+".png")

def plot2c(c,s,u,clip=None,perc=None,png=None):
  n,nlag = len(c[0]),len(c)
  s1 = Sampling(n,1.0,0.0)
  slag = Sampling(nlag,1.0,-(nlag-1)/2)
  panel = PlotPanel(2,1,PlotPanel.Orientation.X1RIGHT_X2UP)
  panel.setHLimits(0,0,n-1)
  panel.setVLimits(0,slag.first,slag.last)
  panel.setVLimits(1,slag.first,slag.last)
  cv0 = panel.addPixels(0,0,s1,slag,c)
  cv1 = panel.addPixels(1,0,s1,slag,c)
  cv0.setInterpolation(PixelsView.Interpolation.NEAREST)
  cv1.setInterpolation(PixelsView.Interpolation.NEAREST)
  cv0.setColorModel(ColorMap.getGray(0.0,0.8))
  cv1.setColorModel(ColorMap.getGray(0.0,0.8))
  if perc:
    cv0.setPercentiles(0,perc)
    cv1.setPercentiles(0,perc)
  elif clip:
    cv0.setClips(0.0,clip)
    cv1.setClips(0.0,clip)
  if s:
    sv = panel.addPoints(1,0,s)
    sv.setLineColor(Color.WHITE)
    sv.setLineStyle(PointsView.Line.DOT)
    sv.setLineWidth(3)
  if u:
    uv = panel.addPoints(1,0,u)
    uv.setLineColor(Color.WHITE)
    uv.setLineWidth(3)
  panel.addColorBar()
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(backgroundColor)
  frame.setFontSizeForPrint(8,240)
  frame.setSize(1000,850)
  frame.setVisible(True)
  if png and pngDir:
    png += "n"+str(int(10*nrms))
    png += "s"+str(int(10*strainMax))
    frame.paintToPng(720,3.33333,pngDir+"/"+png+".png")

def plot3c(c,s,u,cmin=0.0,cmax=0.0,png=None):
  print "c0: min =",min(c[0])," max =",max(c[0])
  print "c1: min =",min(c[1])," max =",max(c[1])
  print "c2: min =",min(c[2])," max =",max(c[2])
  n,nlag = len(c[0][0]),len(c[0])
  s1 = Sampling(n,1.0,0.0)
  slag = Sampling(nlag,1.0,-(nlag-1)/2)
  panel = PlotPanel(3,1,PlotPanel.Orientation.X1RIGHT_X2UP)
  panel.setHLimits(0,0,n-1)
  panel.setVLimits(0,slag.first,slag.last)
  panel.setVLimits(1,slag.first,slag.last)
  panel.setVLimits(2,slag.first,slag.last)
  cv0 = panel.addPixels(0,0,s1,slag,c[0])
  cv1 = panel.addPixels(1,0,s1,slag,c[1])
  cv2 = panel.addPixels(2,0,s1,slag,c[2])
  cv0.setInterpolation(PixelsView.Interpolation.NEAREST)
  cv1.setInterpolation(PixelsView.Interpolation.NEAREST)
  cv2.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmin<cmax:
    cv0.setClips(cmin,cmax)
    cv1.setClips(cmin,cmax)
    cv2.setClips(cmin,cmax)
  if s:
    sv = panel.addPoints(0,0,s[0])
    sv.setLineColor(Color.WHITE)
    sv.setLineStyle(PointsView.Line.DOT)
    sv.setLineWidth(3)
    sv = panel.addPoints(1,0,s[1])
    sv.setLineColor(Color.WHITE)
    sv.setLineStyle(PointsView.Line.DOT)
    sv.setLineWidth(3)
    sv = panel.addPoints(2,0,s[2])
    sv.setLineColor(Color.WHITE)
    sv.setLineStyle(PointsView.Line.DOT)
    sv.setLineWidth(3)
  if u:
    uv = panel.addPoints(0,0,u[0])
    uv.setLineColor(Color.WHITE)
    uv.setLineWidth(3)
    uv = panel.addPoints(1,0,u[1])
    uv.setLineColor(Color.WHITE)
    uv.setLineWidth(3)
    uv = panel.addPoints(2,0,u[2])
    uv.setLineColor(Color.WHITE)
    uv.setLineWidth(3)
  panel.addColorBar()
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(backgroundColor)
  frame.setFontSizeForPrint(8,240)
  frame.setSize(1000,1050)
  frame.setVisible(True)
  if png and pngDir:
    png += "n"+str(int(10*nrms))
    #png += "s"+str(int(10*strainMax))
    frame.paintToPng(720,3.33333,pngDir+"/"+png+".png")

def plot3(f,g=None,k2=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          xyz=None,cells=None,skins=None,fbs=None,surf=None,smax=0.0,
          links=False,curve=False,trace=False,png=None):
  n1 = len(f[0][0])
  n2 = len(f[0])
  n3 = len(f)
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
  if surf:
    tg = TriangleGroup(True, s3, s2, surf)
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
  ipg.setSlices(232,k2,0)
  if cbar:
    sf.setSize(987,700)
  else:
    sf.setSize(850,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  zscale = 0.3*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.5)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(0.0,0.05,-0.08))
  ov.setAzimuthAndElevation(25,45.0)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(137,1,pngDir+png+"cbar.png")
#############################################################################
# Run the function main on the Swing thread
import sys
class _RunMain(Runnable):
  def __init__(self,main):
    self.main = main
  def run(self):
    self.main(sys.argv)
def run(main):
  SwingUtilities.invokeLater(_RunMain(main)) 
run(main)
