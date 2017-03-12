#############################################################################
"""
Demo of dynamic warping for automatic picking
Author: Xinming Wu, University of Texas at Austin
Version: 2016.06.01
"""


from utils import * 
#setupForSubset("semblance")
#setupForSubset("channel")
#setupForSubset("surface")
setupForSubset("env")
#setupForSubset("semblance3d")
#setupForSubset("fd2")
#setupForSubset("dgb")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count
f1,f2,f3 = s1.getFirst(),s2.getFirst(),s3.getFirst()
d1,d2,d3 = s1.getDelta(),s2.getDelta(),s3.getDelta()


#############################################################################
semfile = "avoscnb" # input semblance image
pikfile = "avopikb" # picked path using Sergey's method
chfile = "channel"  # input seismic horizon slice with channels
slfile = "sl" # 3D salt likelihood image
fxfile = "f3d75s" # 
fxfile = "gx" # 
fxfile = "fx126" # 
fxsfile = "fxs" # 
gxfile = "gx" # 
ls1file = "ls1" # 
ls2file = "ls2" # 
p2file = "p2" # seismic slopes
p3file = "p3" # seismic slopes
pngDir = getPngDir()
pngDir = None
plotOnly = False

def main(args):
  #goCorrelation()
  #goHorizonPik()
  #timeMark()
  goAbSemblancePik()
def timeMark():
  fx = readImage(fxfile)
  fx = gain(fx)
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
  dpp = DynamicProgrammingPicker(k,-5,5)
  dpp.setWeights(1.0,0.0,1.0)
  lsf = LocalSlopeFinder(4,2,5) 
  p = zerofloat(n1,n2)
  el= zerofloat(n1,n2)
  lsf.findSlopes(fx,p,el);
  #k1 = [138]; k2 = [0] #fd
  #k1 = [48]; k2 = [n2-1] #fd
  k1 = [319]; k2 = [0]
  #k1 = [124,112]; k2 = [650,750]
  k1 = [382,380,328,358]; k2 = [0,900,700,580]
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
  plot2(s1,s2,el,cmin=0.0,cmax=1)
  plot2(s1,s2,fx,u=u,vint=20,hint=50,cmin=min(fx),cmax=max(fx),color=Color.YELLOW)
  plot2(s1,s2,fx,u=ut,vint=20,hint=50,cmin=min(fx),cmax=max(fx))
  '''
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
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(1.0)
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

def plot2(s1,s2,c,u=None,us=None,ss=None,cps=None,css=None,vint=1,hint=1,
          cmin=0.0,cmax=0.0,cmap=ColorMap.GRAY,color=Color.RED,title=None,perc=None,png=None):
  panel = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT)
          #PlotPanel.AxesPlacement.NONE)
  panel.setHLimits(0,s2.first,s2.last)
  panel.setVLimits(0,s1.first,s1.last)
  panel.setVInterval(0,vint)
  panel.setHInterval(0,hint)
  if title:
    panel.addTitle(title)
  cv = panel.addPixels(0,0,s1,s2,c)
  cv.setInterpolation(PixelsView.Interpolation.LINEAR)
  cv.setColorModel(cmap)
  if perc:
    cv.setPercentiles(100-perc,perc)
  elif cmin<cmax:
    cv.setClips(cmin,cmax)
  if u:
    x2 = rampfloat(0,1,n2)
    uv = panel.addPoints(0,0,u,x2)
    uv.setLineColor(color)
    uv.setLineWidth(2)
  if us:
    colors = [Color.RED,Color.GREEN,Color.BLUE]
    for k in range(len(us)):
      sk = ss[k]
      uk = copy(ss[k].getCount(),0,us[k])
      uv = panel.addPoints(0,0,sk,uk)
      uv.setLineColor(colors[k])
      #uv.setLineColor(Color.WHITE)
      uv.setLineStyle(PointsView.Line.DASH)
      uv.setLineWidth(3)
  if cps and css:
    colors = [Color.RED,Color.GREEN,Color.BLUE]
    for k in range(len(cps)):
      pv = panel.addPoints(0,0,cps[k],css[k])
      pv.setMarkColor(colors[k])
      #uv.setLineColor(Color.WHITE)
      pv.setLineStyle(PointsView.Line.NONE)
      pv.setMarkStyle(PointsView.Mark.FILLED_SQUARE)
      pv.setMarkSize(8)
  #panel.addColorBar()
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(backgroundColor)
  #frame.setFontSizeForPrint(8,240)
  #frame.setSize(470,1000)
  frame.setFontSize(12)
  frame.setSize(n2*1,round(n1*1.5))
  frame.setVisible(True)
  if png and pngDir:
    png += "n"+str(int(10*nrms))
    png += "s"+str(int(10*strainMax))
    frame.paintToPng(720,3.33333,pngDir+"/"+png+".png")

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

def plot3(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
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
  ipg.setSlices(232,63,0)
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
