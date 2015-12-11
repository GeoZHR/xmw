from utils import *

setupForSubset("2d")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

pngDir = None
pngDir = "../../../png/slt/2d/"
gxfile  = "gx" # input image (maybe after bilateral filtering)


def main(args):
  goTest()
  #goSemblance()
def goTest():
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  gx = readImage2d(gxfile)
  lof = LocalOrientFilter(4,4)
  lof.applyForNormal(gx,u1,u2)
  ss = SaltScanner(10,40)
  cx,ax = ss.scan(gx,u1,u2)
  cd = abs(sub(cx,ax))
  plot(gx)
  plot(cx)
  plot(ax)
  plot(cd)
  mul(cd,sub(1,ax),cd)
  plot(cd)
  mul(cd,sub(1,cx),cd)
  plot(cd)
def goSemblance():
  gx = readImage2d(gxfile)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lof = LocalOrientFilter(4,4)
  lof.applyForNormalLinear(gx,u1,u2,el)
  lof = LocalOrientFilter(8,2)
  ets = lof.applyForTensors(gx)
  lsf = LocalSemblanceFilter(4,64)
  sm = lsf.semblance(LocalSemblanceFilter.Direction2.U,ets,el)

  plot(gx)
  plot(el)
  plot(sub(1,sm))
  '''
  plot2(gx,s1,s2,sub(1,el),cmin=0.1,cmax=1,label="saultB")
  plot2(gx,s1,s2,sub(1,sl),cmin=0.1,cmax=1,label="saultA")
  plot2(gx,s1,s2,sub(1,sbd),cmin=0.01,cmax=0.2,label="saultBound")
  '''


def goLineDraw():
  ti = TensorInterp()
  ti.setParameters(20,5,0.5)
  fx = readImage(fxfile)
  sod = StructureOrientDraw()
  sod.setSmoothings(1,10)
  sod.setThreshold(0.6)
  #hx,u1,u2=sod.applyDraw(fx,u1,u2)
  hx,u1,u2=sod.applyDraw(fx)
  #rx,r1,r2=sod.randSample(n1*n2/60,999888,hx,u1,u2) #dave
  #rx,r1,r2=sod.randSample(n1*n2/20,988,hx,u1,u2) #cwp
  rx,r1,r2=sod.randSample(n1*n2/20,988,hx,u1,u2) #sergey
  sm,v1,v2= ti.apply(rx,r1,r2)
  #sm= ti.applyX(rx,r1,r2)
  plot(fx,png="dave")
  plot(hx,png="daveL")
  dx = sub(1,rx)
  plot(dx,png="daveS")
  sm = sub(1,sm)
  plot(sm,cmin=min(sm),cmax=max(sm),png="daveR")

  '''
  au = fillfloat(0.05,n1,n2)
  av = fillfloat(1.0000,n1,n2)
  ets = EigenTensors2(v1,v2,au,av)
  fed = FastExplicitDiffusion()
  sig = 6
  t = 0.5*sig*sig
  wp = fillfloat(1,n1,n2)
  fed.setParameters(t,5,0.1)
  lsf = LocalSmoothingFilter()
  lsf.applySmoothS(rx,rx)
  sm = fed.apply(ets,wp,rx)
  sm = sub(1,sm)
  plot(sm,cmin=min(sm),cmax=max(sm),png="daveRs")
  '''

def goTensorInterpX(fx,u1,u2):
  ti = TensorInterp()
  ti.setParameters(20,5,0.4)
  sm,u1,u2 = ti.apply(fs,v1,v2)
  plot(fe)
  plot(fs)
  plot(sm)


def goTensorInterp():
  fx = readImage(fxfile)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  ti = TensorInterp()
  ti.setParameters(20,5,0.4)
  fe = ti.applyForEdge(fx)
  fe = sub(fe,min(fe))
  fe = div(fe,max(fe))
  lof = LocalOrientFilter(4,4)
  lof.applyForNormal(fx,u1,u2)
  fs = zerofloat(n1,n2)
  v1 = zerofloat(n1,n2)
  v2 = zerofloat(n1,n2)
  for i2 in range(0,n2,3):
    for i1 in range(0,n1,3):
      fs[i2][i1] = fx[i2][i1];
      v1[i2][i1] = u1[i2][i1];
      v2[i2][i1] = u2[i2][i1];
  sm,u1,u2 = ti.apply(fs,v1,v2)
  plot(fe)
  plot(fs)
  plot(sm)

def goFed():
  fx = readImage(fxfile)
  wp = fillfloat(1,n1,n2)
  fs = fillfloat(1,n1,n2)
  u1 = fillfloat(1,n1,n2)
  u2 = fillfloat(1,n1,n2)
  lof = LocalOrientFilter(4,4)
  lof.applyForNormalLinear(fx,u1,u2,wp)
  wp = pow(wp,10.0)
  ets = lof.applyForTensors(fx)
  ets.setEigenvalues(0.0001,1.0)
  fed = FastExplicitDiffusion()
  sig = 4
  t = 0.5*sig*sig
  fed.setParameters(t,3,0.5)
  gx = fed.apply(ets,wp,fx)
  plot(fx,label="fx")
  plot(gx,label="fed")
  lsf = LocalSmoothingFilter()
  lsf.apply(ets,10,wp,fx,fs)
  plot(fs,label="lsf")
  plot(sub(gx,fx),label="gfx")
  plot(sub(fs,fx),label="ffx")


def goTensorInterpFake():
  fx = zerofloat(n1,n2)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  g11 = zerofloat(n1,n2)
  g12 = zerofloat(n1,n2)
  g22 = zerofloat(n1,n2)
  sm = zerofloat(n1,n2)
  for i2 in range(20,60,3):
    if (i2<30) or(i2>40):
      fx[i2][i2] = 1
      u1[i2][i2] = -sin(Math.PI/4)
      u2[i2][i2] = cos(Math.PI/4)
  k = 0
  for i2 in range(20,60,3):
    if (i2<30) or(i2>40):
      fx[i2][50-k] = 1
      u1[i2][50-k] = sin(Math.PI/4)
      u2[i2][50-k] = cos(Math.PI/4)
    k = k+3

  fx[70][60] = 1
  u1[70][60] = 1
  u2[70][60] = 0
  fx[70][80] = 1
  u1[70][80] = 1
  u2[70][80] = 0

  fx[60][70] = 1
  u1[60][70] = 0
  u2[60][70] = 1
  fx[80][70] = 1
  u1[80][70] = 0
  u2[80][70] = 1
  d1 = round(5*sqrt(2))
  d2 = round(5*sqrt(2))
  fx[70-d1][70-d1] = 1
  u1[70-d1][70-d1] = sin(Math.PI/4)
  u2[70-d1][70-d1] = cos(Math.PI/4)
  fx[70+d1][70+d1] = 1
  u1[70+d1][70+d1] = sin(Math.PI/4)
  u2[70+d1][70+d1] = cos(Math.PI/4)
  fx[70-d1][70+d1] = 1
  u1[70-d1][70+d1] = -sin(Math.PI/4)
  u2[70-d1][70+d1] = cos(Math.PI/4)
  fx[70+d1][70-d1] = 1
  u1[70+d1][70-d1] = -sin(Math.PI/4)
  u2[70+d1][70-d1] = cos(Math.PI/4)


  ti = TensorInterp()
  ti.setIterations(100)
  ti.initialTensors(sm,fx,u1,u2,g11,g12,g22)
  sm = ti.apply(fx,u1,u2,g11,g12,g22)
  plot(fx,label="fx")
  plot(sm,label="sm")
def goShockFilter():
  fx = readImage(fxfile)
  plot(fx)
  csf = CoherenceShockFilter()
  csf.setIterations(5)
  csf.setSmoothings(10,10)
  fs = csf.apply(fx)
  fe = csf.applyForEdge(fs)
  fo = csf.applyForEdge(fx)
  plot(fs)
  plot(fe)
  plot(fo)

def goTestXXX(): 
  fx = readImage(fxfile)
  #fx = div(fx,10000)
  u1 = fillfloat(0,n1,n2)
  u2 = fillfloat(0,n1,n2)
  g1 = fillfloat(0,n1,n2)
  g2 = fillfloat(0,n1,n2)
  gx = fillfloat(0,n1,n2)
  lof = LocalOrientFilter(4,4)
  lof.applyForNormal(fx,u1,u2)
  ets = lof.applyForTensors(fx)
  ets.setEigenvalues(0.0001,1.0)
  lsf = LocalSmoothingFilter()
  lsf.apply(ets,10,fx,gx)
  hgf = HermiteGaussFilter2(20,2)
  fs = hgf.applyFilter10(1,fx,u1,u2)
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply10(gx,g1)
  rgf.apply01(gx,g2)
  g1 = mul(g1,u1)
  g2 = mul(g2,u2)
  g2 = add(g1,g2)
  g2 = abs(g2)
  plot(fx)
  plot(abs(fs))
  plot(g2)
  plot(gx)

def goPnz():
  fx = readImage(fxfile)
  tv = TensorVoting2X(6,6)
  gx,os=tv.initialTensorField(0.1,0.9,2,2,fx)
  normalize(gx)
  ss,ps = tv.applyVoting(gx,os)
  normalize(ss)
  ss = tv.findRidges(ss)
  ss,ps = tv.applyVoting(ss,ps)
  normalize(ss)
  ss = tv.findRidges(ss)
  ss,ps = tv.applyVoting(ss,ps)
  ss = copy(n1,n2,ss)
  os = copy(n1,n2,ps)
  normalize(ss)
  #normalize(gx)
  se,fs = tv.smoothEdge(2,2,fx)
  plot(fx)
  plot(fs)
  plot(ss,cmin=0.05,cmax=0.6)
  plot(gx,cmin=0.05,cmax=0.6)
  plot(se,cmin=0.05,cmax=0.6)

  #plot(os,cmap=ColorMap.JET)
  #plot2(fx,s1,s2,g=ss,cmin=0.3,cmax=1.0)
  #plot2(fx,s1,s2,g=gx,cmin=0.3,cmax=1.0)
def goChannel():
  fx = readImage(fxfile)
  fx = div(fx,10000)
  tv = TensorVoting2X(4,8)
  ss,os=tv.initialTensorField(0.1,1.0,2,2,fx)
  normalize(ss)
  plot(ss)
  ss,ps = tv.applyVoting(ss,os)
  ss = tv.findRidges(ss)
  ss,ps = tv.applyVoting(ss,ps)
  ss = tv.findRidges(ss)
  ss,ps = tv.applyVoting(ss,ps)
  ss = tv.findRidges(ss)
  ss,ps = tv.applyVoting(ss,ps)
  ss = tv.findRidges(ss)
  ss,ps = tv.applyVoting(ss,ps)
  normalize(ss)
  ss = copy(n1,n2,ss)
  os = copy(n1,n2,ps)
  plot(fx)
  plot(ss)
  #plot(os,cmap=ColorMap.JET)
  '''
  lof = LocalOrientFilter(4,4)
  ets = lof.applyForTensors(fx)
  ets.setEigenvalues(0.0001,1.0)
  lsf = LocalSmoothingFilter()
  lsf.apply(ets,40,ss,ss)
  plot(ss)
  '''



def goTestXX(): 
  fx = zerofloat(n1,n2)
  u1i = sin(Math.PI/2)
  u2i = -cos(Math.PI/2)
  u1 = fillfloat(u1i,n1,n2)
  u2 = fillfloat(u2i,n1,n2)
  au = fillfloat(0.0001,n1,n2)
  av = fillfloat(1.0000,n1,n2)
  for i2 in range (10,150,5):
    fx[i2][i2] = 1
  et = EigenTensors2(u1,u2,au,av)
  fk = zerofloat(n1,n2)
  lsf = LocalSmoothingFilter()
  lsf.apply(et,40,fx,fk)
  sf = SteerableFilter2(8,10.0)
  fs = sf.applyFilterX(fx,u1,u2)
  plot(fx)
  plot(fk)
  plot(fs)

def goTestX():
  fx = zerofloat(n1,n2)
  ox = fillfloat(0,n1,n2)
  ph = -Math.PI/4
  for i2 in range(100,140,10):
    for i1 in range(100,140,10):
      if (i1==i2):
        fx[i2][i1] = 10
        ox[i2][i1] = ph

  for i2 in range(160,190,10):
    for i1 in range(160,190,10):
      if (i1==i2):
        fx[i2][i1] = 10
        ox[i2][i1] = ph
  '''
  k = 1
  for i2 in range(120,190,1):
    fx[i2][190-k] = 10
    ox[i2][190-k] = -ph
    k=k+1
  '''

  sf = SteerableFilter2(20)
  fs = sf.applyFilter(fx,ox)
  fs = copy(n1,n2,fs)
  plot(fx)
  plot(fs)

def goGaussTest():
  fx = fillfloat(1.0,n1,n2)
  fx[50][50] = 1.0
  #fx = readImage(fxfile)
  g1 = fillfloat(1.0,n1,n2)
  g2 = fillfloat(1.0,n1,n2)
  rgf1 = RecursiveGaussianFilter(10.0,RecursiveGaussianFilter.Method.VAN_VLIET)
  rgf2 = RecursiveGaussianFilterX(10.0)
  rgf1.apply00(fx,g1)
  rgf2.apply00(fx,g2)
  plot(g1,cmin=min(g1),cmax=max(g1))
  plot(g2,cmin=min(g2),cmax=max(g2))
  #plot(abs(sub(fx,g1)),cmin=0,cmax=0.3)
  #plot(abs(sub(fx,g2)),cmin=0,cmax=0.3)
  print g1[50][50]
  print g2[50][50]

  print g1[50][0]
  print g1[50][1]
  print g1[50][2]
  print g1[50][n1-3]
  print g1[50][n1-2]
  print g1[50][n1-1]

  print

  print g2[50][0]
  print g2[50][1]
  print g2[50][2]
  print g2[50][n1-3]
  print g2[50][n1-2]
  print g2[50][n1-1]
def goChannelX():
  fx = readImage(fxfile)
  fx = div(fx,10000)
  lof = LocalOrientFilter(1.0,1.0)
  ets = lof.applyForTensors(fx)
  tv = TensorVoting2(6)
  ss,os=tv.initialTensorField(1,2,2,fx)
  normalize(ss)
  plot(ss)
  ets = tv.applyTensor(8,ss,ets)
  plot(ss)

def goChannelX():
  fx = readImage(fxfile)
  fx = div(fx,10000)
  tv = TensorVoting2(8)
  ss,os=tv.initialTensorField(1,3,4,fx)
  normalize(ss)
  plot(ss)
  ss,ps = tv.applyVoting(ss,os)
  ss = tv.findRidges(ss)
  ss,ps = tv.applyVoting(ss,ps)
  ss = tv.findRidges(ss)
  ss,ps = tv.applyVoting(ss,ps)
  ss = tv.findRidges(ss)
  ss,ps = tv.applyVoting(ss,ps)
  ss = tv.findRidges(ss)
  ss,ps = tv.applyVoting(ss,ps)
  normalize(ss)
  ss = copy(n1,n2,ss)
  os = copy(n1,n2,ps)
  plot(fx)
  plot(ss)
  plot(os,cmap=ColorMap.JET)
  '''
  lof = LocalOrientFilter(4,4)
  ets = lof.applyForTensors(fx)
  ets.setEigenvalues(0.0001,1.0)
  lsf = LocalSmoothingFilter()
  lsf.apply(ets,40,ss,ss)
  plot(ss)
  '''
def normalize(ss):
  sub(ss,min(ss),ss)
  div(ss,max(ss),ss)
def frequencyResponse(x):
  n1 = len(x[0])
  n2 = len(x)
  n1 = FftComplex.nfftSmall(n1)
  n2 = FftComplex.nfftSmall(n2)
  xr = copy(x)
  xi = zerofloat(n1,n2)
  cx = cmplx(xr,xi)
  fft1 = FftComplex(n1)
  fft2 = FftComplex(n2)
  fft1.complexToComplex1(1,n2,cx,cx)
  fft2.complexToComplex2(1,n1,cx,cx)
  ax = cabs(cx)
  a = zerofloat(n1,n2)
  j1 = n1/2
  j2 = n2/2
  copy(n1-j1,n2-j2,0,0,ax,j1,j2,a)
  copy(j1,j2,n1-j1,n2-j2,ax,0,0,a)
  copy(n1-j1,j2,0,n2-j2,ax,j1,0,a)
  copy(j1,n2-j2,n1-j1,0,ax,0,j2,a)
  return a
  
def goExForceF3d():
  sig,sig1,sig2=1,4,4
  fx = readImage(fxfile)
  ws = zerofloat(n1,n2)
  gs = zerofloat(n1,n2,2)
  gvf = GradientVectorFlow()
  ed = gvf.applyForEdge(sig,sig1,sig2,fx)
  et = gvf.applyForGVX(sig,sig1,sig2,ed,gs,ws) #channel 1
  #gvf.setScale(0.01)
  gvf.setScale(0.0)
  wp = pow(ed,1.0)
  ws = pow(ws,1.0)
  wm = pow(wp,1)
  wm = div(wm,max(wm))
  wm = sub(1,wm)
  fs = gvf.applyForGVF(None,gs[0],gs[1],ws,wm)
  writeImage(edfile,ed)
  writeImage(f1file,fs[0])
  writeImage(f2file,fs[1])
  plot(fx)
  plot(wp)
  plot(ws)
  plot(wm)
  plotVectors(10,5,gs[0],gs[1],fx)
  plotVectors(100,5,fs[0],fs[1],fx)

def goExForce():
  sig,sig1,sig2=2,2,2
  fx = readImage(fxfile)
  fx = div(fx,2000)
  ws = zerofloat(n1,n2)
  gs = zerofloat(n1,n2,2)
  gvf = GradientVectorFlow()
  ed = gvf.applyForEdge(sig,sig1,sig2,fx)
  et = gvf.applyForGV(sig,sig1,sig2,ed,gs)
  #et = gvf.applyForGVX(sig,sig1,sig2,ed,gs,ws) #channel 1
  #gvf.setScale(0.01)
  gvf.setScale(0.01)
  wp = pow(ed,1.0)
  fs = gvf.applyForGVF(et,gs[0],gs[1],wp,None)
  writeImage(edfile,ed)
  writeImage(f1file,fs[0])
  writeImage(f2file,fs[1])
  plot(fx)
  plot(wp)
  plotVectors(50,5,gs[0],gs[1],fx)
  plotVectors(50,5,fs[0],fs[1],fx)

def goSnakeReal():
  fx = readImage(fxfile)
  f1 = readImage(f1file)
  f2 = readImage(f2file)
  fx = div(fx,2000)
  #c1,c2=418,345
  c1,c2 = 248,165
  kamma,gamma = 2.0,0.6
  alpha,beta  = 0.1,0.1
  fp = zerofloat(n1,n2,2)
  fn = zerofloat(n1,n2,2)
  ac = ActiveContour(s1,s2,c1,c2,5)
  ac.setParameters(kamma,gamma,alpha,beta)
  ac.exForce(10,25,fx,fp,fn)
  fn = mul(fn,0.5)
  ac.firstUpdate(fp,fn)
  snake = ac.getSnake()
  x1 = snake.getArrayX1()
  x2 = snake.getArrayX2()
  pngName = "2ndChannel"
  plot(fx,ap=[x1,x2],png=None)#pngName+"Initial")
  
  ac.setIterations(10)
  fnn = [f1,f2]
  fnn=mul(fnn,10)
  fpn = mul(fp,0.5)
  for k in range(25):
    ac.releaseSnake(fpn,fnn,False)
    x1 = snake.getArrayX1()
    x2 = snake.getArrayX2()
    k1,k2=[],[]
    for j in range(len(x1)):
      x1i = x1[j]
      x2i = x2[j]
      if(x1i>1 and x1i<n1-6):
       k1.append(x1i)
       k2.append(x2i)
    plot(fx,ap=[k1,k2],png=None)#pngName+str(k))


def goForce():
  #c1,c2 = 248,165 # channel 1
  c1,c2 = 418,345 # channel 2
  #c1,c2 = 220,110
  #c1,c2 = 345,340
  #c1,c2 = 250,250
  g = readImage(ffile)
  g = div(g,2000)
  fp = zerofloat(n1,n2,2)
  fn = zerofloat(n1,n2,2)

  fpn = zerofloat(n1,n2,2)
  fnn = zerofloat(n1,n2,2)
  ups = zerofloat(n1,n2,3)
  gvf = GradientVectorFlow()
  gvf.setScale(0.01)
  et = gvf.applyForGV(2,4,4,g,ups) #channel 1
  wp = pow(ups[2],1)
  fnn = gvf.applyForGVF(et,ups[0],ups[1],wp)
  plot(ups[2])
  plotVectors(100,5,fnn[0],fnn[1],g)

  ac = ActiveContour(s1,s2,c1,c2,5)
  ac.setIterations(1)
  ac.exForce(10,25,g,fp,fn)
  #fn = mul(fn,0.04)
  fn = mul(fn,0.5)
  kamma,gamma = 1.0,0.2
  alpha,beta  = 0.1,0.1
  snake = ac.getSnake()
  for k in range(1):
    ac.releaseSnake(s1,s2,kamma,gamma,alpha,beta,fp,fn,True)
    x1 = snake.getArrayX1()
    x2 = snake.getArrayX2()
  pngName = "2ndChannel"
  plot(g,ap=[x1,x2],png=pngName+"Initial")
  g1 = snake.getArrayU1()
  g2 = snake.getArrayU2()
  [fps,fns]=ac.exForceOnSnake(fp,fn)
  
  acn = ActiveContour(s1,s2,c1,c2,x1,x2)
  acn.setIterations(10)
  fnn=mul(fnn,10)
  fpn = mul(fp,0.5)
  snake = acn.getSnake()
  kamma,gamma = 2.0,0.6
  alpha,beta  = 0.1,5.1
  for k in range(20):
    acn.releaseSnake(s1,s2,kamma,gamma,alpha,beta,fpn,fnn,False)
    x1 = snake.getArrayX1()
    x2 = snake.getArrayX2()
    k1,k2=[],[]
    for j in range(len(x1)):
      x1i = x1[j]
      x2i = x2[j]
      if(x1i>1 and x1i<n1-6):
       k1.append(x1i)
       k2.append(x2i)
    plot(g,ap=[k1,k2],png=pngName+str(k))

def goPaint():
  f,t,x = [1],[250],[165]
  s = readImage(ffile)
  smooth = 1.00 # smoothness for blended gridder
  fnull = -1.0
  sg = SimpleGridder2(f,t,x)
  sg.setNullValue(fnull)
  p = sg.grid(s1,s2)
  tensors = makeImageTensors(s)
  plot2(f,t,x,s,s1,s2,None,"Known value","tp2f",et=None)
  bg = makeBlendedGridder(f,t,x,smooth=smooth)
  bg.setTensors(tensors)
  d = bg.gridNearest(fnull,p)
  plot2(f,t,x,s,s1,s2,d,"Time (samples)","tp2t")

def goSnakeFake():
  sig = 1
  sig1,sig2=2,2
  f = goFakeImage()

  up = zerofloat(n1,n2,3)
  gvf = GradientVectorFlow()
  gvf.setScale(0.5) # fake image
  et = gvf.applyForGV(sig,sig1,sig2,f,up)
  us = gvf.applyForGVF(None,up[0],up[1],up[2])
  plot(f)
  #plot(up[2])
  plotVectors(5,3,up[0],up[1],f)   #fake image
  plotVectors(5,3,us[0],us[1],f)   #fake image
 
  us = mul(us,-1.0)
  x1 = [30,25,25,25,25,25,35,45,45,45,45,45,45,45,30]
  x2 = [20,30,50,60,80,85,85,85,70,65,50,40,30,20,20]
  '''
  acs = ActiveSnakeModel(x1,x2)
  snake = acs.releaseSnake(s1,s2,1.0,0.1,0.1,0.1,us[0],us[1])
  '''
  acs = ActiveContour(x1,x2)
  acs.setIterations(1)
  kamma,gamma = 0.0,0.5
  alpha,beta  = 0.1,0.1
  snake = acs.releaseSnake(s1,s2,kamma,gamma,alpha,beta,us[0],us[1])
  x1s = snake.getArrayX1()
  x2s = snake.getArrayX2()
  plot(f,ap=[x1,x2])
  plot(f,ap=[x1s,x2s])

def goFakeImage():
  f = zerofloat(n1,n2)
  for i2 in range(30,50):
    f[i2][30] = 1
    f[i2][40] = 1
  for i2 in range(55,76):
    f[i2][30] = 1
    f[i2][40] = 1
  for i1 in range(30,40):
    f[30][i1] = 1
    f[75][i1] = 1
  return f

def makeBlendedGridder(f,x,y,smooth=0.5,tmax=FLT_MAX,tmx=False):
  bg = BlendedGridder2(f,x,y)
  bg.setSmoothness(smooth)
  #bg.setBlendingKernel(
  #  LocalDiffusionKernel(LocalDiffusionKernel.Stencil.D21))
  if tmax<FLT_MAX:
    bg.setTimeMax(tmax)
  if tmx:
    bg.setTimeMarkerX(tmx)
  return bg

def makeImageTensors(s):
  """ 
  Returns tensors for guiding along features in specified image.
  """
  sig1,sig2 = 10,30
  n1,n2 = len(s[0]),len(s)
  lof = LocalOrientFilter(sig1,sig2)
  t = lof.applyForTensors(s) # structure tensors
  t.setEigenvalues(0.001,0.999)
  '''
  c = coherence(sig1,t,s) # structure-oriented coherence c
  c = clip(0.0,0.99,c) # c clipped to range [0,1)
  t.scale(sub(1.0,c)) # scale structure tensors by 1-c
  t.invertStructure(1.0,1.0) # invert and normalize
  '''
  return t

def coherence(sigma,t,s):
  lsf = LocalSemblanceFilter(sigma,4*sigma)
  return lsf.semblance(LocalSemblanceFilter.Direction2.V,t,s)

def plotVectors1(dh,b1,b2,u1,u2,f):
  np = len(b1)
  p1 = zerofloat(np)
  p2 = zerofloat(np)
  x1 = zerofloat(2,np)
  x2 = zerofloat(2,np)
  for ip in range(0,np,10):
    x1[ip][0] = b1[ip]
    x2[ip][0] = b2[ip]
    x1[ip][1] = b1[ip]+dh*u1[ip]
    x2[ip][1] = b2[ip]+dh*u2[ip]
    p1[ip] = x1[ip][1]
    p2[ip] = x2[ip][1]
  plot(f,xp=[x1,x2],pp=[p1,p2],ap=[b1,b2])


def plotVectors(dl,dh,u1,u2,f):
  np = int(n1*n2)
  p1 = zerofloat(np)
  p2 = zerofloat(np)
  x1 = zerofloat(2,np)
  x2 = zerofloat(2,np)
  i = 0
  for i2 in range(dh*3,n2-dh*2,dh):
    for i1 in range(dh*2,n1-dh*2,dh):
      x1[i][0] = i1*d1+f1
      x1[i][1] = (i1+dl*u1[i2][i1])*d1+f1
      x2[i][0] = i2*d2+f2
      x2[i][1] = (i2+dl*u2[i2][i1])*d2+f2
      p1[i] = x1[i][1]
      p2[i] = x2[i][1]
      i = i+1
  x1 = copy(2,i,0,0,x1)
  x2 = copy(2,i,0,0,x2)
  p1 = copy(i,0,p1)
  p2 = copy(i,0,p2)

  plot(f,xp=[x1,x2],pp=[p1,p2])
 
def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(10.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2)
  div(x,sqrt(g),y)
  return y

#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET
def jetFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.JET,a)
def jetRamp(alpha):
  return ColorMap.setAlpha(ColorMap.HUE_BLUE_TO_RED,rampfloat(0.0,alpha/256,256))
def grayRamp(alpha):
  return ColorMap.setAlpha(ColorMap.GRAY,rampfloat(0.0,alpha/256,256))

def plot(f,cmap=ColorMap.GRAY,cmin=None,cmax=None,clab=None,
        xp=None,pp=None,ap=None,png=None): 
  orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT;
  panel = PlotPanel(1,1,orientation)#PlotPanel.AxesPlacement.NONE)
  #panel.setVInterval(0.2)
  pxv = panel.addPixels(0,0,s1,s2,f);
  pxv.setColorModel(cmap)
  pxv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmin and cmax:
    pxv.setClips(cmin,cmax)
  if xp:
    ptv = panel.addPoints(xp[0],xp[1])
    ptv.setLineColor(Color.RED)
    ptv.setLineWidth(4.0)
  if pp:
    ptv = panel.addPoints(pp[0],pp[1])
    ptv.setLineStyle(PointsView.Line.NONE)
    ptv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    ptv.setMarkColor(Color.BLUE)
    ptv.setMarkSize(2.0)
  if ap:
    ptv = panel.addPoints(ap[0],ap[1])
    ptv.setLineStyle(PointsView.Line.NONE)
    ptv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    ptv.setMarkColor(Color.YELLOW)
    ptv.setMarkSize(2.0)
  cb = panel.addColorBar();
  #cb.setInterval(0.4)
  if clab:
    cb.setLabel(clab)
  panel.setColorBarWidthMinimum(130)
  moc = panel.getMosaic();
  frame = PlotFrame(panel);
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  #frame.setTitle("normal vectors")
  frame.setVisible(True);
  frame.setSize(890,760)
  #frame.setSize(n2/2,n1/2)
  frame.setSize(1080,650)
  #frame.setSize(1190,760)
  frame.setFontSize(36)
  if pngDir and png:
    frame.paintToPng(300,3.333,pngDir+png+".png")

def plot2(f,s1,s2,g=None,cmin=None,cmax=None,label=None,png=None,et=None):
  n2 = len(f)
  n1 = len(f[0])
  panel = panel2Teapot()
  #panel.setHInterval(2.0)
  #panel.setVInterval(0.2)
  #panel.setHLabel("Distance (km)")
  #panel.setVLabel("Time (s)")
  #panel.setHInterval(100.0)
  #panel.setVInterval(100.0)
  #panel.setHLabel("Pixel")
  #panel.setVLabel("Pixel")
  if label:
    panel.addColorBar(label)
  else:
    panel.addColorBar()
  panel.setColorBarWidthMinimum(130)
  pv = panel.addPixels(s1,s2,f)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setColorModel(ColorMap.GRAY)
  #pv.setClips(min(f),max(f))
  if g:
    alpha = 0.8
  else:
    g = zerofloat(s1.count,s2.count)
    alpha = 0.0
  pv = panel.addPixels(s1,s2,g)
  #pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  #pv.setColorModel(jetFillExceptMin(1.0))
  pv.setColorModel(jetRamp(1.0))
  if cmin and cmax:
    pv.setClips(cmin,cmax)
  frame2Teapot(panel,png)
def panel2Teapot():
  panel = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.NONE)
  return panel
def frame2Teapot(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  #frame.setFontSizeForPrint(8,240)
  #frame.setSize(1240,774)
  frame.setFontSizeForSlide(1.0,0.8)
  frame.setSize(1080,650)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(400,3.2,pngDir+"/"+png+".png")
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

#############################################################################
# utilities

def readImage(name):
  fileName = seismicDir+name+".dat"
  n1,n2 = s1.count,s2.count
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def writeImage(name,image):
  fileName = seismicDir+name+".dat"
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()
  return image

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
