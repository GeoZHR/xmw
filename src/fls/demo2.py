from utils2 import *

#setupForSubset("2dSub1")
#setupForSubset("bag2dTwo")
setupForSubset("bag2dOne")
#setupForSubset("seam2d")
s1,s2 = getSamplings()
n1,n2 = s1.count,s2.count

pngDir = None
pngDir = getPngDir()

fxfile = "st"
fxfile = "gx366"
fxfile = "gs" #bag2dTwo
fxfile = "gs4370" #seam2dsub2
fxfile = "gx1046" #seam2dsub2
fxfile = "st"
fxfile = "gx" #seam2d
fxfile = "gx2538s" #bag2dTwo
fxfile = "fx67cs" #bag2dOne
phfile = "phi"
dpfile = "damp"
#fxfile = "f3d267"
#fxfile = "f3d267Sub"
edfile = "edge"
f1file = "f1"
f2file = "f2"
u1file = "u1"
u2file = "u2"
#ffile = "tp73"

def main(args):
  #goGVF()
  #goSnake()
  #goSnakeReal()
  #goPaint()
  #goForce()
  #goContour()
  #goExForce()
  #goSnakeReal()
  #goExForceF3d()
  #goChannel()
  #goFlsTest()
  #goSeamTest()
  #goSeam2dP()
  #goTensor()
  #goFls()
  #goFlsP()
  #goDLSP()
  #goBagDlsp()
  #goSaltPicker()
  #goBand()
  #goTF()
  #goSaltPickerSub2()
  goBag2dOne()
  goBag2dTwo()
def goBag2dOne():
  gx = readImage(fxfile)
  p1 = []
  p2 = []
  c1 = 150
  c2 = 135
  r = 110
  for k in range(37):
    alpha = k*10*Math.PI/180
    x1 = c1+r*sin(alpha)
    x2 = c2+r*cos(alpha)
    p1.append(x1)
    p2.append(x2)
  sp = SaltPicker2()
  pa = sp.applyForInsAmp(gx)
  #xu = sp.initialBoundary(1,p1,p2)
  xu = sp.regridBoundary(1,p1,p2)
  w1 = round(n1*1.5)
  w2 = round(n2*1.5)
  plot(gx,cmin=-2,cmax=2,w1=w1,w2=w2,png="seis")
  plot(pa,cmin=0,cmax=2,xp=[p1,p2],w1=w1,w2=w2,png="initial")
  plot(pa,cmin=0,cmax=2,xp=[xu[0],xu[1]],w1=w1,w2=w2,)
  plot(pa,cmin=0,cmax=2,xu=xu,w1=w1,w2=w2,png="band")
  bs = sp.refine(50,1,10,2,xu,pa)
  xu = sp.pickNext(10,1,5,2,xu[0],xu[1],pa)
  plot(gx,cmin=-2,cmax=2,xp=[xu[0],xu[1]],w1=w1,w2=w2,png="final")
  opp = OptimalPathPicker(10,2)
  ft = opp.applyTransform(bs)
  m2,m1 = len(bs),len(bs[0])
  wht = opp.applyForWeight(ft)
  tms1 = zerofloat(m2,m1)
  tms2 = zerofloat(m2,m1)
  pik1 = opp.forwardPick(100,wht,tms1)
  pik2 = opp.backwardPick(round(pik1[m2-1]),wht,tms2)
  x2 = zerofloat(m2)
  for i2 in range(m2):
    x2[i2]=i2
  plot(bs,cmin=0.01,cmax=0.5,w1=180,w2=690,png="bandMap")
  plot(bs,cmin=0.01,cmax=0.5,xp=[pik2,x2],w1=180,w2=690,png="pik")

def goBag2dTwo():
  gx = readImage(fxfile)
  gx = gain(gx)
  c1 = [280,250,380,270,360,250,530,150, 70,600,570,790,640,720,600,280]
  c2 = [  0, 80,210,260,310,390,470,600,950,930,660,420,300,200,  0,  0]
  sp = SaltPicker2()
  pa = sp.applyForInsAmp(gx)
  pm = max(pa)/4
  for i1 in range(n1):
    pa[0   ][i1] = pm
    pa[n2-1][i1] = pm
  xu = sp.initialBoundary(1,c1,c2)
  #xu = sp.regridBoundary(1,[c1,c2])
  plot(gx,cmin=-2,cmax=2,png="seis")
  plot(gx,cmin=-2,cmax=2,pp=[c1,c2],png="initial")
  plot(pa,cmin=0,cmax=2,xp=[xu[0],xu[1]])
  plot(gx,cmin=-2,cmax=2,xu=xu,png="band")
  bs = sp.refine(95,1,40,1.5,xu,pa)
  xu = sp.pickNext(10,1,5,1,xu[0],xu[1],pa)
  plot(gx,cmin=-2,cmax=2,xp=[xu[0],xu[1]],png="final")
  opp = OptimalPathPicker(40,2)
  ft = opp.applyTransform(bs)
  m2,m1 = len(bs),len(bs[0])
  wht = opp.applyForWeight(ft)
  tms1 = zerofloat(m2,m1)
  tms2 = zerofloat(m2,m1)
  pik1 = opp.forwardPick(100,wht,tms1)
  pik2 = opp.backwardPick(round(pik1[m2-1]),wht,tms2)
  x2 = zerofloat(m2)
  for i2 in range(m2):
    x2[i2]=i2
  plot(bs,cmin=0.01,cmax=0.5,w1=200,w2=2400)
  plot(bs,cmin=0.01,cmax=0.5,xp=[pik2,x2],w1=200,w2=2400)
def goBagDlsp():
  fx = readImage(fxfile)
  fx = gain(fx)
  c1 = [280,250,380,270,360,250,530,150, 70,600,570,790,640,720,600,280]
  c2 = [  0, 80,210,260,310,390,470,600,950,930,660,420,300,200,  0,  0]
  sp = SaltPicker2()
  xu = sp.initialBoundary(1,c1,c2)
  gx = zerofloat(n1,n2)
  p1 = 0
  p2 = 0
  np = len(xu[0])
  for pk in range(np):
    i1 = round(xu[0][pk])
    i2 = round(xu[1][pk])
    i1 = min(i1,n1-1)
    i2 = min(i2,n2-1)
    i1 = max(i1,0)
    i2 = max(i2,0)
    gx[i2][i1] = 1
    p1 += i1
    p2 += i2
  c1 = [round(p1/np)+150]
  c2 = [round(p2/np)]
  r1 = [10]
  r2 = [10]
  '''
  for i2 in range(200,220,1):
    gx[i2][400] = 1
    gx[i2][420] = 1
  for i1 in range(400,420,1):
    gx[200][i1] = 1
    gx[220][i1] = 1
  c1 = [410]
  c2 = [210]
  r1 = [40]
  r2 = [40]
  '''
  mu,lamda,alpha = 0.2,5,-100
  ls = LevelSet2(mu,lamda,alpha,3,1,100)
  gx = ls.toGrayFloats(gx)
  rgf = RecursiveGaussianFilter(2)
  rgf.apply00(gx,gx)
  gs = div(1,add(1,mul(gx,gx)))
  plot(gx)
  ph = ls.initialLevelSet(n1,n2,c1,c2,r1,r2,2)
  ph0 = copy(ph)
  ls.updateLevelSetP(1.5,gs,ph)
  plot(gx,phi=ph0)
  plot(gx,phi=ph)
  plot(ph,cmin=-2,cmax=2)

def goBag2dThree():
  gx = readImage(fxfile)
  gx = gain(gx)
  c1 = [145,130,250,220,410,250,450, 330, 720,600,812,520,620,400,145]
  c2 = [  0,110,380,500,660,780,880,1128,1080,900,731,285,240,  0,  0]
  sp = SaltPicker2()
  pa = sp.applyForInsAmp(gx)
  pm = max(pa)/2
  for i1 in range(n1):
    pa[0   ][i1] = pm
    pa[n2-1][i1] = pm
  #xu = sp.initialBoundary(1,c1,c2)
  xu = sp.regridBoundary(1,[c1,c2])
  plot(gx,cmin=-2,cmax=2,png="seis")
  plot(gx,cmin=-2,cmax=2,pp=[c1,c2],png="initial")
  plot(pa,cmin=0,cmax=2,xp=[xu[0],xu[1]])
  plot(pa,cmin=0,cmax=2,xu=xu)
  bs = sp.refine(95,1,40,2,xu,pa)
  plot(gx,cmin=-2,cmax=2,xp=[xu[0],xu[1]],png="final")
  opp = OptimalPathPicker(40,2)
  ft = opp.applyTransform(bs)
  m2,m1 = len(bs),len(bs[0])
  wht = opp.applyForWeight(ft)
  tms1 = zerofloat(m2,m1)
  tms2 = zerofloat(m2,m1)
  pik1 = opp.forwardPick(100,wht,tms1)
  pik2 = opp.backwardPick(round(pik1[m2-1]),wht,tms2)
  x2 = zerofloat(m2)
  for i2 in range(m2):
    x2[i2]=i2
  plot(bs,cmin=0.01,cmax=0.5,w1=200,w2=2400)
  plot(bs,cmin=0.01,cmax=0.5,xp=[pik2,x2],w1=200,w2=2400)


def goSaltPickerSub2():
  gx = readImage(fxfile)
  c1 = [390,210,260, 30,285,270,460,280,350,315,400,450]
  c2 = [  0,260,600,820,945,860,720,440,300,180,135,  0]
  sp = SaltPicker2()
  pa = sp.applyForInsAmp(gx)
  xu = sp.initialBoundary(1,c1,c2,pa)
  plot(gx,cmin=-2,cmax=2,xp=[c1,c2],pp=[c1,c2])
  plot(pa,cmin=0,cmax=2,xp=[xu[0],xu[1]])
  plot(pa,cmin=0,cmax=2,xu=xu)
  #bs = sp.bandSample(100,1,xu,pa)
  bs = sp.refine(65,1,xu,pa)
  plot(gx,cmin=-2,cmax=2,xp=[xu[0],xu[1]])
  opp = OptimalPathPicker(40,2.0)
  ft = opp.applyTransform(bs)
  m2,m1 = len(bs),len(bs[0])
  wht = opp.applyForWeight(ft)
  tms1 = zerofloat(m2,m1)
  tms2 = zerofloat(m2,m1)
  pik1 = opp.forwardPick(100,wht,tms1)
  pik2 = opp.backwardPick(round(pik1[m2-1]),wht,tms2)
  x2 = zerofloat(m2)
  for i2 in range(m2):
    x2[i2]=i2
  plot(bs,cmin=0.01,cmax=0.5,w1=200,w2=2400)
  plot(bs,cmin=0.01,cmax=0.5,xp=[pik2,x2],w1=200,w2=2400)

def goTF():
  gx = readImage(fxfile)
  gs = gain(gx)
  plot(gx,cmin=min(gx)/2,cmax=max(gx)/2)
  plot(gs,cmin=min(gs)/2,cmax=max(gs)/2)
  '''

  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  es = zerofloat(n1,n2)
  lof = LocalOrientFilter(4,2)
  lof.applyForNormalLinear(pr,u1,u2,el)
  ets = lof.applyForTensors(pr)
  ets.setEigenvalues(0.05,1.0)
  ss = SaltScanner()
  es = ss.applyForLinear(4,ets,pr)
  plot(es,cmin=0.5,cmax=0.6)
  plot(el,cmin=0.5,cmax=0.6)
  '''

def lowpass(f3db,f):
  bf = ButterworthFilter(f3db,6,ButterworthFilter.Type.LOW_PASS)
  bf.applyForwardReverse(f,f)

def goSaltPicker():
  gx = readImage(fxfile)
  c1 = [65,270, 88, 60,134,200,282,168, 30, 250,386,260,360]
  c2 = [ 0,160,300,530,310,350,815,858,770,1000,705,485,  0]
  u1 = [-1, -1, -1, -1,  1, -1, -1, -1,  1,  -1,  1,  1,  1]
  sp = SaltPicker2()
  pa = sp.applyForInsAmp(gx)
  xu = sp.initialBoundary(1,c1,c2,u1,pa)
  plot(gx,cmin=-5,cmax=5,xp=[c1,c2],pp=[c1,c2])
  plot(pa,cmin=0,cmax=2,xp=[xu[0],xu[1]])
  plot(pa,cmin=0,cmax=2,xu=xu)
  #bs = sp.bandSample(100,1,xu,pa)
  bs = sp.refine(100,1,xu,pa)
  plot(gx,cmin=-5,cmax=5,xp=[xu[0],xu[1]])
  plot(bs,cmin=0.1,cmax=0.8)
  opp = OptimalPathPicker(40,2.5)
  ft = opp.applyTransform(bs)
  m2 = len(bs)
  m1 = len(bs[0])
  wht = opp.applyForWeight(ft)
  tms1 = zerofloat(m2,m1)
  tms2 = zerofloat(m2,m1)
  pik1 = opp.forwardPick(100,wht,tms1)
  pik2 = opp.backwardPick(round(pik1[m2-1]),wht,tms2)
  x2 = zerofloat(m2)
  for i2 in range(m2):
    x2[i2]=i2
  plot(bs,cmin=0.1,cmax=0.8,xp=[pik2,x2])


  '''
  sb = zerofloat(n1,n2)
  np = len(xu[0])
  for ip in range(np):
    i1 = round(xu[0][ip])
    i2 = round(xu[1][ip])
    sb[i2][i1] = 1
  c1 = [250]
  c2 = [300]
  r1 = [2]
  r2 = [2]
  mu,lamda,alpha = 0.2,5,-3
  ls = LevelSet2(mu,lamda,alpha,3,1,2500)
  sb = ls.toGrayFloats(sb)
  rgf = RecursiveGaussianFilter(4)
  rgf.apply00(sb,sb)
  gs = div(1,add(1,mul(sb,sb)))
  plot(sb)
  ph = ls.initialLevelSet(n1,n2,c1,c2,r1,r2,2)
  ph0 = copy(ph)
  ls.updateLevelSetP(1.5,gs,ph)
  plot(gx,phi=ph0)
  plot(gx,phi=ph)
  plot(ph,cmin=-2,cmax=2)
  '''


def goBand():
  gx = readImage(fxfile)
  dp = readImage(dpfile)
  ph = readImage(phfile)
  pa = zerofloat(n1,n2)
  ls = LevelSet2(0.2,2,10,3,1,1500)
  #gx = gain(gx)
  ls.applyForInsAmp(gx,pa)
  plot(pa,cmin=0,cmax=2)
  #fb = ls.bandSample(50,ph,dp)
  r,d=20,1.0
  ps,bs = ls.refine(r,d,ph,pa)
  plot(gx,phi=ph,png="salt")
  plot(gx,xp=ps[0],png="saltRefined")
  opp = OptimalPathPicker(20,1.0)
  ft = opp.applyTransform(bs[0])
  m2 = len(bs[0])
  m1 = len(bs[0][0])
  wht = opp.applyForWeight(ft)
  tms1 = zerofloat(m2,m1)
  tms2 = zerofloat(m2,m1)
  pik1 = opp.forwardPick(r,wht,tms1)
  pik2 = opp.backwardPick(round(pik1[m2-1]),wht,tms2)
  x2 = zerofloat(m2)
  for i2 in range(m2):
    x2[i2]=i2
  plot(bs[0],cmin=0.1,cmax=0.8,xp=[pik2,x2])

def goDLSP():
  gx = readImage(fxfile)
  c1 = [140]
  c2 = [200]
  r1 = [2]
  r2 = [2]
  mu,lamda,alpha = 0.2,5,-3
  ls = LevelSet2(mu,lamda,alpha,3,1,450)
  gx = ls.toGrayFloats(gx)
  rgf = RecursiveGaussianFilter(2)
  rgf.apply00(gx,gx)
  gs = div(10,add(1,mul(gx,gx)))
  plot(gx)
  ph = ls.initialLevelSet(n1,n2,c1,c2,r1,r2,2)
  ph0 = copy(ph)
  ls.updateLevelSetP(1.5,gs,ph)
  plot(gx,phi=ph0)
  plot(gx,phi=ph)
  plot(ph,cmin=-2,cmax=2)

def goDLS():
  gx = readImage(fxfile)
  p2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  #ls = LevelSet2(0.2,5,-3,3,1,2500)
  ls = LevelSet2(0.2,2,50,3,1,1500)
  lsf = LocalSlopeFinder(8,2,5)
  lsf.findSlopes(gx,p2,el)
  lsf = LocalSlopeFinder(2,1,5)
  lsf.findSlopes(gx,p2)
  d1 = fillfloat(0,256)
  d2 = fillfloat(0,256)
  gg = ls.toGrayIntegers(gx)
  pg = ls.toGrayIntegers(p2)
  ls.density(0.6,el,gg,d1,d2)
  print sum(d1)
  print sum(d2)
  c1 = Sampling(256)
  plot1(c1,d1,d2)
  '''
  c1 = [250,250]
  c2 = [250,910]
  r1 = [20,20]
  r2 = [30,30]
  ph = ls.initialLevelSet(n1,n2,c1,c2,r1,r2,2)
  ph0 = copy(ph)
  ls.updateLevelSetPK(1.5,el,gx,p2,ph)
  mu,lamda,alpha = 0.2,10,1
  ls = LevelSet2(mu,lamda,alpha,3,1,500)
  ls.updateLevelSetPK(1.5,el,gx,p2,ph)
  writeImage(phfile,ph)
  ph = readImage(phfile)
  dp = ls.delta(ph)
  ddip = ls.density(0.6,el,p2)
  damp = ls.density(0.6,el,gx)
  writeImage(dpfile,damp)
  dp = readImage(dpfile)
  mk = ls.getMark()
  plot(mk)
  print dp[400][220]
  plot(gx,phi=ph0)
  plot(gx,phi=ph)
  plot(ph,cmin=-2,cmax=2)
  plot(dp)
  plot(ddip)
  plot(damp,cmin=-1,cmax=2)
  print min(ddip)
  print min(damp)
  '''
def goTensor():
  gx = readImage(fxfile)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lof = LocalOrientFilter(1,1)
  st = NonlinearStructureTensor(0,1,0.1)
  #st.applyForNormalLinear(gx,u1,u2,el)
  lof.applyForNormalLinear(gx,u1,u2,el)
  lsf = LocalSlopeFinder(1,1,5)
  lsf.findSlopes(gx,u1)
  #plot(gx)
  plot(u1)
  print min(gx)
  print max(gx)

def goFlsP():
  gx = readImage(fxfile)
  p2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  ft = zerofloat(n1,n2)
  lsf = LocalSlopeFinder(8,2,5)
  lsf.findSlopes(gx,p2,el)
  lsf = LocalSlopeFinder(2,1,5)
  lsf.findSlopes(gx,p2)
  x1 = [220]
  x2 = [900]
  rs = [ 5]
  fls = FastLevelSet2P(n1,n2,x1,x2,rs)
  fls.setIterations(550,6,2)
  dp = fls.applyForDensity(0.4,el,[gx,p2])
  plot(dp,cmin=-1,cmax=1)
  xss = fls.applySegments(9,6,dp,ft)
  plot(gx,xs=[xss[0]],png=None)#pngName+str(k))
  plot(gx,xs=[xss[1]],png=None)#pngName+str(k))
  plot(gx,xs=[xss[1]],png=None)#pngName+str(k))

def goFls():
  gx = readImage(fxfile)
  p2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lsf = LocalSlopeFinder(8,2,5)
  lsf.findSlopes(gx,p2,el)
  lsf = LocalSlopeFinder(2,1,5)
  lsf.findSlopes(gx,p2)
  x1 = [220]
  x2 = [900]
  rs = [ 5]
  x1 = [110]
  x2 = [450]
  rs = [  5]
  gxs = FastLevelSets2.downSample(1,1,gx) 
  els = FastLevelSets2.downSample(1,1,el) 
  p2s = FastLevelSets2.downSample(1,1,p2) 
  m2 = len(gxs)
  m1 = len(gxs[0])
  fls = FastLevelSets2(m1,m2,x1,x2,rs)
  phi = zerofloat(m1,m2)
  plot(gxs)
  fls.setIterations(650,6,2)
  xss = fls.applySegments(9,5,els,gxs,p2s,phi)
  writeImage(phfile,phi)
  rgf = RecursiveGaussianFilter(4)
  rgf.apply00(phi,phi)
  plot(phi)
  plot(gxs,pp=xss[1][0],png=None)#pngName+str(k))
  plot(gxs,phi=phi,png=None)#pngName+str(k))

def goSeamTest():
  gx = readImage(fxfile)
  print min(gx)
  print max(gx)
  rgf = RecursiveGaussianFilter(1)
  rgf.apply00(gx,gx)
  fx = abs(gx)
  ft = copy(gx)
  x1 = [ 80,100,250,400,600]
  x2 = [780,400,200,700,700]
  rs = [ 5 ,  5,  5,  5,  5]
  fls = FastLevelSets2(n1,n2,x1,x2,rs)
  fls.setIterations(500,4,1)
  xss = fls.applySegments(9,3,fx,ft)
  plot(fx)
  plot(ft)
  plot(gx,xs=xss[0],png=None)#pngName+str(k))
  plot(gx,xs=xss[1],png=None)#pngName+str(k))
def goFlsTest():
  fx = readImage(fxfile)
  gx = zerofloat(n1,n2)
  for i2 in range (n2):
    fx[i2][n1-1] = 0.8
  rgf = RecursiveGaussianFilter(5)
  rgf.apply00(fx,gx)
  #gvf = GradientVectorFlow()
  #gvf.setSmoothing(6)
  #gvf.setScale(0.1)
  #g1,g2,gs = gvf.applyForGradient(1,fx)
  #u1,u2 = gvf.applyForGVF(g1,g2,gs)
  #fls = FastLevelSet2(120,210,10,u1,u2)
  fls = FastLevelSet2(120,100,20,gx)
  fls.setIterations(200,10,2)
  ph1 = fls.getPhi()
  fls.updateLevelSet(9,3)
  k1,k2=fls.getLout()
  ph2 = fls.getPhi()
  plot(fx)
  #plot(fx,v1=u1,v2=u2)
  plot(ph1)
  plot(ph2)
  plot(fx,ap=[k1,k2],png=None)#pngName+str(k))
def goTest():
  fx = readImage(fxfile)
  gvf = GradientVectorFlow()
  gvf.setSmoothing(6)
  gvf.setScale(0.1)
  for i2 in range (n2):
    fx[i2][n1-1] = 1
  g1,g2,gs = gvf.applyForGradient(1,fx)
  u1,u2 = gvf.applyForGVF(g1,g2,gs)
  ac2 = ActiveContour2(n1,n2,140,205,20)
  snake = ac2.getSnake()
  for k in range(20):
    k1,k2=[],[]
    x1 = snake.getArrayX1()
    x2 = snake.getArrayX2()
    for j in range(len(x1)):
      x1i = x1[j]
      x2i = x2[j]
      if(x1i>1 and x1i<n1-2):
       k1.append(x1i)
       k2.append(x2i)
    ac2.updateSnake(25,u1,u2)
    plot(fx,ap=[k1,k2],png=None)#pngName+str(k))
  plot(fx)
  plot(fx,v1=u1,v2=u2)

def goChannel():
  fx = readImage(fxfile)
  fx = div(fx,10000)
  el = zerofloat(n1,n2)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  g1 = zerofloat(n1,n2)
  g2 = zerofloat(n1,n2)
  lof = ChannelEnhanceFilter(4,4)
  ets = lof.applyForTensors(fx)
  lof.applyForNormal(fx,u1,u2)
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply10(fx,g1)
  rgf.apply01(fx,g2)
  g1 = mul(g1,g1)
  g2 = mul(g2,g2)
  sc = sqrt(add(g1,g2))
  sc = sub(sc,min(sc))
  sc = div(sc,max(sc))
  ets.setEigenvalues(0.01,1.0)
  etp = lof.applyForLinear(8,ets,sc,fx,el)
  plot(fx)
  plot(sc)
  plot(el)

  
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

'''
def goTest():
  f = goFakeImage()
  x = [20,30,50,60,80,85,85,85,85,70,65,50,40,30,20,20,20]
  y = [20,20,20,20,20,20,25,35,50,55,50,50,50,50,50,30,20]
  acs = ActiveSnake(x,y)
  snake = acs.releaseSnake(s1,s2)
  xs = snake.getArrayX1()
  ys = snake.getArrayX2()
  plot(f,ap=[y,x])
  plot(f,ap=[ys,xs])
'''
  

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
  ref = RecursiveExponentialFilter(50.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2)
  div(x,sqrt(g),y)
  return y


def smooth(sig1,sig2,sigs,x):
  lof = LocalOrientFilter(sig1,sig2)
  ets = lof.applyForTensors(x)
  ets.setEigenvalues(0.001,1.0)
  lsf = LocalSmoothingFilter()
  g = zerofloat(n1,n2)
  lsf.apply(ets,sigs,x,g)
  return g

#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET
def plot1(s1,y1,y2,hlabel="Values",vlabel="Probability",png=None):
  sp = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
  pv1 = sp.addPoints(s1,y1)
  pv1.setLineColor(Color.RED)
  pv1 = sp.addPoints(s1,y2)
  pv1.setLineColor(Color.BLUE)

  #sp.setVLimits(0.1,1.1)
  sp.setSize(800,800)
  sp.setHLabel(hlabel)
  sp.setVLabel(vlabel)
  if png and pngDir:
    sp.paintToPng(300,7.0,pngDir+png+".png")


def plot(f,xp=None,pp=None,xs=None,xu=None,phi=None,v1=None,v2=None,
        cmin=None,cmax=None,w1=None,w2=None,clab=None,png=None): 
  orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT;
  panel = PlotPanel(1,1,orientation);
  #panel.setVInterval(0.2)
  n2 = len(f)
  n1 = len(f[0])
  s2 = Sampling(n2)
  s1 = Sampling(n1)
  panel.setHLabel("Inline (traces)")
  panel.setVLabel("Time (samples)")

  panel.setHLimits(0,0,n2-1)
  panel.setVLimits(0,0,n1-1)
  pxv = panel.addPixels(0,0,s1,s2,f);
  #pxv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pxv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pxv.setColorModel(ColorMap.GRAY)
  if cmin and cmax:
    pxv.setClips(cmin,cmax)
  else:
    pxv.setClips(min(f),max(f))
  #panel.setTitle("normal vectors")
  if phi:
    cv = panel.addContours(phi)
    cv.setContours([0])
    cv.setLineColor(Color.RED)
    cv.setLineWidth(1.0)
  if xp:
    ptv = panel.addPoints(0,0,xp[0],xp[1])
    ptv.setLineColor(Color.RED)
    ptv.setLineWidth(3.0)
  if xu:
    br = 50
    np = len(xu[0])
    for ip in range(np):
      x1c = xu[0][ip]
      x2c = xu[1][ip]
      u1c = xu[2][ip]
      u2c = xu[3][ip]
      x1m = x1c-u1c*br
      x2m = x2c-u2c*br
      x1p = x1c+u1c*br
      x2p = x2c+u2c*br
      x1s = [x1m,x1c,x1p]
      x2s = [x2m,x2c,x2p]
      ptv2 = panel.addPoints(x1s,x2s)
      ptv2.setLineColor(Color.YELLOW)
      ptv2.setLineWidth(1.0)
    ptv1 = panel.addPoints(xu[0],xu[1])
    ptv1.setLineColor(Color.RED)
    ptv1.setLineWidth(3.0)

  if pp:
    ptvl = panel.addPoints(0,0,pp[0],pp[1])
    ptvl.setLineColor(Color.RED)
    ptvl.setLineWidth(2.0)
    ptvp = panel.addPoints(0,0,pp[0],pp[1])
    ptvp.setLineStyle(PointsView.Line.NONE)
    #ptv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    ptvp.setMarkStyle(PointsView.Mark.CROSS)
    ptvp.setMarkColor(Color.RED)
    ptvp.setMarkSize(6.0)
    ptvp.setLineWidth(3.0)
  if xs:
    for ip in range(len(xs)):
      ptv = panel.addPoints(xs[ip][0],xs[ip][1])
      ptv.setLineStyle(PointsView.Line.NONE)
      ptv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
      ptv.setMarkColor(Color.RED)
      ptv.setMarkSize(2.0)
  if(clab):
    cb = panel.addColorBar();
    cb.setLabel(clab)
  panel.setColorBarWidthMinimum(130)
  moc = panel.getMosaic();
  frame = PlotFrame(panel);
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setVisible(True);
  if w1 and w2:
    frame.setSize(w2,w1)
  else:
    frame.setSize(round(n2*0.8),round(n1*0.8))
  #frame.setSize(1190,760)
  frame.setFontSize(14)
  if pngDir and png:
    frame.paintToPng(300,3.333,pngDir+png+".png")

def plot2(f,x1,x2,s,s1,s2,g=None,label=None,png=None,et=None):
  n1 = len(s[0])
  n2 = len(s)
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
  panel.setColorBarWidthMinimum(180)
  pv = panel.addPixels(s1,s2,s)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setColorModel(ColorMap.GRAY)
  if cmin and cmax:
    pv.setClips(cmin,cmax)
  else:
    pv.setClips(min(s),max(s))
  if g:
    alpha = 0.5
  else:
    g = zerofloat(s1.count,s2.count)
    alpha = 0.0
  pv = panel.addPixels(s1,s2,g)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setColorModel(ColorMap.getJet(alpha))
  if label and label[0]=="T":
    #pv.setClips(0.0,1000.0)
    pv.setClips(0.0,400.0)
  else:
    pv.setClips(0.0,1.0)
  cmap = pv.getColorMap()
  if et:
    tv = TensorsView(s1,s2,et)
    tv.setOrientation(TensorsView.Orientation.X1DOWN_X2RIGHT)
    tv.setLineColor(Color.YELLOW)
    tv.setLineWidth(3.0)
    tv.setScale(2.0)
    panel.getTile(0,0).addTiledView(tv)
  else:
    fs,x1s,x2s = makePointSets(cmap,f,x1,x2)
    for i in range(len(fs)):
      color = cmap.getColor(fs[i][0])
      color = Color(color.red,color.green,color.blue)
      pv = panel.addPoints(x1s[i],x2s[i])
      pv.setLineStyle(PointsView.Line.NONE)
      pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
      pv.setMarkSize(10)
      pv.setMarkColor(color)
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
  frame.setSize(880,700)
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
