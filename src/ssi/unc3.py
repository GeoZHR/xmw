import sys

from java.awt import *
from java.awt.image import *
from java.io import *
from java.nio import *
from java.lang import *
from javax.swing import *
from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.ogl.Gl import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *


from ssi import *
#from unc import *

#############################################################################
# real image
n1,n2,n3= 120,200,260
n3h,n2h,n1h = round(n3/2),round(n2/2),round(n1/2)
d1,d2,d3= 0.004,0.025,0.025
f1,f2,f3= 110.0*d1,750.0*d2,10.0*d3
k1f,k2f,k3f=92,166,88
gmin,gmax=-1.5,1.5
d1,d2,d3= 1.0,1.0,1.0
f1,f2,f3= 0.0,0.0,0.0
background = Color.WHITE
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)
s3 = Sampling(n3,d3,f3)
#azimuth,elevation=150,35
azimuth,elevation=146,24
azimuth,elevation=230,14

threshold = 0.8
sigma1s,sigma2s=1.0,2.0
sigma1u,sigma2u = 0.75,max(n2h,n3h)*2
niter = sigma2u+200

pngDir = "../../../png/ssi/"
seismicDir = "../../../data/ssi/"
gxfile = "f3dunc" # input image
fgfile = "fg" # flattened image
ulfile = "ul" # unconformity likelihood image
ulhfile = "ulh" # unconformity likelihood image
ulhsfile = "ulhs" # unconformity likelihood image
ultfile = "ult" # thinned unconformity likelihood image
ulthfile = "ulth" # thinned unconformity likelihood image
sfsfile = "sfs" # thinned unconformity likelihood image
p2file = "p2"
p3file = "p3"
rgtfile = "p3"
#############################################################################

def main(args):
  #displayUnc3D()
  #goUnc()
  #goSurf()
  goFlatten()
def goUnc():
  ul   = zerofloat(n1,n2,n3)
  ulh  = zerofloat(n1h,n2h,n3h)
  ulhs = zerofloat(n1h,n2h,n3h)
  ulth = zerofloat(n1h,n2h,n3h)
  unc = Unconformity()
  gx = readImage(gxfile)
  gx = gain(30,gx)
  '''
  unc.setForLof(sigma1s,sigma2s)
  unc.setForLofu(sigma1u,sigma2u,niter)
  unc.likelihood(gx,ulh)
  unc.interp(2,2,ulh,ul)
  writeImage(ulfile,ul)
  writeImage(ulhfile,ulh)
  '''
  ulh=readImage3D(n1h,n2h,n3h,ulhfile)
  ref = RecursiveExponentialFilter(2.0)
  ref.apply(ulh,ulhs)
  writeImage(ulhsfile,ulhs)
  unc.thin(threshold,gx,ulhs,ulth)
  writeImage(ulthfile,ulth)
  ul = readImage(ulfile)
  displayUnc3D()
def goSurf():
  unc = Unconformity()
  gx  = readImage(gxfile)
  ul  = readImage(ulfile)
  '''
  n3h = round(n3/2)
  n2h = round(n2/2)
  n1h = round(n1/2)
  ulhs = readImage3D(n1h,n2h,n3h,ulhsfile)
  ulth = readImage3D(n1h,n2h,n3h,ulthfile)
  sfs = unc.surfer(n2,n3,threshold,300,ulth,ulhs)
  unc.surfaceUpdate(2.0,2.0,gx,sfs)
  writeImage(sfsfile,sfs)
  ult = unc.thinFromSurf(sfs,sub(1.0,ul))
  writeImage(ultfile,ult)
  '''
  sfs = readImage3D(n2,n3,2,sfsfile)
  add(sfs[0],-1.0,sfs[0]);
  add(sfs[1], 1.0,sfs[1]);
  sfs = unc.trianglesForSurface(ul,sfs)
  world = World()
  gx = gain(20,gx)
  addImageToWorld(s1,s2,s3,world,gx,cmin=-3.2,cmax=3.7)
  ns = len(sfs[0])
  for i in range(ns):
    np = len(sfs[0][i])
    if np>3:
      tg = TriangleGroup(True,sfs[0][i],sfs[1][i])
      world.addChild(tg)
  makeFrame(world,png="uncSurfaces")

  #plot3f(gx,gmap=gray,a=ult,alab="Unconformity likelihood",aint=0.2,png="uncThin")
  #plot3f(gx,gmap=gray,a=sub(1.0,ul),alab="Unconformity likelihood",aint=0.2,png="unc")

def goFlatten():
  '''
  ns = 2
  sfs = readImage3D(n2,n3,ns,sfsfile)
  fh = FlattenHelper(s1)
  gx = readImage(gxfile)
  p2 = zerofloat(n1,n2,n3)
  p3 = zerofloat(n1,n2,n3)
  ep = zerofloat(n1,n2,n3)
  u1 = zerofloat(n1,n2,n3)
  u2 = zerofloat(n1,n2,n3)
  u3 = zerofloat(n1,n2,n3)
  lsf = LocalSlopeFinder(2.0,2.0,4.0)
  lsf.findSlopes(gx,p2,p3,ep)
  ep = pow(ep,8.0)
  ucc = fh.constraintsFromUnconformities(p2,p3,ep,sfs)
  goSlopes(ucc,2.0,8.0)
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  fl3 = Flattener3()
  fl3.setSmoothings(24.0,24.0)
  fl3.setIterations(0.01,200)
  s = zerofloat(n1,n2,n3)
  fl3.computeShifts(p2,p3,ep,ucc,s)
  fh.shiftsToRgt(s)
  gx = gain(20,gx)
  fg = fh.flattenByRgt(s,gx)
  writeImage(fgfile,fg)
  writeImage(rgtfile,s)
  '''
  fg = readImage(fgfile)
  rgt = readImage(rgtfile)
  d1= 0.004
  f1= 110.0*d1
  rgt = mul(rgt,d1)
  rgt = add(rgt,f1)
  '''
  plot3f(fg,gmap=gray,gla1="Relative geologic time",
          gmin=-3.2,gmax=3.7,glab="Amplitude",aint=0.2,png="flatten")
  plot3f(rgt,gmap=jet,glab="Relative geologic time",aint=0.2,gmin=0.43,gmax=0.98,png="rgt")
  '''
  world = World()
  addImageToWorld(s1,s2,s3,world,fg,gmin=-3.2,gmax=3.7)
  makeFrame(world,png="flatten3D") 
  world = World()
  addImageToWorld(s1,s2,s3,world,rgt,cmap=jet,gmin=0.43,gmax=0.98)
  makeFrame(world,png="rgt3D") 


def goSlopes(uc,sigma1,sigma2):
  sg  = fillfloat(1.0,n1,n2,n3)
  sgp = fillfloat(1.0,n1,n2,n3)
  u1 = zerofloat(n1,n2,n3)
  u2 = zerofloat(n1,n2,n3)
  u3 = zerofloat(n1,n2,n3)
  p2 = zerofloat(n1,n2,n3)
  p3 = zerofloat(n1,n2,n3)
  ep = zerofloat(n1,n2,n3)
  fh = FlattenHelper(s1)
  fi,fn=1.0,1.0
  fh.setWeightsFromUnconformities(sgp,uc,fi,fn)
  wlof = WeightedLocalOrientFilter(sigma1,sigma2)
  wlof.setScaleFactor(sg,sgp)
  wlof.setGradientSmoothing(2,2)
  gx = readImage(gxfile)
  gx = gain(20,gx)
  wlof.applyForNormalPlanar(gx,u1,u2,u3,ep)
  nsp = NormalsToSlopes()
  nsp.setLimits(4.0)
  nsp.apply(u1,u2,u3,p2,p3)
  writeImage(p2file,p2)
  writeImage(p3file,p3)
 
def realImage():
  f  = readImage3D(n1,n2,n3,"f3dunc")
  f  = gain(f)
  u  = zerofloat(n1,n2,n3)
  us = zerofloat(n1,n2,n3)
  ut1= fillfloat(0.0,n1,n2,n3)
  ut2= fillfloat(0.0,n1,n2,n3)
  ud = Unconformity()
  ud.setForLofu(0.85,200.0,0.0)
  lof = LocalOrientFilter(2.0,6.0)
  et = lof.applyForTensors(f)
  #lsf = LaplacianSmoother(2.0)
  #lsf.apply(LaplacianSmoother.Direction3.VW,et,f,f)
  ud.uncLikelihood(et,f,u)
  ref = RecursiveExponentialFilter(1.0)
  ref.apply(u,us)
  u = sub(u,min(u))
  u = div(u,max(u))
  us= sub(us,min(us))
  us= div(us,max(us))
  us = readImage3D(n1,n2,n3,"us")
  ud.thin(1,us,ut1,ut2)
  writeImage("us",us)
  writeImage("ut1",ut1)
  writeImage("ut2",ut2)
  ut1 = readImage3D(n1,n2,n3,"ut1")
  us= readImage3D(n1,n2,n3,"us")
  #ud.thin(1,us,ut1,ut2)
  #writeImage("ut1",ut1)
  world = World()
  addImageToWorld(s1,s2,s3,world,f)
  #addImageToWorld(s1,s2,s3,world,u,cmap=jet)
  addImageToWorld(s1,s2,s3,world,ut1,cmap=jet)
  addImageToWorld(s1,s2,s3,world,ut2,cmap=jet)
  makeFrame(world)
  #displayUnc3D(f,us)
  displaySurfaces(f,us,ut1)
def displayUnc3D():
  ul = readImage(ulfile)
  gx = readImage(gxfile)
  gx = gain(20,gx)
  world = World()
  addImageToWorld(s1,s2,s3,world,gx,u=sub(1.0,ul))
  makeFrame(world,png="unc3D")
def activeSurface(f,u,ut1):
  ud = UnconformityDetection()
  surf = ud.applyForSurface(f,ut1,u)
  world = World()
  addImageToWorld(s1,s2,s3,world,f,cmin=-2.5,cmax=2.5)
  ns = len(surf[0])
  for i in range(ns):
    np = len(surf[0][i])
    if np>3:
      tg = TriangleGroup(True,surf[0][i],surf[1][i])
      world.addChild(tg)
  makeFrame(world,png="uncSurfacesO")
 
def displaySurfaces(f,u,ut1):
  ud = Unconformity()
  surf = ud.uncSurfer(0.85,ut1,u)
  world = World()
  '''
  dn = 30
  f = copy(n1-dn,n2,n3,dn,0,0,f)
  d1,d2,d3=1.0,1.0,1.0
  f1,f2,f3=dn,0.0,0.0
  s1 = Sampling(n1-dn,d1,f1)
  s2 = Sampling(n2,d2,f2)
  s3 = Sampling(n3,d3,f3)
  '''
  addImageToWorld(s1,s2,s3,world,f,cmin=-2.5,cmax=2.5)
  ns = len(surf[0])
  for i in range(ns):
    np = len(surf[0][i])
    if np>3:
      tg = TriangleGroup(True,surf[0][i],surf[1][i])
    #qg = QuadGroup(True,surf[0][i],surf[1][i])
      world.addChild(tg)
  makeFrame(world,png="uncSurfaces")

def unweightedLof(f):
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lof = LocalOrientFilter(4.0,2.0)
  lof.applyForNormalLinear(f,u1,u2,el)
  '''
  el = pow(el,8.0)
  flatten(f,el,u1,u2)
  '''
  plot(f,"input",cmap=gray)
  plot(u1,"u1: use RGF (sigma1=4.0, sigma2=2.0)",cmap=jet,cmin=0.0,cmax=1.0)
  plot(u2,"u2: use RGF (sigma1=4.0, sigma2=2.0)",cmap=jet,cmin=0.0,cmax=1.0)
  plotVectors(f,u1,u2,"normal vectors: use RGF (sigma1=4.0, sigma2=2.0)")

def weightedLof():
  f = readImage3D(n1,n2,n3,"f3dunc")
  ut1 = readImage3D(n1,n2,n3,"ut1")
  ut2 = readImage3D(n1,n2,n3,"ut2")
  u1 = zerofloat(n1,n2,n3)
  u2 = zerofloat(n1,n2,n3)
  u3 = zerofloat(n1,n2,n3)
  p2 = zerofloat(n1,n2,n3)
  p3 = zerofloat(n1,n2,n3)
  ep = zerofloat(n1,n2,n3)
  sg = pow(ut2,0)
  sgp= pow(ut2,100)
  wlof = WeightedLocalOrientFilter(12.0,6.0)
  wlof.setScaleFactor(sg,sgp)
  wlof.setGradientSmoothing(2,2)
  wlof.applyForNormalPlanar(f,u1,u2,u3,ep)
  nsp = NormalsToSlopes()
  nsp.setLimits(10.0)
  nsp.apply(u1,u2,u3,p2,p3)
  display(u1)
  display(u2)
  display(p2)
  display(p3)
  writeImage("p2",p2)
  writeImage("p3",p3)
  writeImage("ep",ep)

def flatten():
  f = readImage3D(n1,n2,n3,"f3dunc")
  f = gain(f)
  ut = readImage3D(n1,n2,n3,"ut1")
  ut2= readImage3D(n1,n2,n3,"ut2")
  us = readImage3D(n1,n2,n3,"us")
  p2 = readImage3D(n1,n2,n3,"p2")
  p3 = readImage3D(n1,n2,n3,"p3")
  ep = readImage3D(n1,n2,n3,"ep")
  '''
  w1 = fillfloat(0.01,n1,n2,n3)
  wp = fillfloat(1.0,n1,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(n1):
        if(ut2[i3][i2][i1]<1.0):
          wp[i3][i2][i1] = 0.0
  wp = mul(wp,ep)
  w1 = mul(wp,w1)
  fl = Flattener3unc()
  fl.setIterations(0.01,100)
  fl.setSmoothings(20,20,20)
  mp = fl.getMappingsFromSlopes(s1,s2,s3,p2,p3,wp,w1)
  g = mp.flatten(f)
  r = mp.u1
  '''
  '''
  writeImage("g",g)
  writeImage("rgt",r)
  world = World()
  addImageToWorld(s1,s2,s3,world,f)
  addImageToWorld(s1,s2,s3,world,g)
  addImageToWorld(s1,s2,s3,world,r,cmap=jet)
  makeFrame(world)
  '''
  g = readImage3D(n1,n2,n3,"g")
  r = readImage3D(n1,n2,n3,"rgt")
  d1= 0.004
  f1= 110.0*d1
  r = mul(r,d1)
  r = add(r,f1)
  world = World()
  addImageToWorld(s1,s2,s3,world,f)
  makeFrame(world,png="unflatten3D") 
  world = World()
  addImageToWorld(s1,s2,s3,world,r,cmap=jet)
  makeFrame(world)#,png="rgt3D") 
  g = gain(g)
  plot3f(g,gmap=gray,gla1="RGT",glab="Amplitude",gint=0.5,png="unflatten")
  plot3f(r,gnea=True,gmap=jet, glab="RGT",gint=0.1,png="rgt")
  plot3f(f,gmap=gray,a=sub(1.0,us),alab="Unconformity likelihood",aint=0.2,png="unc")
  plot3f(f,gmap=gray,a=sub(1.0,ut),alab="Unconformity likelihood",aint=0.2,png="uncThin")
  #display(sub(1.0,ut))

def gain(sigma,x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(sigma)
  ref.apply1(g,g)
  y = like(x)
  y=div(x,sqrt(g))
  return y

def like(x):
  n3 = len(x)
  n2 = len(x[0])
  n1 = len(x[0][0])
  return zerofloat(n1,n2,n3)
##################################################################
# plots
jet = ColorMap.JET
gray = ColorMap.GRAY
def plot(x,title,cbar=None,cmap=jet,px1=None,px2=None,cmin=0,cmax=0):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(s1,s2,x)
  #pv.setInterpolation(PixelsView.Interpolation.NEAREST);
  pv.setColorModel(cmap);
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if px1:
    ptv = sp.addPoints(px1,px2);
    ptv.setLineStyle(PointsView.Line.NONE);
    ptv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
    #ptv.setMarkColor(Color.red)
    ptv.setMarkSize(3.0);
  sp.setSize(n2*5+80,n1*5);
  sp.addTitle(title)
  cb=sp.addColorBar(cbar)
  cb.setWidthMinimum(80)
  sp.setHLabel("x2")
  sp.setVLabel("x1")
  sp.setFontSize(24)
  sp.paintToPng(720,3.3,pngDir+title+".png")

def plot3f(g,gnea=False,gmap=gray,gla1="Time (s)",glab=None,gint=None,
           a=None,alab=None,aint=None,gmin=None,gmax=None,png=None):
  pp = PlotPanelPixels3(
    PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT,
    PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM,
    s1,s2,s3,g)
  pp.setSlices(k1f,k2f,k3f)
  pp.setLabel1(gla1)
  pp.setLabel2("Inline (km)")
  pp.setLabel3("Crossline (km)")
  pp.mosaic.setHeightElastic(1,100)
  pp.setColorModel(gmap)
  sg12 = slice12(k3f,g)
  sg13 = slice13(k2f,g)
  sg23 = slice23(k1f,g)
  if gmin:
    pp.setClips(gmin,gmax)
  else :
    gmin = min([sg12,sg13,sg23])
    gmax = max([sg12,sg13,sg23])
    pp.setClips(gmin,gmax)
  pp.setClips(gmin,gmax)
  if (gnea):
    pp.setInterpolation(PixelsView.Interpolation.NEAREST)
  if a:
    pp.setLineColor(Color.BLUE)
    cb = pp.addColorBar(alab)
    if aint:
      cb.setInterval(aint)
  else:
    pp.setLineColor(Color.BLUE)
    cb = pp.addColorBar(glab)
    if gint:
      cb.setInterval(gint)
  pp.setInterval1(0.1)
  pp.setInterval2(2.0)
  pp.setInterval3(2.0)
  if a:
    sa12 = slice12(k3f,a)
    sa13 = slice13(k2f,a)
    sa23 = slice23(k1f,a)
    pv12 = PixelsView(s1,s2,sa12)
    pv12.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv12.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv13 = PixelsView(s1,s3,sa13)
    pv13.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv13.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv23 = PixelsView(s2,s3,sa23)
    pv23.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
    pv23.setInterpolation(PixelsView.Interpolation.NEAREST)
    amin = min([sa12,sa13,sa23])
    amax = max([sa12,sa13,sa23])
    for pv in [pv12,pv13,pv23]:
      pv.setColorModel(ColorMap.getJet(0.8))
      #if amin!=amax:
      pv.setClips(amin,amax)
      updateColorModel(pv,1.0)
    pp.pixelsView12.tile.addTiledView(pv12)
    pp.pixelsView13.tile.addTiledView(pv13)
    pp.pixelsView23.tile.addTiledView(pv23)
  pf = PlotFrame(pp)
  pf.setBackground(background)
  pp.setColorBarWidthMinimum(60)
  #pf.setFontSizeForSlide(1.0,0.8)
  pf.setFontSize(16)
  pf.setSize(850+60,670)
  pf.setVisible(True)
  if png and pngDir:
    png = pngDir+png+str(k1f)
    pf.paintToPng(720,3.3,png+".png")
def updateColorModel(pv,alpha):
    n = 256
    r = zerobyte(n)
    g = zerobyte(n)
    b = zerobyte(n)
    a = zerobyte(n)
    icm = pv.getColorModel()
    icm.getReds(r)
    icm.getGreens(g)
    icm.getBlues(b)
    for i in range(n):
      ai = int(255.0*alpha*i/n)
      #if ai>127:
      if ai>127:
        ai -= 256
      a[i] = ai
    #if alpha<1.0:
      #r[n/2] = r[n/2-1] = -1
      #g[n/2] = g[n/2-1] = -1
      #b[n/2] = b[n/2-1] = -1
      #a[n/2  ] = a[n/2-1] = 0
      #a[n/2+1] = a[n/2-2] = 0
      #a[n/2+2] = a[n/2-3] = 0
      #a[0] = a[1] = 0
    icm = IndexColorModel(8,n,r,g,b,a)
    pv.setColorModel(icm)


def display(u):
  world = World()
  f1,f2,f3=0.0,0.0,0.0
  d1,d2,d3=1.0,1.0,1.0
  s1 = Sampling(n1,d1,f1)
  s2 = Sampling(n2,d2,f2)
  s3 = Sampling(n3,d3,f3)
  ipg = addImageToWorld(s1,s2,s3,world,u,cmap=jet,cmin=min(u),cmax=max(u))
  f = readImage3D(n1,n2,n3,"f3dunc")
  ipg = addImageToWorld(s1,s2,s3,world,f,cmap=gray,cmin=min(f),cmax=max(f))
  makeFrame(world)

def addImageToWorld(s1,s2,s3,world,image,cmap=gray,gmin=None,gmax=None,u=None,cmin=0.0,cmax=0.0):
  ipg = ImagePanelGroup(s1,s2,s3,image)
  ipg.setColorModel(cmap)
  #ipg.setSlices(92,187,243)
  ipg.setSlices(114,122,174)
  #ipg.setSlices(91,187,206)
  sg12 = slice12(k3f,image)
  sg13 = slice13(k2f,image)
  sg23 = slice23(k1f,image)
  if gmin:
    ipg.setClips(gmin,gmax)
  else :
    gmin = min([sg12,sg13,sg23])
    gmax = max([sg12,sg13,sg23])
    ipg.setClips(gmin,gmax)
  world.addChild(ipg)
  if u:
    ipg = ImagePanelGroup(s1,s2,s3,u)
    ipg.setColorModel(ColorMap.getJet(0.8))
    #ipg.setSlices(92,187,244)
    ipg.setSlices(90,186,207)
    su12 = slice12(k3f,u)
    su13 = slice13(k2f,u)
    su23 = slice23(k1f,u)
    umin = min([su12,su13,su23])
    umax = max([su12,su13,su23])
    ipg.setClips(umin,umax)
    updateColorModel(ipg,1.0)
    world.addChild(ipg)

  return ipg

def addImage2ToWorld(world,image1,image2):
  ipg = ImagePanelGroup2(s1,s2,s3,image1,image2)
  ipg.setColorModel1(ColorMap.getGray())
  #ipg.setColorModel2(ColorMap.getJet(0.3))
  ipg.setColorModel2(ColorMap.getHue(0.0,20.0,0.3))
  world.addChild(ipg)
  return ipg

def makeFrame(world,png=None):
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  frame = SimpleFrame(world)
  view = frame.getOrbitView()
  zscale = 0.75*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.1)
  view.setAzimuth(azimuth)
  view.setElevation(elevation)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  #frame.viewCanvas.setBackground(frame.getBackground())
  frame.setSize(1000,900)
  frame.setVisible(True)
  if png:
    frame.paintToFile(pngDir+png+".png")
  return frame

#############################################################################
# read image
def readImage(name):
  fileName = seismicDir+name+".dat"
  image = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image
def readImage3D(n1,n2,n3,name):
  fileName = seismicDir+name+".dat"
  image = zerofloat(n1,n2,n3)
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

def slice12(k3,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n1,n2)
  SimpleFloat3(f).get12(n1,n2,0,0,k3,s)
  return s

def slice13(k2,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n1,n3)
  SimpleFloat3(f).get13(n1,n3,0,k2,0,s)
  return s

def slice23(k1,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n2,n3)
  st= zerofloat(n3,n2)
  SimpleFloat3(f).get23(n2,n3,k1,0,0,s)
  '''
  for i3 in range(n3):
    for i2 in range(n2):
      st[i2][i3] = s[i3][i2]
  return st
  '''
  return s

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
