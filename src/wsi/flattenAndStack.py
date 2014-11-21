import sys

from java.awt import *
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

from wsi import *

#pngDir = None
pngDir = "../../../png/wsi/"
seismicDir = "../../../data/sigsbee/"
#n1,n2,n3=90,420,335
n1,n2,n3=1100,1024,20
#n1,n2,n3=1100,800,100
f1,f2,f3=1.85928,3.048,0
d1,d2,d3=0.00762,0.02286,1
#f1,f2,f3=0,0,0
#d1,d2,d3=1,1,1
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)
s3 = Sampling(n3,d3,f3)
#k1,k2,k3 = 69,419,0; azimuth=130; elevation=40 # for 3D view of horizon 1
#k1,k2,k3 = 82,419,0; azimuth=100; elevation=40 # for 3D view of horizon 2
#k1,k2,k3 = 88,60,160; azimuth=285; elevation=11 # for 3D view of all horizons
k1,k2,k3 = 69,415,8; azimuth=130; elevation=30 # for 3D view of strips
fmin,fmax = -5.5,5.5
k1f,k2f,k3f = 65,406,114
k1f,k2f,k3f = 48,406,114
k1f,k2f,k3f = 48,406,0
gmin,gmax,gint,glab = -2.0,2.0,0.5,"Amplitude"
background = Color.WHITE

fcfile = "sigImgCS"
flfile = "sigImgLS"
fhfile = "sigImgHS"
wcfile = "sigImgCW"
wlfile = "sigImgLW"
whfile = "sigImgHW"
fgcfile = "imgCS"
fglfile = "imgLS"
fghfile = "imgHS"
wgcfile = "imgCW"
wglfile = "imgLW"
wghfile = "imgHW"

plotOnly = True

def main(args):
  goLowVel()
  #goHighVel()
  #goCorrectVel()
def goSmooth():
  f = readImage(flfile,dat=True)
  f2 = zerofloat(n1,n3) 
  w2 = zerofloat(n1,n3) 
  for i3 in range(n3):
    f2[i3] = f[i3][920] 
  ws = WarpAndStack()
  w2 = ws.smooth(8.0,f2)
  plot(s1,s3,f2,clab="Amplitude",vlabel="z (km)",hlabel="x (km)", 
       cmin=min(f2), cmax=max(f2),wide=False,png="gather")
  plot(s1,s3,w2,clab="Warped",vlabel="z (km)",hlabel="x (km)", 
       cmin=min(w2), cmax=max(w2),wide=False,png="gatherS")

def goLowVel():
  '''
  goWarp(low=True)
  goStack(low=True)
  goImageShow(low=True)
  goGatherShow(low=True)
  '''
  goImageShowSub()

def goHighVel():
  goWarp(high=True)
  goStack(high=True)
  goImageShow(high=True)
  goGatherShow(high=True)

def goCorrectVel():
  goWarp(correct=True)
  goStack(correct=True)
  goImageShow(correct=True)
  goGatherShow(correct=True)

def goWarp(low=False,correct=False,high=False):
  if not plotOnly:
    if low:
      ffile,wfile = flfile,wlfile
      minlagSc,maxlagSc = 15,6
    if high:
      ffile,wfile = fhfile,whfile
      minlagSc,maxlagSc = 10,20
    if correct:
      ffile,wfile = fcfile,wcfile
      minlagSc,maxlagSc = 8,8
    esmooth,usmooth=2,1.0
    strainMax1,strainMax2=0.5,1.0
    #strainMax1,strainMax2=0.5,0.5
    f = readImage(ffile,dat=True)
    ws = WarpAndStack()
    ws.setForWarp(minlagSc,maxlagSc,esmooth,usmooth,strainMax1,strainMax2)
    ws.applyWarp(f)
    writeImage(wfile,f)

def goStack(low=False,correct=False,high=False):
  if not plotOnly:
    if low:
      ffile,wfile = flfile,wlfile
      fgfile,wgfile = fglfile,wglfile
    if high:
      ffile,wfile = fhfile,whfile
      fgfile,wgfile = fghfile,wghfile
    if correct:
      ffile,wfile = fcfile,wcfile
      fgfile,wgfile = fgcfile,wgcfile
    f  = readImage(ffile,dat=True)
    w  = readImage(wfile,dat=True)
    fg = zerofloat(n1,n2)
    wg = zerofloat(n1,n2)
    for i3 in range(n3):
      for i2 in range(n2):
        add(f[i3][i2],fg[i2],fg[i2])
        add(w[i3][i2],wg[i2],wg[i2])
    writeImage(fgfile,fg)
    writeImage(wgfile,wg)

def goImageShow(low=False,correct=False,high=False):
  if low:
    pngS = "L"
    fgfile,wgfile = fglfile,wglfile
  if high:
    pngS = "H"
    fgfile,wgfile = fghfile,wghfile
  if correct:
    pngS = "C"
    fgfile,wgfile = fgcfile,wgcfile
  fg = readImage2(fgfile)
  wg = readImage2(wgfile)
  fg = gain(fg,100)
  wg = gain(wg,100)
  fMin = max(min(fg),min(wg))
  fMax = min(max(fg),max(wg))
  plot(s1,s2,fg,clab="Amplitude",vlabel="z (km)",hlabel="x (km)", 
       cmin=fMin, cmax=fMax,wide=True,png="image"+pngS)
  plot(s1,s2,wg,clab="Amplitude",vlabel="z (km)",hlabel="x (km)", 
       cmin=fMin, cmax=fMax,wide=True,png="fmage"+pngS)
  lines=[10,110,210,310,410,510,610,710,810,910,1010]
  plot(s1,s2,fg,i2s=lines,lineColor=Color.blue,clab="Amplitude",vlabel="z (km)",
       hlabel="x (km)",cmin=fMin, cmax=fMax,wide=True,png="image"+pngS+"marked")
  lines=[136,200,267]
  plot(s1,s2,fg,i2s=lines,lineColor=Color.blue,clab="Amplitude",vlabel="z (km)",
       hlabel="x (km)",cmin=fMin, cmax=fMax,wide=True,png="image"+pngS+"defractors")

def goImageShowSub():
  fg = readImage2(fglfile)
  wg = readImage2(wglfile)
  fg = gain(fg,100)
  wg = gain(wg,100)
  nm1,nm2 = 400,420
  nf1,nf2 = 0,600
  fm1,fm2 = nf1*d1+f1,nf2*d2+f2
  sm1 = Sampling(nm1,d1,fm1)
  sm2 = Sampling(nm2,d2,fm2)
  fg = copy(nm1,nm2,nf1,nf2,fg)
  wg = copy(nm1,nm2,nf1,nf2,wg)
  fMin = max(min(fg),min(wg))
  fMax = min(max(fg),max(wg))
  plot(sm1,sm2,fg,clab="Amplitude",vlabel="z (km)",hlabel="x (km)", 
       cmin=fMin, cmax=fMax,wide=True,png="imageLSub")
  plot(sm1,sm2,wg,clab="Amplitude",vlabel="z (km)",hlabel="x (km)", 
       cmin=fMin, cmax=fMax,wide=True,png="fmageLSub")

def goGatherShow(low=False,correct=False,high=False):
  if low:
    pngS = "L"
    ffile,wfile = flfile,wlfile
  if high:
    pngS = "H"
    ffile,wfile = fhfile,whfile
  if correct:
    pngS = "C"
    ffile,wfile = fcfile,wcfile
  f = readImage(ffile,dat=True)
  w = readImage(wfile,dat=True)
  f2 = zerofloat(n1,n3) 
  w2 = zerofloat(n1,n3) 
  for i2 in range(10,n2,100):
    for i3 in range(n3):
      f2[i3] = f[i3][i2] 
      w2[i3] = w[i3][i2] 
    plot(s1,s3,f2,clab="Amplitude",vlabel="z (km)",hlabel="x (km)", 
         cmin=min(f2), cmax=max(f2),wide=False,png="gather"+pngS+str(i2))
    plot(s1,s3,w2,clab="Warped",vlabel="z (km)",hlabel="x (km)", 
         cmin=min(w2), cmax=max(w2),wide=False,png="gatherW"+pngS+str(i2))
  defractors=[136,200,267]
  for i2 in defractors:
    for i3 in range(n3):
      f2[i3] = f[i3][i2] 
      w2[i3] = w[i3][i2] 
    plot(s1,s3,f2,clab="Amplitude",vlabel="z (km)",hlabel="x (km)", 
       cmin=min(f2), cmax=max(f2),wide=False,png="gather"+pngS+str(i2))
    plot(s1,s3,w2,clab="Warped",vlabel="z (km)",hlabel="x (km)", 
       cmin=min(w2), cmax=max(w2),wide=False,png="gatherW"+pngS+str(i2))

def goFlattenTest():
  x = readImage(flfile,dat=True)
  x2 = zerofloat(n1,n3)
  for i3 in range(n3):
    x2[i3] = x[i3][350]
  p2 = zerofloat(n1,n3)
  el = zerofloat(n1,n3)
  lsf = LocalSlopeFinder(20.0,1.0,40.0)
  lsf.findSlopes(x2,p2,el)
  fl2 = Flattener2C()
  fl2.setWeight1(0.5)
  c = zeroint(1)
  fmp = fl2.getMappingsFromSlopes(s1,s3,p2,pow(el,20),c)
  gx = fmp.flatten(x2)
  plot(s1,s3,x2,clab="Amplitude",vlabel="z (km)",hlabel="x (km)", 
     cmin=min(x2), cmax=max(x2),wide=False,png=None)#"gather"+str(i2))
  plot(s1,s3,gx,clab="Flatten",vlabel="z (km)",hlabel="x (km)", 
     cmin=min(gx), cmax=max(gx),wide=False,png=None)

def goWarpTest(x):
  for i2 in range(n2-1):
    dw = DynamicWarping(-1,1)
    dw.setStrainMax(0.5)
    dw.setErrorSmoothing(2)
    dw.setShiftSmoothing(1.0)
    u = dw.findShifts(x[i2],x[i2+1])
    y = dw.applyShifts(u,x[i2+1])
    copy(y,x[i2+1])


def goFirstLook(f,g):
  f3 = zerofloat(n1,n3) 
  g3 = zerofloat(n1,n3) 
  for i2 in range(10,n2,100):
    for i3 in range(n3):
      f3[i3] = f[i3][i2] 
      g3[i3] = g[i3][i2] 
    plot(s1,s3,f3,clab="Amplitude",vlabel="z (km)",hlabel="x (km)", 
       cmin=min(f3), cmax=max(f3),wide=False,png=None)#"gather"+str(i2))
    plot(s1,s3,g3,clab="Amplitude",vlabel="z (km)",hlabel="x (km)", 
       cmin=min(g3), cmax=max(g3),wide=False,png=None)

def goFlattenAndStack(f):
  s1 = Sampling(n1,1.0,0.0)
  s2 = Sampling(n2,1.0,0.0)
  s3 = Sampling(n3,1.0,0.0)
  fas = FlattenAndStack()
  g = fas.apply(s1,s2,s3,f)
  return g

def goFlattenAndStack1(f):
  s1 = Sampling(n1,1.0,0.0)
  s2 = Sampling(n2,1.0,0.0)
  s3 = Sampling(n3,1.0,0.0)
  p2 = zerofloat(n1,n2,n3)
  p3 = zerofloat(n1,n2,n3)
  ep = zerofloat(n1,n2,n3)
  u1 = zerofloat(n1,n2,n3)
  lsf = LocalSlopeFinder(4.0,1.0,2.0,8.0)
  lsf.findSlopes(f,p2,p3,ep,u1)
  p2 = fillfloat(-0.3,n1,n2,n3)
  ep = pow(ep,4.0)
  #fl3 = Flattener3C()
  fl3 = Flattener3()
  fl3.setSmoothings(6.0,6.0,6.0)
  fl3.setWeight1(0.6)
  fl3.setIterations(0.01,200)
  #mp = fl3.getMappingsFromSlopes(s1,s2,s3,p2,p3,ep,u1)
  mp = fl3.getMappingsFromSlopes(s1,s2,s3,p2,p3,ep)
  g = mp.flatten(f)
  return goStack(g)

def goStack2(f):
  img = zerofloat(n1)
  for i3 in range(n3):
    add(f[i3],img,img)
  return img



def slopes():
  f = readImage("gm")
  p2 = copy(f)
  p3 = copy(f)
  ep = copy(f)
  lsf = LocalSlopeFinder(2.0,2.0)
  lsf.findSlopes(f,p2,p3,ep);
  writeImage("gmp2",p2)
  writeImage("gmp3",p3)
  #writeImage("gmep",ep)
  for g in [p2,p3,ep]:
    world = World()
    addImage2ToWorld(world,f,g)
    makeFrame(s1,s2,s3,azimuth,elevation,world)

#############################################################################
# read/write files
def readImage2(name):
  fileName = seismicDir+name+".dat"
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image
  
def readImage(name,dat=False):
  n1,n2,n3 = s1.count,s2.count,s3.count
  n1,n2,n3=1100,1024,20
  image = zerofloat(n1,n2,n3)
  if dat:
    fileName = seismicDir+name+".dat"
    ais = ArrayInputStream(fileName,ByteOrder.BIG_ENDIAN)
  else:
    fileName = seismicDir+name+".bin"
    byteType = "ByteOrder.LITTLE_ENDIAN"
    ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  ais.readFloats(image)
  ais.close()
  return image

def readImageT(name1,name2):
  n1,n2,n3=1100,1024,100
  imageT = zerofloat(n2,n1,n3)
  fileName = seismicDir+name1+".bin"
  byteType = "ByteOrder.LITTLE_ENDIAN"
  ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  ais.readFloats(imageT)
  ais.close()
  image = zerofloat(n1,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(n1):
        image[i3][i2][i1] = imageT[i3][i1][i2]
  writeImage(name2,image)
  return image

def readImageS(name1,name2):
  n1,n2,n3=1100,1024,100
  image = zerofloat(n1,n2,n3)
  fileName = seismicDir+name1+".dat"
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  imageS = zerofloat(n1,n2,20)
  k3 = 0
  for i3 in range(0,n3,5):
    for i2 in range(n2):
      for i1 in range(n1):
        imageS[k3][i2][i1] = image[i3][i2][i1]
    k3 += 1
  writeImage(name2,imageS)
  return image

def writeImage(name,image):
  fileName = seismicDir+name+".dat"
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()
  return image

def readSlice3(name):
  fileName = seismicDir+name+".dat"
  n1,n2 = s1.count,s2.count
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
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
  SimpleFloat3(f).get23(n2,n3,k1,0,0,s)
  return s

def gain(x,e):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(e)
  ref.apply1(g,g)
  n2 = len(x)
  n1 = len(x[0])
  y = zerofloat(n1,n2)
  for i2 in range(n2):
    div(x[i2],sqrt(g[i2]),y[i2])
  return y
def setBounds(kk,wp,w1):
  wpt = zerofloat(n1,n2,n3)
  wpt = copy(wp)
  for i3 in range(n3):
    for i2 in range(n2):
      wt1= wpt[i3][i2][n1-1]
      wt2= wpt[i3][i2][n1-2]
      wpt[i3][i2][n1-1] = 0.5*wt1
      wpt[i3][i2][n1-2] = 0.5*wt2
  np = len(kk[0])
  for ip in range(np):
    i2 = int(kk[1][ip])
    i3 = int(kk[2][ip])
    wpt[i3][i2][n1-1]= wp[i3][i2][n1-1]
    wpt[i3][i2][n1-2]= wp[i3][i2][n1-2]
  copy(wpt,wp)

#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET

def plot(s1,s2,x,i2s=None,i3=None,lineColor=Color.red,
        cmap=ColorMap.GRAY,clab=None,vlabel=None,hlabel=None,
        cmin=0,cmax=0,wide=False,title=None,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  #if title:
  #  sp.setTitle(title)
  cb=sp.addColorBar(clab)
  pv = sp.addPixels(s1,s2,x)
  pv.setColorModel(cmap)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if i2s:
    for i2 in i2s:
      x = zerofloat(n1-20)
      y = zerofloat(n1-20)
      for i1 in range(10,n1-10,1):
        x[i1-10] = i2*d2+f2
        y[i1-10] = i1*d1+f1
      ps = sp.addPoints(y,x)
      ps.setLineWidth(2.0)
      ps.setLineColor(lineColor)
  if i3:
    print i3
    x = zerofloat(n1-4)
    y = zerofloat(n1-4)
    for i1 in range(2,n1-2,1):
      x[i1-2] = i3*d3+f3
      y[i1-2] = i1*d1+f1
    ps = sp.addPoints(y,x)
    ps.setLineWidth(2.0)
    ps.setMarkColor(Color.red)

  if cmin<cmax:
    pv.setClips(cmin,cmax)
  sp.setVLabel(vlabel)
  sp.setHLabel(hlabel)
  if wide:
    sp.setSize(1400,750)
    sp.setFontSizeForPrint(6.0,200)
    sp.plotPanel.setColorBarWidthMinimum(140)
  else:
    sp.setSize(500,1200)
    sp.setFontSizeForPrint(8.0,200)
    sp.plotPanel.setColorBarWidthMinimum(140)
  #sp.setFontSizeForSlide(0.5,0.9,16.0/9.0)
  #sp.plotPanel.setColorBarWidthMinimum(120)
  if pngDir and png:
    sp.paintToPng(720,2.2222,pngDir+png+".png")


def addImageToWorld(s1,s2,s3,world,image,cmap=gray,cmin=0,cmax=0):
  ipg = ImagePanelGroup(s1,s2,s3,image)
  ipg.setColorModel(cmap)
  if cmin<cmax:
    ipg.setClips(cmin,cmax)
  world.addChild(ipg)
  return ipg

def addImage2ToWorld(world,image1,image2):
  ipg = ImagePanelGroup2(s1,s2,s3,image1,image2)
  ipg.setColorModel1(ColorMap.getGray())
  ipg.setColorModel2(ColorMap.getJet())
  #ipg.setColorModel2(ColorMap.getHue(0.0,20.0,0.3))
  world.addChild(ipg)
  return ipg

def makeFrame(s1,s2,s3,azimuth,elevation,world,png=None):
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  lightPosition=[-0.18,-0.4,0.8,0.0] # good for horizons 1 and 2
  lightPosition=[0.,0.,1.0,0.0] #default position
  frame = SimpleFrame(world)
  view = frame.getOrbitView()
  #view.setLightPosition(lightPosition)
  zscale = 0.6*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.1)
  view.setAzimuth(azimuth)
  view.setElevation(elevation)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  #frame.viewCanvas.setBackground(frame.getBackground())
  frame.setSize(1000,900)
  frame.setVisible(True)
  if png and pngDir:
    png = pngDir+png
    frame.paintToFile(png+".png")
  return frame
"""
  frame = SimpleFrame(world)
  view = frame.getOrbitView()
  zscale = 0.75*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.3)
  #view.setAzimuth(75.0)
  #view.setAzimuth(-75.0)
  view.setAzimuth(-65.0)
"""
def plot3f(k1,k2,k3,g,cbar,label1,gmap=gray,cint=0.5,gmin=None,gmax=None,png=None):
  pp = PlotPanelPixels3(
    PlotPanelPixels3.Orientation.X1DOWN_X3RIGHT,
    PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM,
    s1,s2,s3,g)
  pp.setSlices(k1,k2,k3)
  pp.setInterpolation(PixelsView.Interpolation.LINEAR)
  pp.setLabel1(label1)
  pp.setLabel2("Inline (km)")
  pp.setLabel3("Crossline (km)")
  pp.mosaic.setHeightElastic(0,180)
  pp.mosaic.setHeightElastic(1, 70)
  if gmin !=gmax:
    pp.setClips(gmin,gmax)
  pp.setColorModel(gmap)
  pp.setLineColor(Color.YELLOW)
  cb = pp.addColorBar(cbar)
  cb.setInterval(cint)
  pp.setInterval1(0.1)
  pp.setInterval2(2.0)
  pp.setInterval3(2.0)
  pf = PlotFrame(pp)
  pf.setBackground(background)
  pp.setColorBarWidthMinimum(70)
  #pf.setFontSizeForSlide(1.0,0.8)
  pf.setSize(1100,700)
  pf.setVisible(True)
  if png and pngDir:
    png = pngDir+png
    pf.paintToPng(360,7.0,png+".png")
 
def plotFrame(s1,s2,f,h,i3t):
  orient = PlotPanel.Orientation.X1DOWN_X2RIGHT
  panel  = PlotPanel(1,1,orient)
  pxv    = panel.addPixels(0,0,s1,s2,f)
  pxv.setColorModel(ColorMap.GRAY)
  ptv = panel.addPoints(0,0,h[0],h[1])
  ptv.setStyle("r-")
  panel.setTitle("section "+i3t)
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setVisible(True)
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
