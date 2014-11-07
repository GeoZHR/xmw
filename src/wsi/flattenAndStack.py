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
n1,n2,n3=1100,1024,100
f1,f2,f3=1.85928,3.048,0
d1,d2,d3=0.00762,0.02286,1
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

def main(args):
  #slopes()
  fL = readImage("junkL")
  fC = readImage("junkC")
  imgL = goStack(fL)
  imgC = goStack(fC)
  fmgL = goFlattenAndStack(fL)
  writeImage("fmgL",fmgL)
  fmgL = readImage2("fmgL")
  imgL = mul(gain(imgL),1)
  imgC = mul(gain(imgC),1)
  fmgL = mul(gain(fmgL),1)
  fMin = min(min(imgL),min(imgC),min(fmgL))
  fMax = max(max(imgL),max(imgC),max(fmgL))
  plot(s1,s2,imgL,clab="Amplitude",vlabel="z (km)",hlabel="x (km)", 
       cmin=fMin, cmax=fMax,png="image")
  plot(s1,s2,imgC,clab="Amplitude",vlabel="z (km)",hlabel="x (km)", 
       cmin=fMin, cmax=fMax,png="image")
  plot(s1,s2,fmgL,clab="Amplitude",vlabel="z (km)",hlabel="x (km)", 
       cmin=fMin, cmax=fMax,png="image")
def goFlattenAndStack(f):
  sg = zerofloat(n1,n2)
  fg = zerofloat(n1,n3)
  for i2 in range(n2):
    print i2
    for i3 in range(n3):
      fg[i3] = f[i3][i2]
    gg = goFlatten(fg)
    sg[i2] = goStack2(gg)
    '''
    plot(s1,s3,fg,clab="Amplitude",vlabel="z (km)",hlabel="x (km)", 
       cmin=min(fg), cmax=max(fg),png="image")
    plot(s1,s3,gg,clab="Flatten",vlabel="z (km)",hlabel="x (km)", 
       cmin=min(gg), cmax=max(gg),png="image")
    '''
  return sg

def goFlatten(f):
  sigma1,sigma2=4.0,1.0 
  pmax = 4.0
  lsf = LocalSlopeFinder(sigma1,sigma2,pmax)
  p2 = zerofloat(n1,n3)
  wp = zerofloat(n1,n3)
  lsf.findSlopes(f,p2,wp)
  p2 = mul(d1/d3,p2)
  wp = pow(wp,8)
  fl = Flattener2C()
  fl.setWeight1(0.02)
  fl.setIterations(0.01,200)
  fl.setSmoothings(4.0,8.0)
  c = fl.referenceTrace(f)
  fm = fl.getMappingsFromSlopes(s1,s3,p2,wp,c)
  g = fm.flatten(f)
  return g

def goStack2(f):
  img = zerofloat(n1)
  for i3 in range(n3):
    add(f[i3],img,img)
  return img


def goStack(f):
  img = zerofloat(n1,n2)
  for i3 in range(n3):
    for i2 in range(n2):
      add(f[i3][i2],img[i2],img[i2])
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
  
def readImage(name):
  fileName = seismicDir+name+".bin"
  n1,n2,n3 = s1.count,s2.count,s3.count
  #n1,n2,n3 = 120,480,430
  image = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
  ais.readFloats(image)
  ais.close()
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

def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(100.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2)
  div(x,sqrt(g),y)
  return y
def setWeights(kk,wp,ws):
  np = len(kk[0])
  for ip in range(np):
    i1 = int(kk[0][ip])
    i2 = int(kk[1][ip])
    i3 = int(kk[2][ip])
    wp[i3][i2][i1]=0.1*wp[i3][i2][i1]
    ws[i3][i2][i1]=0.1*ws[i3][i2][i1]
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

def plot(s1,s2,x,u=None,c=None,cmap=ColorMap.GRAY,clab=None,vlabel=None,hlabel=None,
  cmin=0,cmax=0,title=None,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  #if title:
  #  sp.setTitle(title)
  cb=sp.addColorBar(clab)
  pv = sp.addPixels(s1,s2,x)
  pv.setColorModel(cmap)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  sp.setVLabel(vlabel)
  sp.setHLabel(hlabel)
  #sp.setSize(1200,600)
  sp.setSize(600,1200)
  sp.setFontSizeForPrint(8.0,200)
  #sp.setFontSizeForSlide(0.5,0.9,16.0/9.0)
  sp.plotPanel.setColorBarWidthMinimum(120)
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
