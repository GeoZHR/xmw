"""
Demonstrate generating an RGT volume or flattening with control points. 
Test using the Aus NW Shelf dataset provided by dGB.
Author: Xinming Wu, Colorado School of Mines
Version: 2015.06.03
"""

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

from hv import *
from poseidon import *
from util import *

# Names and descriptions of image files used below.
gxfile  = "gx" # input image
gtfile  = "gt" # RGT volume
grfile  = "gr" # resampled RGT volume
ghfile  = "gh" # horizon volume
gufile  = "gu" # flattened image 
p2file  = "p2" # inline slopes
p3file  = "p3" # crossline slopes
epfile  = "ep" # eigenvalue-derived planarity

#s1 = Sampling(301,1.40,0.004)
#s2 = Sampling(351,1700,2.000)
#s3 = Sampling(300,1201,1.000)
s1 = Sampling(241,1,0)
s2 = Sampling(1399,1,0)
s3 = Sampling(923,1,0)


n1,n2,n3 = s1.count,s2.count,s3.count
d1,d2,d3 = s1.delta,s2.delta,s3.delta

pngDir = "../../../png/dgb/rgtc/"
seismicDir = "../../../data/seis/poseidon/"
#pngDir = "../../../png/dgb/subset/"
#seismicDir = "../../../data/seis/dgb/subset/"
pngDir = None
plotOnly = True

# Three sets of control points, each set 
# (k11 k12 k13 or k21 k22 k23 or k31 k32 k33) 
# belongs to one seismic horizon

k11 = [ 86, 93, 48, 46, 78, 72] #1st coordinates of the control points
k12 = [300,300,300, 28, 77, 63] #2nd coordinates of the control points
k13 = [ 32, 65,242,271, 85, 31] #3rd coordinates of the control points

k21 = [182,177,175,177,176,175,183,182,181,177]
k22 = [275,220,200,153,178,146, 94,300,300, 33]
k23 = [268,268,268,268,271,271,271,202, 94,289]

k31 = [244,244,257,256,254,254,245]
k32 = [300,300,300,300,244,132, 24]
k33 = [ 18, 73,161,241,271,271,271]

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goSlopes()
  #goFlattenC()
  gx = readImage(gxfile)
  plot3(gx)


def goSlopes():
  print "goSlopes ..."
  if not plotOnly:
    # set half-width of smoother for computing structure tensors
    sig1 = 2.0 #half-width in vertical direction
    sig2 = 2.0 #half-width in one literal direction
    sig3 = 2.0 #half-width in another literal direction
    pmax = 5.0 #maximum slope returned by this slope finder
    gx = readImage(gxfile)
    p2 = zerofloat(n1,n2,n3)
    p3 = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    lsf = LocalSlopeFinder(sig1,sig2,sig3,pmax)
    lsf.findSlopes(gx,p2,p3,ep);
    ep = pow(ep,10.0)
    #control points for extracting the water bottom surface
    c1=[ 31, 68,56]
    c2=[226,275,35]
    c3=[263, 53,35]
    zm = ZeroMask(c1,c2,c3,p2,p3,ep,gx)
    zero = 0.00;
    tiny = 0.01;
    zm.setValue(zero,p2)#set inline slopes for samples above water bottom
    zm.setValue(zero,p3)#set crossline slopes for samples above water bottom
    zm.setValue(tiny,ep)#set planarities for samples above water bottom
    writeImage(p2file,p2)
    writeImage(p3file,p3)
    writeImage(epfile,ep)
    print "p2  min =",min(p2)," max =",max(p2)
    print "p3  min =",min(p3)," max =",max(p3)
    print "ep  min =",min(ep)," max =",max(ep)
  else:
    gx = readImage(gxfile)
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    ep = readImage(epfile)
  plot3(gx)
  plot3(gx,p2, cmin=-1,cmax=1,cmap=jetRamp(1.0),
      clab="Inline slope (sample/sample)",png="p2")
  plot3(gx,p3, cmin=-1,cmax=1,cmap=jetRamp(1.0),
      clab="Crossline slope (sample/sample)",png="p3")
  plot3(gx,pow(ep,4.0),cmin=0,cmax=1,cmap=jetRamp(1.0),
      clab="Planarity")

def goFlattenC():
  print "Flatten with control points..."
  if not plotOnly:
    gx = readImage(gxfile)
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    ep = readImage(epfile)

    sc = SetupConstraints()
    kk1 = sc.extend(k11,k12,k13,n2,n3,p2,p3,ep)
    kk2 = sc.extend(k21,k22,k23,n2,n3,p2,p3,ep)
    kk3 = sc.extend(k31,k32,k33,n2,n3,p2,p3,ep)
    kk4 = sc.extend(k41,k42,k43,n2,n3,p2,p3,ep)
    kk5 = sc.extend(k51,k52,k53,n2,n3,p2,p3,ep)
    kk6 = sc.extend(k61,k62,k63,n2,n3,p2,p3,ep)
    kk7 = sc.extend(k71,k72,k73,n2,n3,p2,p3,ep)
    kk8 = sc.extend(k81,k82,k83,n2,n3,p2,p3,ep)
    kk9 = sc.extend(k91,k92,k93,n2,n3,p2,p3,ep)
    kk10 = sc.extend(k101,k102,k103,n2,n3,p2,p3,ep)
    kk11 = sc.extend(k111,k112,k113,n2,n3,p2,p3,ep)
    kk12 = sc.extend(k121,k122,k123,n2,n3,p2,p3,ep)

    k1 = [kk1[0],kk2[0],kk3[0],kk4[0],kk5[0],kk6[0],kk7[0],
          kk8[0],kk9[0],kk10[0],kk11[0],kk12[0]]
    k2 = [kk1[1],kk2[1],kk3[1],kk4[1],kk5[1],kk6[1],kk7[1],
          kk8[1],kk9[1],kk10[1],kk11[1],kk12[1]]
    k3 = [kk1[2],kk2[2],kk3[2],kk4[2],kk5[2],kk6[2],kk7[2],
          kk8[2],kk9[2],kk10[2],kk11[2],kk12[2]]
    k4 = [kk1[3],kk2[3],kk3[3],kk4[3],kk5[3],kk6[3],kk7[3],
          kk8[3],kk9[3],kk10[3],kk11[3],kk12[3]]

    #p2 = mul(d1/d2,p2)
    #p3 = mul(d1/d3,p3)
    fl = Flattener3C()
    fl.setIterations(0.01,200)
    fl.setSmoothings(6.0,6.0)
    fl.setWeight1(0.05)  
    fl.setScale(0.001)
    fm = fl.getMappingsFromSlopes(s1,s2,s3,p2,p3,ep,k4,k1,k2,k3)
    gu = fm.flatten(gx) # flattened image
    gt = fm.u1 # rgt volume
    gh = fm.x1 # horizon volume
    gr = fl.resampleRgt(s1,gt)
    writeImage(gufile,gu)
    writeImage(gtfile,gt)
    writeImage(ghfile,gh)
    writeImage(grfile,gr)
  else:
    gx = readImage(gxfile)
    gu = readImage(gufile)
    gt = readImage(gtfile)
    gh = readImage(ghfile)
    gr = readImage(grfile)
    '''
  plot3(gx,png="seismic")
  plot3(gu,png="flattened")
  plot3(gx,gt,cmin=min(gt)+20,cmax=max(gt)-10,cmap=jetRamp(1.0),
        clab="Relative geologic time",png="rgt")
  plot3(gx,p2, cmin=-1,cmax=1,cmap=jetRamp(1.0),
      clab="Inline slope (sample/sample)",png="p2")
  plot3(gx,p3, cmin=-1,cmax=1,cmap=jetRamp(1.0),
      clab="Crossline slope (sample/sample)",png="p3")
  #fl = Flattener3C()
  #gr = fl.resampleRgt(s1,gt)
    '''
  f1 = min(gr)+50
  d1 = 5*s1.getDelta()
  n1 = round((max(gr)-f1)/d1)-1
  st = Sampling(n1,d1,f1)
  hfr = HorizonExtraction(s1,s2,s3,None,gr)
  k2 = 320
  k3s = [11,61,91,151,181,251,291]
  for k3 in k3s:
    #hls  = hfr.horizonCurves(st,k2,k3)
    plot3X(gx,k2=k2,k3=k3,png="seismic"+str(k3)+"new")
    #plot3X(gx,k2=k2,k3=k3,curve=True,hs=hls,png="horizonLines"+str(k3)+"new")

#############################################################################
# read/write files
def readImage(name):
  fileName = seismicDir+name+".dat"
  n1,n2,n3 = s1.count,s2.count,s3.count
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

def plot3(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          hs=None,surf=None, png=None):
  n1,n2,n3 = s1.count,s2.count,s3.count
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
      ipg.setClips(-1.0,1.0)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-1.0,1.0)
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
  if hs:
    x1 = readImage(ghfile)
    u1 = readImage(gtfile)
    hfr = HorizonFromRgt(s1,s2,s3,x1,u1)
    for hi in hs:
      [xyz,rgb] = hfr.singleHorizon(hi)
      tg = TriangleGroup(True,xyz,rgb)
      sf.world.addChild(tg)
  if surf:
    tgs = Triangle()
    xyz = tgs.trianglesForSurface(surf,0,n1-1)
    tg  = TriangleGroup(True,xyz)
    sf.world.addChild(tg)
  #ipg.setSlices(290,300,268)
  #ipg.setSlices(290,300,271)
  ipg.setSlices(290,300,111)
  if cbar:
    sf.setSize(937,700)
  else:
    sf.setSize(800,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  ov = sf.getOrbitView()
  zscale = 0.6*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.5)
  #ov.setAzimuthAndElevation(220,25)
  ov.setAzimuthAndElevation(220,45)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(0.0,0.1,0.05))
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(137,1,pngDir+png+"cbar.png")

def plot3X(f,g=None,k1=300,k2=100,k3=100,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          smax=0.0,slices=None,curve=False,hs=None,png=None):
  n1,n2,n3 = s1.count,s2.count,s3.count
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
    cbar.setWidthMinimum(120)
  if hs:
    for hi in hs:
      if not curve:
        tg = TriangleGroup(True,hi[0],hi[1])
        sf.world.addChild(tg)
      else:
        lg = LineGroup(hi[0],hi[1])
        ss = StateSet()
        lg.setStates(ss)
        ls = LineState()
        ls.setWidth(1)
        ls.setSmooth(False)
        ss.add(ls)
        sf.world.addChild(lg)
  ipg.setSlices(290,300,271)
  ipg.setSlices(k1,k2,k3)
  if cbar:
    sf.setSize(937,700)
  else:
    sf.setSize(800,700)
    sf.setSize(1400,900)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  ov = sf.getOrbitView()
  zscale = 0.6*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(2.5)
  ov.setAzimuthAndElevation(230,10)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3-60,f2+40,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(0.0,0.1,0.05))
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")

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
