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

from poseidon import *

# Names and descriptions of image files used below.
gxfile  = "gx" # input image
hzfile  = "hz" # a single horizon
hvfile  = "hv" # horizon volume
p2file  = "p2" # inline slopes
p3file  = "p3" # crossline slopes
epfile  = "ep" # eigenvalue-derived planarity

s1 = Sampling(101,1,0) #copy(120,n2,n3,40,0,0,gx)
s2 = Sampling(102,1,0)
s3 = Sampling(103,1,0)


n1,n2,n3 = s1.count,s2.count,s3.count
d1,d2,d3 = s1.delta,s2.delta,s3.delta

seismicDir = "../../../data/seis/poseidon/fake/"
pngDir = "../../../png/poseidon/fake/"
pngDir = None
plotOnly = False


# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goFakeData()
  #goSlopes()
  #goSingleHorizon()
  goHorizonVolume()

def goFakeData():
  #sequence = 'A' # 1 episode of faulting only
  sequence = 'OA' # 1 episode of folding, followed by one episode of faulting
  #sequence = 'OOOOOAAAAA' # 5 episodes of folding, then 5 of faulting
  #sequence = 'OAOAOAOAOA' # 5 interleaved episodes of folding and faulting
  nplanar = 0 # number of planar faults
  conjugate = False # if True, two large planar faults will intersect
  conical = False # if True, may want to set nplanar to 0 (or not!)
  impedance = False # if True, data = impedance model
  wavelet = True # if False, no wavelet will be used
  noise = 0.8 # (rms noise)/(rms signal) ratio
  gx,p2,p3 = FakeData.seismicAndSlopes3d2014A(
      sequence,nplanar,conjugate,conical,impedance,wavelet,noise)
  writeImage(gxfile,gx)
  print "gx min =",min(gx)," max =",max(gx)
  gmin,gmax,gmap = -3.0,3.0,ColorMap.GRAY
  if impedance:
    gmin,gmax,gmap = 0.0,1.4,ColorMap.JET
  plot3(gx)
  #plot3(gx,cmin=gmin,cmax=gmax,cmap=gmap,clab="Amplitude",png="gx")

def goSlopes():
  print "goSlopes ..."
  if not plotOnly:
    # set half-width of smoother for computing structure tensors
    sig1 = 3.0 #half-width in vertical direction
    sig2 = 2.0 #half-width in one literal direction
    sig3 = 2.0 #half-width in another literal direction
    pmax = 4.0 #maximum slope returned by this slope finder
    gx = readImage(gxfile)
    p2 = zerofloat(n1,n2,n3)
    p3 = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    lsf = LocalSlopeFinder(sig1,sig2,sig3,pmax)
    lsf.findSlopes(gx,p2,p3,ep);
    '''
    zm = ZeroMask(0.1,4,1,1,gx)
    zero,tiny=0.0,0.01
    zm.setValue(zero,p2)
    zm.setValue(zero,p3)
    zm.setValue(tiny,ep)
    '''
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

def goSingleHorizon():
  gx = readImage(gxfile)
  # set a control point
  k1 = [50]
  k2 = [50] 
  k3 = [50]
  if not plotOnly:
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    ep = readImage(epfile)
    ep = pow(ep,8)
    lmt = n1-1.0
    se = SurfaceExtractorM()
    se.setWeights(0.5)
    se.setExternalIterations(12)
    sf = se.surfaceInitialization(n2,n3,lmt,k1,k2,k3)
    se.surfaceUpdateFromSlopes(ep,p2,p3,None,k1,k2,k3,sf)
    writeImage(hzfile,sf)
  else:
    sf = readImage2(n2,n3,hzfile)
  plot3(gx)
  plot3(gx,sfs=[sf])

def goHorizonVolume():
  gx = readImage(gxfile)
  # set a reference trace
  k2,k3 = 50,50
  k1 = zerofloat(n1)
  for i1 in range(n1):
    k1[i1] = i1
  if not plotOnly:
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    ep = readImage(epfile)
    ep = pow(ep,8)
    hv = HorizonVolume()
    hv.setWeights(0.0)
    hv.setExternalIterations(10)
    hs = hv.applyForHorizonVolume(k1,k2,k3,ep,p2,p3)
    writeImage(hvfile,hs)
  else:
    hs = readImageX(n2,n3,len(k1),hvfile)
  plot3(gx,sfs=hs)
  # show only a subset of horizons 
  hss = []
  for k in range(0,len(hs),len(hs)/10):
    hss.append(hs[k])
  plot3(gx,sfs=hss)

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

def readImageX(n1,n2,n3,name):
  fileName = seismicDir+name+".dat"
  image = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def readImage2(n1,n2,name):
  fileName = seismicDir+name+".dat"
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
          sfs=None, png=None):
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
  if sfs:
    for ik in range(len(sfs)):
      tg  = TriangleGroup(True,s3,s2,sfs[ik])
      tg.setColor(Color.CYAN)
      sf.world.addChild(tg)
    '''
    sd = SurfaceDisplay()
    ts = sd.horizonWithAmplitude([-1.0,1.0],surf,f)
    tg = TriangleGroup(True,ts[0],ts[1])
    sf.world.addChild(tg)
    '''
  #ipg.setSlices(290,300,268)
  #ipg.setSlices(290,300,271)
  ipg.setSlices(290,1009,768)
  if cbar:
    sf.setSize(937,700)
  else:
    sf.setSize(1100,900)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  ov = sf.getOrbitView()
  zscale = 0.8*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.4)
  #ov.setAzimuthAndElevation(220,25)
  ov.setAzimuthAndElevation(220,25)
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
      ipg.setClips(-1.0,1.0) # use for subset plots
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
    sf.setSize(1200,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  ov = sf.getOrbitView()
  zscale = 0.4*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(2.0)
  '''
  ov.setAzimuthAndElevation(310,15)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(0.2,0.08,-0.20))
  '''
  ov.setAzimuthAndElevation(390,15)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(0.0,0.08,-0.20))
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
