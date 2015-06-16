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
from dgb import *
from util import *

# Names and descriptions of image files used below.
gxfile  = "gx" # input image
gtfile  = "gt" # RGT volume
ghfile  = "gh" # horizon volume
gufile  = "gu" # flattened image 
p2file  = "p2" # inline slopes
p3file  = "p3" # crossline slopes
epfile  = "ep" # eigenvalue-derived planarity

#s1 = Sampling(301,1.40,0.004)
#s2 = Sampling(351,1700,2.000)
#s3 = Sampling(300,1201,1.000)
s1 = Sampling(301,1,0)
s2 = Sampling(351,1,0)
s3 = Sampling(300,1,0)

n1,n2,n3 = s1.count,s2.count,s3.count
d1,d2,d3 = s1.delta,s2.delta,s3.delta

#pngDir = None
pngDir = "../../../png/dgb/rgtc/"
seismicDir = "../../../data/seis/dgb/rgtc/"
plotOnly = False

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

k11 = [ 56, 49, 50, 52, 47, 55, 49, 56, 70, 73, 70, 72, 72, 54, 61, 49, 29, 30, 19, 21, 53, 54, 63, 43, 28, 23, 16, 25, 32, 50, 43, 66, 56, 19, 18, 19, 35, 58]
k12 = [ 40, 82, 32, 32, 145, 173, 175, 193, 246, 300, 307, 333, 32, 32, 23, 32, 32, 32, 20, 32, 111, 175, 175, 175, 175, 175, 163, 87, 57, 169, 219, 307, 307, 292, 307, 307, 307, 307]
k13 = [ 50, 50, 27, 42, 50, 50, 34, 50, 50, 50, 33, 50, 119, 60, 150, 167, 211, 248, 280, 275, 150, 133, 85, 169, 211, 266, 280, 280, 280, 150, 150, 75, 108, 280, 266, 239, 209, 181]

k21 = [ 41, 35, 43, 37, 46, 42, 40, 41, 57, 60, 52, 51, 57, 64, 68, 69, 72, 75, 73, 75, 39, 50, 34, 70, 69, 72, 69, 56, 71, 70, 66, 59, 58, 56, 63, 65, 73, 56, 67, 48, 39, 44, 38, 33, 36, 71, 75]
k22 = [ 141, 170, 216, 251, 282, 307, 307, 307, 307, 307, 307, 307, 307, 258, 206, 307, 307, 307, 307, 307, 175, 175, 175, 175, 130, 53, 12, 296, 175, 175, 175, 157, 108, 68, 188, 211, 265, 32, 32, 32, 32, 110, 43, 11, 32, 32, 32]
k23 = [ 280, 280, 280, 280, 280, 269, 224, 211, 192, 181, 170, 163, 138, 150, 150, 105, 72, 55, 21, 4, 264, 206, 296, 165, 150, 150, 150, 150, 133, 92, 71, 50, 50, 50, 50, 50, 50, 192, 157, 225, 240, 280, 280, 280, 267, 135, 105]

k31 = [ 60, 60, 53, 51, 51, 50, 50, 53, 50, 54, 46, 53, 59, 48, 48, 54, 67, 76, 77, 83, 79, 84, 86, 78, 73, 60, 83, 82, 72, 63, 60, 61, 58, 83]
k32 = [ 307, 298, 264, 212, 178, 175, 165, 126, 100, 49, 6, 175, 175, 32, 32, 32, 32, 32, 22, 32, 32, 60, 135, 174, 175, 175, 175, 175, 222, 283, 307, 307, 307, 175]
k33 = [ 288, 280, 280, 280, 280, 295, 280, 280, 280, 280, 280, 249, 203, 267, 235, 220, 169, 152, 150, 140, 109, 150, 150, 150, 170, 199, 108, 92, 150, 150, 139, 244, 204, 127]

k41 = [ 56, 66, 66, 70, 72, 75, 76, 72, 72, 77, 78, 79, 84, 86, 91, 94, 93, 90, 95, 81, 69, 84, 87, 90, 100, 93, 89, 87, 83, 83, 82, 57, 60, 66, 81]
k42 = [ 6, 81, 96, 175, 151, 175, 175, 267, 307, 307, 307, 307, 307, 307, 307, 307, 307, 307, 307, 300, 201, 277, 243, 201, 176, 123, 86, 41, 19, 6, 32, 32, 32, 32, 32]
k43 = [ 280, 280, 280, 279, 280, 235, 208, 280, 213, 184, 154, 140, 112, 96, 68, 51, 32, 2, 18, 150, 280, 150, 150, 150, 150, 150, 150, 150, 150, 150, 163, 262, 217, 199, 185]
# FQ E5
k51 = [ 72, 72, 73, 71, 76, 84, 82, 83, 89, 92, 99, 106, 89, 98, 89, 104, 104, 99, 86, 77, 84, 98, 100, 91, 95, 88, 91, 87, 114, 119, 122, 108, 93, 86, 95, 85, 71, 69, 69, 83, 85, 93, 88, 104, 106, 99, 90, 78, 76, 72, 73, 75, 74, 66, 72, 71, 73, 78, 77, 82, 95, 104, 92, 89, 87, 90, 97, 89, 93, 82, 89, 83, 90, 85, 90, 92, 88, 97, 95, 94, 97, 95, 99, 100, 96, 105, 110, 118, 109, 105, 102, 107, 101, 102, 100, 93, 90]
k52 = [ 4, 27, 61, 11, 94, 329, 307, 298, 279, 261, 237, 210, 174, 144, 117, 175, 175, 175, 175, 175, 175, 175, 156, 82, 130, 17, 60, 5, 175, 175, 175, 175, 175, 175, 175, 157, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 21, 32, 32, 32, 32, 32, 37, 66, 109, 130, 151, 173, 179, 200, 233, 262, 273, 286, 307, 307, 307, 307, 307, 307, 307, 307, 307, 307, 307, 307, 307, 307, 307, 307, 306, 184, 200, 208, 218, 228, 234, 243, 270, 291, 300, 175, 175, 175]
k53 = [ 280, 280, 280, 280, 280, 280, 287, 280, 280, 280, 280, 280, 280, 280, 280, 264, 253, 240, 214, 194, 178, 160, 150, 150, 150, 150, 150, 150, 129, 115, 104, 89, 71, 59, 51, 50, 247, 223, 200, 185, 169, 160, 146, 129, 117, 98, 85, 67, 52, 50, 42, 33, 26, 20, 1, 50, 50, 50, 50, 50, 50, 150, 150, 150, 150, 150, 150, 151, 163, 276, 268, 229, 258, 208, 181, 169, 148, 114, 96, 90, 85, 74, 54, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 39, 9, 0]



# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  goSlopes()
  goFlattenC()

def goSlopes():
  print "goSlopes ..."
  if not plotOnly:
    # set half-width of smoother for computing structure tensors
    sig1 = 4.0 #half-width in vertical direction
    sig2 = 1.0 #half-width in one literal direction
    sig3 = 1.0 #half-width in another literal direction
    pmax = 5.0 #maximum slope returned by this slope finder
    gx = readImage(gxfile)
    p2 = zerofloat(n1,n2,n3)
    p3 = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    lsf = LocalSlopeFinder(sig1,sig2,sig3,pmax)
    lsf.findSlopes(gx,p2,p3,ep);
    ep = pow(ep,6.0)
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
    kk1 = sc.extend(k11,k12,k13,n2,n3,p2,p3,ep,gx)
    kk2 = sc.extend(k21,k22,k23,n2,n3,p2,p3,ep,gx)
    kk3 = sc.extend(k31,k32,k33,n2,n3,p2,p3,ep,gx)
    kk4 = sc.extend(k41,k42,k43,n2,n3,p2,p3,ep,gx)
    kk5 = sc.extend(k51,k52,k53,n2,n3,p2,p3,ep,gx)
    k1 = [kk1[0],kk2[0],kk3[0],kk4[0],kk5[0]]
    k2 = [kk1[1],kk2[1],kk3[1],kk4[1],kk5[1]]
    k3 = [kk1[2],kk2[2],kk3[2],kk4[2],kk5[2]]
    k4 = [kk1[3],kk2[3],kk3[3],kk4[3],kk5[3]]

    p2 = mul(d1/d2,p2)
    p3 = mul(d1/d3,p3)
    fl = Flattener3C()
    fl.setIterations(0.01,500)
    fl.setSmoothings(6.0,6.0)
    #fl.setWeight1(0.05)  
    fl.setScale(0.0001)
    fm = fl.getMappingsFromSlopes(s1,s2,s3,p2,p3,ep,k4,k1,k2,k3)
    gu = fm.flatten(gx) # flattened image
    gt = fm.u1 # rgt volume
    gh = fm.x1 # horizon volume
    writeImage(gufile,gu)
    writeImage(gtfile,gt)
    writeImage(ghfile,gh)
  else:
    gx = readImage(gxfile)
    gu = readImage(gufile)
    gt = readImage(gtfile)
    gh = readImage(ghfile)
  plot3(gx,png="seismic")
  plot3(gu,png="flattened")
  plot3(gx,gt,cmin=min(gt)+20,cmax=max(gt)-10,cmap=jetRamp(1.0),
        clab="Relative geologic time",png="rgt")
  ha = []
  hs = [280,260,240,220,200,180,160,140,120,100,80,60,40]
  for ih, h in enumerate(hs):
    ha.append(h)
    plot3(gx,hs=ha,png="horizon"+str(ih))

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
  ipg.setSlices(290,300,271)
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
  ov.setAzimuthAndElevation(220,25)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(0.0,0.1,0.05))
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
