"""
Demonstrate structure-oriented processing
Author: Xinming Wu, University of Texas at Austin
Version: 2016.09.16
"""
import sys

from java.awt import *
from java.io import *
from java.lang import *
from javax.swing import *
from java.util import *


from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.interp import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from util import *

pngDir = "../../../png/mda/fed2/"
pngDir = None

seismicDir = "../../../data/seis/mda/"
gxfile = "fingerPrint"
g1file = "g1"
g2file = "g2"
f1,f2 = 0,0
d1,d2 = 1,1
n1,n2 = 1312,1316
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goGradients()
  #goNormals()
  goTensors()
  #goGaussianSmoothing()
  #goAnisotropicDiffusion()

def goGradients():
  print "goGradients ..."
  gx = readImage(n1,n2,gxfile)
  g1 = readImage(n1,n2,g1file)
  g2 = readImage(n1,n2,g2file)
  plot2(s1,s2,gx,cmin=0.3,cmax=1.0)
  plot2(s1,s2,gx,u1=g1,u2=g2,dx=40,scale=500,cmin=0.3,cmax=1)

def goNormals():
  print "goNormals ..."
  gx = readImage(n1,n2,gxfile)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lof = LocalOrientFilterP(20,20)
  lof.applyForNormalLinear(gx,u1,u2,el)
  plot2(s1,s2,gx,u1=u1,u2=u2,dx=40,scale=20,cmin=0.3,cmax=1)

def goTensors():
  print "goTensors ..."
  gx = readImage(n1,n2,gxfile)
  el = zerofloat(n1,n2)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  lof = LocalOrientFilterP(20,20)
  et = lof.applyForTensors(gx)
  lof.applyForNormalLinear(gx,u1,u2,el)
  eu = zerofloat(n1,n2)
  ev = zerofloat(n1,n2)
  et.getEigenvalues(eu,ev)
  eu=clip(0.00005,max(eu),eu)
  ev=clip(0.00005,max(ev),ev)
  et.setEigenvalues(eu,ev)
  et.invertStructure(1.0,1.0)
  et.getEigenvalues(eu,ev)
  #plotTensors(gx,s1,s2,d=et,dscale=20,mk=mk,cmin=-2,cmax=2,png="tensors")
  #plot2(s1,s2,gx)
  plotTensors(gx,s1,s2,d=et,dscale=1,ne=20,cmin=0.3,cmax=1,png="tensors")
  #plot2(s1,s2,el,cmin=0.1,cmax=1.0)

def goGaussianSmoothing():
  print "goGaussianSmoothing ..."
  gx = readImage(n1,n2,gxfile)
  gs = zerofloat(n1,n2)
  rgf = RecursiveGaussianFilterP(30)
  rgf.apply00(gx,gs)
  plot2(s1,s2,gx,cmin=0.3,cmax=1)
  plot2(s1,s2,gs,cmin=0.3,cmax=1)

def goIsotropicDiffusion():
  print "goIsotropicDiffusion ..."
  gx = readImage(n1,n2,gxfile)
  gs = zerofloat(n1,n2)
  lof = LocalOrientFilterP(4,2)
  et = lof.applyForTensors(gx)
  et.setEigenvalues(1.0,1.0)
  lsf = LocalSmoothingFilter()
  lsf.apply(et,400,gx,gs)
  plot2(s1,s2,gx)
  plot2(s1,s2,gs)

def goAnisotropicDiffusion():
  print "goAnisotropicDiffusion ..."
  gx = readImage(n1,n2,gxfile)
  g1 = zerofloat(n1,n2)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lof = LocalOrientFilterP(20,20)
  lof.applyForNormalLinear(gx,u1,u2,el)
  et = lof.applyForTensors(gx)
  et.setEigenvalues(0.001,1.0)
  lsf = LocalSmoothingFilter(0.01,200)
  lsf.apply(et,1000,gx,g1)
  plot2(s1,s2,gx,cmin=0.3,cmax=1)
  plot2(s1,s2,g1,cmin=0.3,cmax=1)

def like(x):
  n2 = len(x)
  n1 = len(x[0])
  return zerofloat(n1,n2)

def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(80.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2)
  div(x,sqrt(g),y)
  return y

def readImage(n1,n2,name):
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

#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET
backgroundColor = Color(0xfd,0xfe,0xff) # easy to make transparent
def jetFill(alpha):
  return ColorMap.setAlpha(ColorMap.JET,alpha)

def jetFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.JET,a)
def bwrNotch(alpha):
  a = zerofloat(256)
  for i in range(len(a)):
    if i<128:
      a[i] = alpha*(128.0-i)/128.0
    else:
      a[i] = alpha*(i-127.0)/128.0
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,a)


def bwrFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,a)

def jetRamp(alpha):
  return ColorMap.setAlpha(ColorMap.JET,rampfloat(0.0,alpha/256,256))

def bwrRamp(alpha):
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,rampfloat(0.0,alpha/256,256))

def grayRamp(alpha):
  return ColorMap.setAlpha(ColorMap.GRAY,rampfloat(0.0,alpha/256,256))

def plotTensors(g,s1,s2,d=None,dscale=1,ne=20,mk=None,cmin=0,cmax=0,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setBackground(backgroundColor)
  sp.setHLabel("Samples")
  sp.setVLabel("Samples")

  sp.setHInterval(200)
  sp.setVInterval(200)
  sp.setFontSize(24)
  #sp.setFontSizeForPrint(8,240)
  #sp.setFontSizeForSlide(1.0,0.9)
  sp.setSize(800,750)
  pv = sp.addPixels(s1,s2,g)
  pv.setColorModel(ColorMap.GRAY)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  else:
    pv.setPercentiles(1,99)
  if d:
    tv = TensorsView(s1,s2,d)
    tv.setOrientation(TensorsView.Orientation.X1DOWN_X2RIGHT)
    tv.setLineColor(Color.RED)
    tv.setLineWidth(3)
    if(mk):
      tv.setEllipsesDisplayed(mk)
    else:
      tv.setEllipsesDisplayed(ne)
    tv.setScale(dscale)
    tile = sp.plotPanel.getTile(0,0)
    tile.addTiledView(tv)
  if pngDir and png:
    sp.paintToPng(360,3.3,pngDir+png+".png")
    #sp.paintToPng(720,3.3,pngDir+png+".png")

def plot2(s1,s2,f,g=None,u1=None,u2=None,dx=10,scale=1,
        cmin=None,cmax=None,cmap=None,label=None,png=None):
  n2 = len(f)
  n1 = len(f[0])
  f1,f2 = s1.getFirst(),s2.getFirst()
  d1,d2 = s1.getDelta(),s2.getDelta()
  panel = panel2Teapot()
  panel.setHInterval(200)
  panel.setVInterval(200)
  panel.setHLabel("Samples")
  panel.setVLabel("Samples")
  #panel.setHInterval(100.0)
  #panel.setVInterval(100.0)
  #panel.setHLabel("Pixel")
  #panel.setVLabel("Pixel")
  #panel.setColorBarWidthMinimum(80)
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
  if (u1 and u2):
    x1 = zerofloat(2)
    x2 = zerofloat(2)
    for i2 in range(dx,n2-dx,dx):
      for i1 in range(dx,n1-dx,dx):
        x2[0] = (i2-u2[i2][i1]*scale)*d2+f2
        x2[1] = (i2+u2[i2][i1]*scale)*d2+f2
        x1[0] = (i1-u1[i2][i1]*scale)*d1+f1
        x1[1] = (i1+u1[i2][i1]*scale)*d1+f1
        pvu = panel.addPoints(x1,x2)
        pvu.setLineWidth(3)
        pvu.setLineColor(Color.RED)
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
  frame.setSize(800,750)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(400,3.2,pngDir+png+".png")
  return frame


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
