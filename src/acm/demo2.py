import sys

from java.awt import *
from java.io import *
from java.lang import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from acm import *

pngDir = None
#pngDir = "./png/"

seismicDir = "../../data/"
ffile = "tp13m"
f1,f2 = 0,0
d1,d2 = 1,1
s1 = Sampling(150,d1,f1)
s2 = Sampling(150,d2,f2)
n1,n2 = s1.count,s2.count

def main(args):
  f = goFakeImage()
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply1X(f,u1)
  rgf.applyX1(f,u2)
  wp = add(pow(u1,2),pow(u2,2))
  '''
  lof = LocalOrientFilter(6.0,6.0)
  lof.applyForNormalLinear(f,u1,u2,el)
  '''
  wp = div(wp,max(wp))
  wp = pow(wp,6)
  gvf = GradientVectorFlow()
  us = gvf.applyForGVF(u1,u2,wp)
  print max(us)
  print min(us)
  plot(f)
  plot(wp)
  plotVectors(20,u1,u2,f)
  plotVectors(40,us[0],us[1],f)
  

def goFakeImage():
  f = zerofloat(n1,n2)
  for i2 in range(30,40):
    f[i2][30] = 1
    f[i2][60] = 1
  for i2 in range(45,56):
    f[i2][30] = 1
    f[i2][60] = 1
  for i1 in range(30,60):
    f[30][i1] = 1
    f[55][i1] = 1
  return f

def plotVectors(dh,u1,u2,f):
  np = int(n1*n2)
  x1 = zerofloat(2,np)
  x2 = zerofloat(2,np)
  i = 0
  for i2 in range(5,n2-2,6):
    for i1 in range(1,n1-5,6):
      x1[i][0] = i1*d1+f1
      x1[i][1] = (i1+dh*u1[i2][i1])*d1+f1
      x2[i][0] = i2*d2+f2
      x2[i][1] = (i2+dh*u2[i2][i1])*d2+f2
      i = i+1
  x1 = copy(2,i,0,0,x1)
  x2 = copy(2,i,0,0,x2)
  plot(f,[x1,x2])
 
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

def plot(f,xp=None,png=None): 
  orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT;
  panel = PlotPanel(1,1,orientation);
  #panel.setVInterval(0.2)
  pxv = panel.addPixels(0,0,s1,s2,f);
  pxv.setColorModel(ColorMap.GRAY)
  #panel.setTitle("normal vectors")
  if xp:
    ptv = panel.addPoints(xp[0],xp[1])
    ptv.setLineColor(Color.RED)
    ptv.setLineWidth(3.0)
  cb = panel.addColorBar();
  cb.setInterval(1.0)
  cb.setLabel("Amplitude")
  panel.setColorBarWidthMinimum(130)
  moc = panel.getMosaic();
  frame = PlotFrame(panel);
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setTitle("normal vectors")
  frame.setVisible(True);
  frame.setSize(790,660)
  frame.setFontSizeForSlide(0.5,1.0,16.0/9.0)
  if pngDir and png:
    frame.paintToPng(300,3.333,pngDir+png+".png")

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
