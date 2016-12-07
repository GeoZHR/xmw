import sys

from java.awt import *
from java.io import *
from java.nio import *
from java.lang import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.util import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util.ArrayMath import *

from util import *
from sso import *

pngDir = None
seismicDir = "../../../data/seis/sso/2d/pwd/"
fxfile = "fx70" #3d image
pkfile = "pk70" #2d image
f1,f2,f3 = 0,0,0
d1,d2,d3 = 1,1,1
n1,n2,n3 = 451,101,153
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)
s3 = Sampling(n3,d3,f3)
plotOnly = False

def main(args):
  goPwd()
def goPwd():
  fx = readImage2d(n1,n2,fxfile)
  pk = readImage2d(n1,n2,pkfile)
  pwd = PlaneWaveDestructor(-5,5)
  pwd.setSmoothness(2,0)
  p2 = pwd.findSlopes(fx)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  lof = LocalOrientFilterP(8,1)
  lof.applyForNormal(fx,u1,u2)
  p2s = mul(-1,div(u2,u1))
  plot2(fx)
  plot2(pk,cmap=ColorMap.JET,cmin=-1.5,cmax=1.5,label="pk")
  plot2(p2,cmap=ColorMap.JET,cmin=-1.5,cmax=1.5,label="p2")
  plot2(p2s,cmap=ColorMap.JET,cmin=-1.5,cmax=1.5,label="p2s")
  plot2(abs(sub(p2,pk)),cmap=ColorMap.JET,cmin=0.0001,cmax=0.2,label="dp")
  plot2(abs(sub(p2s,pk)),cmap=ColorMap.JET,cmin=0.0001,cmax=0.2,label="dps")

def goSmooth2d():
  fx = readImage2d(n1,n2,fxfile)     #input seismic
  sigma1,sigma2=4,2 #smoothing parameters in computing structure tensors
  lof = LocalOrientFilter(sigma1,sigma2)
  ets = lof.applyForTensors(fx) #compute structure tensors
  ets.setEigenvalues(0.001,1.0) #set smoothing to be parallel not perpendicular to reflections
  small,niter=0.001,200 #set CG method in structure-oriented smoothing filter
  lsf = LocalSmoothingFilter(small,niter)#structure-oriented smoothing filter 
  scale = 20 #smoothing extent
  fs = zerofloat(n1,n2) #initialize an output array
  lsf.apply(ets,scale,fx,fs)
  plot2(fx)
  plot2(fs)

def goSmooth3d():
  fx = readImage(fxfile)     #input seismic
  sigma1,sigma2,sigma3=4,2,2 #smoothing parameters in computing structure tensors
  lof = LocalOrientFilter(sigma1,sigma2,sigma3)
  ets = lof.applyForTensors(fx) #compute structure tensors
  ets.setEigenvalues(0.001,1.0,1.0) #set smoothing to be parallel not perpendicular to reflections
  small,niter=0.001,200 #set CG method in structure-oriented smoothing filter
  lsf = LocalSmoothingFilter(small,niter)#structure-oriented smoothing filter 
  scale = 10 #smoothing extent
  fs = zerofloat(n1,n2,n3) #initialize an output array
  lsf.apply(ets,scale,fx,fs)
  plot3(fx)
  plot3(fs)

def readImage(basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  image = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def readImage2d(n1,n2,basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image


def writeImage(basename,image):
  """ 
  Writes an image to a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()
  return image

#############################################################################
# graphics

def addColorBar(frame,clab=None,cint=None):
  cbar = ColorBar(clab)
  if cint:
    cbar.setInterval(cint)
  cbar.setFont(Font("Arial",Font.PLAIN,32)) # size by experimenting
  cbar.setWidthMinimum
  cbar.setBackground(Color.WHITE)
  frame.add(cbar,BorderLayout.EAST)
  return cbar
def plot2(f,cmap=ColorMap.GRAY,cmin=None,cmax=None,label=None,png=None): 
  orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT # for pnz data
  n2 = len(f)
  n1 = len(f[0])
  s1,s2=Sampling(n1),Sampling(n2)
  panel = PlotPanel(1,1,orientation)
  panel.setVInterval(50)
  panel.setHInterval(50)
  pxv = panel.addPixels(0,0,s1,s2,f);
  pxv.setColorModel(cmap)
  pxv.setInterpolation(PixelsView.Interpolation.LINEAR)
  if cmin and cmax:
    pxv.setClips(cmin,cmax)
  cb = panel.addColorBar();
  cb.setInterval(1.0)
  if label:
    cb.setLabel(label)
  panel.setColorBarWidthMinimum(50)
  moc = panel.getMosaic();
  frame = PlotFrame(panel);
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  #frame.setTitle("normal vectors")
  frame.setVisible(True);
  frame.setSize(round(n2*3)+50,round(n1*3))
  frame.setFontSize(24)
  if pngDir and png:
    frame.paintToPng(720,3.333,pngDir+png+".png")

def plot3(f,cmin=None,cmax=None,cmap=None,clab=None,cint=None,png=None):
  n3 = len(f)
  n2 = len(f[0])
  n1 = len(f[0][0])
  s1,s2,s3=Sampling(n1),Sampling(n2),Sampling(n3)
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  sf = SimpleFrame(AxesOrientation.XRIGHT_YOUT_ZDOWN)
  cbar = None
  ipg = sf.addImagePanels(s1,s2,s3,f)
  if cmap!=None:
    ipg.setColorModel(cmap)
  if cmin!=None and cmax!=None:
    ipg.setClips(cmin,cmax)
  else:
    ipg.setClips(-2.0,2.0)
  if clab:
    cbar = addColorBar(sf,clab,cint)
    ipg.addColorMapListener(cbar)
  if cbar:
    cbar.setWidthMinimum(120)
  ipg.setSlices(401,24,24)
  if cbar:
    sf.setSize(887,700)
  else:
    sf.setSize(750,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.48*sqrt(n1*n1+n2*n2+n3*n3)
  zscale = 0.80*max(n2*d2,n3*d3)/(n1*d1)
  ov = sf.getOrbitView()
  ov.setAxesScale(1.0,1.0,zscale)
  #ov.setWorldSphere(BoundingSphere(0.5*n1,0.4*n2,0.4*n3,radius))
  #ov.setAzimuthAndElevation(140.0,40.0)
  #ov.setTranslate(Vector3(-0.06,0.12,-0.27))
  #ov.setScale(1.25)
  ov.setScale(2.25)
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
