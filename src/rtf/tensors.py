import sys

from java.awt import *
from java.lang import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from rtf import *

def main(args):
  goTensors()
  #goGradients()

def goGradients():
  f,s1,s2 = readAwImage()
  plotAw(f,s1,s2,png="awef")
  rgf = RecursiveGaussianFilter(1.0)
  g1 = copy(f); rgf.apply1X(f,g1); rgf.applyX0(g1,g1)
  g2 = copy(f); rgf.applyX1(f,g2); rgf.apply0X(g2,g2)
  plotAw(g1,s1,s2,png="aweg1")
  plotAw(g2,s1,s2,png="aweg2")
  g11 = mul(g1,g1)
  g12 = mul(g1,g2)
  g22 = mul(g2,g2)
  ggmax = max(max(g11),max(g12),max(g22))
  cmin,cmax = -0.1*ggmax,0.1*ggmax
  plotAw(g11,s1,s2,cmin=cmin,cmax=cmax,png="aweg11")
  plotAw(g12,s1,s2,cmin=cmin,cmax=cmax,png="aweg12")
  plotAw(g22,s1,s2,cmin=cmin,cmax=cmax,png="aweg22")
  rgf = RecursiveGaussianFilter(4.0)
  rgf.apply00(g11,g11)
  rgf.apply00(g12,g12)
  rgf.apply00(g22,g22)
  ggmax = max(max(g11),max(g12),max(g22))
  cmin,cmax = -0.2*ggmax,0.2*ggmax
  plotAw(g11,s1,s2,cmin=cmin,cmax=cmax,png="awes11")
  plotAw(g12,s1,s2,cmin=cmin,cmax=cmax,png="awes12")
  plotAw(g22,s1,s2,cmin=cmin,cmax=cmax,png="awes22")

def goTensors():
  g,s1,s2 = readTpImage()
  n1=s1.getCount()
  n2=s2.getCount()
  u1=zerofloat(n1,n2)
  u2=zerofloat(n1,n2)
  el=zerofloat(n1,n2)
  lof=LocalOrientFilter(4,2)
  lof.applyForNormalLinear(g,u1,u2,el)
  et=makeImageTensors(el)
  plotTp(g,s1,s2)
  plotTp(el,s1,s2,et,dscale=1,png="11")

def makeImageTensors(s):
  """ 
  Returns tensors for guiding along features in specified image.
  """
  sigma = 4
  n1,n2 = len(s[0]),len(s)
  lof = LocalOrientFilter(sigma)
  t = lof.applyForTensors(s) # structure tensors
  ed=edge(s)
  c = coherence(sigma,t,ed) # structure-oriented coherence c
  c = pow(c,10)
  c = sub(c,min(c))
  c = div(c,max(c))
  for i2 in range(n2):
    for i1 in range(n1):
      if c[i2][i1]<0.2:
        c[i2][i1] = 0.001
  t.scale(c) # scale structure tensors by 1-c
  t.invertStructure(1.0,1.0)
  return t

def coherence(sigma,t,s):
  #lsf = LocalSemblanceFilter(sigma,4*sigma)
  lsf = LocalSemblanceFilter(32,4)
  return lsf.semblance(LocalSemblanceFilter.Direction2.V,t,s)
def edge(g):
  g = gain(g)
  n2 = len(g)
  n1 = len(g[0])
  print n1
  print n2
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  g1 = zerofloat(n1,n2)
  g2 = zerofloat(n1,n2)
  rgf = RecursiveGaussianFilter(2.0)
  lof = LocalOrientFilter(4.0)
  lof.applyForNormal(g,u1,u2)
  rgf.apply10(g,g1)
  rgf.apply01(g,g2)
  gu = add(mul(u1,g1),mul(u2,g2))
  gu = abs(gu)
  gu = sub(gu,min(gu))
  gu = div(gu,max(gu))
  #gu = sub(1.0,gu)
  #gu = clip(0.0001,1.0,gu)
  return gu

def gain(x):
  n2 = len(x)
  n1 = len(x[0])
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(100.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2)
  div(x,sqrt(g),y)
  return y



#############################################################################
# plotting

#pngDir = "../../png/" # where to put PNG images of plots
pngDir = None # for no PNG images

backgroundColor = Color(0xfd,0xfe,0xff) # easy to make transparent

def plotAw(g,s1,s2,d=None,dscale=1,cmin=0,cmax=0,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setBackground(backgroundColor)
  sp.setHLabel("Inline (km)")
  sp.setVLabel("Crossline (km)")
  sp.setHInterval(2.0)
  sp.setVInterval(2.0)
  #sp.setFontSizeForPrint(8,240)
  sp.setFontSizeForSlide(1.0,0.9)
  sp.setSize(910,950)
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
    tv.setLineColor(Color.YELLOW)
    tv.setLineWidth(2.0)
    #tv.setScale(16.0)
    tv.setScale(1.0)
    sp.plotPanel.getTile(0,0).addTiledView(tv)
    '''
    tv = TensorsView(s1,s2,d)
    tv.setOrientation(TensorsView.Orientation.X1DOWN_X2RIGHT)
    tv.setLineColor(Color.YELLOW)
    tv.setLineWidth(3)
    tv.setScale(2.0)
    #tv.setEllipsesDisplayed(20)
    #tv.setScale(dscale)
    tile = sp.plotPanel.getTile(0,0)
    tile.addTiledView(tv)
    '''
  if pngDir and png:
    sp.paintToPng(360,3.3,pngDir+png+".png")
    #sp.paintToPng(720,3.3,pngDir+png+".png")

def plotTp(g,s1,s2,d=None,dscale=1,cmin=0,cmax=0,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setBackground(backgroundColor)
  sp.setHLabel("Distance (km)")
  sp.setVLabel("Time (s)")
  sp.setHInterval(2.0)
  sp.setVInterval(0.2)
  #sp.setFontSizeForPrint(8,240)
  sp.setFontSizeForSlide(1.0,0.9)
  sp.setSize(910,670)
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
    tv.setLineColor(Color.YELLOW)
    tv.setLineWidth(3)
    tv.setEllipsesDisplayed(20)
    tv.setScale(dscale)
    tile = sp.plotPanel.getTile(0,0)
    tile.addTiledView(tv)
  if pngDir and png:
    sp.paintToPng(360,3.3,pngDir+png+".png")
    #sp.paintToPng(720,3.3,pngDir+png+".png")

#############################################################################
# data input/output

# John Mathewson's subsets, resampled to 20 x 20 m trace spacing
# atwj1.dat = John's area1 file
# atwj1s.dat = horizontal slice of John's atwj1 for i1=40
# atwj3.dat = John's area3 file
#n1= 129; d1=0.0040; f1=0.0000
#n2= 500; d2=0.0200; f2=0.0000
#n3= 500; d3=0.0200; f3=0.0000

def readAwImage():
  s1 = Sampling(500,0.02,0.0)
  s2 = Sampling(500,0.02,0.0)
  g = readImage("../../../data/seis/acm/atwj1s",s1,s2)
  return g,s1,s2

def readFakeImage():
  s1 = Sampling(152,0.02,0.0)
  s2 = Sampling(153,0.02,0.0)
  g1 = readImage("../../../data/seis/swt/fake/fake110",s1,s2)
  g2 = readImage("../../../data/seis/swt/fake/fake111",s1,s2)
  gs = add(g1,g2)
  gs = div(gs,2)
  return gs,s1,s2


def readTpImage():
  s1 = Sampling(251,0.004,0.500)
  s2 = Sampling(357,0.025,0.000)
  g  = readImage("../../../data/seis/rtf/tp73",s1,s2)
  return g,s1,s2

def readImage(fileName,s1,s2):
  n1,n2 = s1.count,s2.count
  ais = ArrayInputStream(fileName+".dat")
  x = zerofloat(n1,n2)
  ais.readFloats(x)
  ais.close()
  return x

def writeImage(fileName,x):
  aos = ArrayOutputStream(fileName+".dat")
  aos.writeFloats(x)
  aos.close()

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
