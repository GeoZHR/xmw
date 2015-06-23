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

from dnp import *
from fault import *
from util import FakeData

seismicDir = "../../data/"
ffile = "tp13m"
ft = 0.700
fx = 0.250
dt = 0.004
dx = 0.025
s1 = Sampling(150,dt,ft)
s2 = Sampling(150,dx,fx)
#s1 = Sampling(251,1.0,0.0)
#s2 = Sampling(161,1.0,0.0)
n1,n2 = s1.count,s2.count
d1,d2 = s1.delta,s2.delta

def main(args):
  plotVectors()
  flatten()
def plotVectors():
  f = readImage(ffile)
  f = gain(f)
  u1  = zerofloat(n1,n2)  
  u2  = zerofloat(n1,n2)  
  lof = LocalOrientFilter(4.0,1.0)
  lof.applyForNormal(f,u1,u2)
  np = int(n1*n2)
  x1 = zerofloat(2,np)
  x2 = zerofloat(2,np)
  i = 0
  dth = 8
  for i2 in range(2,n2,10):
    for i1 in range(1,n1-12,12):
      x1[i][0] = i1*dt+ft
      x1[i][1] = (i1+dth*u1[i2][i1])*dt+ft
      x2[i][0] = i2*dx+fx
      x2[i][1] = (i2+dth*u2[i2][i1])*dx+fx
      i = i+1
  x1 = copy(2,i,0,0,x1)
  x2 = copy(2,i,0,0,x2)
   
  orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT;
  panel = PlotPanel(1,1,orientation);
  panel.setHLabel("Inline (km)")
  panel.setVLabel("Time (s)")
  panel.setVInterval(0.2)
  pxv = panel.addPixels(0,0,s1,s2,f);
  pxv.setColorModel(ColorMap.GRAY)
  #panel.setTitle("normal vectors")
  ptv = panel.addPoints(x1,x2)
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
  frame.setSize(660,650)
  frame.setFontSizeForSlide(0.5,1.0,16.0/9.0)
  png = "vectors"
  frame.paintToPng(300,3.333,pngDir+png+".png")
 

 
def flatten():
  #f = FakeData.seismic2d2011A(n1,n2,30)
  f = readImage(ffile)
  f = gain(f) 
  #sigma = 1.0 # good for fake data
  sigma = 1.0 # good for Teapot Dome image tp73
  pmax = 10.0
  lsf = LocalSlopeFinder(4.0,1.0,pmax)
  p2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  lsf.findSlopes(f,p2,el)
  p2 = mul(d1/d2,p2)
  el = pow(el,8)
  #plot(s1,s2,el,cmap=jet,title="Weights",png="w")
  plot(s1,s2,f,clab="Amplitude",cint=1.0,title="Input",png="f")
  #plot(s1,s2,p2,cmap=jet,cmin=-0.1,cmax=0.1,)
  fl = Flattener2()
  fl.setWeight1(0.01)
  fl.setIterations(0.01,1000)
  #fl.setSmoothings(4.0,8.0)
  fm = fl.getMappingsFromSlopes(s1,s2,p2,el)
  g = fm.flatten(f)
  h = fm.unflatten(g)
  s = fm.getShiftsS()
  u0 = zerofloat(n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      u0[i2][i1] = ft+i1*dt
  #hv = fm.hv
  #hv = mul(hv,-1.0)
  plot(s1,s2,fm.u1,cmap=jet,clab="Relative geologic time",cint=0.2,title="rgt",png="rgt")
  plot(s1,s2,s,cmap=jet,clab="Shift",cint=0.04,title="shifts",png="shifts")
  plot(s1,s2,u0,cmap=jet,clab="Relative geologic time",cint=0.2,title="initial",png="initial")
  #plot(s1,s2,hv,cmap=jet,title="Horizon volume",png="hv")
  plot(s1,s2,f,u=fm.u1,clab="Relative geologic time",cint=0.2,title="Horizons",png="fu")
  plot(s1,s2,g,clab="Amplitude",cint=1.0,vlabel="Relative geologic time",title="Flattened",png="g")
  plotm(s1,s2,f,u=fm.x1,u1=fm.u1,clab="Relative geologic time",cint=0.2,title="Horizons",png="fu")
  plotm(s1,s2,g,u=fm.x1,u1=fm.u1,flat=True,clab="Relative geologic time",
        cint=0.2,vlabel="Relative geologic time",title="Horizons",png="slices")
    #plot(s1,s2,s,cmap=jet,title="Shifts"+tsuffix,png="s"+psuffix)
    
def faults(f):
  fs = FaultSemblance()
  p2 = fs.slopes(f)
  snd = fs.semblanceNumDen(p2,f)
  fs = FaultScanner2(20,snd)
  (fl,ft) = fs.scan(-20,20)
  #(fl,ft) = fs.thin((fl,ft))
  return fl
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
#pngDir = None
pngDir = "./png/"
def plotm(s1,s2,x,u=None,u1=None,flat=False,cmap=ColorMap.GRAY,clab=None,cint=None,
          vlabel=None,cmin=0,cmax=0,title=None,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  #if title:
    #sp.setTitle(title)
  cb = sp.addColorBar(clab)
  cb.setLabel(clab)
  cb.setInterval(cint)
  sp.setSize(660,650)
  sp.plotPanel.setColorBarWidthMinimum(130)
  pv = sp.addPixels(s1,s2,x)
  pv.setColorModel(cmap)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if u and not flat:
    n2t = n2-6
    fxt = 3*dx+fx
    px = zerofloat(n2t)
    pz = zerofloat(n2t)
    mp = ColorMap(ft,ft+n1*dt,ColorMap.JET)
    for hi in range(19):
      tau = 5 + hi*8
      if(hi==18):
        tau -=3 
      for k in range(n2t):
        px[k] = fxt+k*dx
        pz[k] = u[k+3][tau]
      pv = sp.addPoints(pz,px)
      pv.setLineColor(mp.getColor(tau*dt+ft))
      pv.setLineWidth(4.0)
      cv = sp.addContours(s1,s2,u1)
      cv.setContours(0)
      cv.setColorModel(jet)
  if u and flat:
    n2t = n2-6
    fxt = 3*dx+fx
    px = zerofloat(n2t)
    pz = zerofloat(n2t)
    mp = ColorMap(ft,ft+n1*dt,ColorMap.JET)
    for hi in range(19):
      tau = 5 + hi*8
      if(hi==18):
        tau -=3 
      for k in range(n2t):
        px[k] = fxt+k*dx
        pz[k] = tau*s1.delta+s1.first#u[k+3][tau]
      pv = sp.addPoints(pz,px)
      pv.setLineColor(mp.getColor(tau*dt+ft))
      pv.setLineWidth(4.0)
      cv = sp.addContours(s1,s2,u1)
      cv.setContours(0)
      cv.setColorModel(jet)

    #cv.setLineColor(Color.YELLOW)
  sp.setFontSizeForSlide(0.5,1.0,16.0/9.0)
  sp.setHLabel("Inline (km)")
  sp.setVLabel("Time (s)")
  if vlabel:
    sp.setVLabel(vlabel)
  sp.setVInterval(0.2)
  if pngDir and png:
    sp.paintToPng(300,3.333,pngDir+png+".png")
def plot(s1,s2,x,u=None,cmap=ColorMap.GRAY,clab=None,cint=None,vlabel=None,cmin=0,cmax=0,
         title=None,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  #if title:
    #sp.setTitle(title)
  cb = sp.addColorBar(clab)
  cb.setLabel(clab)
  cb.setInterval(cint)
  sp.setSize(660,650)
  sp.plotPanel.setColorBarWidthMinimum(130)
  pv = sp.addPixels(s1,s2,x)
  pv.setColorModel(cmap)
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if u:
    n1t = n1-2
    n2t = n2-2
    fxt = 1*dx+fx
    ftt = 1*dt+ft
    s1t = Sampling(n1t,dt,ftt)
    s2t = Sampling(n2t,dx,fxt)
    ut = copy(n1t,n2t,1,1,u)
    cv = sp.addContours(s1t,s2t,ut)
    cv.setContours(50)
    cv.setColorModel(jet)
    cv.setLineWidth(4.0)
    #cv.setLineColor(Color.YELLOW)
  sp.setFontSizeForSlide(0.5,1.0,16.0/9.0)
  sp.setHLabel("Inline (km)")
  sp.setVLabel("Time (s)")
  if vlabel:
    sp.setVLabel(vlabel)
  sp.setVInterval(0.2)
  if pngDir and png:
    sp.paintToPng(300,3.333,pngDir+png+".png")

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
