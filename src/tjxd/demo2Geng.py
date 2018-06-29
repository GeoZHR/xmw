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
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *
from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.interp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.ogl.Gl import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from he import *
from mhe import *
from tjxd import *
_dataDir = "../../../data/seis/tjxd/2d/interp/"
_pngDir = "../../../png/tjxd/2d/interp/"
_pntDir = None
plotOnly = True
plotOnly = False
uxfile = "ux" #RGT volume
fxfile = "seis"
global s1,s2
global n1,n2
global seismicDir,pngDir
pngDir = _pngDir
seismicDir = _dataDir
d1,f1 = 0.004,0.0
s1 = Sampling(1500,d1,f1)
s2 = Sampling(7000,1,37001)
#s1 = Sampling(1500,1,0)
#s2 = Sampling(7000,1,0)

n1,n2 = s1.count,s2.count


def main(args):
  #goHorizonVolume()
  goInterpolation()

def goHorizonVolume():
  fx = readImage(fxfile)
  if not plotOnly:
    p2 = zerofloat(n1,n2)
    wp = zerofloat(n1,n2)
    sigma1,sigma2=8.0,2.0
    pmax = 5
    lsf = LocalSlopeFinder(sigma1,sigma2,pmax)
    lsf.findSlopes(fx,p2,wp) # estimate slopes and linearity
    wp = pow(wp,4)
    k1 = [200,197,198, 217, 225, 225,212, 209]
    k2 = [190,720,1220,1330,3900,4600,5660,6400]

    k10 = [226,221,225,252, 266, 255, 250]
    k20 = [190,720,1220,1330,3900,5660,6400]

    k11 = [283,283,320,336,322]
    k21 = [360,1290,1400,3900,5550]

    k12 = [320,321, 364, 411, 392,397, 379]
    k22 = [360,1320,1410,3000,3720,3900,5500]


    k13 = [403,453, 486, 490, 498, 450]
    k23 = [280,1750,2240,2940,4750,6350]
    k14 = [406,530, 547, 502, 500, 453]
    k24 = [280,1750,2240,2940,4750,6350]
    k15 = [409,570,690, 603,565,455]
    k25 = [280,1750,2290,2800,4750,6350]
    k16 = [412,580, 825,764, 530, 460]
    k26 = [280,1750,2250,2750,3650,6350]

    k17 = [414,590, 850, 910, 677, 463,459]
    k27 = [280,1750,2250,2700,3650,6350,6810]

    k18 = [415,593, 851, 911, 758, 465,460]
    k28 = [280,1750,2250,2700,3650,6350,6810]

    k19 = [416,595, 852, 912, 1010,1020,570, 465, 462]
    k29 = [280,1750,2250,2700,3300,3650,5350,6350,6810]

    k1s = [k1,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19]
    k2s = [k2,k20,k21,k22,k23,k24,k25,k26,k27,k28,k29]
    hv2 = HorizonVolume2()
    k1 = rampfloat(0,10,100)
    hv2.setSmoothings(12)
    hv2.setWeights(0.001)
    hv2.setExternalIterations(20)
    hvs = hv2.applyForHorizonVolume(k1s,k2s,wp,p2)
    ux = hv2.rgtFromHorizonVolume(n1,10,hvs)
    writeImage(uxfile,ux)
    plot(s1,s2,fx,hs=hvs,k2=k2s,k1=k1s,w1=900,w2=1400)
  else:
    ux = readImage(uxfile)
  #plot2(None,None,None,fx,s1,s2,label="Amplitude",png="fx")
  plot2(None,None,None,fx,s1,s2,g=ux,vmin=0,vmax=100,label="RGT",png="ux")

def goInterpolation():
  ffile = "seis"
  curve,vifile,label,vmin,vmax = "density", "deni","Density (g/cc)", 1.8,2.8
  curve,vifile,label,vmin,vmax = "vs", "vsi","Vs (km/s)", 0.6,2.5
  curve,vifile,label,vmin,vmax = "vp", "vpi","Vp (km/s)", 1.5,4.5
  fx = readImage(ffile)
  ux = readImage(uxfile)
  ri = RgtInterpolator(s1,s2,s2,0.001)
  t,x,f=getLogSamples(curve)
  s1i = Sampling(3000,0.002,0.0)
  f,t,x=ri.regridLogSamples(s1,f,t,x)
  vi = ri.apply([f],[t],[x],ux)
  writeImage(vifile,vi)
  m1,b1 = 660,250
  fx = copy(m1,n2,b1,0,fx)
  vi = copy(m1,n2,b1,0,vi)
  ux = copy(m1,n2,b1,0,ux)
  c1 = Sampling(m1,d1,f1+b1*d1)
  plot2(None,t,x,fx,c1,s2,label="Amplitude",png="fx")
  plot2(None,t,x,fx,c1,s2,g=ux,vmin=0,vmax=110,label="RGT",png="ux")
  plot2(f,t,x,fx,c1,s2,vmin=vmin,vmax=vmax,label=label,png="fx"+curve)
  plot2(f,t,x,fx,c1,s2,g=vi,vmin=vmin,vmax=vmax,label=label,png=vifile)


def getLogSamples(curve):
  #depth density time vp vs
  nt = 33668
  wy = 40679 #km
  NULL = -999.25
  log = readImageX(nt,5,"log")
  t,x,f=[],[],[]
  if (curve=="vp"):
    for it in range(0,nt,1):
      tx = log[2][it] #time
      fx = log[3][it] #vp
      if (tx!=NULL and fx!=NULL):
        x.append(wy)
        t.append(tx*0.001)
        f.append(304.8/fx) #convert from us/ft to km/s
  elif (curve=="vs"):
    for it in range(0,nt,1):
      tx = log[2][it] #time
      fx = log[4][it] #vs
      if (tx!=NULL and fx!=NULL):
        x.append(wy)
        t.append(tx*0.001)
        f.append(fx*0.001) #convert from m/s to km/s
  elif (curve=="density"):
    for it in range(0,nt,1):
      tx = log[2][it] #time
      fx = log[1][it] #density
      if (tx!=NULL and fx!=NULL):
        x.append(wy)
        t.append(tx*0.001)
        f.append(fx)
  else:
    return null
  return t,x,f


def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(20.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2)
  div(x,sqrt(g),y)
  return y

#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET

# plotting
backgroundColor = Color.WHITE
cjet = ColorMap.JET
alpha = fillfloat(1.0,256); alpha[0] = 0.0
ajet = ColorMap.setAlpha(cjet,alpha)

def plot2(f,x1,x2,s,s1,s2,g=None,vmin=1.2,vmax=2.0,
        label=None,png=None,et=None):
  n1 = len(s[0])
  n2 = len(s)
  panel = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT)
    #PlotPanel.AxesPlacement.NONE)
  #panel.setHInterval(100.0)
  #panel.setVInterval(100.0)
  panel.setVLimits(s1.getFirst(),s1.getLast())
  panel.setHLabel("Inline (trace)")
  panel.setVLabel("Time (s)")
  if label:
    panel.addColorBar(label)
  else:
    panel.addColorBar()
  panel.setColorBarWidthMinimum(60)
  pv = panel.addPixels(s1,s2,s)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ColorMap.GRAY)
  pv.setClips(-2,2)
  if g:
    alpha = 0.8
  else:
    g = zerofloat(n1,n2)
    alpha = 0.0
  pv = panel.addPixels(s1,s2,g)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ColorMap.getJet(alpha))
  pv.setClips(vmin,vmax)
  if et:
    tv = TensorsView(s1,s2,et)
    tv.setOrientation(TensorsView.Orientation.X1DOWN_X2RIGHT)
    tv.setLineColor(Color.YELLOW)
    tv.setLineWidth(3.0)
    tv.setScale(2.0)
    panel.getTile(0,0).addTiledView(tv)
  else:
    cmap = ColorMap(vmin,vmax,ColorMap.JET)
    if f and x1 and x2:
      fs,x1s,x2s = makePointSets(cmap,f,x1,x2)
      for i in range(len(fs)):
        color = cmap.getColor(fs[i][0])
        pv = panel.addPoints(x1s[i],x2s[i])
        pv.setLineStyle(PointsView.Line.NONE)
        pv.setMarkStyle(PointsView.Mark.FILLED_SQUARE)
        pv.setMarkSize(5)
        pv.setMarkColor(color)
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSize(18)
  frame.setSize(1600,700)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(720,3.33,pngDir+"/"+png+".png")
  return frame


def makePointSets(cmap,f,x1,x2):
  sets = {}
  for i in range(len(f)):
    if f[i] in sets:
      points = sets[f[i]]
      points[0].append(f[i])
      points[1].append(x1[i])
      points[2].append(x2[i])
    else:
      points = [[f[i]],[x1[i]],[x2[i]]] # lists of f, x1, x2
      sets[f[i]] = points
  ns = len(sets)
  fs = zerofloat(1,ns)
  x1s = zerofloat(1,ns)
  x2s = zerofloat(1,ns)
  il = 0
  for points in sets:
    fl = sets[points][0]
    x1l = sets[points][1]
    x2l = sets[points][2]
    nl = len(fl)
    fs[il] = zerofloat(nl)
    x1s[il] = zerofloat(nl)
    x2s[il] = zerofloat(nl)
    copy(fl,fs[il])
    copy(x1l,x1s[il])
    copy(x2l,x2s[il])
    il += 1
  return fs,x1s,x2s

def plot(s1,s2,x,h=None,hs=None,k2=None,k1=None,v1=None,v2=None,
         w1=1000,w2=500,cmap=ColorMap.GRAY,cmin=0,cmax=0,color=None,title=None,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  #if title:
  #  sp.setTitle(title)
  #cb=sp.addColorBar(clab)
  pv = sp.addPixels(s1,s2,x)
  sp.setHLimits(0,n2-1)
  sp.setVLimits(0,n1-1)
  sp.setHLabel("Inline (trace number)")
  sp.setVLabel("Time (sample)")
  pv.setColorModel(cmap)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmax>0:
    pv.setClips(cmin,cmax)
  else:
    pv.setClips(-2,2)
  if (v1 and v2):
    x1 = zerofloat(2)
    x2 = zerofloat(2)
    dx1,dx2 = 5, 12
    scale = 4
    d1,d2=1,1
    f1,f2=0,0
    for i2 in range(0,n2,dx2):
      for i1 in range(0,n1,dx1):
        x2[0] = (i2-v2[i2][i1]*scale)*d2+f2
        x2[1] = (i2+v2[i2][i1]*scale)*d2+f2
        x1[0] = (i1-v1[i2][i1]*scale)*d1+f1
        x1[1] = (i1+v1[i2][i1]*scale)*d1+f1
        pvu = sp.addPoints(x1,x2)
        pvu.setLineWidth(4)
        pvu.setLineColor(Color.CYAN)

  if h:
    for hi in h:
      x = rampfloat(0,1,n2)
      pv = sp.addPoints(hi,x)
      pv.setLineColor(Color.YELLOW)
      pv.setLineWidth(3.0)
      sp.add(pv)

  if hs:
    nh = len(hs)
    cp = ColorMap(0,nh,ColorMap.JET)
    for k in range(nh):
      x = rampfloat(0,1,n2)
      pv = sp.addPoints(hs[k],x)
      pv.setLineColor(cp.getColor(k))
      pv.setLineWidth(4.0)
      if(k<5):
        pv.setLineStyle(PointsView.Line.SOLID)
      sp.add(pv)
      k = k+1
  if k2 and k1:
    np = len(k2)
    for ip in range(np):
      pv = PointsView([k1[ip]],[k2[ip]])
      pv.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT)
      pv.setLineStyle(PointsView.Line.NONE)
      pv.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE)
      pv.setMarkColor(Color.YELLOW)
      pv.setMarkSize(12)
      pv.setLineWidth(6)
      pv.setMarkSize(8)
      pv.setLineWidth(4)
      sp.add(pv)
  sp.setSize(w2,w1)
  sp.setFontSize(14)
  sp.plotPanel.setColorBarWidthMinimum(45)
  if pngDir and png:
    sp.paintToPng(720,2.2222,pngDir+png+".png")

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
def readImageX(n1,n2,name):
  """ 
  Reads an image from a file with specified name.
  name: base name of image file; e.g., "seis"
  """
  fileName = seismicDir+name+".dat"
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image


def readImageL(name):
  fileName = seismicDir+name+".dat"
  n1,n2 = s1.count,s2.count
  image = zerofloat(n1,n2)
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
