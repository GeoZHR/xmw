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

from mhe import *
_dataDir = "../../../data/seis/mhe/2d/"
_pngDir = "../../../png/mhe/2d/"
_pntDir = None


def main(args):
  #goF3dPick()
  #goF3dPickPredict()
  #goSlopeDispF3d()
  goDtw()
def goDtw():
  ffile = "f3d178Sub2"
  global s1,s2
  global n1,n2
  global seismicDir,pngDir
  pngDir = _pngDir+"f3d/"
  seismicDir = _dataDir+"f3d/"
  s1 = Sampling(242,1.0,0.0)
  s2 = Sampling(951,1.0,0.0)
  n1,n2 = s1.count,s2.count
  f = readImage(ffile)
  f = gain(f)
  f1 = f[540]
  f2 = f[650]
  us = zerofloat(n1)
  dw = DynamicWarping(-35,0)
  dw.setStrainMax(0.4)
  dw.findShifts(f1,f2,us)
  g2 = dw.applyShifts(us,f2)
  plot1s(s1,[f1,f2])
  plot1s(s1,[f1,f2],us=us,png="dtw")
  plotc(s1,s2,f,k21=540,k22=650,us=us,w1=400,w2=650,png="f3dc")
def goF3dPick():
  ffile = "f3d178Sub2"
  global s1,s2
  global n1,n2
  global seismicDir,pngDir
  pngDir = _pngDir+"f3d/"
  seismicDir = _dataDir+"f3d/"
  s1 = Sampling(242,1.0,0.0)
  s2 = Sampling(951,1.0,0.0)
  n1,n2 = s1.count,s2.count
  f = readImage(ffile)
  f = gain(f)
  fs = zerofloat(n1,n2)
  p2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  sigma1,sigma2=8.0,2.0
  pmax = 5
  lsf = LocalSlopeFinder(sigma1,sigma2,pmax)
  lsf.findSlopes(f,p2,el) # estimate slopes and linearity
  wp = pow(el,8)
  he2 = MultigridHorizon2()
  he2.setSmoothing(12)
  he2.setWeight(0.001)
  he2.setExternalIterations(20)
  # pick a horizon with 2 control points
  k11,k12=[189],[806]  
  k21,k22=[170],[880]  
  k31,k32=[128,96],[442,805]  
  k41,k42=[42 ],[330]  
  k51,k52=[65],[756]  
  k51,k52=[59],[753]  
  k1s = [k11,k21,k31,k41,k51]
  k2s = [k12,k22,k32,k42,k52]
  cv1s,cv2s,cv3s = [],[],[]
  p1s,p2s = [],[]
  np = len(k1s)
  for ip in range(np):
    k1 = k1s[ip]
    k2 = k2s[ip]
    cv1i = he2.curveInitialization(n2,n1-1,k1,k2)
    cv2i = copy(cv1i)
    cv1i = he2.curveUpdateFromSlopes(18,35,60,f,wp,p2,k1,k2,cv1i);
    cv2i = he2.curveUpdateFromSlopes(wp,p2,k1,k2,cv2i);
    cv3i = he2.dipPredict1(k1[0],k2[0],p2)
    cv1s.append(cv1i)
    cv2s.append(cv2i)
    cv3s.append(cv3i)
    for ik in range(len(k1)):
      p1s.append(k1[ik])
      p2s.append(k2[ik])
  plot(s1,s2,pow(el,4),w1=400,w2=650,cmin=0.01,cmax=1.0,png="f3dEp")
  #plots(s1,s2,f,hs=cv1s,k2=p2s,k1=p1s,w1=400,w2=650,png="f3dm")
  #plots(s1,s2,f,hs=cv2s,k2=p2s,k1=p1s,w1=400,w2=650,png="f3ds")
  ''' 
  plot(s1,s2,f,w1=400,w2=650,png="f3d")
  plot(s1,s2,f,h=cv1s,k2=p2s,k1=p1s,w1=400,w2=650,png="f3dm")
  plot(s1,s2,f,h=cv2s,k2=p2s,k1=p1s,w1=400,w2=650,png="f3ds")
  p1s[3] = -100
  v1 = zerofloat(n1,n2)
  v2 = zerofloat(n1,n2)
  lof = LocalOrientFilter(8,2)
  lof.apply(f,None,None,None,v1,v2,None,None,None)
  plot(s1,s2,f,h=cv3s,k2=p2s,k1=p1s,w1=400,w2=650,png="f3dp")
  ''' 
  #plot(s1,s2,f,v1=v1,v2=v2,w1=400,w2=650,png="f3dSlope")

def goSlopeDispF3d():
  ffile = "f3d178Sub2"
  global s1,s2
  global n1,n2
  global seismicDir,pngDir
  pngDir = _pngDir+"f3d/"
  seismicDir = _dataDir+"f3d/"
  s1 = Sampling(242,1.0,0.0)
  s2 = Sampling(951,1.0,0.0)
  n1,n2 = s1.count,s2.count
  f = readImage(ffile)
  f = gain(f)
  '''
  v1 = zerofloat(n1,n2)
  v2 = zerofloat(n1,n2)
  lof = LocalOrientFilter(4,2)
  lof.apply(f,None,None,None,v1,v2,None,None,None)
  '''
  p2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  pt = fillfloat(1,n1,n2)
  sigma1,sigma2=8.0,2.0
  pmax = 5
  lsf = LocalSlopeFinder(sigma1,sigma2,pmax)
  lsf.findSlopes(f,p2,el) # estimate slopes and linearity
  plot(s1,s2,f,v1=p2,v2=pt,w1=400,w2=650,png="f3dSlope")

def goSlopeDispCrf():
  ffile = "gx596Sub"
  global s1,s2
  global n1,n2
  global seismicDir,pngDir
  pngDir = _pngDir+"crf/"
  seismicDir = _dataDir+"crf/"
  s1 = Sampling(130,1.0,0.0)
  s2 = Sampling(400,1.0,0.0)
  n1,n2 = s1.count,s2.count
  f = readImage(ffile)
  f = gain(f)
  '''
  v1 = zerofloat(n1,n2)
  v2 = zerofloat(n1,n2)
  lof = LocalOrientFilter(4,2)
  lof.apply(f,None,None,None,v1,v2,None,None,None)
  '''
  p2 = zerofloat(n1,n2)
  el = zerofloat(n1,n2)
  pt = fillfloat(1,n1,n2)
  sigma1,sigma2=8.0,2.0
  pmax = 5
  lsf = LocalSlopeFinder(sigma1,sigma2,pmax)
  lsf.findSlopes(f,p2,el) # estimate slopes and linearity
  plot(s1,s2,f,v1=p2,v2=pt,w1=400,w2=650,png="crfSlope")
def goF3dPickPredict():
  ffile = "f3d178Sub2"
  global s1,s2
  global n1,n2
  global seismicDir,pngDir
  pngDir = _pngDir+"f3d/"
  seismicDir = _dataDir+"f3d/"
  s1 = Sampling(242,1.0,0.0)
  s2 = Sampling(951,1.0,0.0)
  n1,n2 = s1.count,s2.count
  f = readImage(ffile)
  f = gain(f)
  fs = zerofloat(n1,n2)
  p2 = zerofloat(n1,n2)
  wp = zerofloat(n1,n2)
  sigma1,sigma2=8.0,2.0
  pmax = 5
  lsf = LocalSlopeFinder(sigma1,sigma2,pmax)
  lsf.findSlopes(f,p2,wp) # estimate slopes and linearity
  wp = pow(wp,8)
  he2 = MultigridHorizon2()
  he2.setSmoothing(12)
  he2.setWeight(0.001)
  he2.setExternalIterations(20)
  # pick a horizon with 2 control points
  k11,k12=[189],[806]  
  k21,k22=[170],[880]  
  k31,k32=[128],[442]  
  k41,k42=[42 ],[330]  
  k51,k52=[65],[756]  
  k51,k52=[59],[753]  
  k1s = [k11,k21,k31,k41,k51]
  k2s = [k12,k22,k32,k42,k52]
  cv1s,cv2s,cv3s = [],[],[]
  p1s,p2s = [],[]
  np = len(k1s)
  for ip in range(np):
    k1 = k1s[ip]
    k2 = k2s[ip]
    cv3i = he2.dipPredict1(k1[0],k2[0],p2)
    cv3s.append(cv3i)
    p1s.append(k1[0])
    p2s.append(k2[0])
  '''
  plot(s1,s2,f,w1=400,w2=650,png="f3d")
  plot(s1,s2,f,h=cv1s,k2=p2s,k1=p1s,w1=400,w2=650,png="f3dm")
  plot(s1,s2,f,h=cv2s,k2=p2s,k1=p1s,w1=400,w2=650,png="f3ds")
  p1s[3] = -100
  v1 = zerofloat(n1,n2)
  v2 = zerofloat(n1,n2)
  lof = LocalOrientFilter(8,2)
  lof.apply(f,None,None,None,v1,v2,None,None,None)
  plot(s1,s2,f,h=cv3s,k2=p2s,k1=p1s,w1=400,w2=650,png="f3dp")
  '''
  plotm(s1,s2,f,h=cv3s,k2=p2s,k1=p1s,w1=400,w2=650,png="f3dp")
  #plot(s1,s2,f,v1=v1,v2=v2,w1=400,w2=650,png="f3dSlope")

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

def plot1s(s1,ys,us=None,vmin=None,vmax=None,color=Color.BLACK,hcolor=Color.CYAN,
  hlabel="Time (sample)",vlabel="Time (sample)",png=None):
  sp = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
  sf = 1.0
  yf = sf
  sp.setVLimits(-10,50)
  if vmin and vmax:
    sp.setVLimits(vmin,vmax)
  sp.setHLimits(0,n1-1)
  sp.setHInterval(20)
  sp.setVInterval(20)
  y = zerofloat(n1)
  x = rampfloat(0,1,n1)
  for ic in range(len(ys)):
    ya = sum(ys[ic])/n1
    ys[ic] = sub(ys[ic],ya)
    ys[ic] = mul(ys[ic],2)
    ys[ic] = add(ys[ic],yf)
    pv1 = sp.addPoints(s1,ys[ic])
    pv1.setLineColor(color)
    pv1.setLineWidth(3.0)
    yf = yf+40
  if us:
    x1 = zerofloat(2)
    x2 = zerofloat(2)
    for i1 in range(23,n1,10):
      x1[0] = i1
      x1[1] = i1+us[i1]
      x2[0] = ys[0][i1]
      x2[1] = ys[1][round(i1+us[i1])]
      pvu = sp.addPoints(x1,x2)
      pvu.setLineWidth(2)
      pvu.setLineColor(Color.RED)
  sp.setSize(650,200)
  sp.setHLabel(hlabel)
  sp.setVLabel(vlabel)
  sp.setFontSize(14) #for print
  #sp.setFontSize(30) #for slides
  if png and pngDir:
    sp.paintToPng(300,7.0,pngDir+png+".png")
def plotc(s1,s2,x,k21=540,k22=650,us=None,
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

  t11 = rampfloat(0,1,n1)
  t12 = fillfloat(k21,n1)
  t21 = rampfloat(0,1,n1)
  t22 = fillfloat(k22,n1)
  pv1 = sp.addPoints(t11,t12)
  pv2 = sp.addPoints(t21,t22)
  pv1.setLineColor(Color.YELLOW)
  pv1.setLineWidth(2.5)
  pv2.setLineColor(Color.YELLOW)
  pv2.setLineWidth(2.5)
  if us:
    x1 = zerofloat(2)
    x2 = zerofloat(2)
    for i1 in range(23,n1,10):
      x1[0] = i1
      x1[1] = i1+us[i1]
      x2[0] = k21
      x2[1] = k22
      pvu = sp.addPoints(x1,x2)
      pvu.setLineWidth(3)
      pvu.setLineColor(Color.RED)
      for ip in range(len(x1)):
        pv = PointsView([x1[ip]],[x2[ip]])
        pv.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT)
        pv.setLineStyle(PointsView.Line.NONE)
        pv.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE)
        pv.setMarkColor(Color.YELLOW)
        pv.setMarkSize(6)
        pv.setLineWidth(3)
        sp.add(pv)
  sp.setSize(w2,w1)
  sp.setFontSize(14)
  sp.plotPanel.setColorBarWidthMinimum(45)
  if pngDir and png:
    sp.paintToPng(720,2.2222,pngDir+png+".png")

def plots(s1,s2,fx,hs=None,k2=None,k1=None,v1=None,v2=None,
         w1=1000,w2=500,cmap=ColorMap.GRAY,cmin=0,cmax=0,title=None,png=None):
  #if title:
  #  sp.setTitle(title)
  #cb=sp.addColorBar(clab)
  np = len(hs[0])
  for ip in range(np):
    sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
    pv = sp.addPixels(s1,s2,fx)
    sp.setHLimits(0,n2-1)
    sp.setVLimits(0,n1-1)
    sp.setHLabel("Inline (trace number)")
    sp.setVLabel("Time (sample)")
    pv.setColorModel(cmap)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setClips(-2,2)
    for ik in range(len(hs)):
      hi = hs[ik][ip]
      x = rampfloat(0,1,n2)
      pvh = sp.addPoints(hi,x)
      pvh.setLineColor(Color.YELLOW)
      pvh.setLineWidth(3.0)
      sp.add(pvh)
    for ik in range(len(k2)):
      pv = PointsView([k1[ik]],[k2[ik]])
      pv.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT)
      pv.setLineStyle(PointsView.Line.NONE)
      pv.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE)
      pv.setMarkColor(Color.RED)
      pv.setMarkSize(8)
      pv.setLineWidth(4)
      sp.add(pv)
    sp.setSize(w2,w1)
    sp.setFontSize(14)
    sp.plotPanel.setColorBarWidthMinimum(45)
    if pngDir and png:
      sp.paintToPng(720,2.2222,pngDir+png+str(ip)+".png")

def plotm(s1,s2,fx,h=None,k2=None,k1=None,v1=None,v2=None,
         w1=1000,w2=500,cmap=ColorMap.GRAY,cmin=0,cmax=0,title=None,png=None):
  #if title:
  #  sp.setTitle(title)
  #cb=sp.addColorBar(clab)
  if k2 and k1:
    np = len(k2)
    for m2 in range(10,900,10):
      sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
      pv = sp.addPixels(s1,s2,fx)
      sp.setHLimits(0,n2-1)
      sp.setVLimits(0,n1-1)
      sp.setHLabel("Inline (trace number)")
      sp.setVLabel("Time (sample)")
      pv.setColorModel(cmap)
      pv.setInterpolation(PixelsView.Interpolation.NEAREST)
      pv.setClips(-2,2)
      for ip in range(np):
        pv = PointsView([k1[ip]],[k2[ip]])
        pv.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT)
        pv.setLineStyle(PointsView.Line.NONE)
        pv.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE)
        pv.setMarkColor(Color.RED)
        pv.setMarkSize(8)
        pv.setLineWidth(4)
        sp.add(pv)
        hi = h[ip]
        b2 = max(k2[ip]-m2,0)
        e2 = min(k2[ip]+m2,n2-1)
        y = copy(e2-b2+1,b2,hi)
        x = rampfloat(b2,1,e2-b2+1)
        pvh = sp.addPoints(y,x)
        pvh.setLineColor(Color.YELLOW)
        pvh.setLineWidth(3.0)
        sp.add(pvh)
      sp.setSize(w2,w1)
      sp.setFontSize(14)
      sp.plotPanel.setColorBarWidthMinimum(45)
      if pngDir and png:
        sp.paintToPng(360,2.2222,pngDir+"movie/"+png+str(m2)+".png")

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
    dx1,dx2 = 5, 15
    scale = 4
    dx1,dx2 = 10, 20
    scale = 6
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
    colors = [Color.BLUE,Color.CYAN,Color.GREEN,
              Color.ORANGE,Color.MAGENTA,Color.RED]
    for k in range(1):
      x = rampfloat(0,1,n2)
      pv = sp.addPoints(hs[k],x)
      #pv.setLineColor(colors[k])
      pv.setLineColor(color)
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
