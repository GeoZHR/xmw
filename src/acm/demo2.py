import sys

from java.awt import *
from java.io import *
from java.lang import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.interp import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from acm import *

pngDir = None
pngDir = "../../../png/acm/"

seismicDir = "../../../data/seis/acm/"
fxfile = "atwj1s"
#fxfile = "f3d267"
#fxfile = "f3d267Sub"
edfile = "edge"
f1file = "f1"
f2file = "f2"
u1file = "u1"
u2file = "u2"
#ffile = "tp73"
f1,f2 = 0,0
d1,d2 = 1,1
n1,n2 = 251,357
n1,n2 = 500,500
#n1,n2 = 462,951
#n1,n2 = 140,350
#n1,n2 = 100,101
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)

def main(args):
  #goGVF()
  #goSnake()
  #goSnakeReal()
  #goPaint()
  #goForce()
  #goContour()
  
  #goExForce()
  #goSnakeReal()
  #goExForceF3d()
  goChannel()
def goChannel():
  fx = readImage(fxfile)
  fx = div(fx,10000)
  el = zerofloat(n1,n2)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  g1 = zerofloat(n1,n2)
  g2 = zerofloat(n1,n2)
  lof = ChannelEnhanceFilter(4,4)
  ets = lof.applyForTensors(fx)
  lof.applyForNormal(fx,u1,u2)
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply10(fx,g1)
  rgf.apply01(fx,g2)
  g1 = mul(g1,g1)
  g2 = mul(g2,g2)
  sc = sqrt(add(g1,g2))
  sc = sub(sc,min(sc))
  sc = div(sc,max(sc))
  ets.setEigenvalues(0.01,1.0)
  etp = lof.applyForLinear(8,ets,sc,fx,el)
  plot(fx)
  plot(sc)
  plot(el)

  
def goExForceF3d():
  sig,sig1,sig2=1,4,4
  fx = readImage(fxfile)
  ws = zerofloat(n1,n2)
  gs = zerofloat(n1,n2,2)
  gvf = GradientVectorFlow()
  ed = gvf.applyForEdge(sig,sig1,sig2,fx)
  et = gvf.applyForGVX(sig,sig1,sig2,ed,gs,ws) #channel 1
  #gvf.setScale(0.01)
  gvf.setScale(0.0)
  wp = pow(ed,1.0)
  ws = pow(ws,1.0)
  wm = pow(wp,1)
  wm = div(wm,max(wm))
  wm = sub(1,wm)
  fs = gvf.applyForGVF(None,gs[0],gs[1],ws,wm)
  writeImage(edfile,ed)
  writeImage(f1file,fs[0])
  writeImage(f2file,fs[1])
  plot(fx)
  plot(wp)
  plot(ws)
  plot(wm)
  plotVectors(10,5,gs[0],gs[1],fx)
  plotVectors(100,5,fs[0],fs[1],fx)

def goExForce():
  sig,sig1,sig2=2,2,2
  fx = readImage(fxfile)
  fx = div(fx,2000)
  ws = zerofloat(n1,n2)
  gs = zerofloat(n1,n2,2)
  gvf = GradientVectorFlow()
  ed = gvf.applyForEdge(sig,sig1,sig2,fx)
  et = gvf.applyForGV(sig,sig1,sig2,ed,gs)
  #et = gvf.applyForGVX(sig,sig1,sig2,ed,gs,ws) #channel 1
  #gvf.setScale(0.01)
  gvf.setScale(0.01)
  wp = pow(ed,1.0)
  fs = gvf.applyForGVF(et,gs[0],gs[1],wp,None)
  writeImage(edfile,ed)
  writeImage(f1file,fs[0])
  writeImage(f2file,fs[1])
  plot(fx)
  plot(wp)
  plotVectors(50,5,gs[0],gs[1],fx)
  plotVectors(50,5,fs[0],fs[1],fx)

def goSnakeReal():
  fx = readImage(fxfile)
  f1 = readImage(f1file)
  f2 = readImage(f2file)
  fx = div(fx,2000)
  #c1,c2=418,345
  c1,c2 = 248,165
  kamma,gamma = 2.0,0.6
  alpha,beta  = 0.1,0.1
  fp = zerofloat(n1,n2,2)
  fn = zerofloat(n1,n2,2)
  ac = ActiveContour(s1,s2,c1,c2,5)
  ac.setParameters(kamma,gamma,alpha,beta)
  ac.exForce(10,25,fx,fp,fn)
  fn = mul(fn,0.5)
  ac.firstUpdate(fp,fn)
  snake = ac.getSnake()
  x1 = snake.getArrayX1()
  x2 = snake.getArrayX2()
  pngName = "2ndChannel"
  plot(fx,ap=[x1,x2],png=None)#pngName+"Initial")
  
  ac.setIterations(10)
  fnn = [f1,f2]
  fnn=mul(fnn,10)
  fpn = mul(fp,0.5)
  for k in range(25):
    ac.releaseSnake(fpn,fnn,False)
    x1 = snake.getArrayX1()
    x2 = snake.getArrayX2()
    k1,k2=[],[]
    for j in range(len(x1)):
      x1i = x1[j]
      x2i = x2[j]
      if(x1i>1 and x1i<n1-6):
       k1.append(x1i)
       k2.append(x2i)
    plot(fx,ap=[k1,k2],png=None)#pngName+str(k))


def goForce():
  #c1,c2 = 248,165 # channel 1
  c1,c2 = 418,345 # channel 2
  #c1,c2 = 220,110
  #c1,c2 = 345,340
  #c1,c2 = 250,250
  g = readImage(ffile)
  g = div(g,2000)
  fp = zerofloat(n1,n2,2)
  fn = zerofloat(n1,n2,2)

  fpn = zerofloat(n1,n2,2)
  fnn = zerofloat(n1,n2,2)
  ups = zerofloat(n1,n2,3)
  gvf = GradientVectorFlow()
  gvf.setScale(0.01)
  et = gvf.applyForGV(2,4,4,g,ups) #channel 1
  wp = pow(ups[2],1)
  fnn = gvf.applyForGVF(et,ups[0],ups[1],wp)
  plot(ups[2])
  plotVectors(100,5,fnn[0],fnn[1],g)

  ac = ActiveContour(s1,s2,c1,c2,5)
  ac.setIterations(1)
  ac.exForce(10,25,g,fp,fn)
  #fn = mul(fn,0.04)
  fn = mul(fn,0.5)
  kamma,gamma = 1.0,0.2
  alpha,beta  = 0.1,0.1
  snake = ac.getSnake()
  for k in range(1):
    ac.releaseSnake(s1,s2,kamma,gamma,alpha,beta,fp,fn,True)
    x1 = snake.getArrayX1()
    x2 = snake.getArrayX2()
  pngName = "2ndChannel"
  plot(g,ap=[x1,x2],png=pngName+"Initial")
  g1 = snake.getArrayU1()
  g2 = snake.getArrayU2()
  [fps,fns]=ac.exForceOnSnake(fp,fn)
  
  acn = ActiveContour(s1,s2,c1,c2,x1,x2)
  acn.setIterations(10)
  fnn=mul(fnn,10)
  fpn = mul(fp,0.5)
  snake = acn.getSnake()
  kamma,gamma = 2.0,0.6
  alpha,beta  = 0.1,5.1
  for k in range(20):
    acn.releaseSnake(s1,s2,kamma,gamma,alpha,beta,fpn,fnn,False)
    x1 = snake.getArrayX1()
    x2 = snake.getArrayX2()
    k1,k2=[],[]
    for j in range(len(x1)):
      x1i = x1[j]
      x2i = x2[j]
      if(x1i>1 and x1i<n1-6):
       k1.append(x1i)
       k2.append(x2i)
    plot(g,ap=[k1,k2],png=pngName+str(k))

def goPaint():
  f,t,x = [1],[250],[165]
  s = readImage(ffile)
  smooth = 1.00 # smoothness for blended gridder
  fnull = -1.0
  sg = SimpleGridder2(f,t,x)
  sg.setNullValue(fnull)
  p = sg.grid(s1,s2)
  tensors = makeImageTensors(s)
  plot2(f,t,x,s,s1,s2,None,"Known value","tp2f",et=None)
  bg = makeBlendedGridder(f,t,x,smooth=smooth)
  bg.setTensors(tensors)
  d = bg.gridNearest(fnull,p)
  plot2(f,t,x,s,s1,s2,d,"Time (samples)","tp2t")

def goSnakeFake():
  sig = 1
  sig1,sig2=2,2
  f = goFakeImage()

  up = zerofloat(n1,n2,3)
  gvf = GradientVectorFlow()
  gvf.setScale(0.5) # fake image
  et = gvf.applyForGV(sig,sig1,sig2,f,up)
  us = gvf.applyForGVF(None,up[0],up[1],up[2])
  plot(f)
  #plot(up[2])
  plotVectors(5,3,up[0],up[1],f)   #fake image
  plotVectors(5,3,us[0],us[1],f)   #fake image
 
  us = mul(us,-1.0)
  x1 = [30,25,25,25,25,25,35,45,45,45,45,45,45,45,30]
  x2 = [20,30,50,60,80,85,85,85,70,65,50,40,30,20,20]
  '''
  acs = ActiveSnakeModel(x1,x2)
  snake = acs.releaseSnake(s1,s2,1.0,0.1,0.1,0.1,us[0],us[1])
  '''
  acs = ActiveContour(x1,x2)
  acs.setIterations(1)
  kamma,gamma = 0.0,0.5
  alpha,beta  = 0.1,0.1
  snake = acs.releaseSnake(s1,s2,kamma,gamma,alpha,beta,us[0],us[1])
  x1s = snake.getArrayX1()
  x2s = snake.getArrayX2()
  plot(f,ap=[x1,x2])
  plot(f,ap=[x1s,x2s])

def goTest():
  f = goFakeImage()
  x = [20,30,50,60,80,85,85,85,85,70,65,50,40,30,20,20,20]
  y = [20,20,20,20,20,20,25,35,50,55,50,50,50,50,50,30,20]
  acs = ActiveSnake(x,y)
  snake = acs.releaseSnake(s1,s2)
  xs = snake.getArrayX1()
  ys = snake.getArrayX2()
  plot(f,ap=[y,x])
  plot(f,ap=[ys,xs])
  

def goFakeImage():
  f = zerofloat(n1,n2)
  for i2 in range(30,50):
    f[i2][30] = 1
    f[i2][40] = 1
  for i2 in range(55,76):
    f[i2][30] = 1
    f[i2][40] = 1
  for i1 in range(30,40):
    f[30][i1] = 1
    f[75][i1] = 1
  return f

def makeBlendedGridder(f,x,y,smooth=0.5,tmax=FLT_MAX,tmx=False):
  bg = BlendedGridder2(f,x,y)
  bg.setSmoothness(smooth)
  #bg.setBlendingKernel(
  #  LocalDiffusionKernel(LocalDiffusionKernel.Stencil.D21))
  if tmax<FLT_MAX:
    bg.setTimeMax(tmax)
  if tmx:
    bg.setTimeMarkerX(tmx)
  return bg

def makeImageTensors(s):
  """ 
  Returns tensors for guiding along features in specified image.
  """
  sig1,sig2 = 10,30
  n1,n2 = len(s[0]),len(s)
  lof = LocalOrientFilter(sig1,sig2)
  t = lof.applyForTensors(s) # structure tensors
  t.setEigenvalues(0.001,0.999)
  '''
  c = coherence(sig1,t,s) # structure-oriented coherence c
  c = clip(0.0,0.99,c) # c clipped to range [0,1)
  t.scale(sub(1.0,c)) # scale structure tensors by 1-c
  t.invertStructure(1.0,1.0) # invert and normalize
  '''
  return t

def coherence(sigma,t,s):
  lsf = LocalSemblanceFilter(sigma,4*sigma)
  return lsf.semblance(LocalSemblanceFilter.Direction2.V,t,s)

def plotVectors1(dh,b1,b2,u1,u2,f):
  np = len(b1)
  p1 = zerofloat(np)
  p2 = zerofloat(np)
  x1 = zerofloat(2,np)
  x2 = zerofloat(2,np)
  for ip in range(0,np,10):
    x1[ip][0] = b1[ip]
    x2[ip][0] = b2[ip]
    x1[ip][1] = b1[ip]+dh*u1[ip]
    x2[ip][1] = b2[ip]+dh*u2[ip]
    p1[ip] = x1[ip][1]
    p2[ip] = x2[ip][1]
  plot(f,xp=[x1,x2],pp=[p1,p2],ap=[b1,b2])


def plotVectors(dl,dh,u1,u2,f):
  np = int(n1*n2)
  p1 = zerofloat(np)
  p2 = zerofloat(np)
  x1 = zerofloat(2,np)
  x2 = zerofloat(2,np)
  i = 0
  for i2 in range(dh*3,n2-dh*2,dh):
    for i1 in range(dh*2,n1-dh*2,dh):
      x1[i][0] = i1*d1+f1
      x1[i][1] = (i1+dl*u1[i2][i1])*d1+f1
      x2[i][0] = i2*d2+f2
      x2[i][1] = (i2+dl*u2[i2][i1])*d2+f2
      p1[i] = x1[i][1]
      p2[i] = x2[i][1]
      i = i+1
  x1 = copy(2,i,0,0,x1)
  x2 = copy(2,i,0,0,x2)
  p1 = copy(i,0,p1)
  p2 = copy(i,0,p2)

  plot(f,xp=[x1,x2],pp=[p1,p2])
 
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

def plot(f,xp=None,pp=None,ap=None,png=None): 
  orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT;
  panel = PlotPanel(1,1,orientation);
  #panel.setVInterval(0.2)
  pxv = panel.addPixels(0,0,s1,s2,f);
  pxv.setColorModel(ColorMap.GRAY)
  #panel.setTitle("normal vectors")
  if xp:
    ptv = panel.addPoints(xp[0],xp[1])
    ptv.setLineColor(Color.RED)
    ptv.setLineWidth(4.0)
  if pp:
    ptv = panel.addPoints(pp[0],pp[1])
    ptv.setLineStyle(PointsView.Line.NONE)
    ptv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    ptv.setMarkColor(Color.BLUE)
    ptv.setMarkSize(2.0)
  if ap:
    ptv = panel.addPoints(ap[0],ap[1])
    ptv.setLineStyle(PointsView.Line.NONE)
    ptv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    ptv.setMarkColor(Color.YELLOW)
    ptv.setMarkSize(2.0)
  cb = panel.addColorBar();
  #cb.setInterval(0.2)
  cb.setLabel("Amplitude")
  panel.setColorBarWidthMinimum(130)
  moc = panel.getMosaic();
  frame = PlotFrame(panel);
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setTitle("normal vectors")
  frame.setVisible(True);
  frame.setSize(890,760)
  #frame.setSize(1190,760)
  frame.setFontSize(36)
  if pngDir and png:
    frame.paintToPng(300,3.333,pngDir+png+".png")

def plot2(f,x1,x2,s,s1,s2,g=None,label=None,png=None,et=None):
  n1 = len(s[0])
  n2 = len(s)
  panel = panel2Teapot()
  #panel.setHInterval(2.0)
  #panel.setVInterval(0.2)
  #panel.setHLabel("Distance (km)")
  #panel.setVLabel("Time (s)")
  #panel.setHInterval(100.0)
  #panel.setVInterval(100.0)
  #panel.setHLabel("Pixel")
  #panel.setVLabel("Pixel")
  if label:
    panel.addColorBar(label)
  else:
    panel.addColorBar()
  panel.setColorBarWidthMinimum(180)
  pv = panel.addPixels(s1,s2,s)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setColorModel(ColorMap.GRAY)
  pv.setClips(min(s),max(s))
  if g:
    alpha = 0.5
  else:
    g = zerofloat(s1.count,s2.count)
    alpha = 0.0
  pv = panel.addPixels(s1,s2,g)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setColorModel(ColorMap.getJet(alpha))
  if label and label[0]=="T":
    #pv.setClips(0.0,1000.0)
    pv.setClips(0.0,400.0)
  else:
    pv.setClips(0.0,1.0)
  cmap = pv.getColorMap()
  if et:
    tv = TensorsView(s1,s2,et)
    tv.setOrientation(TensorsView.Orientation.X1DOWN_X2RIGHT)
    tv.setLineColor(Color.YELLOW)
    tv.setLineWidth(3.0)
    tv.setScale(2.0)
    panel.getTile(0,0).addTiledView(tv)
  else:
    fs,x1s,x2s = makePointSets(cmap,f,x1,x2)
    for i in range(len(fs)):
      color = cmap.getColor(fs[i][0])
      color = Color(color.red,color.green,color.blue)
      pv = panel.addPoints(x1s[i],x2s[i])
      pv.setLineStyle(PointsView.Line.NONE)
      pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
      pv.setMarkSize(10)
      pv.setMarkColor(color)
  frame2Teapot(panel,png)
def panel2Teapot():
  panel = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.NONE)
  return panel
def frame2Teapot(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  #frame.setFontSizeForPrint(8,240)
  #frame.setSize(1240,774)
  frame.setFontSizeForSlide(1.0,0.8)
  frame.setSize(880,700)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(400,3.2,pngDir+"/"+png+".png")
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
