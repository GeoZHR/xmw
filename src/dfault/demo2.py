"""
Demo of generating avo models
Author: Xinming Wu, University of Texas at Austin
Version: 2017.05.05
"""

from utils import *
setupForSubset("fake2d")
#setupForSubset("tp")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

gxfile = "gx"


# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
pngDir = None
pngDir = "../../../png/dfault/fake/2d/"
plotOnly = False
# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goTrainData()
  #goTestImage(140)
  #goFaultStack(440)
  goFaultStackX(20)
  #goNormalize()
  #goMarkSmooth()
  #goLabelCheck(440)
  #goLabel()
  #goF3dExample()
def goF3dExample():
  #gx = readImage2D(101,102,"gx50")
  el = zerofloat(222,440)
  u1 = zerofloat(222,440)
  u2 = zerofloat(222,440)
  gx = readImage2D(222,440,"f3d75s")
  lof = LocalOrientFilter(4,1)
  est = lof.applyForTensors(gx)
  dst = DstCoherence(est,8)
  dst.setEigenvalues(1,0.2)
  dst.applyForLinear(gx,el)
  lof.applyForNormalLinear(gx,u1,u2,el)
  el = pow(el,8)
  el = sub(1,el)
  hp = Helper()
  et = hp.thin(el)
  sd = hp.pickSeeds(4,0.6,et)
  print len(sd)
  print len(sd[0])
  writeImageL("f3dSeeds",sd)
  tx = hp.setTestImage(32,32,sd,gx)
  print len(tx)
  print len(tx[0])
  print len(tx[0][0])
  writeImageL("seisF3d",tx)
  plot(gx,et,cmin=0.01,cmax=0.8,cmap=jetRamp(1.0),png="Fault probability")
  '''
  st = fd.setTestImage(32,32,gx)
  writeImageL("seisF3d",st)
  plot(gx)
  plot(gx)
  '''

def goTestExample():
  fd = FakeData2(5,2,25)
  #gx = readImage2D(101,102,"gx50")
  gx = readImage2D(222,440,"f3d75s")
  st = fd.setTestImage(32,32,gx)
  writeImageL("seisF3d",st)
  plot(gx)
def goFaultStackX(ip):
  fd = FakeData2(5,2,25)
  fb = fd.getFilterBank(32,32)
  gx = readImage2D(222,440,"f3d75s")
  sd = readImage2DL(2,818,"f3dSeeds")
  fx = readImage2DL(818,1,"faultF3d"+str(ip))
  ft = readImage2DL(818,1,"dipF3d"+str(ip))
  fs,tt = fd.faultStack(222,440,sd,ft[0],fx[0],fb)
  fs = sub(fs,min(fs))
  fs = div(fs,max(fs))
  plot(gx,png="seisF3d")
  plot(gx,tt,cmin=min(tt),cmax=max(tt),cmap=jetFill(1.0),png="fdF3d")
  plot(gx,fs,cmin=0.01,cmax=0.8,cmap=jetRamp(1.0),png="sfpF3d")

def goFaultStack(ip):
  fd = FakeData2(5,2,25)
  fb = fd.getFilterBank(32,32)
  nb = len(fb)
  for ib in range(nb):
    bmax = max(fb[ib])
    print bmax
    plot(fb[ib],cmin=0.0,cmax=bmax,png="bank"+str(ib))
  gx = readImage2D(101,102,"gx50")
  fx = readImage2DL(101,102,"fault"+str(ip))
  ft = readImage2DL(101,102,"dip"+str(ip))
  fs,tt = fd.faultStack(ft,fx,fb)
  fs = sub(fs,min(fs))
  fs = div(fs,max(fs))
  plot(gx,png="seis")
  plot(gx,fx,cmin=0.01,cmax=0.8,cmap=jetRamp(1.0),png="Fault probability")
  plot(gx,tt,cmin=min(tt),cmax=max(tt),cmap=jetFill(1.0),png="Fault dip")
  plot(gx,fs,cmin=0.01,cmax=0.8,cmap=jetRamp(1.0),png="Stacked fault probability")

def goTestImage(ip):
  gx = readImage2DL(101,102,"gx56")
  fx = readImage2DL(101,102,"fault"+str(ip))
  plot(gx)
  plot(gx,fx,cmin=0.01,cmax=0.8,cmap=jetRamp(1.0))
def goLabel():
  np = 1000
  for ip in range(np):
    gb = zerofloat(100)
    gx = readImage2DL(22,100,"mark"+str(ip))
    for k in range(100):
      gb[k] = max(gx[k])
    writeImageL("label"+str(ip),gb)

def goNormalize():
  np = 1000
  for ip in range(np):
    gx = readImage3DL(28,28,100,"seis"+str(ip))
    gx = sub(gx,min(gx))
    gx = div(gx,max(gx))
    writeImageL("datan"+str(ip),gx)
def goMarkSmooth():
  np = 101
  rgf = RecursiveGaussianFilterP(1)
  for ip in range(np):
    gs = zerofloat(22,100)
    gx = readImage2DL(22,100,"labeld"+str(ip))
    for k in range(100):
      if(max(gx[k])>0):
        gk = zerofloat(22)
        rgf.apply0(gx[k],gk)
        gs[k] = div(gk,sum(gk))
    writeImageL("labelds"+str(ip),gs)
def goLabelCheck(ip):
  sx = readImage3DL(32,32,100,"seisd"+str(1))
  gx = readImage2DL(23,100,"labeld"+str(1))
  gy = readImage2DL(23,100,"test"+str(ip))
  fd = FakeData2(5,2,25)
  fx = fd.dipToFaultImage(32,32,gx)
  fy = fd.dipToFaultImage(32,32,gy)
  for k in range(0,100,3):
    plot(sx[k],fx[k],cmin=0.1,cmax=0.5,cmap=jetFillExceptMin(1.0),png="true"+str(k))
    plot(sx[k],fy[k],cmin=0.1,cmax=0.5,cmap=jetFillExceptMin(1.0),png="predict"+str(k))
  #plot(gx)
  #plot(gy)

def goTrainData():
  n3 = 100
  np = 1000
  fd = FakeData2(5,2,25)
  for ip in range(np):
    print ip
    gx = zerofloat(32,32,n3)
    #el = zerofloat(32,32,n3)
    mark = fd.getTrainDataAndLabels(gx)
    writeImageL("seisd"+str(ip),gx)
    #writeImageL("attr"+str(ip),el)
    writeImageL("labeld"+str(ip),mark)

def like(x):
  n2 = len(x)
  n1 = len(x[0])
  return zerofloat(n1,n2)

def like(x):
  n3 = len(x)
  n2 = len(x[0])
  n1 = len(x[0][0])
  return zerofloat(n1,n2,n3)

def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(10.0)
  ref.apply1(g,g)
  y = like(x)
  div(x,sqrt(g),y)
  return y

#############################################################################
# graphics
#############################################################################
# plotting
backgroundColor = Color.WHITE
cjet = ColorMap.JET
alpha = fillfloat(1.0,256); alpha[0] = 0.0
ajet = ColorMap.setAlpha(cjet,alpha)

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

def plot(f,g=None,ps=None,t=None,cmap=None,cmin=None,cmax=None,cint=None,
        clab=None,neareast=False,png=None): 
  orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT;
  n1,n2=len(f[0]),len(f)
  s1,s2=Sampling(n1),Sampling(n2)
  panel = PlotPanel(1,1,orientation)#,PlotPanel.AxesPlacement.NONE)
  panel.setVInterval(5)
  panel.setHInterval(5)
  '''
  panel.setHLabel("Image index in a batch")
  panel.setVLabel("Fault dip")
  '''
  pxv = panel.addPixels(0,0,s1,s2,f);
  pxv.setColorModel(ColorMap.GRAY)
  pxv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if g:
    pxv.setClips(-2,2)
  else:
    if cmin and cmax:
      pxv.setClips(cmin,cmax)
  if g:
    pv = panel.addPixels(s1,s2,g)
    pv.setInterpolation(PixelsView.Interpolation.LINEAR)
    pv.setColorModel(cmap)
    if cmin and cmax:
      pv.setClips(cmin,cmax)
  if ps:
    uv = panel.addPoints(0,0,ps[0],ps[1])
    uv.setLineColor(Color.YELLOW)
    uv.setLineWidth(2)
  if clab:
    panel.addColorBar(clab)
  moc = panel.getMosaic();
  frame = PlotFrame(panel);
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  #frame.setTitle("normal vectors")
  frame.setVisible(True);
  #frame.setSize(1400,700)
  frame.setSize(500,460)
  frame.setFontSize(12)
  if pngDir and png:
    frame.paintToPng(720,3.333,pngDir+png+".png")


#############################################################################
run(main)
