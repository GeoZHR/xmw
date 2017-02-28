#############################################################################
"""
Demo of dynamic warping for well-log flattening
Author: Xinming Wu, University of Texas at Austin
Version: 2016.06.01
"""
from common import * 

#############################################################################
fxfile = "gx" # 
wxfile = "deepSub"
wufile = "deepSubFlatten"
curves  = ["v", "p", "d"]

logDir = "../../../data/seis/hdw/logs/"
pngDir = None
pngDir =  "../../../png/hdw/logs/"

plotOnly = False

lmin,lmax=-350,350

def main(args):
  #goLogs()
  #goDeepLogSub()
  #goDeepLogs()
  #goDeepLogSubFlatten()
  #goFlatten()
  #goDeepPorosityLogs()
  #goPorisityFlatten()
  #goDensityFlatten()
  #goVelocityFlatten()
  #goGammaFlatten()
  goGammaFlattenError()
def goGammaFlatten():
  lmin,lmax=-250,250
  sz,dp = goDeepGammaLogs()
  wh = WellHelper()
  d = 3
  dc = wh.resample(d,dp)
  fk = 0.0003048  #1 ft = 0.0003048 km
  nz = len(dc[0])
  sz = Sampling(nz,d*fk,sz.getFirst())
  nw = len(dc)
  sw = Sampling(nw,1,1)
  ww = WellFlattener(lmin,lmax)
  ww.setStrainMax(0.1)
  ww.setErrorExponent(1)
  ww.setErrorExtrapolation(ww.ErrorExtrapolation.AVERAGE)
  clab4 = "Gamma ray (API units)"
  vlab = "Relative geologic time"
  wh = round(15.5*nw)
  plot(sz,sw,dc,wh=wh,cmin=0.001,cmax=200,hint=5,cbar=clab4,png="sg")
  for w in [1,nw]:
    ww.setGate(w,0.85)
    df = ww.flatten(dc)
    plot(sz,sw,df,wh=wh,cmin=0.001,cmax=200,hint=5,
        vlab=vlab,cbar=clab4,png="sgf"+str(w))

def goGammaFlattenError():
  lmin,lmax=-250,250
  sz,dp = goDeepGammaLogs()
  wh = WellHelper()
  d = 3
  dc = wh.resample(d,dp)
  nz = len(dc[0])
  dc = copy(nz,16,0,0,dc)
  fk = 0.0003048  #1 ft = 0.0003048 km
  sz = Sampling(nz,d*fk,sz.getFirst())
  nw = len(dc)
  sw = Sampling(nw,1,1)
  ww = WellFlattener(lmin,lmax)
  ww.setStrainMax(0.1)
  ww.setErrorExponent(1)
  ww.setErrorExtrapolation(ww.ErrorExtrapolation.AVERAGE)
  clab4 = "Gamma ray (API units)"
  vlab = "Relative geologic time"
  wh = round(15.5*nw)
  plot(sz,sw,dc,wh=wh,cmin=0.001,cmax=200,hint=5,cbar=clab4)
  for w in [1,nw]:
    ww.setGate(w,0.85)
    nl = lmax-lmin+1
    et = zerofloat(nl,nz,2)
    df = ww.flatten(dc,et)
    ut = zerofloat(nz)
    ww.backtrackReverse(10,lmin,et[1],et[0],ut)
    mz = 1531-160+1
    et = copy(nl,mz,2,0,160,0,et)
    ut = copy(mz,160,ut)
    if(w==nw):
      plotc(dtran(et[0]),s=ut,perc=95,clab="Error",png="ep"+str(w))
    plotc(dtran(et[1]),u=ut,perc=90,clab="Accumulated error",png="d"+str(w))
    plotc(dtran(et[0]),perc=95,clab="Error",png="e"+str(w))
    plot(sz,sw,df,wh=wh,cmin=0.001,cmax=200,hint=5,
        vlab=vlab,cbar=clab4)

def goVelocityFlatten():
  lmin,lmax=-200,350
  sz,dp = goDeepVelosityLogs()
  wh = WellHelper()
  d = 3
  dc = wh.resample(d,dp)
  fk = 0.0003048  #1 ft = 0.0003048 km
  nz = len(dc[0])
  sz = Sampling(nz,d*fk,0)#sz.getFirst())
  nw = len(dc)
  sw = Sampling(nw,1,1)
  ww = WellFlattener(lmin,lmax)
  ww.setStrainMax(0.1)
  df = ww.flatten(dc)
  clab4 = "Velocity (km/s)"
  wh = 20*nw
  plot(sz,sw,dc,wh=wh,cmin=2.0,cmax=6.0,hint=2,cbar=clab4,png="sv")
  vlab = "Relative geologic time"
  for w in [1,nw]:
    ww.setGate(w,0.85)
    df = ww.flatten(dc)
    plot(sz,sw,df,wh=wh,cmin=2.0,cmax=6.0,hint=2,vlab=vlab,
        cbar=clab4,png="svf"+str(w))

def goDensityFlatten():
  lmin,lmax=-150,350
  sz,dp = goDeepDensityLogs()
  wh = WellHelper()
  d = 3
  dc = wh.resample(d,dp)
  fk = 0.0003048  #1 ft = 0.0003048 km
  nz = len(dc[0])
  sz = Sampling(nz,d*fk,0)
  nw = len(dc)
  sw = Sampling(nw,1,1)
  ww = WellFlattener(lmin,lmax)
  ww.setGate(nw,0.95)
  ww.setStrainMax(0.1)
  #ww.setErrorExtrapolation(ww.ErrorExtrapolation.AVERAGE)
  ww.setErrorExponent(0.5)
  df = ww.flatten(dc)
  clab = "Density (g/cc)"
  wh = 12*nw
  plot(sz,sw,dc,wh=wh,cmin=2.0,cmax=2.8,hint=5,cbar=clab,png="sd")
  vlab = "Relative geologic time"
  plot(sz,sw,df,wh=wh,cmin=2.0,cmax=2.8,hint=5,vlab=vlab,cbar=clab,png="sdf")

def goPorisityFlatten():
  lmin,lmax=-350,350
  sz,dp = goDeepPorosityLogs()
  wh = WellHelper()
  d = 3
  dc = wh.resample(d,dp)
  fk = 0.0003048  #1 ft = 0.0003048 km
  nz = len(dc[0])
  sz = Sampling(nz,d*fk,0)#sz.getFirst())
  nw = len(dc)
  sw = Sampling(nw,1,1)
  ww = WellFlattener(lmin,lmax)
  ww.setStrainMax(0.1)
  ww.setGate(nw,0.90)
  #ww.setGate(1,0.95)
  df = ww.flatten(dc)
  clab = "Porosity (%)"
  wh = 12*nw
  plot(sz,sw,mul(dc,100),wh=wh,cmin=0.001,cmax=45.,hint=5,cbar=clab,png="sp")
  vlab = "Relative geologic time"
  plot(sz,sw,mul(df,100),wh=wh,cmin=0.001,cmax=45.,hint=5,
          vlab=vlab,cbar=clab,png="spf")

def goDeepLogSubFlatten():
  sz,logs = goDeepLogSub()
  nw = len(logs[0])
  nz = len(logs[0][0])
  #sz = Sampling(nz)
  sw = Sampling(nw,1,1)
  ww = WellFlattener(lmin,lmax)
  ww.setStrainMax(0.1)
  nl = lmax-lmin+1
  flogs = ww.flatten(logs)
  writeImage(wufile,flogs)
  sl = Sampling(nl)
  s1 = Sampling(nz)
  clab1 = "Velocity (km/s)"
  clab2 = "Density (g/cc)"
  clab3 = "Gamma ray (API units)"
  clab4 = "Porosity (%)"
  plot(sz,sw,logs[0],cmin=2.0,cmax=6.0,cbar=clab1,png="v")
  plot(sz,sw,logs[1],cmin=2.0,cmax=2.8,cbar=clab2,png="d")
  plot(sz,sw,logs[2],cmin=0.001,cmax=200,cbar=clab3,png="g")
  plot(sz,sw,mul(logs[3],100),cmin=0.001,cmax=45.,cbar=clab4,png="p")
  vlab = "Relative geologic time"
  plot(sz,sw,flogs[0],cmin=2.0,cmax=6.0,vlab=vlab,cbar=clab1,png="vf")
  plot(sz,sw,flogs[1],cmin=2.0,cmax=2.8,vlab=vlab,cbar=clab2,png="df")
  plot(sz,sw,flogs[2],cmin=0.001,cmax=200,cbar=clab3,png="gf")
  plot(sz,sw,mul(flogs[3],100),cmin=0.001,cmax=45.,vlab=vlab,cbar=clab4,png="pf")

def goDeepLogSub():
  ws = [161,1,6,27,29,33,46,48,53,62,66,85,90,91,101,131,132,155,160,161,166]
  #ws = [1,6,27,29,33,46,48,62,66,85,91,101,131,132,155,160,161,166]
  wd = getDeepLogs()
  wh = WellHelper()
  ndf = zerodouble(3)
  dd = wh.toArray(wd,ndf)
  nw = len(ws)
  nc = len(dd)
  nz = len(dd[0][0])
  ds = zerofloat(nz,nw,nc)
  iw = 0
  for k in ws:
    ds[0][iw] = dd[0][k]
    ds[1][iw] = dd[1][k]
    ds[2][iw] = dd[2][k]
    ds[3][iw] = dd[3][k]
    iw = iw+1
  d = 3
  dc = wh.resample(d,ds)
  fk = 0.0003048  #1 ft = 0.0003048 km
  nz = len(dc[0][0])
  sz = Sampling(nz,d*fk,ndf[2])
  #sz = Sampling(nz)
  sw = Sampling(nw)
  writeImage(wxfile,dc)
  print min(dc)
  '''
  plot(sz,sw,dc[0],cmin=2.0,cmax=6.0,cbar="Velocity (km/s)")
  plot(sz,sw,dc[1],cmin=2.0,cmax=2.8,cbar="Density (g/cc)")
  plot(sz,sw,dc[2],cmin=0.001,cmax=200,cbar="Gamma ray (API units)")
  plot(sz,sw,dc[3],cmin=0.001,cmax=0.45,cbar="Porosity (%)")
  '''
  return sz,dc

def goDeepLogs():
  wd = getDeepLogs()
  wh = WellHelper()
  ndf = zerodouble(3)
  dd = wh.toArray(wd,ndf)
  print ndf[0]
  print ndf[1]
  print ndf[2]
  sz = Sampling(int(ndf[0]),ndf[1],ndf[2])
  sz = Sampling(len(dd[0][0]))
  sw = Sampling(len(dd[0]))
  plot(sz,sw,dd[0],cmin=2.0,cmax=6.0,cbar="Velocity (km/s)")
  plot(sz,sw,dd[1],cmin=2.0,cmax=2.8,cbar="Density (g/cc)")
  plot(sz,sw,dd[2],cmin=0.001,cmax=200,cbar="Gamma ray (API units)")
  plot(sz,sw,dd[3],cmin=0.001,cmax=0.45,cbar="Porosity (%)")

def goDeepPorosityLogs():
  wd = getDeepLogs()
  wh = WellHelper()
  ndf = zerodouble(3)
  dd = wh.toArray(wd,ndf)
  dp = []
  for wp in dd[3]:
    if(max(wp)>0):
      dp.append(wp)
  sw = Sampling(len(dp))
  dr = wh.sortLogs(dp)
  dt = wh.trim(ndf,dr)
  fk = 0.0003048  #1 ft = 0.0003048 km
  sz = Sampling(int(ndf[0]),ndf[1]*fk*2,ndf[2]*fk*2)
  plot(sz,sw,dp,cmin=0.001,cmax=0.45,cbar="Porosity (%)")
  #plot(sz,sw,dr,cmin=0.001,cmax=0.45,cbar="Porosity (%)")
  return sz,dt


def goDeepGammaLogs():
  wd = getDeepLogs()
  wh = WellHelper()
  ndf = zerodouble(3)
  dd = wh.toArray(wd,ndf)
  dp = []
  for wp in dd[2]:
    if(max(wp)>0):
      dp.append(wp)
  sw = Sampling(len(dp))
  dr = wh.sortLogs(dp)
  dt = wh.trim(ndf,dr)
  sz = Sampling(int(ndf[0]),ndf[1],ndf[2])
  #plot(sz,sw,dp,cmin=0.001,cmax=0.45,cbar="Porosity (%)")
  nz = len(dt[0])
  dc = copy(nz,113,0,0,dt)
  #plot(sz,sw,dr,cmin=0.001,cmax=0.45,cbar="Porosity (%)")
  return sz,dc

def goDeepDensityLogs():
  wd = getDeepLogs()
  wh = WellHelper()
  ndf = zerodouble(3)
  dd = wh.toArray(wd,ndf)
  dp = []
  for wp in dd[1]:
    if(max(wp)>0):
      dp.append(wp)
  sw = Sampling(len(dp))
  dr = wh.sortLogs(dp)
  dr0 = dr[0]
  dr[0] = dr[1]
  dr[1] = dr0
  dt = wh.trim(ndf,dr)
  sz = Sampling(int(ndf[0]),ndf[1],ndf[2])
  #plot(sz,sw,dp,cmin=0.001,cmax=0.45,cbar="Porosity (%)")
  #plot(sz,sw,dr,cmin=0.001,cmax=0.45,cbar="Porosity (%)")
  return sz,dt

def goDeepVelosityLogs():
  wd = getDeepLogs()
  wh = WellHelper()
  ndf = zerodouble(3)
  dd = wh.toArray(wd,ndf)
  dp = []
  for wp in dd[0]:
    if(max(wp)>0):
      dp.append(wp)
  sw = Sampling(len(dp))
  dr = wh.sortLogs(dp)
  dt = wh.trim(ndf,dr)
  sz = Sampling(int(ndf[0]),ndf[1],ndf[2])
  #plot(sz,sw,dp,cmin=0.001,cmax=0.45,cbar="Porosity (%)")
  #plot(sz,sw,dr,cmin=0.001,cmax=0.45,cbar="Porosity (%)")
  return sz,dt


def goLogs():
  sz,logs = getLogs(curves)
  print sz.getCount()
  nw = len(logs[0])
  sw = Sampling(nw,1,1)
  plot(sz,sw,logs[0],cmin=2.0,cmax=6.0,cbar="Velocity (km/s)")
  plot(sz,sw,mul(logs[1],100),cmin=0.001,cmax=45.,cbar="Porosity (%)")
  plot(sz,sw,logs[2],cmin=2.0,cmax=2.8,cbar="Density (g/cc)")

def goFlatten():
  sz,logs = getLogs(curves)
  nw = len(logs[0])
  nz = len(logs[0][0])
  sw = Sampling(nw,1,1)
  ww = WellFlattener(lmin,lmax)
  ww.setStrainMax(0.2)
  nl = lmax-lmin+1
  flogs = ww.flatten(logs)
  sl = Sampling(nl)
  s1 = Sampling(nz)
  clab1 = "Velocity (km/s)"
  clab2 = "Porosity (%)"
  clab3 = "Density (g/cc)"
  plot(sz,sw,logs[0],cmin=2.0,cmax=6.0,cbar=clab1,png="v")
  plot(sz,sw,mul(logs[1],100),cmin=0.001,cmax=45.,cbar=clab2,png="p")
  plot(sz,sw,logs[2],cmin=2.0,cmax=2.8,cbar=clab3,png="d")
  vlab = "Relative geologic time"
  plot(sz,sw,flogs[0],cmin=2.0,cmax=6.0,vlab=vlab,cbar=clab1,png="vf")
  plot(sz,sw,mul(flogs[1],100),cmin=0.001,cmax=45.,vlab=vlab,cbar=clab2,png="pf")
  plot(sz,sw,flogs[2],cmin=2.0,cmax=2.8,vlab=vlab,cbar=clab3,png="df")

def getDeepLogs():
  wlName = logDir+"tpwd.dat"
  wldata = WellLog.Data.readBinary(wlName)
  logs = []
  for log in wldata.getAll():
    logs.append(log)
  return logs

def arrangeOrder(x):
  ks = [0,10,6,1,9,5,8,12,11,2,7,4,3]
  y = like(x)
  nw = len(x)
  for iw in range(nw):
    ip = 0
    for ik in range(len(ks)):
      y[iw][ip] = x[iw][ks[ip]]
      ip = ip+1
  return y
  

def like(x):
  n3 = len(x)
  n2 = len(x[0])
  n1 = len(x[0][0])
  return zerofloat(n1,n2,n3)

def gain(x):
  n2 = len(x)
  n1 = len(x[0])
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(5.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2)
  div(x,sqrt(g),y)
  return y

def smooth(sig,u):
  v = copy(u)
  rgf = RecursiveGaussianFilterP(sig)
  rgf.apply0(u,v)
  return v

def smooth2(sig1,sig2,u):
  v = copy(u)
  rgf1 = RecursiveGaussianFilterP(sig1)
  rgf2 = RecursiveGaussianFilterP(sig2)
  rgf1.apply0X(u,v)
  rgf2.applyX0(v,v)
  return v


def normalize(e):
  emin = min(e)
  emax = max(e)
  return mul(sub(e,emin),1.0/(emax-emin))
def dtran(d):
  return transpose(d)



def getLogs(curves):
  nc = len(curves)

  fileName = logDir+curves[0]+"logs.txt"
  ifile = open(fileName,'r+')
  lines = ifile.readlines()
  nz = len(lines)
  ifile.close()
  ifile = open(fileName,'r+')
  c = 0
  start = False
  while (start!=True):
    line = ifile.readline()
    c += 1
    if line == "":
      print 'End of file'
      break
    elif line[0] == '~':
      start = True

  line  = ifile.readline()
  wdata = line.split('\t')
  nl = len(wdata)-1
  logs = zerofloat(nz-c,nl-1,nc)
  depth = zerodouble(nz-c)
  for ic,cv in enumerate(curves):
    fileName = logDir+cv+"logs.txt"
    ifile = open(fileName,'r+')
    c = 0
    start = False
    while (start!=True):
      line = ifile.readline()
      c += 1
      if line == "":
        print 'End of file'
        break
      elif line[0] == '~':
        start = True

    line = ifile.readline()
    wdata = line.split('\t')
    depth[0] = float(wdata[0])
    for l in range(1,nl):
      logs[ic][l-1][0] = float(wdata[l])
    i = 1
    start = False
    while (start!=True):
      line = ifile.readline()
      if line == "":
        start = True
      else:
        wdata = line.split('\t')
        depth[i] = float(wdata[0])
        for l in range(1,nl):
          logs[ic][l-1][i] = float(wdata[l])
        i += 1
    ifile.close()
  zs = Sampling(depth)
  sz = Sampling(zs.count,zs.delta*0.001,zs.first)
  return sz,logs

# read/write images
def readImage(n1,n2,n3,basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = logDir+basename+".dat"
  image = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def writeImage(basename,image):
  """ 
  Writes an image to a file with specified basename
  """
  fileName = logDir+basename+".dat"
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()
  return image


#############################################################################
# plotting
backgroundColor = Color.WHITE
cjet = ColorMap.JET
alpha = fillfloat(1.0,256); alpha[0] = 0.0
ajet = ColorMap.setAlpha(cjet,alpha)
def plot(sz,sw,fx,wh=500,wv=900,cmin=None,cmax=None,
    hint=2,vlab="Depth (km)",cbar=None,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(wh,wv)
  sp.setHInterval(hint)
  sp.setVLabel(vlab)
  sp.setHLabel("Log index")
  sp.addColorBar(cbar)
  sp.plotPanel.setColorBarWidthMinimum(72)
  pv = sp.addPixels(sz,sw,fx)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ajet)
  sp.setFontSize(20)
  if cmin and cmax:
    pv.setClips(cmin,cmax)
  if png and pngDir:
    sp.paintToPng(1080,3.33333,pngDir+png+".png")

def plotc(c,s=None,u=None,cmin=0.0,cmax=0.0,perc=None,clab="Error", png=None):
  n,nlag = len(c[0]),len(c)
  fk = 0.0003048  #1 ft = 0.0003048 km
  f1 = 160*fk
  d1 = fk
  s1 = Sampling(n,d1,f1)
  slag = Sampling(nlag,d1,-d1*(nlag-1)/2)
  panel = PlotPanel(1,1,PlotPanel.Orientation.X1RIGHT_X2UP)
  panel.setHLimits(0,s1.first,s1.last)
  panel.setVLimits(0,slag.first,slag.last)
  panel.setVLabel("Shifts (km)")
  panel.setHLabel("Depth (km)")
  cv = panel.addPixels(0,0,s1,slag,c)
  cv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if perc:
    cv.setPercentiles(100-perc,perc)
  elif cmin<cmax:
    cv.setClips(cmin,cmax)
  if s:
    sv = panel.addPoints(0,0,s1,mul(s,d1))
    sv.setLineColor(Color.WHITE)
    sv.setLineStyle(PointsView.Line.DOT)
    sv.setLineWidth(3)
  if u:
    uv = panel.addPoints(0,0,s1,mul(u,d1))
    uv.setLineColor(Color.WHITE)
    uv.setLineWidth(3)
  cbar = panel.addColorBar(clab)
  panel.setColorBarWidthMinimum(72)
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(backgroundColor)
  frame.setSize(round(len(c[0])*0.8)+72,round(len(c)*0.8))
  frame.setVisible(True)
  frame.setFontSize(20)
  if png and pngDir:
    frame.paintToPng(720,3.33333,pngDir+"/"+png+".png")



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
