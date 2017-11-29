from utils import *
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

# Names and descriptions of image files used below.
gxfile  = "gx" # input image (maybe after bilateral filtering)
gsxfile = "gsx" # image after lsf with fault likelihoods
epfile  = "ep" # eigenvalue-derived planarity
elfile  = "el" # eigenvalue-derived planarity
p2file  = "p2" # inline slopes
p3file  = "p3" # crossline slopes
p2kfile = "p2k" # inline slopes (known)
p3kfile = "p3k" # crossline slopes (known)
flfile  = "fl" # fault likelihood
fpfile  = "fp" # fault strike (phi)
ftfile  = "ft" # fault dip (theta)
klfile  = "kl" # fault likelihood
kpfile  = "kp" # fault strike (phi)
ktfile  = "kt" # fault dip (theta)
fltfile = "flt" # fault likelihood thinned
fptfile = "fpt" # fault strike thinned
fttfile = "ftt" # fault dip thinned
kltfile = "klt" # fault likelihood thinned
kptfile = "kpt" # fault strike thinned
kttfile = "ktt" # fault dip thinned
fs1file = "fs1" # fault slip (1st component)
fs2file = "fs2" # fault slip (2nd component)
fs3file = "fs3" # fault slip (3rd component)
fskbase = "fsk" # fault skin (basename only)
v1file = "v1"
v2file = "v2"
v3file = "v3"
w1file = "w1"
w2file = "w2"
w3file = "w3"
h70file = "shb5_hor70"
h74file = "shb5_hor74"
h76file = "shb5_hor76"
h80file = "shb5_hor80"
h83file = "shb5_hor83"
smfile = "sm"
ddfile = "dd"
pdfile = "pd"
pcfile = "pc"
ncfile = "nc"

# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
minPhi,maxPhi = 0,360
minTheta,maxTheta = 70,89
sigmaPhi,sigmaTheta = 4,30
minTheta,maxTheta = 75,90
sigmaPhi,sigmaTheta = 6,20

# These parameters control the construction of fault skins.
# See the class FaultSkinner for more information.
lowerLikelihood = 0.2
upperLikelihood = 0.5
minSkinSize = 3000

# These parameters control the computation of fault dip slips.
# See the class FaultSlipper for more information.
minThrow = 0.01
maxThrow = 15.0

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
pngDir = None
plotOnly = True
pngDir = "../../../png/gu/"

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goSeisData()
  #goHorizons()
  #goHorizons2D()
  #goInlineSlices()
  #goTimeSlices()
  #goSlopeVectors()
  #goKaustEdge()
  #goSlopes()
  #goAmplitudeCurvature()
  #goPlaneWaveDestruction()
  #goKaustEdgeEnhance()
  goThin()
  #goSemblance()

def goSeisData():
  gx = readImage(gxfile)
  gmin,gmax = -2,2
  gmap = ColorMap.GRAY
  plot3(gx,cmin=gmin,cmax=gmax,cmap=gmap,clab="Amplitude",png="gx")

def goSlopeVectors():
  gx = readImage(gxfile)
  v1 = zerofloat(n1,n2,n3)
  v2 = zerofloat(n1,n2,n3)
  v3 = zerofloat(n1,n2,n3)
  w1 = zerofloat(n1,n2,n3)
  w2 = zerofloat(n1,n2,n3)
  w3 = zerofloat(n1,n2,n3)
  lof = LocalOrientFilter(4,1)
  lof.apply(gx,None,None,None,None,None,v1,v2,v3,w1,w2,w3,None,None,None,None,None)
  writeImage(v1file,v1)
  writeImage(v2file,v2)
  writeImage(v3file,v3)
  writeImage(w1file,w1)
  writeImage(w2file,w2)
  writeImage(w3file,w3)

def goPlaneWaveDestruction():
  gx = readImage(gxfile)
  if not plotOnly:
    ke = KaustEdge()
    pd = ke.planeWaveDestruction(8,12,gx)
    pd = sub(pd,min(pd));
    pd = div(pd,max(pd));
    ke = KaustEdge()
    pd = ke.scale(pd)
    writeImage(pdfile+"12",pd)
  else:
    pd = readImage(pdfile)
  '''
  plot3(gx,pd,cmin=0.2,cmax=1.0,cmap=jetRamp(1.0),png="pd")
  plot3(gx,pd,cmin=0.05,cmax=0.2,cmap=jetRamp(1.0),
        clab="Fault likelihood",png="kl")
  '''

def goKaustEdge():
  gx = readImage(gxfile)
  if not plotOnly:
    v1 = readImage(v1file)
    v2 = readImage(v2file)
    v3 = readImage(v3file)
    w1 = readImage(w1file)
    w2 = readImage(w2file)
    w3 = readImage(w3file)
    dd = zerofloat(n1,n2,n3)
    ke = KaustEdge()
    dd = ke.directionalDifference(gx,v1,v2,v3,w1,w2,w3)
    dd = sub(dd,min(dd))
    dd = div(dd,max(dd))
    writeImage(ddfile,dd)
  else:
    dd = readImage(ddfile)
  plot3(dd,cmin=0.01,cmax=0.5,clab="Kaust Edge",cint=0.1,png="dd")
def goKaustEdgeEnhance():
  gx = readImage(gxfile)
  if not plotOnly:
    pd = readImage(pdfile)
    pd = KaustEdgeScanner.taper(10,0,0,pd)
    ks = KaustEdgeScanner(sigmaPhi,sigmaTheta)
    kl,kp,kt = ks.scan(minPhi,maxPhi,minTheta,maxTheta,pd)
    print "kl min =",min(kl)," max =",max(kl)
    print "kp min =",min(kp)," max =",max(kp)
    print "kt min =",min(kt)," max =",max(kt)
    writeImage(klfile+"620",kl)
    writeImage(kpfile+"620",kp)
    writeImage(ktfile+"620",kt)
  else:
    kl = readImage(klfile)
  plot3(gx,kl,cmin=0.01,cmax=0.2,cmap=jetRamp(1.0),
        clab="Fault likelihood",png="kl")
def goTimeSlices():
  gx = readImage(gxfile)
  kl = readImage(kltfile+"20")
  hd = HorizonDisplay()
  ke = KaustEdge()
  for i1 in range(640,740,20):
    sa = hd.getTimeSlice(i1,gx)
    sk = hd.getTimeSlice(i1,kl)
    sk = ke.scale(sk)
    sa = mul(sa,-1)
    sa = hd.flip1(sa)
    sk = hd.flip1(sk)
    plot2(sa,neareast=True,cint=0.5,label="Amplitude",png=str(i1)+"amplitude")
    plot2(sa,sk,cmap=jetFillExceptMin(1.0),cmin=0.2,cmax=1,
          neareast=True,cint=0.1,label="Fracture probability",png=str(i1)+"fracture")
def goInlineSlices():
  gx = readImage(gxfile)
  kl = readImage(kltfile+"20")
  hd = HorizonDisplay()
  ke = KaustEdge()
  k2s = [638,675,933]
  for k2 in k2s:
    sa = hd.getInlineSlice(k2,gx)
    sk = hd.getInlineSlice(k2,kl)
    sk = ke.scale(sk)
    sa = mul(sa,-1)
    plot2(sa,c1=s1,c2=s3,neareast=True,vlabel="Time (s)",
           vint=0.1,cint=0.5,label="Amplitude",png="ampInline"+str(k2))
    plot2(sa,sk,c1=s1,c2=s3,cmap=jetFillExceptMin(1.0),cmin=0.2,cmax=1,
          neareast=True,vlabel="Time (s)",vint=0.1, 
          cint=0.1,label="Fracture probability",
          png="fractureInline"+str(k2))

def goHorizons2D():
  #h70file = "shb5_hor70"
  h74file = "shb5_hor74"
  h76file = "shb5_hor76"
  h80file = "shb5_hor80"
  h83file = "shb5_hor83"
  h74file = "ShunB5_Hor_T74"
  h76file = "ShunB5_Hor_T76"
  h80file = "ShunB5_Hor_T80"
  h83file = "ShunB5_Hor_T83"

  #h90file = "shb5_hor90"
  hs = [h74file,h76file,h80file,h83file]
  if not plotOnly:
    gx = readImage(gxfile)
    kl = readImage(kltfile+"20")
    hd = HorizonDisplay()
    for hi in hs:
      hs = readHorizon(hi)
      sa = hd.amplitudeOnHorizon(gx,hs)
      sk = hd.amplitudeOnHorizon(kl,hs)
      hd = HorizonDisplay()
      sa = mul(sa,-1)
      ke = KaustEdge()
      sk = ke.scale(sk)
      sa = hd.flip1(sa)
      sk = hd.flip1(sk)
      writeHorizon(hi+"a",sa)
      writeHorizon(hi+"k",sk)
      plot2(sa,neareast=True,png=hi+"amplitude")
      plot2(sa,sk,cmap=jetFillExceptMin(1.0),cmin=0.2,cmax=1,
            neareast=True,png=hi+"fracture")
  else:
    for hi in hs:
      sa = readHorizon(hi+"a")
      sk = readHorizon(hi+"k")
      plot2(sa,neareast=True,cint=0.5,label="Amplitude",png=hi+"amplitude")
      plot2(sa,sk,cmap=jetFillExceptMin(1.0),cmin=0.2,cmax=1,
            neareast=True,cint=0.1,label="Fracture probability",png=hi+"fracture")
def goHorizons():
  gx = readImage(gxfile)
  kl = readImage(kltfile)
  gx = mul(gx,-1)
  h70file = "shb5_hor70"
  h74file = "shb5_hor74"
  h76file = "shb5_hor76"
  h80file = "shb5_hor80"
  h83file = "shb5_hor83"
  h90file = "shb5_hor90"
  hs = [h70file,h74file,h76file,h80file]
  hd = HorizonDisplay()
  hs = readHorizon(h76file)
  ha = readHorizon(h76file+'a')
  hk = readHorizon(h76file+'k')
  xyza,rgba=hd.buildTrigs(s3,s2,n1,-3,3,hs,ha)
  xyzk,rgbk=hd.buildTrigs(s3,s2,n1,0.0,1,hs,hk)
  tga=TriangleGroup(True,xyza,rgba)
  tgk=TriangleGroup(True,xyzk,rgbk)
  plot3(gx,tgs=tga,k1=n1-1,cmin=-3,cmax=3)
  plot3(kl,tgs=tgk,k1=n1-1,cmin=0.0,cmax=1)

  '''
  for hi in hs:
    hs = readHorizon(hi)
    sa = hd.amplitudeOnHorizon(gx,hs)
    sk = hd.amplitudeOnHorizon(kl,hs)
    #hd.fillHoles(hs)
    #mp1 = ColorMap(-3.0,3.0,ColorMap.GRAY)
    #mp2 = ColorMap(0.2,1.0,ColorMap.JET)
    #r,g,b = hd.amplitudeRgb(mp1,gx,hs) 
    #hz = [hs,r,g,b]
    #plot3(gx,hz=hz,k1=n1-1,cmin=-3,cmax=3)
    writeImage(hi+"a",sa)
    writeImage(hi+"a",sk)
  plot3(gx,kl,hz=hz,k1=n1-1,cmin=0.2,cmax=1,cmap=jetFillExceptMin(1.0),
        png=hfile+"k")
  '''

def goSta():
  gx = readImage(gxfile)
  if not plotOnly:
    sig1,sig2=8,2
    ep = zerofloat(n1,n2,n3)
    el = zerofloat(n1,n2,n3)
    lof = LocalOrientFilter(sig1,sig2)
    et = lof.applyForTensors(gx)
    sta = StructureTensorAttribute(et,40)
    sta.setEigenvalues(1.0,0.01,0.2)
    sta.applyForPlanarLinear(gx,ep,el)
    writeImage(epfile,ep)
    writeImage(elfile,el)
  else:
    ep = readImage(epfile)
  plot3(ep,cmin=0.2,cmax=1.0,clab="Planarity",cint=0.1,png="ep")

def goSemblance():
  gx = readImage(gxfile)
  if not plotOnly:
    sig1,sig2 = 4,2
    p2 = zerofloat(n1,n2,n3)
    p3 = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    lsf = LocalSlopeFinder(sig1,sig2,sig2,5)
    lsf.findSlopes(gx,p2,p3,ep)
    cov = Covariance()
    em,es = cov.covarianceEigen(12,p2,p3,gx)
    sem = div(em,es)
    writeImage(smfile,sem)
  else:
    sem = readImage(smfile)
  #sem = pow(sem,2)
  #sem = sub(sem,min(sem))
  #sem = div(sem,max(sem))
  plot3(ep,cmin=0.2,cmax=1.0,clab="Semblance",cint=0.1,png="sm")

def goSlopes():
  print "goSlopes ..."
  gx = readImage(gxfile)
  ke = KaustEdge()
  p2,p3 = ke.findSlopes(8,2,gx)
  writeImage(p2file,p2)
  writeImage(p3file,p3)
  print "p2  min =",min(p2)," max =",max(p2)
  print "p3  min =",min(p3)," max =",max(p3)
def goAmplitudeCurvature():
  gx = readImage(gxfile)
  if not plotOnly:
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    ke = KaustEdge()
    pc,nc = ke.amplitudeCurvature(gx,p2,p3)
    writeImage(pcfile,pc)
    writeImage(ncfile,nc)
  else:
    pc = readImage(pcfile)
    nc = readImage(ncfile)
  plot3(gx,pc,cmin=min(pc)*0.2,cmax=max(pc)*0.2,cmap=jetFill(1.0),
        clab="Most positive",png="pc")
  plot3(gx,nc,cmin=min(pc)*0.2,cmax=max(pc)*0.2,cmap=jetFill(1.0),
        clab="Most negative",png="nc")


def goScan():
  print "goScan ..."
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  gx = readImage(gxfile)
  gx = FaultScanner.taper(10,0,0,gx)
  fs = FaultScanner(sigmaPhi,sigmaTheta)
  fl,fp,ft = fs.scan(minPhi,maxPhi,minTheta,maxTheta,p2,p3,gx)
  print "fl min =",min(fl)," max =",max(fl)
  print "fp min =",min(fp)," max =",max(fp)
  print "ft min =",min(ft)," max =",max(ft)
  writeImage(flfile,fl)
  writeImage(fpfile,fp)
  writeImage(ftfile,ft)
  plot3(gx,clab="Amplitude")
  plot3(gx,fl,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
        clab="Fault likelihood",png="fl")
  plot3(gx,fp,cmin=0,cmax=360,cmap=hueFill(1.0),
        clab="Fault strike (degrees)",cint=45,png="fp")
  plot3(gx,convertDips(ft),cmin=25,cmax=65,cmap=jetFill(1.0),
        clab="Fault dip (degrees)",png="ft")

def goThin():
  print "goThin ..."
  gx = readImage(gxfile)
  if not plotOnly:
    kl = readImage(klfile+"420")
    kp = readImage(kpfile+"420")
    kt = readImage(ktfile+"420")
    klt,kpt,ktt = FaultScanner.thin([kl,kp,kt])
    writeImage(kltfile+"420",klt)
  else:
    klt = readImage(kltfile+"20")
    ke = KaustEdge()
    klt = ke.scale(klt)
    gx = mul(gx,-1)
    #writeImage(kltfile,klt)
  plot3(gx,k1=658,clab="Amplitude",png="s"+str(658))
  plot3(gx,klt,k1=658,cmin=0.2,cmax=1,cmap=jetFillExceptMin(1.0),
        clab="Fracture probability",png="f"+str(658))
  '''
  for k1 in range(705,905,10):
    plot3(gx,klt,k1=k1,cmin=0.2,cmax=1,cmap=jetFillExceptMin(1.0),
          png="f"+str(k1))
  '''

def goStat():
  def plotStat(s,f,slabel=None):
    sp = SimplePlot.asPoints(s,f)
    sp.setVLimits(0.0,max(f))
    sp.setVLabel("Frequency")
    if slabel:
      sp.setHLabel(slabel)
  fl = readImage(fltfile)
  fp = readImage(fptfile)
  ft = readImage(fttfile)
  fs = FaultScanner(sigmaPhi,sigmaTheta)
  sp = fs.getPhiSampling(minPhi,maxPhi)
  st = fs.getThetaSampling(minTheta,maxTheta)
  pfl = fs.getFrequencies(sp,fp,fl); pfl[-1] = pfl[0] # 360 deg = 0 deg
  tfl = fs.getFrequencies(st,ft,fl)
  plotStat(sp,pfl,"Fault strike (degrees)")
  plotStat(st,tfl,"Fault dip (degrees)")

def goSmooth():
  print "goSmooth ..."
  flstop = 0.1
  fsigma = 8.0
  gx = readImage(gxfile)
  flt = readImage(fltfile)
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  gsx = FaultScanner.smooth(flstop,fsigma,p2,p3,flt,gx)
  writeImage(gsxfile,gsx)
  plot3(gx,clab="Amplitude")
  plot3(gsx,clab="Amplitude",png="gsx")

def goSkin():
  print "goSkin ..."
  gx = readImage(gxfile)
  gsx = readImage(gsxfile)
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  fs = FaultSkinner()
  fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  fs.setMinSkinSize(minSkinSize)
  cells = fs.findCells([fl,fp,ft])
  skins = fs.findSkins(cells)
  for skin in skins:
    skin.smoothCellNormals(4)
  print "total number of cells =",len(cells)
  print "total number of skins =",len(skins)
  print "number of cells in skins =",FaultSkin.countCells(skins)
  removeAllSkinFiles(fskbase)
  writeSkins(fskbase,skins)
  plot3(gx,cells=cells,png="cells")
  plot3(gx,skins=skins,png="skins")
  for iskin,skin in enumerate(skins):
    plot3(gx,skins=[skin],links=True,png="skin"+str(iskin))

def goSlip():
  print "goSlip ..."
  gx = readImage(gxfile)
  gsx = readImage(gsxfile)
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  skins = readSkins(fskbase)
  fsl = FaultSlipper(gsx,p2,p3)
  fsl.setOffset(2.0) # the default is 2.0 samples
  fsl.setZeroSlope(False) # True only if we want to show the error
  fsl.computeDipSlips(skins,minThrow,maxThrow)
  print "  dip slips computed, now reskinning ..."
  print "  number of skins before =",len(skins),
  fsk = FaultSkinner() # as in goSkin
  fsk.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  fsk.setMinSkinSize(minSkinSize)
  fsk.setMinMaxThrow(minThrow,maxThrow)
  skins = fsk.reskin(skins)
  print ", after =",len(skins)
  removeAllSkinFiles(fskbase)
  writeSkins(fskbase,skins)
  smark = -999.999
  s1,s2,s3 = fsl.getDipSlips(skins,smark)
  writeImage(fs1file,s1)
  writeImage(fs2file,s2)
  writeImage(fs3file,s3)
  plot3(gx,skins=skins,smax=10.0,png="skinss1")
  plot3(gx,s1,cmin=-0.01,cmax=10.0,cmap=jetFillExceptMin(1.0),
        clab="Fault throw (samples)",png="gxs1")
  s1,s2,s3 = fsl.interpolateDipSlips([s1,s2,s3],smark)
  plot3(gx,s1,cmin=0.0,cmax=10.0,cmap=jetFill(0.3),
        clab="Vertical shift (samples)",png="gxs1i")
  plot3(gx,s2,cmin=-2.0,cmax=2.0,cmap=jetFill(0.3),
        clab="Inline shift (samples)",png="gxs2i")
  plot3(gx,s3,cmin=-1.0,cmax=1.0,cmap=jetFill(0.3),
        clab="Crossline shift (samples)",png="gxs3i")
  gw = fsl.unfault([s1,s2,s3],gx)
  plot3(gx)
  plot3(gw,clab="Amplitude",png="gw")

def goScanOneStrikeDip():
  print "goScanOneStrikeDip ..."
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  gx = readImage(gxfile)
  gx = FaultScanner.taper(10,0,0,gx)
  fs = FaultScanner(sigmaPhi,sigmaTheta)
  def cd(theta):
    return toDegrees(atan(tan(toRadians(theta))*5.0))
  a1,a2,a3,a4,a5,a6 = 330,350,10,30,50,190
  t1,t2,t3,t4,t5 = 20,25,30,35,40
  angles = [(a3,t1),(a3,t2),(a3,t3),(a3,t4),(a3,t5),(a1,t5),
            (a1,t3),(a2,t3),(a3,t3),(a4,t3),(a5,t3),(a6,t3)]
  for phi,theta in angles:
    suffix = "_"+str(phi)+"_"+str(theta)
    theta = cd(theta)
    fl,fp,ft = fs.scan(phi,phi,theta,theta,p2,p3,gx)
    plot3(gx,fl,cmin=0.25,cmax=1,cmap=jetRamp(1.0),
          clab="Fault likelihood",png="fl"+suffix)

def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(100.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2,n3)
  div(x,sqrt(g),y)
  return y

#############################################################################
# graphics

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

def addColorBar(frame,clab=None,cint=None):
  cbar = ColorBar(clab)
  if cint:
    cbar.setInterval(cint)
  cbar.setFont(Font("Arial",Font.PLAIN,32)) # size by experimenting
  cbar.setWidthMinimum
  cbar.setBackground(Color.WHITE)
  frame.add(cbar,BorderLayout.EAST)
  return cbar

def convertDips(ft):
  return FaultScanner.convertDips(0.2,ft) # 5:1 vertical exaggeration


def plot2(f,g=None,c1=s2,c2=s3,cmap=None,cmin=None,cmax=None,cint=None,
        label=None,hlabel="Crossline (traces)", vlabel="Inline (traces)", 
        hint=100, vint=100, neareast=True,png=None): 
  orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT;
  panel = PlotPanel(1,1,orientation)#,PlotPanel.AxesPlacement.NONE)
  panel.setVInterval(vint)
  panel.setHInterval(hint)
  panel.setHLabel(hlabel)
  panel.setVLabel(vlabel)
  pxv = panel.addPixels(0,0,c1,c2,f);
  pxv.setColorModel(ColorMap.GRAY)
  pxv.setInterpolation(PixelsView.Interpolation.LINEAR)
  if g:
    pxv.setClips(-3,3)
  else:
    if cmin and cmax:
      pxv.setClips(cmin,cmax)
    else:
      pxv.setClips(-3,3)
  if g:
    pv = panel.addPixels(c1,c2,g)
    if neareast:
      pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    else:
      pv.setInterpolation(PixelsView.Interpolation.LINEAR)
    pv.setColorModel(cmap)
    if cmin and cmax:
      pv.setClips(cmin,cmax)
  if label:
    cbar = panel.addColorBar(label)
    if cint:
      cbar.setInterval(cint)
  panel.setColorBarWidthMinimum(55)
  moc = panel.getMosaic();
  frame = PlotFrame(panel);
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setVisible(True);
  #frame.setSize(1400,700)
  frame.setSize(len(f)+55,len(f[0]))
  frame.setFontSize(16)
  if pngDir and png:
    frame.paintToPng(720,3.333,pngDir+png+".png")

def plot3(f,g=None,hz=None,tgs=None,k1=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          xyz=None,cells=None,skins=None,smax=0.0,
          links=False,curve=False,trace=False,png=None):
  n1 = len(f[0][0])
  n2 = len(f[0])
  n3 = len(f)
  s1 = Sampling(n1)
  s2 = Sampling(n2)
  s3 = Sampling(n3)
  sf = SimpleFrame(AxesOrientation.XRIGHT_YOUT_ZDOWN)
  cbar = None
  if g==None:
    ipg = sf.addImagePanels(s1,s2,s3,f)
    if cmap!=None:
      ipg.setColorModel(cmap)
    if cmin!=None and cmax!=None:
      ipg.setClips(cmin,cmax)
    else:
      ipg.setClips(-3.0,3.0)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-3.0,3.0)
    if cmin!=None and cmax!=None:
      ipg.setClips2(cmin,cmax)
    if cmap==None:
      cmap = jetFill(0.8)
    ipg.setColorModel2(cmap)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMap2Listener(cbar)
    sf.world.addChild(ipg)
  if cbar:
    cbar.setWidthMinimum(120)
  if hz:
    tg=TriangleGroup(True,s3,s2,hz[0],hz[1],hz[2],hz[3])
    sf.world.addChild(tg)
  if tgs:
    sf.world.addChild(tgs)

  if xyz:
    pg = PointGroup(0.2,xyz)
    ss = StateSet()
    cs = ColorState()
    cs.setColor(Color.YELLOW)
    ss.add(cs)
    pg.setStates(ss)
    #ss = StateSet()
    #ps = PointState()
    #ps.setSize(5.0)
    #ss.add(ps)
    #pg.setStates(ss)
    sf.world.addChild(pg)
  if cells:
    ss = StateSet()
    lms = LightModelState()
    lms.setTwoSide(True)
    ss.add(lms)
    ms = MaterialState()
    ms.setSpecular(Color.GRAY)
    ms.setShininess(100.0)
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE)
    ms.setEmissiveBack(Color(0.0,0.0,0.5))
    ss.add(ms)
    cmap = ColorMap(0.0,1.0,ColorMap.JET)
    xyz,uvw,rgb = FaultCell.getXyzUvwRgbForLikelihood(0.5,cmap,cells,False)
    qg = QuadGroup(xyz,uvw,rgb)
    qg.setStates(ss)
    sf.world.addChild(qg)
  if skins:
    sg = Group()
    ss = StateSet()
    lms = LightModelState()
    lms.setTwoSide(True)
    ss.add(lms)
    ms = MaterialState()
    ms.setSpecular(Color.GRAY)
    ms.setShininess(100.0)
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE)
    if not smax:
      ms.setEmissiveBack(Color(0.0,0.0,0.5))
    ss.add(ms)
    sg.setStates(ss)
    size = 2.0
    if links:
      size = 0.5 
    for skin in skins:
      if smax>0.0: # show fault throws
        cmap = ColorMap(0.0,smax,ColorMap.JET)
        xyz,uvw,rgb = skin.getCellXyzUvwRgbForThrow(size,cmap,False)
      else: # show fault likelihood
        cmap = ColorMap(0.0,1.0,ColorMap.JET)
        xyz,uvw,rgb = skin.getCellXyzUvwRgbForLikelihood(size,cmap,False)
      qg = QuadGroup(xyz,uvw,rgb)
      qg.setStates(None)
      sg.addChild(qg)
      if curve or trace:
        cell = skin.getCellNearestCentroid()
        if curve:
          xyz = cell.getFaultCurveXyz()
          pg = PointGroup(0.5,xyz)
          sg.addChild(pg)
        if trace:
          xyz = cell.getFaultTraceXyz()
          pg = PointGroup(0.5,xyz)
          sg.addChild(pg)
      if links:
        xyz = skin.getCellLinksXyz()
        lg = LineGroup(xyz)
        sg.addChild(lg)
    sf.world.addChild(sg)
  ipg.setSlices(k1,5,676)
  #ipg.setSlices(95,5,95)
  if cbar:
    sf.setSize(1337,900)
  else:
    sf.setSize(1200,900)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  ov.setWorldSphere(BoundingSphere(0.5*n1,0.5*n2,0.5*n3,radius))
  ov.setAzimuthAndElevation(-65.0,45.0)
  ov.setTranslate(Vector3(-0.08,0.5,0.25))
  ov.setAxesScale(1.0,1.0,0.6)
  ov.setScale(1.75)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(137,1,pngDir+png+"cbar.png")

#############################################################################
run(main)
