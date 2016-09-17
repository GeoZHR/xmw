from fakeutils import *
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

gxfile  = "gx" # input image (maybe after bilateral filtering)
vxfile  = "vx" # input image (maybe after bilateral filtering)
dxfile  = "dx" # input image (maybe after bilateral filtering)
gsxfile = "gsx" # image after lsf with fault likelihoods
epfile  = "ep" # eigenvalue-derived planarity
p2file  = "p2" # inline slopes
p3file  = "p3" # crossline slopes
p2kfile = "p2k" # inline slopes (known)
p3kfile = "p3k" # crossline slopes (known)
flfile  = "fl" # fault likelihood
fpfile  = "fp" # fault strike (phi)
ftfile  = "ft" # fault dip (theta)
fltfile = "flt" # fault likelihood thinned
fptfile = "fpt" # fault strike thinned
fttfile = "ftt" # fault dip thinned
fs1file = "fs1" # fault slip (1st component)
fs2file = "fs2" # fault slip (2nd component)
fs3file = "fs3" # fault slip (3rd component)
fskbase = "fsk" # fault skin (basename only)
fslbase = "fsl" # fault skin (basename only)
fskgood = "fsg" # fault skin (basename only)
fwsfile = "fws" # unfaulted image
sw1file = "sw1" # 1st component of unfaulting shifts
sw2file = "sw2" # 2nd component of unfaulting shifts
sw3file = "sw3" # 3rd component of unfaulting shifts
gufile = "gu" # flattened image
gsfile = "gs" # flattening shifts
x1file = "x1" # horizon volume
u1file = "u1" # relateive geologic time volume
ulfile = "ul"
uncfile = "unc"
fgfile = "fg"
rgtfile = "rgt"
vqxfile = "vqx"
dqxfile = "dqx"
vqwfile = "vqw"
dqwfile = "dqw"
vqufile = "vqu"
dqufile = "dqu"

logType = "vel"
logType = "den"

pngDir = "../../../png/swt/fake/"
plotOnly = True

def main(args):
  goTensorsDisplay()
  #goGradients()
  #goInterp()

def goInterp():
  #p=getImageWithLogPoints(logType)
  #samplesW=getLogPointsW(logType)
  #samplesX=getLogPointsX(logType)
  su1 = Sampling(n1,1,1)
  gx = readImage(gxfile)
  gw = readImage(fwsfile)
  gt = readImage(rgtfile)
  fl = Flattener3Unc()
  gu = fl.flatten(s1,su1,gt,gw)
  #gu = readImage(gufile)
  sw1 = readImage(sw1file)
  sw2 = readImage(sw2file)
  sw3 = readImage(sw3file)
  skins = readSkins(fskgood)
  logs = getLogs(logType)
  cp = ConvertPoints()
  ps = cp.setUnfaultCoord(logs,skins,sw1,sw2,sw3)
  ps = cp.setFlattenedCoord(s1,s2,s3,gt,ps)
  fw,w1,w2,w3=cp.getSamplesW(ps)
  fx,x1,x2,x3=cp.getSamplesX(ps)
  fu,u1,u2,u3=cp.getSamplesU(ps)
  samplesW = fw,w1,w2,w3
  samplesX = fx,x1,x2,x3
  samplesU = fu,u1,u2,u3
  if not plotOnly:
    ri = RgtInterp3(ps)
    ri.setScales(0.001,1.0)
    ri.setRgt(gt)
    uf  = UnfaultS(4.0,2.0)
    fqu,fqw = ri.gridX(s1,s2,s3,gw)
    fqu = copy(n1,n2,n3,fqu)
    fqx = zerofloat(n1,n2,n3)
    uf.applyShiftsX([sw1,sw2,sw3],fqw,fqx)
    if logType=="vel":
      writeImage(vqxfile,fqx)
      writeImage(vqwfile,fqw)
      writeImage(vqufile,fqu)
      fqx = div(fqx,1000)
      fqw = div(fqw,1000)
      fqu = div(fqu,1000)
    if logType=="den":
      writeImage(dqxfile,fqx)
      writeImage(dqwfile,fqw)
      writeImage(dqufile,fqu)
  else:
    if logType=="vel":
      fqu = readImage(vqufile)
      fqw = readImage(vqwfile)
      fqx = readImage(vqxfile)
      fqx = div(fqx,1000)
      fqw = div(fqw,1000)
    if logType=="den":
      fqu = readImage(dqufile)
      fqw = readImage(dqwfile)
      fqx = readImage(dqxfile)
  if logType=="vel":
    fx = readImage(vxfile)
    vmin,vmax,vmap= 3.0,5.5,jetFill(0.5)
    clab = "Velocity (km/s)"
    pngx = "vqx"
    pngw = "vqw"
    pngu = "vqu"
    pngf = "vgx"
    pngxw = "gxvw"
    pngww = "gwvw"
    pnguw = "guvw"
    pngi = "vxi"
  if logType=="den":
    fx = readImage(dxfile)
    vmin,vmax,vmap= 2.1,3.2,jetFill(0.5)
    clab = "Density (g/cc)"
    pngx = "dqx"
    pngw = "dqw"
    pngu = "dqu"
    pngf = "dgx"
    pngxw = "gxdw"
    pngww = "gwdw"
    pnguw = "gudw"
    pngi = "dxi"
  plot3(gx,samples=samplesX,png=pngxw)
  plot3(gw,samples=samplesW,png=pngww)
  plot3(gu,samples=samplesU,png=pnguw)
  plot3(gu,fqu,cmin=vmin,cmax=vmax,cmap=vmap,clab=clab,samples=samplesU,png=pngu)
  plot3(gw,fqw,cmin=vmin,cmax=vmax,cmap=vmap,clab=clab,samples=samplesW,png=pngw)
  plot3(gx,fqx,cmin=vmin,cmax=vmax,cmap=vmap,clab=clab,samples=samplesX,png=pngx)
  plot3(gx,fqx,cmin=vmin,cmax=vmax,cmap=vmap,clab=clab,png=pngi)
  plot3(gx,fx, cmin=vmin,cmax=vmax,cmap=vmap,clab=clab,png=pngf)

def goTensorsDisplay():
  su1 = Sampling(n1,1,1)
  gw = readImage(fwsfile)
  gt = readImage(rgtfile)
  qu = readImage(dqufile)
  fl = Flattener3Unc()
  gu = fl.flatten(s1,su1,gt,gw)
  sw1 = readImage(sw1file)
  sw2 = readImage(sw2file)
  sw3 = readImage(sw3file)
  skins = readSkins(fskgood)
  logs = getLogs("den")
  cp = ConvertPoints()
  ps = cp.setUnfaultCoord(logs,skins,sw1,sw2,sw3)
  ps = cp.setFlattenedCoord(s1,s2,s3,gt,ps)
  fu,u1,u2,u3=cp.getSamplesU(ps)
  samplesU = fu,u1,u2,u3
  fs = zerofloat(n2,n3)
  for ip in range(len(u1)):
    c1 = round(u1[ip])
    c2 = round(u2[ip])
    c3 = round(u3[ip])
    if (c1==109):
     fs[c3][c2] = fu[ip]
     fs[c3+1][c2] = fu[ip]
     fs[c3-1][c2] = fu[ip]
     fs[c3][c2+1] = fu[ip]
     fs[c3][c2-1] = fu[ip]
     fs[c3-1][c2+1] = fu[ip]
     fs[c3-1][c2-1] = fu[ip]
     fs[c3+1][c2-1] = fu[ip]
     fs[c3+1][c2+1] = fu[ip]
  ri = RgtInterp3(ps)
  ri.setScales(0.001,1.0)
  ri.setRgt(gt)
  fa = zerofloat(n2,n3)
  v1 = zerofloat(n2,n3)
  v2 = zerofloat(n2,n3)
  el = zerofloat(n2,n3)
  et = ri.gridTest(s1,s2,s3,gw,fa)
  lof = LocalOrientFilter(4)
  gs = zerofloat(n2,n3)
  qs = zerofloat(n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      gs[i3][i2] = gu[i3][i2][109]
      qs[i3][i2] = qu[i3][i2][109]
  et = ri.makeImageTensorsX(fa)
  '''
  lof.applyForNormalLinear(fa,v1,v2,el)
  ets = lof.applyForTensors(fa)
  ets.setEigenvalues(0.001,1.0)
  lsf = LocalSmoothingFilter()
  lsf.apply(ets,16,el,el)
  for i3 in range(n3):
    for i2 in range(n2):
      if(el[i3][i2]<0.5):
        el[i3][i2] = 0
  '''
  plot2(gs,s2,s3,d=et,fmin=-2.0,fmax=1.5,dscale=2.5,png="tensorsNew")
  plot2s(s2,s3,gs,g=qs,fmin=-2.0,
          fmax=1.5,gmin=2.1,gmax=3.2,cmap=jetFill(0.5))#,png="interp2d")
  '''
  plot2s(s2,s3,gs,g=fs,fmin=-2.0,
        fmax=1.5,gmin=2.1,gmax=3.2,nearest=True,
        cmap=jetFillExceptMin(1.0),png="samples")
  plot2(gs,s2,s3,d=et,g=fs,dscale=2.5,fmin=-2,fmax=1.5,
          gmin=2.1,gmax=3.2,cmap=jetFillExceptMin(0.5))
  '''
def getLogs(logType):
  k1u = [105,105,105,105,105, 21, 21, 21, 21, 21, 21, 21]
  k2u = [127,100, 81, 56, 20,142, 70, 20,  8, 96, 50,140]
  k3u = [ 90, 78, 90, 93, 65,145,130,135,  8, 20, 30, 10]
  np = len(k1u)
  r1 = readImage(sw1file)
  r2 = readImage(sw2file)
  r3 = readImage(sw3file)
  f = zerofloat(n1,n2,n3)
  if logType == "vel":
    f = readImage(vxfile)
  if logType == "den":
    f = readImage(dxfile)
  fx = zerofloat(n1,np,4)
  for ip in range(np):
    i1 = round(k1u[ip])
    i2 = round(k2u[ip])
    i3 = round(k3u[ip])
    k2 = round(k2u[ip]+r2[i3][i2][i1])
    k3 = round(k3u[ip]+r3[i3][i2][i1])
    for k1 in range(n1):
      fx[1][ip][k1]=k1
      fx[2][ip][k1]=k2
      fx[3][ip][k1]=k3
      fx[0][ip][k1]=f[k3][k2][k1]
  return fx

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

backgroundColor = Color(0xfd,0xfe,0xff) # easy to make transparent
def jetFill(alpha):
  return ColorMap.setAlpha(ColorMap.JET,alpha)

def jetFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.JET,a)
def jetRamp(alpha):
  return ColorMap.setAlpha(ColorMap.JET,rampfloat(0.0,alpha/256,256))

def addColorBar(frame,clab=None,cint=None):
  cbar = ColorBar(clab)
  if cint:
    cbar.setInterval(cint)
  cbar.setFont(Font("Arial",Font.PLAIN,32)) # size by experimenting
  cbar.setWidthMinimum
  cbar.setBackground(Color.WHITE)
  frame.add(cbar,BorderLayout.EAST)
  return cbar

def updateColorModel(pv,alpha):
    n = 256
    r = zerobyte(n)
    g = zerobyte(n)
    b = zerobyte(n)
    a = zerobyte(n)
    icm = pv.getColorModel()
    icm.getReds(r)
    icm.getGreens(g)
    icm.getBlues(b)
    for i in range(n):
      ai = int(255.0*alpha*i/n)
      if ai>127:
        ai -= 256
      a[i] = ai
    icm = IndexColorModel(8,n,r,g,b,a)
    pv.setColorModel(icm)

def plot2(f,s1,s2,d=None,g=None,dscale=1,
        fmin=0,fmax=0,gmin=0,gmax=0,cmap=None,png=None):
  pp = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.NONE)
  pv = pp.addPixels(s1,s2,f)
  if fmin<fmax:
    pv.setClips(fmin,fmax)
  else:
    pv.setPercentiles(1,99)
  if d:
    tv = TensorsView(s1,s2,d)
    tv.setOrientation(TensorsView.Orientation.X1DOWN_X2RIGHT)
    tv.setLineColor(Color.YELLOW)
    tv.setLineWidth(4.0)
    tv.setScale(dscale)
    tv.setEllipsesDisplayed(15)
    pp.getTile(0,0).addTiledView(tv)
  if g:
    alpha = 0.0
    pv = pp.addPixels(s1,s2,g)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    if gmin<gmax:
      pv.setClips(gmin,gmax)
    else:
      pv.setPercentiles(1,99)
    pv.setColorModel(cmap)
  pf = PlotFrame(pp)
  pf.setSize(900,900)
  pf.setVisible(True)
  if png and pngDir:
    pf.paintToPng(720,3.333,pngDir+png+".png")

def plot2s(s1,s2,f,g=None,fmin=0,fmax=0,gmin=0,gmax=0,
           cmap=None,nearest=False,png=None):
  n1 = len(f[0])
  n2 = len(f)
  panel = panel2()
  pv = panel.addPixels(s1,s2,f)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setColorModel(ColorMap.GRAY)
  if fmin<fmax:
    pv.setClips(fmin,fmax)
  else:
    pv.setPercentiles(1,99)
  if g:
    alpha = 0.0
    pv = panel.addPixels(s1,s2,g)
    if nearest:
      pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    else:
      pv.setInterpolation(PixelsView.Interpolation.LINEAR)
    if gmin<gmax:
      pv.setClips(gmin,gmax)
    else:
      pv.setPercentiles(1,99)
    pv.setColorModel(cmap)
    #updateColorModel(pv,1.0)
  frame2(panel,png)


def panel2():
  #panel = PlotPanel(1,1,
  #  PlotPanel.Orientation.X1DOWN_X2RIGHT,
  #  PlotPanel.AxesPlacement.NONE)
  panel = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.NONE)
  return panel

def frame2(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setSize(900,900);
  frame.setVisible(True)
  if png:
    frame.paintToPng(720,3.333,pngDir+png+".png")
  return frame

def plot3(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          xyz=None,cells=None,skins=None,smax=0.0,slices=None,
          links=False,curve=False,trace=False,htgs=None,
          uncs=None,samples=None,png=None):
  n3 = len(f)
  n2 = len(f[0])
  n1 = len(f[0][0])
  s1,s2,s3=Sampling(n1),Sampling(n2),Sampling(n3)
  d1,d2,d3=s1.getDelta(),s2.getDelta(),s3.getDelta()
  sf = SimpleFrame(AxesOrientation.XRIGHT_YOUT_ZDOWN)
  cbar = None
  if g==None:
    ipg = sf.addImagePanels(s1,s2,s3,f)
    if cmap!=None:
      ipg.setColorModel(cmap)
    if cmin!=None and cmax!=None:
      ipg.setClips(cmin,cmax)
    else:
      #ipg.setClips(-2.0,2.0)
      ipg.setClips(-2.0,1.5) # use for subset plots
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-2.0,1.5)
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
    xyz,uvw,rgb = FaultCell.getXyzUvwRgbForLikelihood(0.7,cmap,cells,False)
    qg = QuadGroup(xyz,uvw,rgb)
    qg.setStates(ss)
    sf.world.addChild(qg)
  if htgs:
    for htg in htgs:
      sf.world.addChild(htg)
  if skins:
    sg = Group()
    ss = StateSet()
    lms = LightModelState()
    lms.setLocalViewer(True)
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
      size = 0.65 
      ls = LineState()
      ls.setWidth(1.5)
      ls.setSmooth(True)
      ss.add(ls)
    ct = 0
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
        if ct==0:
          r,g,b=0,0,0
        if ct==1:
          r,g,b=0,0,1
        if ct==2:
          r,g,b=0,1,1
        if ct==3:
          #r,g,b=0.627451,0.12549,0.941176
          r,g,b=1,1,1
        r,g,b=0,0,1
        xyz = skin.getCellLinksXyz()
        rgb = skin.getCellLinksRgb(r,g,b,xyz)
        lg = LineGroup(xyz,rgb)
        #lg = LineGroup(xyz)
        sg.addChild(lg)
        #ct = ct+1
    sf.world.addChild(sg)
  ipg.setSlices(109,138,59)
  #ipg.setSlices(92,140,59)
  if uncs:
    sg = Group()
    ss = StateSet()
    lms = LightModelState()
    lms.setLocalViewer(True)
    lms.setTwoSide(True)
    ss.add(lms)
    ms = MaterialState()
    ms.setSpecular(Color.GRAY)
    ms.setShininess(100.0)
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE)
    ss.add(ms)
    sg.setStates(ss)
    us = UncSurfer()
    ul=readImage(ulfile)
    #ul = div(exp(ul),exp(1.0))
    for unc in uncs:
      [xyz,rgb]=us.buildTrigs(n1,s3,s2,-0.1,unc,ul)
      #[xyz,rgb]=us.buildTrigs(n1,s3,s2,0.01,unc,ul)
      tg  = TriangleGroup(True,xyz,rgb)
      sg.addChild(tg)
    sf.world.addChild(sg)
  if samples:
    fx,x1,x2,x3 = samples
    if logType=="vel":
      fx = div(fx,1000)
      vmin,vmax,vmap= 3.0,5.5,ColorMap.JET
    if logType=="den":
      vmin,vmax,vmap= 2.1,3.2,ColorMap.JET
    pg = makePointGroup(fx,x1,x2,x3,vmin,vmax,None)
    sf.world.addChild(pg)
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
  ov.setWorldSphere(BoundingSphere(0.5*n1,0.4*n2,0.4*n3,radius))
  ov.setAzimuthAndElevation(120.0,25.0)
  ov.setTranslate(Vector3(0.02,0.16,-0.27))
  ov.setScale(1.25)
  # for subset plots
  #ov.setWorldSphere(BoundingSphere(0.5*n1,0.5*n2,0.5*n3,radius))
  #ov.setAzimuthAndElevation(-40.0,25.0)
  #ov.setTranslate(Vector3(0.0241,-0.0400,0.0103))
  #ov.setScale(1.3) #use only for subset plots
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")

#############################################################################
# data input/output

# John Mathewson's subsets, resampled to 20 x 20 m trace spacing
# atwj1.dat = John's area1 file
# atwj1s.dat = horizontal slice of John's atwj1 for i1=40
# atwj3.dat = John's area3 file
#n1= 129; d1=0.0040; f1=0.0000
#n2= 500; d2=0.0200; f2=0.0000
#n3= 500; d3=0.0200; f3=0.0000

def writeImage(fileName,x):
  aos = ArrayOutputStream(fileName+".dat")
  aos.writeFloats(x)
  aos.close()

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
