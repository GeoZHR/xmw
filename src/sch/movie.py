"""
Plots for unfaulting and unfolding movie
Author: Xinming Wu, Colorado School of Mines
Version: 2015.03.29
"""
from uff import *
from schutils import *
setupForSubset("suf")
s1,s2,s3=getSamplings()
n1,n2,n3=s1.count,s2.count,s3.count

# Names and descriptions of files.
gxfile = "gx"  # raw input seismic image
t1file = "ft1" # unfaulting shifts (1st component)
t2file = "ft2" # unfaulting shifts (2st component)
t3file = "ft3" # unfaulting shifts (3st component)
r1file = "fr1" # unfolding shifts (1st component)
r2file = "fr2" # unfolding shifts (2st component)
r3file = "fr3" # unfolding shifts (3st component)
fwfile = "fws" # unfaulted image 

pngDir = "../../../png/sch/movie/ua/"

def main(args):
  gx = readImage(gxfile) 
  fw = readImage(fwfile) 
  t1 = readImage(t1file)
  t2 = readImage(t2file)
  t3 = readImage(t3file)
  r1 = readImage(r1file)
  r2 = readImage(r2file)
  r3 = readImage(r3file)
  uf = UnfaultS(4.0,2.0)
  fr = FlattenerRT(4.0,4.0)
  fu = zerofloat(n1,n2,n3)
  he = HorizonExtractor()
  [c1,c2,c3]=he.computeCompositeShifts([t1,t2,t3],[r1,r2,r3])
  st = [0.05,0.10,0.15,0.20,0.25,0.3,0.35,0.40,0.45,0.50,
        0.55,0.60,0.65,0.70,0.75,0.8,0.85,0.90,0.95,1.00]
  slices=(111,185,174) # most plots use these
  plot3(gx,png="u0")
  '''
  slices=(111,185,174) # most plots use these
  for it,si in enumerate(st):
    fw = zerofloat(n1,n2,n3)
    t1i = mul(t1,si)
    t2i = mul(t2,si)
    t3i = mul(t3,si)
    uf.applyShifts([t1i,t2i,t3i],gx,fw)
    if it<10:
      pngName="ua0"+str(it)
    else:
      pngName="ua"+str(it)
    plot3(fw,png=pngName)
  slices=(93,185,174) # most plots use these
  for it,si in enumerate(st):
    it = it+20
    r1i = mul(r1,si)
    r2i = mul(r2,si)
    r3i = mul(r3,si)
    fr.applyShifts([r1i,r2i,r3i],fw,fu)
    pngName="uo"+str(it)
    plot3(fu,slices=slices,png=pngName)
  slices=(93,185,174) # most plots use these
  for it,si in enumerate(st):
    it = it+20
    c1i = mul(c1,si)
    c2i = mul(c2,si)
    c3i = mul(c3,si)
    fr.applyShifts([c1i,c2i,c3i],gx,fu)
    pngName="uao"+str(it)
    plot3(fu,slices=slices,png=pngName)
  '''
def plot3(f,g=None,qg=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          slices=None,png=None):
  n3 = len(f)
  n2 = len(f[0])
  n1 = len(f[0][0])
  s1,s2,s3=Sampling(n1,1.0,0),Sampling(n2),Sampling(n3)
  #sf = SimpleFrame(AxesOrientation.XRIGHT_YOUT_ZDOWN)
  sf = SimpleFrame(AxesOrientation.XRIGHT_YIN_ZDOWN)
  cbar = None
  if g==None:
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
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-2.0,2.0)
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
  if qg:
    ss = StateSet()
    lms = LightModelState()
    #lms.setTwoSide(True)
    ss.add(lms)
    ms = MaterialState()
    ms.setSpecular(Color.GRAY)
    ms.setShininess(100.0)
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE)
    ms.setEmissiveBack(Color(0.0,0.0,0.5))
    ss.add(ms)
    for qgi in qg:
      qgi.setStates(ss)
      sf.world.addChild(qgi)
  if slices:
    k1,k2,k3 = slices
  else:
    k1,k2,k3 = (111,185,174) # most plots use these
  ipg.setSlices(k1,k2,k3)
  if cbar:
    sf.setSize(985,700) # for sch data
  else:
    sf.setSize(848,700) # for sch data
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  ov.setEyeToScreenDistance(3018.87) # for consistency with brooks
  ov.setWorldSphere(BoundingSphere(0.5*n1,0.5*n2,0.5*n3,radius))
  ov.setAzimuthAndElevation(-35.0,40.0)
  ov.setScale(1.4)
  ov.setTranslate(Vector3(-0.240,-0.258,-0.00))
  #ov.setTranslate(Vector3(-0.15,-0.15,-0.23))
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(137,1,pngDir+png+"cbar.png")

#############################################################################
run(main)
