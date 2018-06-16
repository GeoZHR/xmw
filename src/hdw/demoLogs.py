#############################################################################
"""
Demo of dynamic warping for well-log flattening
Author: Xinming Wu, University of Texas at Austin
Version: 2016.06.01
"""
from common import * 

#############################################################################
fxfile = "tpsz" # 
wxfile = "deepSub"
wufile = "deepSubFlatten"
curves  = ["v", "p", "d"]
logType = "v"; logLabel = "Velocity (km/s)"; vmin,vmax = 2.4,5.6
logType = "p"; logLabel = "Porosity"; vmin,vmax = 0.0,0.4
logType = "g"; logLabel = "Gamma ray (API units)"; vmin,vmax = 0.0,200.0
logType = "d"; logLabel = "Density (g/cc)"; vmin,vmax = 2.0,2.8
smin,smax = -5.5,5.5
k1,k2,k3 = 366,15,96 # 3D displays


n1,n2,n3=401,357,161
d1,d2,d3=0.004,0.025,0.025
f1,f2,f3=0.6,0.00,0.00
s1,s2,s3=Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)

logDir = "../../../data/seis/hdw/logs/"
seismicDir = "../../../data/seis/hdw/logs/"
pngDir = None
pngDir =  "../../../png/hdw/logs/seg/"

plotOnly = True
plotOnly = False

lmin,lmax=-350,350

def main(args):
  #goLogs()
  #goDeepLogSub()
  #goDeepLogs()
  #goDeepLogSubFlatten()
  #goFlatten()
  #goDeepDensityLogs()
  #goDeepPorosityLogs()
  #goPorisityFlatten()
  goDensityFlatten()
  #goVelocityFlatten()
  #goVelocityFlattenSlides()
  #goGammaFlatten()
  #goGammaFlattenError()
  #goDisplay()
  #goDeepDensityLogs()
def goDisplay():
  sz,dp,lgs = goDeepGammaLogs()
  fx = readImage(n1,n2,n3,fxfile)
  skins = readSkins("fsg")
  fl = zerofloat(n1,n2,n3)
  fsk = FaultSkinnerX()
  fsk.getFl(skins,fl)
  fls = zerofloat(n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      fls[i3][i2] = sum(fl[i3][i2])
  x2s,x3s = getLogPositions(lgs,logType)
  wh = WellHelper()
  ks = wh.sortLogIds(dp)
  lm = wh.lengthMap(s2,s3,lgs,dp)
  nw = len(dp)
  c2,c3 = getLogPosition(lgs[ks[nw-2]],logType)
  rgf1 = RecursiveGaussianFilterP(1)
  rgf2 = RecursiveGaussianFilterP(2)
  rgf1.apply00(fls,fls)
  rgf2.apply00(lm,lm)
  lm = normalize(lm)
  fls = normalize(fls)
  pv = mul(lm,sub(1,fls))
  pv = normalize(pv)
  pv = clip(0.001,1,pv)
  plot2(s2,s3,fls,label="Fault likelihood")
  plot2(s2,s3,lm,label="Normalized well-log length")
  ids = zeroint(2)
  max(pv,ids)
  c2 = ids[0]*d2
  c3 = ids[1]*d3
  plot2(s2,s3,pv,pp=[x2s,x3s],source=[[c2],[c3]],label="Pseudo velocity")
  #plot3(fx,fl,cmin=0.25,cmax=1.0,cmap=jetFillExceptMin(1.0)),
  #display(fx,logs=lgs,cmin=vmin,cmax=vmax,png=None)
  u1 = fillfloat(1,n2,n3)
  u2 = fillfloat(1,n2,n3)
  tt = fillfloat(1,n2,n3)
  mk = fillint(0,n2,n3)
  tt[ids[1]][ids[0]] = 0
  mk[ids[1]][ids[0]] = 1
  tensors = EigenTensors2(u1,u2,pv,pv)
  tm = TimeMarker2(n2,n3,tensors)
  tm.apply(tt,mk)
  print min(tt)
  print max(tt)
  tt = clip(0,400,tt)
  plot2(s2,s3,tt,label="Travel time")

def goTimeVel(dp,lgs):
  wh = WellHelper()
  nw = len(dp)
  lm = wh.lengthMap(s2,s3,lgs,dp)
  rgf1 = RecursiveGaussianFilterP(1)
  rgf2 = RecursiveGaussianFilterP(10)
  lm = normalize(lm)
  lm = pow(lm,8)
  rgf2.apply00(lm,lm)
  pv = normalize(lm)
  ids = zeroint(2)
  pmax = max(pv,ids)
  c2 = ids[0]*d2
  c3 = ids[1]*d3
  pv = clip(0.001,1,pv)
  plot2(s2,s3,pv,label="Pseudo velocity")
  '''
  va = wh.varianceMap(s2,s3,lgs,dp)
  va = normalize(va)
  va = sub(1,va)
  va = pow(va,10)
  va = normalize(va)
  pv = mul(pv,va)
  pv = normalize(pv)
  '''

  u1 = fillfloat(1,n2,n3)
  u2 = fillfloat(1,n2,n3)
  tt = fillfloat(1,n2,n3)
  mk = fillint(0,n2,n3)
  tt[ids[1]][ids[0]] = 0
  mk[ids[1]][ids[0]] = 1

  tensors = EigenTensors2(u1,u2,pv,pv)
  tm = TimeMarker2(n2,n3,tensors)
  tm.apply(tt,mk)
  print min(tt)
  print max(tt)
  plot2(s2,s3,clip(0,2000,tt),cmap=ColorMap.JET,label="Travel time")
  return tt

def goTimeGamma(dp,lgs):
  fx = readImage(n1,n2,n3,fxfile)
  skins = readSkins("fsg")
  fl = zerofloat(n1,n2,n3)
  fsk = FaultSkinnerX()
  fsk.getFl(skins,fl)
  fls = zerofloat(n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      fls[i3][i2] = max(fl[i3][i2])
  x2s,x3s = getLogPositions(lgs,logType)
  wh = WellHelper()
  nw = len(dp)
  lm = wh.lengthMap(s2,s3,lgs,dp)
  rgf1 = RecursiveGaussianFilterP(1)
  rgf2 = RecursiveGaussianFilterP(20) #gamma
  fx1 = zerofloat(n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      fx1[i3][i2] = fx[i3][i2][100]
  fls = pow(fls,2)
  fls = sub(1,fls)
  fls = normalize(fls)
  lm = normalize(lm)
  lm = pow(lm,4)
  rgf1.apply00(fls,fls)
  rgf2.apply00(lm,lm)
  lm = normalize(lm)
  fls = normalize(fls)
  pv = mul(lm,fls)
  pv = normalize(pv)
  ids = zeroint(2)
  pmax = max(pv,ids)
  c2 = ids[0]*d2
  c3 = ids[1]*d3
  pv = clip(0.00001,1,pv)
  '''
  va = wh.varianceMap(s2,s3,lgs,dp)
  va = normalize(va)
  va = sub(1,va)
  va = pow(va,10)
  va = normalize(va)
  pv = mul(pv,va)
  pv = normalize(pv)
  '''

  u1 = fillfloat(1,n2,n3)
  u2 = fillfloat(1,n2,n3)
  tt = fillfloat(1,n2,n3)
  mk = fillint(0,n2,n3)
  tt[ids[1]][ids[0]] = 0
  mk[ids[1]][ids[0]] = 1

  #pv = fillfloat(1,n2,n3)
  #pv = fillfloat(1,n2,n3)
  tensors = EigenTensors2(u1,u2,pv,pv)
  tm = TimeMarker2(n2,n3,tensors)
  tm.apply(tt,mk)
  return tt

def goTime(dp,lgs):
  fx = readImage(n1,n2,n3,fxfile)
  skins = readSkins("fsg")
  fl = zerofloat(n1,n2,n3)
  fsk = FaultSkinnerX()
  fsk.getFl(skins,fl)
  fls = zerofloat(n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      fls[i3][i2] = max(fl[i3][i2])
  x2s,x3s = getLogPositions(lgs,logType)
  print len(x2s)
  wh = WellHelper()
  nw = len(dp)
  lm = wh.lengthMap(s2,s3,lgs,dp)
  '''
  va = wh.varianceMap(s2,s3,lgs,dp)
  va = normalize(va)
  va = sub(1,va)
  va = pow(va,1)
  va = normalize(va)
  lm = mul(lm,va)
  '''
  rgf1 = RecursiveGaussianFilterP(1)
  #rgf2 = RecursiveGaussianFilterP(20) #gamma
  rgf2 = RecursiveGaussianFilterP(10) # density
  plot2(s2,s3,fls,label="Fault likelihood",png="fault")
  fx1 = zerofloat(n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      fx1[i3][i2] = fx[i3][i2][100]
  plot2(s2,s3,fx1,pp=[x2s,x3s],label="Amplitude",png="wells")
  plot2(s2,s3,fls,pp=[x2s,x3s],label="Fault likelihood",png="faultw")
  fls = pow(fls,2)
  fls = sub(1,fls)
  fls = normalize(fls)
  lm = normalize(lm)
  plot2(s2,s3,lm,nearest=True,label="Well-log length",png="length")
  lm = pow(lm,4)
  rgf1.apply00(fls,fls)
  rgf2.apply00(lm,lm)
  lm = normalize(lm)
  fls = normalize(fls)
  pv = mul(lm,fls)
  pv = normalize(pv)
  ids = zeroint(2)
  pmax = max(pv,ids)
  c2 = ids[0]*d2
  c3 = ids[1]*d3
  pv = clip(0.00001,1,pv)
  plot2(s2,s3,lm,label="Smoothed well-log length",png="lengthSmoothed")
  plot2(s2,s3,pv,source=[[c2],[c3]],cmin=0.001,cmax=1,label="Pseudo velocity",
          png="pvel")
  plot2(s2,s3,pv,pp=[x2s,x3s],source=[[c2],[c3]],label="Pseudo velocity",png="pvelw")
  '''
  va = wh.varianceMap(s2,s3,lgs,dp)
  va = normalize(va)
  va = sub(1,va)
  va = pow(va,10)
  va = normalize(va)
  pv = mul(pv,va)
  pv = normalize(pv)
  '''

  u1 = fillfloat(1,n2,n3)
  u2 = fillfloat(1,n2,n3)
  tt = fillfloat(1,n2,n3)
  mk = fillint(0,n2,n3)
  tt[ids[1]][ids[0]] = 0
  mk[ids[1]][ids[0]] = 1

  #pv = fillfloat(1,n2,n3)
  #pv = fillfloat(1,n2,n3)
  tensors = EigenTensors2(u1,u2,pv,pv)
  tm = TimeMarker2(n2,n3,tensors)
  tm.apply(tt,mk)
  print min(tt)
  print max(tt)
  plot2(s2,s3,clip(0,500,tt),source=[[c2],[c3]],
          cmap=ColorMap.JET,contour=True,label="Travel time",png="time")
  plot2(s2,s3,clip(0,500,tt),pp=[x2s,x3s],source=[[c2],[c3]],
          cmap=ColorMap.JET,contour=True,label="Travel time",png="timew")
  return tt

def goGammaFlatten():
  lmin,lmax=-100,300
  sz,dp,lgs = goDeepGammaLogs()
  wh = WellHelper()
  d = 3
  dc = wh.resample(d,dp)
  fk = 0.0003048  #1 ft = 0.0003048 km
  nz = len(dc[0])
  sz = Sampling(nz,d*fk,sz.getFirst())
  nw = len(dc)
  sw = Sampling(nw,1,1)
  ww = WellFlattener(lmin,lmax)
  ww.setStrainMax(0.01)
  ww.setErrorExponent(1)
  ww.setErrorExtrapolation(ww.ErrorExtrapolation.AVERAGE)
  clab = "Gamma ray (API units)"
  vlab = "Relative geologic time"
  wh = round(15.5*nw)
  #plot(sz,sw,dc,wh=wh,cmin=0.001,cmax=200,hint=5,cbar=clab4,png="sg")
  for w in [nw]:
    sm = zerofloat(nw)
    ww.setGate(w,0.85)
    df = ww.flatten(dc,sm)
    hp = WellHelper()
    dr = hp.sortLogs(df)
    '''
    plotCurve(sm,wh=wh-72,png="cg"+str(w))
    plot(sz,sw,df,wh=wh,cmin=0.001,cmax=200,hint=5,
        vlab=vlab,cbar=clab4,png="sgf"+str(w))
    '''
    fz = sz.first
    dz = sz.delta
    sz = Sampling(nz+100,dz,fz)
    dfe = fillfloat(-20,nz+100,nw)
    dce = fillfloat(-20,nz+100,nw)
    for iw in range(nw):
      for iz in range(nz):
        dce[iw][iz] = dc[iw][iz]
        dfe[iw][iz] = df[iw][iz]
    plot(sz,sw,dce,wh=wh,cmin=0.001,cmax=200,hint=5,cbar=clab)#,png="sgt")
    plot(sz,sw,dfe,wh=wh,cmin=0.001,cmax=200,hint=5,vlab=vlab,cbar=clab)#,png="sgf")
    cm = ww.confidence(df)
    '''
    for kw in range(nw):
      sm = []
      for k2 in range(kw+1):
        dce[k2] = dfe[k2]
        sm.append(sz.last-cm[k2]*0.3+0.1)
      plot(sz,sw,dce,wh=wh,k=kw,sm=sm,cmin=0.001,cmax=200,hint=5,
            vlab=vlab,cbar=clab,png="sgf"+str(kw))
    '''

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

def goVelocityFlattenSlides():
  lmin,lmax=-150,250
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
  ww.setErrorExponent(0.25)
  #ww.setErrorExtrapolation(ww.ErrorExtrapolation.AVERAGE)
  ww.setStrainMax(0.1)
  df = ww.flatten(dc)
  clab4 = "Velocity (km/s)"
  m = 3
  mw = nw*m
  dcp = fillfloat(-999.25,nz,mw)
  for ik in range(0,nw,1):
    dcp[ik*m] = dc[ik]
  wh = 15*mw
  sw = Sampling(mw,1.0/m,1)
  plot(sz,sw,dcp,wh=wh,cmin=2.0,cmax=6.0,hint=1,cbar=clab4,png="svSlide")
  vlab = "Relative geologic time"
  for w in [nw]:
    ww.setGate(w,0.85)
    df,us = ww.flattenX(dc)
    dfp = fillfloat(-999.25,nz,mw)
    it = int(1.348/(d*fk))
    it = int(0.8/(d*fk))
    its = [0.22,0.55,0.84,1.02,1.38,1.673]
    hz = []
    hw = []
    ht = []
    for kt in range(len(its)):
      it = int(its[kt]/(d*fk))
      hzt,hwt = ww.pickTops(it,sz,Sampling(nw,1,1),us,dc)
      htt = fillfloat(its[kt],len(hzt))
      hz.append(hzt)
      hw.append(hwt)
      ht.append(htt)
    for ik in range(0,nw,1):
      dfp[ik*m] = df[ik]
    plot(sz,sw,dfp,wh=wh,cmin=2.0,cmax=6.0,hint=1,vlab=vlab,
        cbar=clab4,png="svfSlide"+str(w))
    if w>1:
      plot(sz,sw,dfp,z=ht,w=hw,wh=wh,cmin=2.0,cmax=6.0,hint=1,vlab=vlab,
        cbar=clab4,png="svfSlidePikt"+str(w))
      plot(sz,sw,dcp,z=hz,w=hw,wh=wh,cmin=2.0,cmax=6.0,hint=1,vlab="Depth (km)",
        cbar=clab4,png="svfSlidePikz"+str(w))


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
  ww = WellFlattener(lmin,lmax)
  if not plotOnly:
    sz,dp = goDeepDensityLogs()
    wh = WellHelper()
    d = 3
    dc = wh.resample(d,dp)
    fk = 0.0003048  #1 ft = 0.0003048 km
    nz = len(dc[0])
    sz = Sampling(nz,d*fk,0)
    nw = len(dc)
    sw = Sampling(nw,1,1)
    ww.setGate(nw,0.95)
    ww.setStrainMax(0.1)
    ww.setErrorExtrapolation(ww.ErrorExtrapolation.AVERAGE)
    ww.setErrorExponent(0.125)
    #clab = "Density (g/cc)"
    #plot(sz,sw,dc,wh=wh,cmin=2.0,cmax=2.8,hint=5,cbar=clab,png="sdr")
    #df = ww.flatten(dc)
    sm = zerofloat(nw)
    df = ww.flatten(dc,sm)
    #print len(df)
    #print len(df[0])
    writeImage("denx",dc)
    writeImage("fdenx",df)
    ndfz = [sz.count,sz.delta,sz.first]
    ndfw = [sw.count,sw.delta,sw.first]
    writeImage("denndfzx",ndfz)
    writeImage("denndfwx",ndfw)
  else:
    ndfz = readImage1D(3,"denndfzx")
    ndfw = readImage1D(3,"denndfwx")
    sz = Sampling(int(ndfz[0]),ndfz[1],ndfz[2])
    sw = Sampling(int(ndfw[0]),ndfw[1],ndfw[2])
    nz = sz.count
    nw = sw.count
    dc = readImage2D(nz,nw,"denx")
    df = readImage2D(nz,nw,"fdenx")
  sz = Sampling(int(ndfz[0]),ndfz[1],ndfz[2])
  sw = Sampling(int(ndfw[0]),ndfw[1],ndfw[2])
  nz = sz.count
  nw = sw.count
  wh = 12*nw
  vlab = "Relative geologic time"
  clab = "Density (g/cc)"
  fz = sz.first
  dz = sz.delta
  sz = Sampling(nz+100,dz,fz)
  dfe = fillfloat(-20,nz+100,nw)
  dce = fillfloat(-20,nz+100,nw)
  for iw in range(nw):
    for iz in range(nz):
      dce[iw][iz] = dc[iw][iz]
      dfe[iw][iz] = df[iw][iz]
  plot(sz,sw,dce,wh=wh,cmin=2.0,cmax=2.8,hint=5,cbar=clab,png="sdtx")
  plot(sz,sw,dfe,wh=wh,cmin=2.0,cmax=2.8,hint=5,vlab=vlab,cbar=clab,png="sdfx")
  cm = ww.confidence(df)
  for kw in range(nw-1,nw,1):
    sm = []
    for k2 in range(kw+1):
      dce[k2] = dfe[k2]
      sm.append(sz.last-cm[k2]*0.3+0.1)
    plot(sz,sw,dce,wh=wh,k=kw,sm=sm,cmin=2.0,cmax=2.8,hint=5,vlab=vlab,cbar=clab,
          png="sdf"+str(kw))

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
  ww.setErrorExponent(0.5)
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
  lgs = []
  for iw in range(len(wd)):
    wg = dd[2][iw]
    if(max(wg)>0):
      dp.append(wg)
      lgs.append(wd[iw])
  sw = Sampling(len(dp))
  tt = goTimeGamma(dp,lgs)
  #dr = wh.sortLogs(dp)
  dr = wh.sortByTravelTime(s2,s3,tt,lgs,dp)
  dt = wh.trim(ndf,dr)
  sz = Sampling(int(ndf[0]),ndf[1],ndf[2])
  #plot(sz,sw,dp,cmin=0.001,cmax=0.45,cbar="Porosity (%)")
  nz = len(dt[0])
  dc = copy(nz,113,0,0,dt)
  #plot(sz,sw,dr,cmin=0.001,cmax=0.45,cbar="Porosity (%)")
  return sz,dc,lgs

def goDeepDensityLogs():
  wd = getDeepLogs()
  wh = WellHelper()
  ndf = zerodouble(3)
  dd = wh.toArray(wd,ndf)
  dp = []
  lgs = []
  for iw in range(len(wd)):
    wg = dd[1][iw]
    if(max(wg)>0):
      dp.append(wg)
      lgs.append(wd[iw])
  sw = Sampling(len(dp))
  #dr = wh.sortLogs(dp)
  #dr0 = dr[0]
  #dr[0] = dr[1]
  #dr[1] = dr0
  tt = goTime(dp,lgs)
  ds = wh.sortLogs(dp)
  dr = []
  for k in range(len(ds)-1,-1,-1):
    dr.append(ds[k])

  #dr = wh.sortByTravelTime(s2,s3,tt,lgs,dp)
  #dr = dp
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
  lgs = []
  for iw in range(len(wd)):
    wg = dd[0][iw]
    if(max(wg)>0):
      dp.append(wg)
      lgs.append(wd[iw])
  '''
  dp = []
  for wp in dd[0]:
    if(max(wp)>0):
      dp.append(wp)
  '''
  sw = Sampling(len(dp))
  #dr = wh.sortLogs(dp)
  tt = goTimeVel(dp,lgs)
  #dr = wh.sortLogs(dp)
  dr = wh.sortByTravelTime(s2,s3,tt,lgs,dp)
  dr2 = dr[2]
  dr[2] = dr[0]
  dr[0] = dr2
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
def readImage1D(n1,basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = logDir+basename+".dat"
  image = zerofloat(n1)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def readImage2D(n1,n2,basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = logDir+basename+".dat"
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

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
def skinName(basename,index):
  return basename+("%03i"%(index))

def skinIndex(basename,fileName):
  assert fileName.startswith(basename)
  i = len(basename)
  return int(fileName[i:i+3])
def readSkin(basename,index):
  """ Reads one skin with specified basename and index. """
  return FaultSkin.readFromFile(seismicDir+skinName(basename,index)+".dat")

def readSkins(basename):
  """ Reads all skins with specified basename. """
  fileNames = []
  for fileName in File(seismicDir).list():
    if fileName.startswith(basename):
      fileNames.append(fileName)
  fileNames.sort()
  skins = []
  for iskin,fileName in enumerate(fileNames):
    index = skinIndex(basename,fileName)
    skin = readSkin(basename,index)
    skins.append(skin)
  return skins

#############################################################################
# plotting
backgroundColor = Color.WHITE
cjet = ColorMap.JET
alpha = fillfloat(1.0,256); alpha[0] = 0.0
ajet = ColorMap.setAlpha(cjet,alpha)
def plot(sz,sw,fx,k=None,sm=None,z=None,w=None,wh=500,wv=800,cmin=None,cmax=None,
    hint=2,vlab="Depth (km)",cbar=None,png=None):
  fz = sz.first
  dz = sz.delta
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(wh,wv)
  sp.setHInterval(hint)
  sp.setVLabel(vlab)
  sp.setHLabel("Log index")
  sp.addColorBar(cbar)
  sp.plotPanel.setColorBarWidthMinimum(72)
  sp.setHLimits(sw.first-sw.delta/2,sw.last)
  sp.setVLimits(sz.first,sz.last)
  pv = sp.addPixels(sz,sw,fx)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ajet)
  n1 = len(fx[0])
  if k:
    #pv1 = sp.addPoints([sz.last-250*dz,sz.last-250*dz],[1,k])
    #pv1.setLineColor(Color.RED)
    #pv1.setLineWidth(3.0)
    e1 = n1-1
    nullValue=-999.25
    for k1 in range(n1-1,0,-1):
      if(fx[k][k1]>0):
        e1 = k1
        break
    pv2 = sp.addPoints([sz.last,(e1-10)*dz+fz],[k+1,k+1])
    pv2.setLineColor(Color.BLUE)
    pv2.setLineStyle(PointsView.Line.DASH)
    pv2.setLineWidth(2.0)
  if sm:
    wm = []
    for i in range(len(sm)):
      wm.append(i+1)
    pv1 = sp.addPoints(sm,wm)
    pv1.setLineColor(Color.RED)
    pv1.setLineWidth(3.0)
  if z and w:
    for k in range(len(z)):
      hpv = sp.addPoints(z[k],w[k])
      hpv.setLineColor(Color.RED)
      hpv.setLineWidth(3)
  sp.setFontSize(20)
  if cmin and cmax:
    pv.setClips(cmin,cmax)
  if png and pngDir:
    sp.paintToPng(1080,3.33333,pngDir+png+".png")

def plotCurve(sm,wh=500,wv=220,hint=5,png=None):
  n = len(sm)
  panel = PlotPanel(1,1,PlotPanel.Orientation.X1RIGHT_X2UP)
  panel.mosaic.setHeightElastic(0,25)
  panel.setHLimits(0,1,n)
  s1 = Sampling(n,1,1)
  sv = panel.addPoints(0,0,s1,sm)
  sv.setLineWidth(2)
  panel.setHInterval(hint)
  panel.setVLabel("Confidence")
  panel.setHLabel("Log index")
  panel.setVLimits(0,0.6,1.0)
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(backgroundColor)
  frame.setFontSize(20)
  frame.setSize(wh,wv)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(720,3.33333,pngDir+"/"+png+".png")

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

def addColorBarX(frame,clab=None,cint=None):
  cbar = ColorBar(clab)
  if cint:
    cbar.setInterval(cint)
  cbar.setFont(Font("Arial",Font.PLAIN,32)) # size by experimenting
  cbar.setWidthMinimum
  cbar.setBackground(Color.WHITE)
  frame.add(cbar,BorderLayout.EAST)
  return cbar

def plot3(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          xyz=None,cells=None,skins=None,smax=0.0,
          links=False,curve=False,trace=False,png=None):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s1,s2,s3 = Sampling(n1),Sampling(n2),Sampling(n3)
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  world = World()
  sf = SimpleFrame(world)
  cbar = None
  if g==None:
    ipg = sf.addImagePanels(s1,s2,s3,f)
    if cmap!=None:
      ipg.setColorModel(cmap)
    if cmin!=None and cmax!=None:
      ipg.setClips(cmin,cmax)
    else:
      ipg.setClips(smin,smax)
    if clab:
      cbar = addColorBarX(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setColorModel1(ColorMap.getGray())
    ipg.setClips1(smin,smax)
    if cmin!=None and cmax!=None:
      ipg.setClips2(cmin,cmax)
    if cmap==None:
      cmap = jetFill(0.8)
    ipg.setColorModel2(cmap)
    if clab:
      cbar = addColorBarX(sf,clab,cint)
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
    xyz,uvw,rgb = FaultCell.getXyzUvwRgbForLikelihood(0.5,cmap,cells,False)
    qg = QuadGroup(xyz,uvw,rgb)
    qg.setStates(ss)
    sf.world.addChild(qg)
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
      size = 0.6 
    ct = 0
    for skin in skins:
      if smax>0.0: # show fault throws
        cmap = ColorMap(-10,smax,ColorMap.JET)
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
          r,g,b=0,0,1
        if ct==1:
          r,g,b=0,1,0
        if ct==2:
          r,g,b=0,1,1
        if ct==3:
          r,g,b=0.627451,0.12549,0.941176
        xyz = skin.getCellLinksXyz()
        rgb = skin.getCellLinksRgb(r,g,b,xyz)
        #lg = LineGroup(xyz,rgb)
        lg = LineGroup(xyz)
        sg.addChild(lg)
        ct = ct+1
    sf.world.addChild(sg)
  ipg.setSlices(498,296,110)
  ipg.setSlices(k1,k2,k3)
  #ipg.setSlices(n1-1,0,n3-1)
  if cbar:
    sf.setSize(837,700)
    sf.setSize(1250,900)
  else:
    sf.setSize(700,700)
    sf.setSize(1250,900)

  view = sf.getOrbitView()
  zscale = 0.75*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.3)
  #view.setAzimuth(75.0)
  #view.setAzimuth(-75.0)
  view.setAzimuth(-65.0)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  sf.viewCanvas.setBackground(sf.getBackground())
  sf.setSize(1250,900)

  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")

def plot2(s1,s2,f,g=None,pp=None,source=None,nearest=False,
        cmin=None,cmax=None,cmap=None,contour=False,label=None,png=None):
  n2 = len(f)
  n1 = len(f[0])
  f1,f2 = s1.getFirst(),s2.getFirst()
  l1,l2 = s1.getLast(),s2.getLast()
  d1,d2 = s1.getDelta(),s2.getDelta()
  panel = panel2Teapot()
  panel.setHLimits(0,f1,l1)
  panel.setVLimits(0,f2,l2)
  panel.setHInterval(1.0)
  panel.setVInterval(1.0)
  panel.setHLabel("Inline (km)")
  panel.setVLabel("Crossline (km)")
  if label:
    panel.addColorBar(label)
  else:
    panel.addColorBar()
  panel.setColorBarWidthMinimum(45)
  if nearest:
    ft = copy(f)
    for i2 in range(1,n2-1,1):
      for i1 in range(1,n1-1,1):
        fti = ft[i2][i1]
        if(fti):
          f[i2  ][i1-1] = fti
          f[i2  ][i1+1] = fti
          f[i2-1][i1-1] = fti
          f[i2-1][i1+1] = fti
          f[i2+1][i1-1] = fti
          f[i2+1][i1+1] = fti
          f[i2-1][i1  ] = fti
          f[i2+1][i1  ] = fti
  pv = panel.addPixels(s1,s2,f)
  if nearest:
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  else:
    pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setColorModel(ColorMap.GRAY)
  if cmap:
    pv.setColorModel(cmap)
  if contour:
    cv = panel.addContours(s1,s2,f)
    cv.setContours(100)
    #cv.setColorModel(ColorMap.JET)
    cv.setLineWidth(2.0)
    cv.setLineColor(Color.BLACK)
  if pp:
    ptvp = panel.addPoints(0,0,pp[0],pp[1])
    ptvp.setLineStyle(PointsView.Line.NONE)
    ptvp.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE)
    #ptvp.setMarkStyle(PointsView.Mark.CROSS)
    ptvp.setMarkColor(Color.YELLOW)
    ptvp.setMarkSize(6.0)
    ptvp.setLineWidth(3.0)
  if source:
    ptvp = panel.addPoints(0,0,source[0],source[1])
    ptvp.setLineStyle(PointsView.Line.NONE)
    ptvp.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE)
    ptvp.setMarkColor(Color.RED)
    ptvp.setMarkSize(6.0)
    ptvp.setLineWidth(3.0)
  if g:
    pv = panel.addPixels(s1,s2,g)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setColorModel(cmap)
    if label:
      panel.addColorBar(label)
    else:
      panel.addColorBar()
  if cmin and cmax:
    pv.setClips(cmin,cmax)
  else:
    pv.setClips(min(f),max(f))
  frame2Teapot(panel,png)
def panel2Teapot():
  panel = PlotPanel(1,1,
    PlotPanel.Orientation.X1RIGHT_X2UP)#,PlotPanel.AxesPlacement.NONE)
  return panel
def frame2Teapot(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  #frame.setFontSizeForPrint(8,240)
  #frame.setSize(1240,774)
  #frame.setFontSizeForSlide(1.0,0.9)
  frame.setFontSize(12)
  frame.setSize(round(n2*1.5)+45,round(n3*1.5))
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(720,3.2,pngDir+png+".png")
  return frame

def display(s,logs=None,cmin=0,cmax=0,png=None):
  world = World()
  ipg = addImageToWorld(world,s)
  ipg.setClips(smin,smax)
  ipg.setSlices(k1,k2,k3)
  frame = makeFrame(world)
  if logs:
    cbar = addColorBar(frame,logLabel)
    addLogsToWorld(world,logs,logType,cmin,cmax,cbar,smooth=smooth)
  if png and pngDir:
    frame.paintToFile(pngDir+png+".png")
def addImageToWorld(world,image):
  ipg = ImagePanelGroup(s1,s2,s3,image)
  world.addChild(ipg)
  return ipg
def makeFrame(world):
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  frame = SimpleFrame(world)
  view = frame.getOrbitView()
  zscale = 0.75*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.3)
  #view.setAzimuth(75.0)
  #view.setAzimuth(-75.0)
  view.setAzimuth(-65.0)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  frame.viewCanvas.setBackground(frame.getBackground())
  frame.setSize(1250,900)
  frame.setVisible(True)
  return frame

def addLogsToWorld(world,logs,type,cmin=0,cmax=0,cbar=None,smooth=0):
  samples = getLogSamples(logs,type,smooth)
  #print "number of logs =",len(samples[0])
  lg = makeLogPoints(samples,type,cmin,cmax,cbar)
  #lg = makeLogLines(samples,type,cmin,cmax)
  states = StateSet()
  cs = ColorState()
  cs.setColor(Color.YELLOW)
  states.add(cs)
  lg.setStates(states)
  world.addChild(lg)

def getLogPositions(logs,type):
  """ 
  Reads log curves from the specified set that have the specified type.
  set: "s" for shallow, "d" for deep, or "a" for all
  type: "v" (velocity), "d" (density), "p" (porosity), or "g" (gamma)
  smooth: half-width of Gaussian smoothing filter
  Returns a tuple (f,x1,x2,x3) of lists of arrays of samples f(x1,x2,x3)
  """
  x2l,x3l = [],[]
  for log in logs:
    samples = log.getSamples(type,s1,s2,s3)
    if samples:
      f,x1,x2,x3 = samples
      x2l.append(x2[0])
      x3l.append(x3[0])
  return x2l,x3l

def getLogPosition(log,type):
  """ 
  Reads log curves from the specified set that have the specified type.
  set: "s" for shallow, "d" for deep, or "a" for all
  type: "v" (velocity), "d" (density), "p" (porosity), or "g" (gamma)
  smooth: half-width of Gaussian smoothing filter
  Returns a tuple (f,x1,x2,x3) of lists of arrays of samples f(x1,x2,x3)
  """
  samples = log.getSamples(type,s1,s2,s3)
  if samples:
    f,x1,x2,x3 = samples
  return [x2[0],x3[0]]

def getLogSamples(logs,type,smooth=0):
  """ 
  Reads log curves from the specified set that have the specified type.
  set: "s" for shallow, "d" for deep, or "a" for all
  type: "v" (velocity), "d" (density), "p" (porosity), or "g" (gamma)
  smooth: half-width of Gaussian smoothing filter
  Returns a tuple (f,x1,x2,x3) of lists of arrays of samples f(x1,x2,x3)
  """
  fl,x1l,x2l,x3l = [],[],[],[]
  for log in logs:
    samples = log.getSamples(type,s1,s2,s3)
    if samples:
      f,x1,x2,x3 = samples
      fl.append(f)
      x1l.append(x1)
      x2l.append(x2)
      x3l.append(x3)
  return fl,x1l,x2l,x3l

def makeLogPoints(samples,type,cmin,cmax,cbar):
  lg = Group()
  fl,x1l,x2l,x3l = samples
  for i,f in enumerate(fl):
    f = fl[i]
    x1 = x1l[i]
    x2 = x2l[i]
    x3 = x3l[i]
    pg = makePointGroup(f,x1,x2,x3,cmin,cmax,cbar)
    lg.addChild(pg)
  return lg

def makePointGroup(f,x1,x2,x3,cmin,cmax,cbar):
  n = len(x1)
  xyz = zerofloat(3*n)
  copy(n,0,1,x3,0,3,xyz)
  copy(n,0,1,x2,1,3,xyz)
  copy(n,0,1,x1,2,3,xyz)
  rgb = None
  if cmin<cmax:
    cmap = ColorMap(cmin,cmax,ColorMap.JET)
    if cbar:
      cmap.addListener(cbar)
    rgb = cmap.getRgbFloats(f)
  pg = PointGroup(xyz,rgb)
  ps = PointState()
  ps.setSize(4)
  ps.setSmooth(False)
  ss = StateSet()
  ss.add(ps)
  pg.setStates(ss)
  return pg

def addColorBar(frame,label):
  cbar = ColorBar(logLabel)
  cbar.setFont(cbar.getFont().deriveFont(36.0))
  frame.add(cbar,BorderLayout.EAST)
  cbar.setWidthMinimum
  #frame.viewCanvas.setBackground(frame.getBackground())
  return cbar

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
