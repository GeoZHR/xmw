"""
Demonstrate simultaneous multiple-well ties
Author: Xinming Wu, Colorado School of Mines
Version: 2015.06.03
"""

from utils import *
setupForSubset("subw")
global sz,sl,sc
global nz,nl,nc
global dz,dl,dc

denWeight = 2.0
velWeight = 1.0
# power of norm for alignment errors
denPower = 0.125
velPower = 0.250
ms = 350 # maximum shift

# Names and descriptions of image files used below.
wdfile  = "wd" # well log data 
wxfile  = "wx" # array of well log curves 
wshifts = "wellShifts"
wshiftd = "wellShiftsDimensions"
wsample = "wellSampling"

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
pngDir = None
pngDir = "../../../png/swt/"

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  goFlatten()

def goFlatten():
  print "goFlatten..."
  logs = getLogs()
  zs = zerofloat(3)
  sw = SeismicWellTie()
  wx = sw.logsToArray(logs,zs)
  sz = Sampling((int)(zs[0]),zs[1],zs[2])
  sl = Sampling(len(wx[0]),1,0)
  wd = wx[0] #array of density logs
  wv = wx[1] #array of velocity logs

  wlw = WellLogWarping()
  wlw.setMaxShift(ms)
  wlw.setPowError([denPower,velPower])
  s = wlw.findShifts([denWeight,velWeight],wx)

  sdm = zerofloat(2)
  sdm[0] = len(s)
  sdm[1] = len(s[0])
  writeImage(wshifts,s)
  writeImage(wsample,zs)
  writeImage(wshiftd,sdm)

  gd = wlw.applyShifts(wd,s)
  gv = wlw.applyShifts(wv,s)

  dcbar = "Density (g/cc)"
  vcbar = "Velocity (km/s)"
  plot2(wd,sz,sl,wmin=2.0,wmax=3.0,cbar=dcbar,png="den")
  plot2(wv,sz,sl,wmin=2.0,wmax=6.0,cbar=vcbar,png="vel")

  vlabel = "Relative geologic time"
  plot2(gd,sz,sl,wmin=2.0,wmax=3.0,vlabel=vlabel,cbar=dcbar,png="denFlatten")
  plot2(gv,sz,sl,wmin=2.0,wmax=6.0,vlabel=vlabel,cbar=vcbar,png="velFlatten")

############################################################################
cjet = ColorMap.JET
alpha = fillfloat(1.0,256); alpha[0] = 0.0
ajet = ColorMap.setAlpha(cjet,alpha)

def plot2(w,sz,sl,wmin=0.0,wmax=0.0,vlabel="Depth (km)",cbar=None,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(500,900)
  sp.setVLabel(vlabel)
  sp.setHLabel("Log index")
  sp.addColorBar(cbar)
  sp.plotPanel.setColorBarWidthMinimum(90)
  pv = sp.addPixels(sz,sl,w)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ajet)
  pv.setClips(wmin,wmax)
  if png and pngDir:
    sp.paintToPng(300,7.0,pngDir+png+".png")


#############################################################################
run(main)
