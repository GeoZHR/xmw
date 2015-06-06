"""
Demonstrate simultaneous multiple-well ties
Author: Xinming Wu, Colorado School of Mines
Version: 2015.06.03
"""

from utils import *
setupForSubset("subw")
sz,sl,sc = getSamplings()
nz,nl,nc = sz.count,sc.count,sl.count
dz,dl,dc = sz.delta,sc.delta,sl.delta

denWeight = 2.0
velWeight = 1.0
# power of norm for alignment errors
denPower = 0.125
velPower = 0.250
ms = 350 # maximum shift

# Names and descriptions of image files used below.
wdfile  = "wd" # well log data 
wxfile  = "wx" # array of well log curves 

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
pngDir = None
pngDir = "../../../png/swt/"

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  #goLogArray()
  goFlatten()

def goLogArray():
  logs = getLogs() #get logs with large depth ranges
  writeLogs(wdfile,logs)
  wldata = readLogData(wdfile)
  wx = writeLogDataToArray(wxfile,wldata)
  wd,wv = wx[0],wx[1]

  dcbar = "Density (g/cc)"
  vcbar = "Velocity (km/s)"
  plot2(wd,wmin=2.0,wmax=3.0,cbar=dcbar,png="den")
  plot2(wv,wmin=2.0,wmax=6.0,cbar=vcbar,png="vel")

def goFlatten():
  wx = readLogArray(wxfile)
  wd,wv = wx[0],wx[1]

  wlw = WellLogWarping()
  wlw.setMaxShift(ms)
  wlw.setPowError([denPower,velPower])
  s = wlw.findShifts([denWeight,velWeight],wx)

  gd = wlw.applyShifts(wd,s)
  gv = wlw.applyShifts(wv,s)

  wd = wlw.replaceNulls(wd,2.0)
  wv = wlw.replaceNulls(wv,2.0)
  gd = wlw.replaceNulls(gd,2.0)
  gv = wlw.replaceNulls(gv,2.0)

  dcbar = "Density (g/cc)"
  vcbar = "Velocity (km/s)"
  plot2(wd,wmin=2.0,wmax=3.0,cbar=dcbar,png="den")
  plot2(wv,wmin=2.0,wmax=6.0,cbar=vcbar,png="vel")

  vlabel = "Relative geologic time"
  plot2(gd,wmin=2.0,wmax=3.0,vlabel=vlabel,cbar=dcbar,png="denFlatten")
  plot2(gv,wmin=2.0,wmax=6.0,vlabel=vlabel,cbar=vcbar,png="velFlatten")

############################################################################
cjet = ColorMap.JET
alpha = fillfloat(1.0,256); alpha[0] = 0.0
ajet = ColorMap.setAlpha(cjet,alpha)

def plot2(w,wmin=0.0,wmax=0.0,vlabel="Depth (km)",cbar=None,png=None):
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
