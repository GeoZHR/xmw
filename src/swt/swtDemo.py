"""
Demonstrate simultaneous multiple-well ties
Author: Xinming Wu, Colorado School of Mines
Version: 2015.06.08
"""

from utils import *
setupForSubset("subt")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count
d1,d2,d3 = s1.delta,s2.delta,s3.delta

# Names and descriptions of image files used below.
gxfile  = "gx" # input seismic image 
gtfile  = "gt" # RGT volume
ghfile  = "gh" # horizon volume
gufile  = "gu" # flattened image 

# Directory for saved png images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
pngDir = None
pngDir = "../../../png/swt/"
plotOnly = False

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  goDisplay()
  #goSynSeis()

def goSynSeis():
  logs = getLogs()

def goDisplay():
  gx = readImage(gxfile)
  gu = readImage(gufile)
  logs = getLogs()
  yxs,yus=[],[]
  yf = 10
  logs = logs
  gxs = zerofloat(n1,12,1)
  for il, log in enumerate(logs):
    model = SynSeis.getModel(log)
    i2 = s2.indexOfNearest(model.x2)
    i3 = s3.indexOfNearest(model.x3)
    gxs[0][il] = gx[i3][i2]
    yx = add(gx[i3][i2],yf)
    yu = add(gu[i3][i2],yf)
    yf = yf+10
    yxs.append(yx)
    yus.append(yu)
  wlw = WellLogWarping()
  wlw.setMaxShift(50)
  wlw.setPowError([2.0])
  s = wlw.findShifts([1.0],gxs)
  gus = wlw.applyShifts(gxs[0],s)
  wus = []
  uf = 10
  for i in range(12):
    for i1 in range(n1):
      if(gus[i][i1]<-10.0):
        gus[i][i1] = 0.0
    wus.append(add(gus[i],uf))
    uf = uf+10
  plot1(s1,yxs,png="originalSeisTraces")
  plot1(s1,yus,vlabel="Relative geologic time",png="flattenedSeisTraces")
  plot1(s1,wus,vlabel="Relative geologic time (dw)",png="flattenedSeisTraces")


def plot1(s1,ys,hlabel="Seismic traces",vlabel="depth (km)",png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  for y in ys:
    pv = sp.addPoints(s1,y)
    pv.setLineColor(Color.BLACK)
  #sp.setVLimits(0.1,1.1)
  sp.setSize(800,800)
  sp.setHLabel(hlabel)
  sp.setVLabel(vlabel)
  if png and pngDir:
    sp.paintToPng(300,7.0,pngDir+png+".png")



#############################################################################
run(main)
