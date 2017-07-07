"""
Demonstrate simultaneous multiple-well ties
Author: Xinming Wu, Colorado School of Mines
Version: 2015.06.08
"""

import sys
import jarray

from java.awt import *
from java.io import *
from java.lang import *
from java.nio import *
from java.util import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.interp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.ogl.Gl import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from csu import *


'''
n1,n2,n3 = s1.count,s2.count,s3.count
d1,d2,d3 = s1.delta,s2.delta,s3.delta
f1,f2,f3 = s1.first,s2.first,s3.first
'''

# Names and descriptions of image files used below.
Ixfile = "BY_2_L100P2200_79A025.dat.txt" # input seismic image 
U1file = "BY_2_L100P2200_79A050.dat.txt" # input seismic image 
U2file = "BY_2_L100P2200_79A191.dat.txt" # input seismic image 
U3file = "BY_2_L100P2200_79A195.dat.txt" # input seismic image 

pngDir = "../../../png/csu/ip/"
seisDir = "../../../data/seis/csu/ip/"
pngDir = None
plotOnly = True

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  ixs = readData(Ixfile)
  ux1 = readData(U1file)
  ux2 = readData(U2file)
  ux3 = readData(U3file)
  print len(ixs[1])
  np = 1024
  #ixs = copy(np,4,ixs)
  #ux1 = copy(np,4,ux1)
  #ux2 = copy(np,4,ux2)
  #ux3 = copy(np,4,ux3)
  Ipr = IpResistivity()
  ixt = Ipr.stack(100,1024,ixs[1])
  uxt = Ipr.stack(100,1024,ux1[1])
  Icx = Ipr.applyForComplex(ixt)
  Ucx = Ipr.applyForComplex(uxt)
  rcx = Ipr.applyForResistivityAndPhase(Icx,Ucx)
  n1 = len(rcx[0])
  d1 = 0.5*64/n1
  s1 = Sampling(n1,d1,0)
  plot1(s1,[Icx[2]])
  plot1(s1,[Ucx[2]])
  plot1(s1,[rcx[0]])
  plot1(s1,[rcx[1]])
  '''
  plot1s(ixt,ux1)
  plot1s(ixs[1],ux2)
  plot1s(ixs[1],ux3)
  '''

def readData(basename):
  tr = TxtReader()
  fileName = seisDir+basename
  ips = tr.readIps(fileName)
  return ips

#############################################################################
# plotting
backgroundColor = Color.WHITE

def plot1(s1,ys,hlabel="Frequency",vlabel="Amplitude",png=None):
  sp = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
  for y in ys:
    pv = sp.addPoints(s1,y)
    pv.setLineColor(Color.BLACK)
  #sp.setVLimits(0.1,1.1)
  sp.setSize(800,800)
  sp.setHLabel(hlabel)
  sp.setVLabel(vlabel)
  if png and pngDir:
    sp.paintToPng(300,7.0,pngDir+png+".png")

def plot1s(I,Us,png=None):
  n = len(I)
  panel = PlotPanel(5,1,PlotPanel.Orientation.X1RIGHT_X2UP)
  panel.mosaic.setHeightElastic(0,25)
  panel.mosaic.setHeightElastic(1,25)
  panel.setHLimits(0,0,n-1)
  Iv = panel.addPoints(0,0,I)
  U1 = panel.addPoints(1,0,Us[0])
  U2 = panel.addPoints(2,0,Us[1])
  U3 = panel.addPoints(3,0,Us[2])
  U4 = panel.addPoints(4,0,Us[3])
  Iv.setLineWidth(2)
  U1.setLineWidth(2)
  U2.setLineWidth(2)
  U3.setLineWidth(2)
  U4.setLineWidth(2)
  '''
  panel.setVLimits(0,-550,550)
  panel.setVLimits(1,-80,80)
  panel.setVLimits(2,-80,80)
  panel.setVLimits(3,-80,80)
  panel.setVLimits(4,-80,80)
  '''
  #panel.setHLabel("sample index i")
  #panel.setVLabel(0,"f")
  #panel.setVLabel(1,"g")
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(backgroundColor)
  #frame.setFontSizeForPrint(8,240)
  frame.setSize(1600,850)
  frame.setVisible(True)
  if png and pngDir:
    png += "n"+str(int(10*nrms))
    frame.paintToPng(720,3.33333,pngDir+"/"+png+".png")

def plot1ss(s1,rs,ss,ys,sm,ym,vmin=None,vmax=None,
    hlabel="Log index",vlabel="Time (s)",png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sf = 1.0
  yf = sf
  sp.setVLimits(0.1,1.0)
  if vmin and vmax:
    sp.setVLimits(vmin,vmax)
  sp.setHLimits(0.5,11.5)
  sp.setHInterval(2)
  for il,y in enumerate(ys):
    ya = sum(y)/len(y)
    y = sub(y,ya)
    y = div(y,10)
    y = add(y,yf)
    pv = sp.addPoints(ss[il],y)
    pv.setLineColor(Color.GREEN)
    pv.setLineWidth(3.0)
    yf = yf+sf
  yf = sf
  for il,y in enumerate(ym):
    ya = sum(y)/len(y)
    y = sub(y,ya)
    y = div(y,10)
    y = add(y,yf)
    pv = sp.addPoints(sm[il],y)
    pv.setLineColor(Color.RED)
    pv.setLineWidth(3.0)
    yf = yf+sf
  rf = sf
  for il,r in enumerate(rs):
    ra = sum(r)/len(r)
    r = sub(r,ra)
    r = div(r,10)
    r = add(r,rf)
    pv = sp.addPoints(s1,r)
    pv.setLineColor(Color.BLACK)
    pv.setLineWidth(3.0)
    rf = rf+sf
  sp.setSize(600,650)
  sp.setHLabel(hlabel)
  sp.setVLabel(vlabel)
  sp.setFontSize(20) #for print
  #sp.setFontSize(30) #for slides
  sp.setVInterval(0.2)
  if png and pngDir:
    sp.paintToPng(300,7.0,pngDir+png+".png")

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
