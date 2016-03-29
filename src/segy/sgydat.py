"""
Converts sgy image files (SEG-Y format) to dat files (3D arrays of floats)
Removes all SEG-Y headers and, if necessary, converts the data format to
IEEE floats. The byte order for floats in dat files is BIG_ENDIAN. 

Author: Dave Hale, Colorado School of Mines
Version: 2012.06.18
"""
from imports import *
global n1,n2,n3

#############################################################################
def main(args):
  #goHongliu()
  #goF3dUnc()
  #goJake()
  #goNathan()
  #goLulia()
  #goCranfield2007()
  #goCranfield2010()
  goSeam()
def goSeam():
  '''
  ****** beginning of SEG-Y file info ******
  file name = ../../../data/seis/seam/SEAM_Interpretation_Challenge_1_Time.sgy
  byte order = BIG_ENDIAN
  number of bytes = 4268359272
  number of traces = 1171338
  format = 1 (4-byte IBM floating point)
  units for spatial coordinates: m (will be converted to km)
  indices and coordinates from trace headers:
    i2min =  1499, i2max =  8507 (inline indices)
    i3min =  1499, i3max =  7505 (crossline indices)
    xmin =    2.490000, xmax =   32.520000 (x coordinates, in km)
    ymin =    2.490000, ymax =   37.530000 (y coordinates, in km)
  grid sampling:
    n1 =   851 (number of samples per trace)
    n2 =  7009 (number of traces in inline direction)
    n3 =  6007 (number of traces in crossline direction)
    d1 = 0.012000 (time sampling interval, in s)
    d2 = 0.005000 (inline sampling interval, in km)
    d3 = 0.005000 (crossline sampling interval, in km)
  grid corner points:
    i2min =  1499, i3min =  1499, x =    2.490000, y =    2.490000
    i2max =  8507, i3min =  1499, x =    2.490000, y =   37.530000
    i2min =  1499, i3max =  7505, x =   32.520000, y =    2.490000
    i2max =  8507, i3max =  7505, x =   32.520000, y =   37.530000
  grid azimuth:  0.00 degrees
  ****** end of SEG-Y file info ******
  '''
  firstLook = False # fast, does not read all trace headers
  secondLook = False # slow, must read all trace headers
  writeImage = False # reads all traces, writes an image
  showImage = True # displays the image
  basedir = "../../../data/seis/seam/"
  sgyfile = basedir+"SEAM_Interpretation_Challenge_1_Time.sgy"
  datfile = basedir+"gx.dat"
  i1min,i1max,i2min,i2max,i3min,i3max = 0,850,1499,8507,1499,7505
  i1min,i1max,i2min,i2max,i3min,i3max = 0,850,0,1168,0,1001
  n1,n2,n3 = 1+i1max-i1min,1+i2max-i2min,1+i3max-i3min
  si = SegyImage(sgyfile)
  if firstLook:
    si.printSummaryInfo();
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
  if secondLook:
    si.printAllInfo()
    plot23(si)
    plotXY(si)
  if writeImage:
    scale = 1.00
    si.writeFloats(datfile,scale,i1min,i1max,i2min,i2max,i3min,i3max,6,6)
  si.close()
  if showImage:
    x = readImage(datfile,n1,n2,n3)
    show3d(x,clip=max(x))

def goCranfield2010():
  '''
  ****** beginning of SEG-Y file info ******
  file name = ../../../data/seis/cfd/2010_Match_dvd.sgy
  byte order = BIG_ENDIAN
  number of bytes = 633336744
  number of traces = 51726
  format = 1 (4-byte IBM floating point)
  units for spatial coordinates: ft (will be converted to km)
  indices and coordinates from trace headers:
    i2min =   114, i2max =   346 (inline indices)
    i3min =     3, i3max =   224 (crossline indices)
    xmin =   70.873839, xmax =   76.432491 (x coordinates, in km)
    ymin =  116.073488, ymax =  121.908684 (y coordinates, in km)
  grid sampling:
    n1 =  3001 (number of samples per trace)
    n2 =   233 (number of traces in inline direction)
    n3 =   222 (number of traces in crossline direction)
    d1 = 0.002000 (time sampling interval, in s)
    d2 = 0.025146 (inline sampling interval, in km)
    d3 = 0.025146 (crossline sampling interval, in km)
  grid corner points:
    i2min =   114, i3min =     3, x =   76.431105, y =  116.073488
    i2max =   346, i3min =     3, x =   76.432491, y =  121.907360
    i2min =   114, i3max =   224, x =   70.873839, y =  116.074812
    i2max =   346, i3max =   224, x =   70.875225, y =  121.908684
  grid azimuth:  0.01 degrees
  
  grid azimuth:  0.01 degrees
  ****** end of SEG-Y file info ******
  '''
  firstLook = False # fast, does not read all trace headers
  secondLook = True # slow, must read all trace headers
  writeImage = True # reads all traces, writes an image
  showImage = True # displays the image
  basedir = "../../../data/seis/cfd/"
  sgyfile = basedir+"2010_Match_dvd.sgy"
  datfile = basedir+"gx2010.dat"
  i1min,i1max,i2min,i2max,i3min,i3max = 0,3000,114,346,3,224
  n1,n2,n3 = 1+i1max-i1min,1+i2max-i2min,1+i3max-i3min
  si = SegyImage(sgyfile)
  if firstLook:
    si.printSummaryInfo();
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
  if secondLook:
    si.printAllInfo()
    plot23(si)
    plotXY(si)
  if writeImage:
    scale = 1.00
    si.writeFloats(datfile,scale,i1min,i1max,i2min,i2max,i3min,i3max,1,1)
  si.close()
  if showImage:
    x = readImage(datfile,n1,n2,n3)
    show3d(x,clip=max(x))

def goCranfield2007():
  '''
  ****** beginning of SEG-Y file info ******
  file name = ../../../data/seis/cfd/2007_baseline_dvd.sgy
  byte order = BIG_ENDIAN
  number of bytes = 660518424
  number of traces = 53946
  format = 1 (4-byte IBM floating point)
  units for spatial coordinates: ft (will be converted to km)
  indices and coordinates from trace headers:
    i2min =   109, i2max =   351 (inline indices)
    i3min =     3, i3max =   224 (crossline indices)
    xmin =   70.873811, xmax =   76.432519 (x coordinates, in km)
    ymin =  115.947758, ymax =  122.034414 (y coordinates, in km)
  grid sampling:
    n1 =  3001 (number of samples per trace)
    n2 =   243 (number of traces in inline direction)
    n3 =   222 (number of traces in crossline direction)
    d1 = 0.002000 (time sampling interval, in s)
    d2 = 0.025146 (inline sampling interval, in km)
    d3 = 0.025146 (crossline sampling interval, in km)
  grid corner points:
    i2min =   109, i3min =     3, x =   76.431077, y =  115.947758
    i2max =   351, i3min =     3, x =   76.432519, y =  122.033090
    i2min =   109, i3max =   224, x =   70.873811, y =  115.949082
    i2max =   351, i3max =   224, x =   70.875253, y =  122.034414
  grid azimuth:  0.01 degrees
  ****** end of SEG-Y file info ******
  '''
  firstLook = False # fast, does not read all trace headers
  secondLook = False # slow, must read all trace headers
  writeImage = False # reads all traces, writes an image
  showImage = True # displays the image
  basedir = "../../../data/seis/cfd/"
  sgyfile = basedir+"2007_baseline_dvd.sgy"
  datfile = basedir+"gx.dat"
  i1min,i1max,i2min,i2max,i3min,i3max = 0,3000,109,351,3,224
  n1,n2,n3 = 1+i1max-i1min,1+i2max-i2min,1+i3max-i3min
  si = SegyImage(sgyfile)
  if firstLook:
    si.printSummaryInfo();
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
  if secondLook:
    si.printAllInfo()
    plot23(si)
    plotXY(si)
  if writeImage:
    scale = 1.00
    si.writeFloats(datfile,scale,i1min,i1max,i2min,i2max,i3min,i3max,1,1)
  si.close()
  if showImage:
    x = readImage(datfile,n1,n2,n3)
    show3d(x,clip=max(x))

def goLulia():
  '''
  ****** beginning of SEG-Y file info ******
  file name = ../../../data/seis/beg/lulia/WFX_PSTM_90trim.sgy
  byte order = BIG_ENDIAN
  number of bytes = 5466696624
  number of traces = 1460121
  format = 1 (4-byte IBM floating point)
  units for spatial coordinates: ft (will be converted to km)
  indices and coordinates from trace headers:
    i2min =  9700, i2max = 10830 (inline indices)
    i3min = 12200, i3max = 13490 (crossline indices)
    xmin =  745.952966, xmax =  803.062275 (x coordinates, in km)
    ymin =   96.613770, ymax =  152.346782 (y coordinates, in km)
  grid sampling:
    n1 =   876 (number of samples per trace)
    n2 =  1131 (number of traces in inline direction)
    n3 =  1291 (number of traces in crossline direction)
    d1 = 0.004000 (time sampling interval, in s)
    d2 = 0.033528 (inline sampling interval, in km)
    d3 = 0.033528 (crossline sampling interval, in km)
  grid corner points:
    i2min =  9700, i3min = 12200, x =  767.438546, y =   96.613770
    i2max = 10830, i3min = 12200, x =  745.952966, y =  127.819050
    i2min =  9700, i3max = 13490, x =  803.062275, y =  121.141501
    i2max = 10830, i3max = 13490, x =  781.576694, y =  152.346782
  grid azimuth: -34.55 degrees
  ****** end of SEG-Y file info ******
  '''
  firstLook = False # fast, does not read all trace headers
  secondLook = False # slow, must read all trace headers
  writeImage = True # reads all traces, writes an image
  showImage = True # displays the image
  basedir = "../../../data/seis/beg/lulia/"
  sgyfile = basedir+"WFX_PSTM_90trim.sgy"
  datfile = basedir+"gx.dat"
  i1min,i1max,i2min,i2max,i3min,i3max = 0,875,9700,10830,12200,13490
  n1,n2,n3 = 1+i1max-i1min,1+i2max-i2min,1+i3max-i3min
  si = SegyImage(sgyfile)
  if firstLook:
    si.printSummaryInfo();
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
  if secondLook:
    si.printAllInfo()
    plot23(si)
    plotXY(si)
  if writeImage:
    scale = 1.00
    si.writeFloats(datfile,scale,i1min,i1max,i2min,i2max,i3min,i3max,1,1)
  si.close()
  if showImage:
    x = readImage(datfile,n1,n2,n3)
    show3d(x,clip=max(x))

def goNathan():
  '''
  ****** beginning of SEG-Y file info ******
  file name = ../../../data/seis/beg/nathan/Costa_Rica_It6-0-3000m.segy
  byte order = BIG_ENDIAN
  number of bytes = 8679832560
  number of traces = 3282840
  format = 1 (4-byte IBM floating point)
  units for spatial coordinates: m (will be converted to km)
  indices and coordinates from trace headers:
    i2min =  2220, i2max =  7193 (inline indices)
    i3min =  2012, i3max =  2671 (crossline indices)
    xmin =    0.000000, xmax =    0.000000 (x coordinates, in km)
    ymin =    0.000000, ymax =    0.000000 (y coordinates, in km)
  grid sampling:
    n1 =   601 (number of samples per trace)
    n2 =  4974 (number of traces in inline direction)
    n3 =   660 (number of traces in crossline direction)
    d1 = 0.005000 (time sampling interval, in s)
    d2 = 0.000000 (inline sampling interval, in km)
    d3 = 0.000000 (crossline sampling interval, in km)
  grid corner points:
    i2min =  2220, i3min =  2012, x =    0.000000, y =    0.000000
    i2max =  7193, i3min =  2012, x =    0.000000, y =    0.000000
    i2min =  2220, i3max =  2671, x =    0.000000, y =    0.000000
    i2max =  7193, i3max =  2671, x =    0.000000, y =    0.000000
  grid azimuth: 90.00 degrees
  ****** end of SEG-Y file info ******
  '''
  firstLook = False # fast, does not read all trace headers
  secondLook = False # slow, must read all trace headers
  writeImage = True # reads all traces, writes an image
  showImage = True # displays the image
  basedir = "../../../data/seis/beg/nathan/"
  sgyfile = basedir+"Costa_Rica_It6-0-3000m.segy"
  datfile = basedir+"gx.dat"
  i1min,i1max,i2min,i2max,i3min,i3max = 0,600,2220,7193,2012,2671
  n1,n2,n3 = 1+i1max-i1min,1+i2max-i2min,1+i3max-i3min
  si = SegyImage(sgyfile)
  if firstLook:
    si.printSummaryInfo();
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
  if secondLook:
    si.printAllInfo()
    plot23(si)
    plotXY(si)
  if writeImage:
    scale = 1.00
    si.writeFloats(datfile,scale,i1min,i1max,i2min,i2max,i3min,i3max,1,1)
  si.close()
  if showImage:
    x = readImage(datfile,n1,n2,n3)
    show3d(x,clip=max(x))

def goJake():
  '''
  ****** beginning of SEG-Y file info ******
  file name = ../../../data/seis/beg/jake/sub2/Subset2.sgy
  byte order = BIG_ENDIAN
  number of bytes = 5449815144
  number of traces = 2803401
  format = 1 (4-byte IBM floating point)
  units for spatial coordinates: m (will be converted to km)
  indices and coordinates from trace headers:
    i2min =  6500, i2max = 10500 (inline indices)
    i3min =  3300, i3max =  4700 (crossline indices)
    xmin = 2559.479500, xmax = 2589.977250 (x coordinates, in km)
    ymin = 6246.370000, ymax = 6274.663000 (y coordinates, in km)
  grid sampling:
    n1 =   426 (number of samples per trace)
    n2 =  4001 (number of traces in inline direction)
    n3 =  1401 (number of traces in crossline direction)
    d1 = 0.004000 (time sampling interval, in s)
    d2 = 0.006250 (inline sampling interval, in km)
    d3 = 0.012500 (crossline sampling interval, in km)
  grid corner points:
    i2min =  6500, i3min =  3300, x = 2569.010750, y = 6246.370000
    i2max = 10500, i3min =  3300, x = 2589.977250, y = 6259.986000
    i2min =  6500, i3max =  4700, x = 2559.479500, y = 6261.046500
    i2max = 10500, i3max =  4700, x = 2580.446000, y = 6274.662500
  grid azimuth: 57.00 degrees
  ****** end of SEG-Y file info ******
  '''
  firstLook = False # fast, does not read all trace headers
  secondLook = False # slow, must read all trace headers
  writeImage = True # reads all traces, writes an image
  showImage = False # displays the image
  basedir = "../../../data/seis/beg/jake/sub2/"
  sgyfile = basedir+"Subset2.sgy"
  datfile = basedir+"gx.dat"
  i1min,i1max,i2min,i2max,i3min,i3max = 0,425,6500,10500,3300,4700
  n1,n2,n3 = 1+i1max-i1min,1+i2max-i2min,1+i3max-i3min
  si = SegyImage(sgyfile)
  if firstLook:
    si.printSummaryInfo();
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
  if secondLook:
    si.printAllInfo()
    plot23(si)
    plotXY(si)
  if writeImage:
    scale = 1.00
    si.writeFloats(datfile,scale,i1min,i1max,i2min,i2max,i3min,i3max,2,1)
  si.close()
  if showImage:
    x = readImage(datfile,n1,n2,n3)
    show3d(x,clip=max(x))

def goF3dUnc():
  '''
  ****** beginning of SEG-Y file info ******
  file name = ../../../data/seis/f3d/unc.sgy
  byte order = BIG_ENDIAN
  number of bytes = 1292686488
  number of traces = 619101
  format = 1 (4-byte IBM floating point)
  units for spatial coordinates: m (will be converted to km)
  indices and coordinates from trace headers:
    i2min =   300, i2max =  1250 (inline indices)
    i3min =   100, i3max =   750 (crossline indices)
    xmin =  605.381000, xmax =  629.576000 (x coordinates, in km)
    ymin = 6073.556000, ymax = 6090.463000 (y coordinates, in km)
  grid sampling:
    n1 =   462 (number of samples per trace)
    n2 =   951 (number of traces in inline direction)
    n3 =   651 (number of traces in crossline direction)
    d1 = 0.004000 (time sampling interval, in s)
    d2 = 0.025000 (inline sampling interval, in km)
    d3 = 0.025001 (crossline sampling interval, in km)
  grid corner points:
    i2min =   300, i3min =   100, x =  605.835000, y = 6073.556000
    i2max =  1250, i3min =   100, x =  629.576000, y = 6074.220000
    i2min =   300, i3max =   750, x =  605.381000, y = 6089.800000
    i2max =  1250, i3max =   750, x =  629.122000, y = 6090.464000
  grid azimuth: 88.40 degrees
  ****** end of SEG-Y file info ******
  '''
  firstLook = False # fast, does not read all trace headers
  secondLook = False # slow, must read all trace headers
  writeImage = True # reads all traces, writes an image
  showImage = True # displays the image
  basedir = "../../../data/seis/f3d/"
  sgyfile = basedir+"unc.sgy"
  datfile = basedir+"unc.dat"
  i1min,i1max,i2min,i2max,i3min,i3max = 100,461,300,1250,100,750
  n1,n2,n3 = 1+i1max-i1min,1+i2max-i2min,1+i3max-i3min
  si = SegyImage(sgyfile)
  if firstLook:
    si.printSummaryInfo();
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
  if secondLook:
    si.printAllInfo()
    plot23(si)
    plotXY(si)
  if writeImage:
    scale = 1.00
    #si.writeFloats(datfile,scale,i1min,i1max,i2min,i2max,i3min,i3max,1,1)
    si.writeFloats(datfile,scale,i1min,i1max,i2min,i2max,i3min,i3max,1,1)
  si.close()
  if showImage:
    x = readImage(datfile,n1,n2,n3)
    show3d(x,clip=max(x))

def goHongliu():
  '''
  ****** beginning of SEG-Y file info ******
  file name = ../../../data/seis/hongliu/sgr.sgy
  byte order = BIG_ENDIAN
  number of bytes = 131064444
  number of traces = 40401
  format = 1 (4-byte IBM floating point)
  units for spatial coordinates: ft (will be converted to km)
  indices and coordinates from trace headers:
    i2min =   400, i2max =   600 (inline indices)
    i3min =  1500, i3max =  1700 (crossline indices)
    xmin =  534.667016, xmax =  544.012641 (x coordinates, in km)
    ymin =   90.798510, ymax =  100.137296 (y coordinates, in km)
  grid sampling:
    n1 =   751 (number of samples per trace)
    n2 =   201 (number of traces in inline direction)
    n3 =   201 (number of traces in crossline direction)
    d1 = 0.004000 (time sampling interval, in s)
    d2 = 0.033527 (inline sampling interval, in km)
    d3 = 0.033546 (crossline sampling interval, in km)
  grid corner points:
    i2min =   400, i3min =  1500, x =  534.667016, y =   94.644686
    i2max =   600, i3min =  1500, x =  540.159778, y =   90.798510
    i2min =   400, i3max =  1700, x =  538.519878, y =  100.137296
    i2max =   600, i3max =  1700, x =  544.012641, y =   96.291121
  grid azimuth: 125.00 degrees
  ****** end of SEG-Y file info ******
  '''
  firstLook = False # fast, does not read all trace headers
  secondLook = False # slow, must read all trace headers
  writeImage = False # reads all traces, writes an image
  showImage = True # displays the image
  basedir = "../../../data/seis/hongliu/"
  sgyfile = basedir+"sgr.sgy"
  datfile = basedir+"gx.dat"
  i1min,i1max,i2min,i2max,i3min,i3max = 0,750,400,600,1500,1700
  n1,n2,n3 = 1+i1max-i1min,1+i2max-i2min,1+i3max-i3min
  si = SegyImage(sgyfile)
  if firstLook:
    si.printSummaryInfo();
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
  if secondLook:
    si.printAllInfo()
    plot23(si)
    plotXY(si)
  if writeImage:
    scale = 1.00
    si.writeFloats(datfile,scale,i1min,i1max,i2min,i2max,i3min,i3max,1,1)
  si.close()
  if showImage:
    x = readImage(datfile,n1,n2,n3)
    show3d(x,clip=max(x))

def goBp():
  '''
  ****** beginning of SEG-Y file info ******
  file name = ../../../data/seis/bp/1_120_angle30_tvsw_ftrim_55Hz.sgy
  byte order = BIG_ENDIAN
  number of bytes = 9002929356
  number of traces = 4011999
  format = 5 (4-byte IEEE floating point)
  units for spatial coordinates: m (will be converted to km)
  indices and coordinates from trace headers:
    i2min =  5400, i2max =  7509 (inline indices)
    i3min = 10744, i3max = 13561 (crossline indices)
    xmin =  422.267000, xmax =  477.717000 (x coordinates, in km)
    ymin = 2377.265000, ymax = 2445.636000 (y coordinates, in km)
  grid sampling:
    n1 =   501 (number of samples per trace)
    n2 =  2110 (number of traces in inline direction)
    n3 =  2818 (number of traces in crossline direction)
    d1 = 0.004000 (time sampling interval, in s)
    d2 = 0.025000 (inline sampling interval, in km)
    d3 = 0.025000 (crossline sampling interval, in km)
  grid corner points:
    i2min =  5400, i3min = 10744, x =  414.668656, y = 2392.604687
    i2max =  7509, i3min = 10744, x =  457.858656, y = 2362.362687
    i2min =  5400, i3max = 13561, x =  455.063545, y = 2450.293444
    i2max =  7509, i3max = 13561, x =  498.253545, y = 2420.051444
  grid azimuth: 125.00 degrees
  ****** end of SEG-Y file info ******
  '''
  firstLook = False # fast, does not read all trace headers
  secondLook = False # slow, must read all trace headers
  writeImage = True # reads all traces, writes an image
  showImage = True # displays the image
  basedir = "../../../data/seis/beg/bp/"
  sgyfile = basedir+"1_120_angle30_tvsw_ftrim_55Hz.sgy"
  #sgyfile = basedir+"bp.sgy"
  datfile = basedir+"gxs.dat"
  i1min,i1max,i2min,i2max,i3min,i3max =260,450,5400,6404,11400,12650
  n1,n2,n3 = 1+i1max-i1min,1+i2max-i2min,1+i3max-i3min
  si = SegyImage(sgyfile)
  if firstLook:
    si.printSummaryInfo();
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
  if secondLook:
    si.printAllInfo()
    plot23(si)
    plotXY(si)
  if writeImage:
    scale = 1.00
    si.writeFloats(datfile,scale,i1min,i1max,i2min,i2max,i3min,i3max,1,1)
  si.close()
  if showImage:
    x = readImage(datfile,n1,n2,n3)
    show3d(x,clip=max(x))

def goClyde():
  firstLook = False # fast, does not read all trace headers
  secondLook = False # slow, must read all trace headers
  writeImage = True # reads all traces, writes an image
  showImage = True # displays the image
  basedir = "../../../data/seis/pdgm/clyde/"
  sgyfile = basedir+"clyde.sgy"
  datfile = basedir+"gx.dat"
  i1min,i1max,i2min,i2max,i3min,i3max = 0,499,0,240299,0,0
  #i1min,i1max,i2min,i2max,i3min,i3max = 0,1500,1,345,2,188
  n1,n2,n3 = 1+i1max-i1min,1+i2max-i2min,1+i3max-i3min
  si = SegyImage(sgyfile)
  if firstLook:
    si.printSummaryInfo();
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
  if secondLook:
    si.printAllInfo()
    plot23(si)
    plotXY(si)
  if writeImage:
    scale = 1.00
    #i1min,i1max,i2min,i2max,i3min,i3max = 0,589,29,28+240,46,80+45
    si.writeFloats(datfile,scale,i1min,i1max,i2min,i2max,i3min,i3max)
  si.close()
  if showImage:
    x = readImage(datfile,n1,n2,n3)
    show3d(x,clip=max(x))

def goAns():
  """
  ***************************************************************************
  ****** beginning of SEG-Y file info ******
  file name = ../../../data/seis/ans/seismic.sgy
  byte order = BIG_ENDIAN
  number of bytes = 2367279210
  number of traces = 1587710
  format = 8 (1-byte two's complement integer)
  units for spatial coordinates: m (will be converted to km)
  indices and coordinates from trace headers:
    i2min =   976, i2max =  5370 (inline indices)
    i3min =  1001, i3max =  1800 (crossline indices)
    xmin =  161.359300, xmax =  217.928500 (x coordinates, in km)
    ymin = 7806.933400, ymax = 7851.231100 (y coordinates, in km)
  grid sampling:
    n1 =  1251 (number of samples per trace)
    n2 =  4395 (number of traces in inline direction)
    n3 =   800 (number of traces in crossline direction)
    d1 = 0.004000 (time sampling interval, in s)
    d2 = 0.012500 (inline sampling interval, in km)
    d3 = 0.025000 (crossline sampling interval, in km)
  grid corner points:
    i2min =   976, i3min =  1001, x =  172.507956, y = 7805.807225
    i2max =  5370, i3min =  1001, x =  218.005273, y = 7836.576535
    i2min =   976, i3max =  1800, x =  161.317772, y = 7822.353566
    i2max =  5370, i3max =  1800, x =  206.815090, y = 7853.122876
  grid azimuth: 55.93 degrees
  ****** end of SEG-Y file info ******
  """
  firstLook = False # fast, does not read all trace headers
  secondLook = False # slow, must read all trace headers
  writeImage = True # reads all traces, writes an image
  showImage = True # displays the image
  basedir = "../../../data/seis/ans/"
  sgyfile = basedir+"seismic.sgy"
  datfile = basedir+"gxSub.dat"
  i1min,i1max,i2min,i2max,i3min,i3max = 455,655,3836,5210,1065,1700
  #i1min,i1max,i2min,i2max,i3min,i3max = 655,955,3836,5210,1065,1700
  n1,n2,n3 = 1+i1max-i1min,1+(i2max-i2min)/2,1+i3max-i3min
  print n1
  print n2
  print n3
  si = SegyImage(sgyfile)
  if firstLook:
    si.printSummaryInfo();
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
  if secondLook:
    si.printAllInfo()
    plot23(si)
    plotXY(si)
  if writeImage:
    scale = 0.02
    si.writeFloats(datfile,scale,i1min,i1max,i2min,i2max,i3min,i3max,2,1)
  si.close()
  if showImage:
    x = readImage(datfile,n1,n2,n3)
    show3d(x,clip=1.0)

def goGbc():
  """
  ***************************************************************************
  ****** beginning of SEG-Y file info ******
  file name = /data/seis/gbc/sgy/pwave_fs_NOCUT_npa.sgy
  byte order = BIG_ENDIAN
  number of bytes = 176949360
  number of traces = 21474
  format = 1 (4-byte IBM floating point)
  units for spatial coordinates: ft (will be converted to km)
  indices from trace headers:
    i2min =  5669, i2max =  5818 (inline indices)
    i3min =  5947, i3max =  6091 (crossline indices)
  grid sampling:
    n1 =  2000 (number of samples per trace)
    n2 =   150 (number of traces in inline direction)
    n3 =   145 (number of traces in crossline direction)
    d1 = 0.002000 (time sampling interval, in s)
    d2 = 0.033531 (inline sampling interval, in km)
    d3 = 0.033529 (crossline sampling interval, in km)
  ****** end of SEG-Y file info ******
  """
  imageType = "p" # which image to process
  firstLook = False # fast, does not read all trace headers
  secondLook = False # slow, must read all trace headers
  writeImage = True # reads all traces, writes an image
  showImage = True # displays the image
  basedir = "/data/seis/gbc/"
  if imageType=="p":
    sgyfile = basedir+"sgy/pwave_fs_NOCUT_npa.sgy"
    datfile = basedir+"dat/p.dat"
  elif imageType=="s1":
    sgyfile = basedir+"sgy/s1_nmut_fs_ps_NOCUT_npa.sgy"
    datfile = basedir+"dat/s1.dat"
  elif imageType=="s2":
    sgyfile = basedir+"sgy/s2_nmut_fs_ps_NOCUT_npa.sgy"
    datfile = basedir+"dat/s2.dat"
  i1min,i1max = 0,1999
  i2min,i2max = 5669,5818
  i3min,i3max = 5947,6091
  n1,n2,n3 = 1+i1max-i1min,1+i2max-i2min,1+i3max-i3min
  si = SegyImage(sgyfile)
  if firstLook:
    si.printSummaryInfo();
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
    plotIbmIeeeFloats(si)
  if secondLook:
    si.printAllInfo()
    plot23(si)
    plotXY(si)
  if writeImage:
    scale = 0.0001
    si.writeFloats(datfile,scale,i1min,i1max,i2min,i2max,i3min,i3max)
  si.close()
  if showImage:
    f = readImage(datfile,n1,n2,n3)
    show3d(f,clip=1.0)

def goMbs():
  """
  ***************************************************************************
  PstmSmall (pstm_fraw.sgy):
  ****** beginning of SEG-Y file info ******
  file name = /data/seis/mbs/PstmSmall/Marathon20070228/pstm_fraw.sgy
  byte order = BIG_ENDIAN
  number of bytes = 959341200
  number of traces = 153740
  format = 1 (4-byte IBM floating point)
  units for spatial coordinates: ft (will be converted to km)
  indices from trace headers:
    i2min =   601, i2max =  1048 (inline indices)
    i3min =  1001, i3max =  1422 (crossline indices)
  grid sampling:
    n1 =  1500 (number of samples per trace)
    n2 =   448 (number of traces in inline direction)
    n3 =   422 (number of traces in crossline direction)
    d1 = 0.002000 (time sampling interval, in s)
    d2 = 0.016764 (inline sampling interval, in km)
    d3 = 0.016764 (crossline sampling interval, in km)
  ****** end of SEG-Y file info ******
  ***************************************************************************
  PstmLarge
  ****** beginning of SEG-Y file info ******
  file name = /data/seis/mbs/PstmLarge/pstm_raw_cut.sgy
  byte order = BIG_ENDIAN
  number of bytes = 4647376320
  number of traces = 795783
  format = 1 (4-byte IBM floating point)
  units for spatial coordinates: ft (will be converted to km)
  indices from trace headers:
    i2min =   350, i2max =  1468 (inline indices)
    i3min =   234, i3max =  1422 (crossline indices)
  grid sampling:
    n1 =  1400 (number of samples per trace)
    n2 =  1119 (number of traces in inline direction)
    n3 =  1189 (number of traces in crossline direction)
    d1 = 0.002000 (time sampling interval, in s)
    d2 = 0.016764 (inline sampling interval, in km)
    d3 = 0.016764 (crossline sampling interval, in km)
  ****** end of SEG-Y file info ******
  good subset:
    i1min,i1max = 150, 650,  m1 = 501
    i2min,i2max = 490,1258,  m2 = 769
    i3min,i3max = 358, 917,  m3 = 560
  """
  imageType = "PstmSmall" # which image to process
  firstLook = False # fast, does not read all trace headers
  secondLook = False # slow, must read all trace headers
  writeImage = False # reads all traces, writes an image
  showImage = True # displays the image
  basedir = "/data/seis/mbs/"
  if imageType=="PstmSmall":
    sgyfile = basedir+"PstmSmall/Marathon20070228/pstm_fraw.sgy"
    datfile = basedir+"dat/pstm_fraw_s1.dat"
    i1min,i1max,i2min,i2max,i3min,i3max = 150,650,601,1048,1001,1422
  elif imageType=="PstmLarge":
    sgyfile = basedir+"PstmLarge/pstm_raw_cut.sgy"
    datfile = basedir+"dat/pstm_raw_s1.dat"
    i1min,i1max,i2min,i2max,i3min,i3max = 150,650,490,1258,358,917
  n1,n2,n3 = 1+i1max-i1min,1+i2max-i2min,1+i3max-i3min
  si = SegyImage(sgyfile)
  if firstLook:
    si.printSummaryInfo();
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
    plotIbmIeeeFloats(si)
  if secondLook:
    si.printAllInfo()
    #plot23(si)
    #plotXY(si)
  if writeImage:
    scale = 0.0001
    si.writeFloats(datfile,scale,i1min,i1max,i2min,i2max,i3min,i3max)
  si.close()
  if showImage:
    f = readImage(datfile,n1,n2,n3)
    show3d(f,clip=1.0)

def goNorne():
  """
  ***************************************************************************
  ****** beginning of SEG-Y file info ******
  file name = /data/seis/norne/sgy/norne4d_2006-full.sgy
  byte order = BIG_ENDIAN
  number of bytes = 1363689924
  number of traces = 321321
  format = 1 (4-byte IBM floating point)
  units for spatial coordinates: m (will be converted to km)
  indices and coordinates from trace headers:
    i2min =  1300, i2max =  2300 (inline indices)
    i3min =   970, i3max =  1290 (crossline indices)
    xmin =  453.320000, xmax =  464.634000 (x coordinates, in km)
    ymin = 7317.354880, ymax = 7329.340160 (y coordinates, in km)
  grid sampling:
    n1 =  1001 (number of samples per trace)
    n2 =  1001 (number of traces in inline direction)
    n3 =   321 (number of traces in crossline direction)
    d1 = 0.004000 (time sampling interval, in s)
    d2 = 0.012501 (inline sampling interval, in km)
    d3 = 0.012500 (crossline sampling interval, in km)
  grid corner points:
    i2min =  1300, i3min =   970, x =  453.320000, y = 7320.021120
    i2max =  2300, i3min =   970, x =  461.652000, y = 7329.340160
    i2min =  1300, i3max =  1290, x =  456.302000, y = 7317.354880
    i2max =  2300, i3max =  1290, x =  464.634000, y = 7326.673920
  grid azimuth: 41.80 degrees
  ****** end of SEG-Y file info ******
  Full Norne corner coordinates (Knut, 16/7-09)
  line	trace	X	Y
  970	2300	461652	7329340.16
  970	1300	453320	7320021.12
  1290	1300	456302	7317355
  1290	2300	464634	7326674
  ***************************************************************************
  """
  firstLook = False # fast, does not read all trace headers
  secondLook = False # slow, must read all trace headers
  writeImage = False # reads all traces, writes an image
  showImage = True # displays the image
  i1min,i1max,i2min,i2max,i3min,i3max = 0,1000,1300,2300,970,1290
  n1,n2,n3 = 1+i1max-i1min,1+i2max-i2min,1+i3max-i3min
  basedir = "/data/seis/norne/"
  sgyfiles = [
    "norne4d_2001-full.sgy",
    "norne4d_2003-full.sgy",
    "norne4d_2004-full.sgy",
    "norne4d_2006-full.sgy",
    "norne4d_2006-near.sgy",
    "norne4d_2006-mid.sgy",
    "norne4d_2006-far.sgy"]
  datfiles = [
    "norne2001full.dat",
    "norne2003full.dat",
    "norne2004full.dat",
    "norne2006full.dat",
    "norne2006near.dat",
    "norne2006mid.dat",
    "norne2006far.dat"]
  nfile = len(sgyfiles)
  for ifile in range(nfile):
    sgyfile = basedir+"sgy/"+sgyfiles[ifile]
    datfile = basedir+"dat/"+datfiles[ifile]
    si = SegyImage(sgyfile)
    if firstLook:
      si.printSummaryInfo();
      #si.printBinaryHeader()
      #si.printTraceHeader(0)
      #si.printTraceHeader(1)
      #plotIbmIeeeFloats(si)
    if secondLook:
      si.printAllInfo()
      #plot23(si)
      #plotXY(si)
    if writeImage:
      scale = 0.001
      si.writeFloats(datfile,scale,i1min,i1max,i2min,i2max,i3min,i3max)
    si.close()
    if showImage:
      x = readImage(datfile,n1,n2,n3)
      show3d(x,clip=1.0)

def goSino():
  """
  Two 2D (x- and z-component) images
  nt = 2001, dt = 0.004 s
  nx =  721, dx = 0.030 km?
  samples [0:1250] in z component ~ samples [0:2000] in x component
  """
  basedir = "/data/seis/sino/"
  for component in ["x","z"]:
    sgyfile = basedir+"260_"+component+"_201-921_stack.segy"
    datfile = basedir+component+"260.dat"
    i1min,i1max,i2min,i2max = 0,2000,201,921
    n1,n2 = 1+i1max-i1min,1+i2max-i2min
    si = SegyImage(sgyfile,ByteOrder.LITTLE_ENDIAN) # non-standard byte order!
    si.printSummaryInfo()
    #si.printBinaryHeader()
    #si.printTraceHeader(0)
    #si.printTraceHeader(1)
    #plotIbmIeeeFloats(si)
    si.setD2(0.015) # a guess, from looking at group coordinates in headers
    if component=="x":
      si.setFormat(5) # formats appear to be IEEE for x and IBM for z!???
    #si.printAllInfo()
    si.writeFloats(datfile)
    si.close()
    f = readImage(datfile,n1,n2)
    if component=="z":
      stretch(1.6,f)
    gain(500,f)
    show2d(f,title=component+" component")

def goF3d():
  """
  ***************************************************************************
  ****** beginning of SEG-Y file info ******
  file name = /data/seis/f3d/f3draw.sgy
  byte order = BIG_ENDIAN
  number of bytes = 699003060
  number of traces = 600515
  format = 3 (2-byte two's complement integer)
  units for spatial coordinates: m (will be converted to km)
  indices and coordinates from trace headers:
    i2min =   300, i2max =  1250 (inline indices)
    i3min =   100, i3max =   750 (crossline indices)
    xmin =  605.416700, xmax =  629.576300 (x coordinates, in km)
    ymin = 6073.556400, ymax = 6090.463200 (y coordinates, in km)
  grid sampling:
    n1 =   462 (number of samples per trace)
    n2 =   951 (number of traces in inline direction)
    n3 =   651 (number of traces in crossline direction)
    d1 = 0.004000 (time sampling interval, in s)
    d2 = 0.025000 (inline sampling interval, in km)
    d3 = 0.024999 (crossline sampling interval, in km)
  grid corner points:
    i2min =   300, i3min =   100, x =  605.835500, y = 6073.556400
    i2max =  1250, i3min =   100, x =  629.576300, y = 6074.219900
    i2min =   300, i3max =   750, x =  605.381800, y = 6089.799700
    i2max =  1250, i3max =   750, x =  629.122600, y = 6090.463200
  grid azimuth: 88.40 degrees
  ****** end of SEG-Y file info ******
  good subset with no dead traces
  i1min,i1max,i2min,i2max,i3min,i3max = 0,461,300,1250,100,690
  n1,n2,n3 = 462,951,591
  ***************************************************************************
  """
  firstLook = False # fast, does not read all trace headers
  secondLook = False # slow, must read all trace headers
  writeImage = True # reads all traces, writes an image
  showImage = True # displays the image
  basedir = "/data/seis/f3d/"
  sgyfile = basedir+"odt/seis/f3d.sgy"
  datfile = basedir+"seta/g.dat"
  #i1min,i1max,i2min,i2max,i3min,i3max = 0,461,300,1250,100,690
  i1min,i1max,i2min,i2max,i3min,i3max = 0,461,300,1250,100,750
  n1,n2,n3 = 1+i1max-i1min,1+i2max-i2min,1+i3max-i3min
  si = SegyImage(sgyfile)
  if firstLook:
    si.printSummaryInfo();
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
  if secondLook:
    si.printAllInfo()
    plot23(si)
    plotXY(si)
  if writeImage:
    scale = 0.0001
    si.writeFloats(datfile,scale,i1min,i1max,i2min,i2max,i3min,i3max)
  si.close()
  if showImage:
    x = readImage(datfile,n1,n2,n3)
    show3d(x,clip=1.0)

def goPnz():
  """
  ***************************************************************************
  ****** beginning of SEG-Y file info ******
  file name = /data/seis/par/sgy/Parihaka3d_raw.sgy
  byte order = BIG_ENDIAN
  number of bytes = 69479795156
  number of traces = 11127449
  format = 1 (4-byte IBM floating point)
  units for spatial coordinates: m (will be converted to km)
  indices and coordinates from trace headers:
    i2min =  2050, i2max = 14026 (inline indices)
    i3min =  1665, i3max =  5599 (crossline indices)
    xmin = 2546.537000, xmax = 2619.591000 (x coordinates, in km)
    ymin = 6226.154000, ymax = 6288.372000 (y coordinates, in km)
  grid sampling:
    n1 =  1501 (number of samples per trace)
    n2 = 11977 (number of traces in inline direction)
    n3 =  3935 (number of traces in crossline direction)
    d1 = 0.004000 (time sampling interval, in s)
    d2 = 0.006250 (inline sampling interval, in km)
    d3 = 0.012500 (crossline sampling interval, in km)
  grid corner points:
    i2min =  2050, i3min =  1665, x = 2556.816611, y = 6214.081487
    i2max = 14026, i3min =  1665, x = 2619.590620, y = 6254.847271
    i2min =  2050, i3max =  5599, x = 2530.034337, y = 6255.322955
    i2max = 14026, i3max =  5599, x = 2592.808346, y = 6296.088740
  grid azimuth: 57.00 degrees
  ****** end of SEG-Y file info ******
  NOTE:
  In the SEG-Y file, inline indices increment by 2, 
  so that the number n2 = 11977 of traces per line
  includes traces not actually present in the file.
  Ignoring those missing traces, d2 = d3 = 0.0125.
  ***************************************************************************
  """
  firstLook = True # fast, does not read all trace headers
  secondLook = False # slow, must read all trace headers
  writeImage = False # reads all traces, writes an image
  showImage = False # displays the image
  basedir = "/data/seis/par/"
  sgyfile = basedir+"sgy/Parihaka3d_raw.sgy"
  #sgyfile = basedir+"sgy/Parihaka3d_full.sgy"
  #datfile = basedir+"dat/pnztest.dat"
  #i1min,i1max,i2min,i2max,i3min,i3max = 0,1500,6500,10500,2100,4100
  #n1,n2,n3 = 1+i1max-i1min,1+i2max-i2min,1+i3max-i3min
  si = SegyImage(sgyfile)
  si.setInlineXlineBytes(197,201)
  if firstLook:
    si.printSummaryInfo();
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
  if secondLook:
    si.printAllInfo()
    plot23(si)
    plotXY(si)
  if writeImage:
    scale = 0.00001
    si.writeFloats(datfile,scale,i1min,i1max,i2min,i2max,i3min,i3max)
  si.close()
  if showImage:
    x = readImage(datfile,n1,n2,n3)
    show3d(x,clip=1.0)

def goSch():
  """
  ***************************************************************************
  ****** beginning of SEG-Y file info ******
  file name = /data/seis/sch/sgy/pack_0.sgy
  byte order = BIG_ENDIAN
  number of bytes = 2143537580
  number of traces = 343295
  format = 1 (4-byte IBM floating point)
  units for spatial coordinates: m (will be converted to km)
  indices and coordinates from trace headers:
    i2min =  1000, i2max =  1622 (inline indices)
    i3min =  1000, i3max =  1646 (crossline indices)
    xmin =  251.044000, xmax =  256.212000 (x coordinates, in km)
    ymin =  518.203000, ymax =  523.179000 (y coordinates, in km)
  grid sampling:
    n1 =  1501 (number of samples per trace)
    n2 =   623 (number of traces in inline direction)
    n3 =   647 (number of traces in crossline direction)
    d1 = 0.002000 (time sampling interval, in s)
    d2 = 0.008000 (inline sampling interval, in km)
    d3 = 0.008000 (crossline sampling interval, in km)
  grid corner points:
    i2min =  1000, i3min =  1000, x =  251.044000, y =  518.203000
    i2max =  1622, i3min =  1000, x =  251.044000, y =  523.179000
    i2min =  1000, i3max =  1646, x =  256.212000, y =  518.203000
    i2max =  1622, i3max =  1646, x =  256.212000, y =  523.179000
  grid azimuth:  0.00 degrees
  ****** end of SEG-Y file info ******
  NOTE:
  Most interesting faults appear to die out after i1max = 900.
  The file pack_0.sgy is larger than the file pack_1.sgy, and both
  have missing traces. Subset of pack_0.sgy with no missing traces:
  i1min,i1max,i2min,i2max,i3min,i3max = 0,900,1000,1550,1220,1646
  Another good subset of pack_0.sgy with no missing traces:
  i1min,i1max,i2min,i2max,i3min,i3max = 0,900,1000,1375,1100,1646
  ***************************************************************************
  """
  firstLook = False # fast, does not read all trace headers
  secondLook = False # slow, must read all trace headers
  writeImage = True # reads all traces, writes an image
  showImage = True # displays the image
  basedir = "/data/dhale/sch/"
  sgyfile = basedir+"sgy/pack_0.sgy"
  #datfile = basedir+"dat/s1/g0.dat" # subset 1
  #i1min,i1max,i2min,i2max,i3min,i3max = 0,900,1000,1550,1220,1646
  datfile = basedir+"dat/s2/g0.dat" # subset 2
  i1min,i1max,i2min,i2max,i3min,i3max = 0,900,1000,1375,1100,1646
  n1,n2,n3 = 1+i1max-i1min,1+i2max-i2min,1+i3max-i3min
  si = SegyImage(sgyfile)
  if firstLook:
    si.printSummaryInfo();
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
  if secondLook:
    si.printAllInfo()
    plot23(si)
    plotXY(si)
  if writeImage:
    scale = 0.0001
    si.writeFloats(datfile,scale,i1min,i1max,i2min,i2max,i3min,i3max)
  si.close()
  if showImage:
    x = readImage(datfile,n1,n2,n3)
    show3d(x,clip=1.0)

#############################################################################

def show2d(f,clip=None,title=None):
  print "show2d: f min =",min(f)," max =",max(f)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(f)
  if clip:
    pv.setClips(-clip,clip)
  else:
    pv.setPercentiles(2,98)
  if title:
    sp.setTitle(title)
  sp.setSize(600,1100)
  
def show3d(f,clip=None):
  print "show3d: f min =",min(f)," max =",max(f)
  frame = SimpleFrame()
  ipg = frame.addImagePanels(f)
  if clip:
    ipg.setClips(-clip,clip)
  frame.orbitView.setScale(2.0)
  frame.setSize(1000,1000)

def plot23(si):
  i2 = si.getI2sAsFloats()
  i3 = si.getI3sAsFloats()
  sp = SimplePlot()
  sp.setHLabel("inline sample index i2")
  sp.setVLabel("crossline sample index i3")
  pv = sp.addPoints(i2,i3)
  pv.setMarkStyle(PointsView.Mark.POINT);
  pv.setLineStyle(PointsView.Line.NONE);
  w,h = goodWidthHeight(i2,i3)
  sp.setSize(w,h)

def plotXY(si):
  x = si.getXs()
  y = si.getYs()
  sp = SimplePlot()
  sp.setHLabel("x coordinate (km)")
  sp.setVLabel("y coordinate (km)")
  pv = sp.addPoints(x,y)
  pv.setMarkStyle(PointsView.Mark.POINT);
  pv.setLineStyle(PointsView.Line.NONE);
  w,h = goodWidthHeight(x,y)
  sp.setSize(w,h)

def goodWidthHeight(x,y):
  xmin,xmax = min(x),max(x)
  ymin,ymax = min(y),max(y)
  w,h = 1000,1000
  if (xmax-xmin)>(ymax-ymin):
    h = int(h*(ymax-ymin)/(xmax-xmin))
  else:
    w = int(w*(xmax-xmin)/(ymax-ymin))
  return w,h

def readImage(datfile,n1,n2,n3=1):
  if n3==1:
    x = zerofloat(n1,n2)
  else:
    x = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(datfile)
  ais.readFloats(x)
  ais.close()
  return x

def writeImage(datfile,x):
  aos = ArrayOutputStream(datfile)
  aos.writeFloats(x)
  aos.close()

def plotIbmIeeeFloats(si):
  ntrace = si.countTraces()
  itrace = ntrace/2
  fmt = si.getFormat()
  si.setFormat(1) # IBM floats
  fibm = si.getTrace(itrace)
  si.setFormat(5) # IEEE floats
  fieee = si.getTrace(itrace)
  si.setFormat(fmt)
  pp = PlotPanel(2,1)
  pp.setTitle("IBM (top) versus IEEE (bottom)")
  pp.setHLabel(0,"Sample index")
  pp.setVLabel(0,"IBM amplitude")
  pp.setVLabel(1,"IEEE amplitude")
  pp.addPoints(0,0,fibm)
  pp.addPoints(1,0,fieee)
  pf = PlotFrame(pp)
  pf.setSize(1000,800)
  pf.setVisible(True)

def lowpass(f3db,f):
  bf = ButterworthFilter(f3db,6,ButterworthFilter.Type.LOW_PASS)
  bf.apply1ForwardReverse(f,f)

def gain(hw,f):
  g = mul(f,f)
  RecursiveExponentialFilter(hw).apply1(g,g)
  div(f,sqrt(g),f)

def stretch(c,f):
  n1,n2 = len(f[0]),len(f)
  t = rampfloat(0.0,1.0/c,n1)
  si = SincInterpolator()
  g = zerofloat(n1)
  for i2 in range(n2):
    si.interpolate(n1,1.0,0.0,f[i2],n1,1.0/c,0.0,g)
    copy(g,f[i2])

def spow(p,f):
  return mul(sgn(f),pow(abs(f),p))

def slog(f):
  return mul(sgn(f),log(add(1.0,abs(f))))

def sexp(f):
  return mul(sgn(f),sub(exp(abs(f)),1.0))

#############################################################################
# Parihaka functions below are outdated; kept here for updating later

def goParihakaOld():
  """
  INLINES  : 1665 - 5599 (INC 1)   CROSSLINES : 2050 - 14026 (INC 2)
  TIME: 1501 samples  (6 sec, interval 4ms)
  nbytes = 69,479,807,644
  ntrace = (nbytes-nhead-nbhed)/(240+4*n1)
  """
  fmt = 1
  n1 = 1501
  datdir = "/data/seis/nz/par/"
  sgyfile = datdir+"Parihaka3d_raw.sgy"
  datfile = datdir+"par11.dat"
  #displayParihaka(datfile)
  #makeSubsetsParihaka(datdir)
  #displaySubsetParihaka(datfile,0,500,500,751,501,501)
  #bigSubsetParihaka(n1,sgyfile,datfile)
  #convert(n1,n2,n3,sgyfile,datfile)

def displayParihaka(datfile):
  n1,n2,n3 = 751,1001,1001
  x = readImage(datfile,n1,n2,n3)
  display3d(x,1.0e5)

def makeSubsetsParihaka(datdir):
  m1,m2,m3 = 751,1001,1001
  x = readParihaka(datdir+"parBig.dat",0,0,0,m1,m2,m3)
  writeImage(datdir+"par00.dat",x)
  x = None
  x = readParihaka(datdir+"parBig.dat",0,1000,0,m1,m2,m3)
  writeImage(datdir+"par01.dat",x)
  x = None
  x = readParihaka(datdir+"parBig.dat",0,0,1000,m1,m2,m3)
  writeImage(datdir+"par10.dat",x)
  x = None
  x = readParihaka(datdir+"parBig.dat",0,1000,1000,m1,m2,m3)
  writeImage(datdir+"par11.dat",x)
  display3d(x,1.0e5)

def displaySubsetParihaka(datfile,j1,j2,j3,m1,m2,m3):
  x = readSubsetParihaka(datfile,j1,j2,j3,m1,m2,m3)
  display3d(x,1.0e5)

def readParihaka(datfile,j1,j2,j3,m1,m2,m3):
  n1,n2,n3 = 1501,2001,2001
  m1 = min(m1,n1)
  m2 = min(m2,n2)
  m3 = min(m3,n3)
  j1 = min(j1,n1-m1)
  j2 = min(j2,n2-m2)
  j3 = min(j3,n3-m3)
  x = zerofloat(m1,m2,m3)
  ais = ArrayInputStream(datfile)
  for i3 in range(j3):
    ais.skipBytes(4*n1*n2)
  ais.skipBytes(4*(j1+n1*j2))
  for i3 in range(m3):
    if i3%10==0: 
      print "i3 =",i3
    for i2 in range(m2):
      ais.readFloats(x[i3][i2])
      ais.skipBytes(4*(n1-m1))
    ais.skipBytes(4*n1*(n2-m2))
  return x

def bigSubsetParihaka(n1,sgyfile,datfile):
  """ big subset 1501x2001x2001 ~ 24 GB
  i1min,i1max =    0, 1500 # time samples
  i2min,i2max = 6500,10500 # increment by 2
  i3min,i3max = 2100, 4100 # increment by 1
  """
  """ A: small subset 501x501x501 ~ 500 MB
  i1min,i1max =    0, 500 # time samples
  i2min,i2max = 7500,8500 # increment by 2
  i3min,i3max = 2500,3000 # increment by 1
  """
  """ B: subset 751x1001x1001 ~ 500 MB
  i1min,i1max =    0, 750 # time samples
  i2min,i2max = 7000,9000 # increment by 2
  i3min,i3max = 3000,4000 # increment by 1
  """
  i1min,i1max =    0, 1500 # time samples
  i2min,i2max = 6500,10500 # increment by 2
  i3min,i3max = 2100, 4100 # increment by 1
  m1 = 1+i1max-i1min
  m2 = 1+(i2max-i2min)/2
  m3 = 1+i3max-i3min
  m23 = m2*m3
  ais = ArrayInputStream(sgyfile,bo)
  aos = ArrayOutputStream(datfile)
  ais.skipBytes(nhead)
  ais.skipBytes(nbhed)
  h = zeroint(nthed/4)
  x = zeroint(n1)
  y = zeroint(m1)
  z = zerofloat(m1)
  nread = 0
  i23 = 0
  while i23<m23:
    nread += 1
    ais.readInts(h)
    i3,i2 = h[49],h[50]
    if i2min<=i2 and i2<=i2max and i3min<=i3 and i3<=i3max:
      i23 += 1
      if i2==i2min:
        print "nread =",nread," i3min =",i3min," i3 =",i3," i3max =",i3max
      ais.readInts(x) # read trace samples
      copy(m1,x,y) # keep only m1 samples
      IbmIeee.ibmToFloat(y,z) # ibm to ieee
      aos.writeFloats(z) # write trace samples
    else:
      ais.skipBytes(4*n1)
  ais.close()
  aos.close()

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
