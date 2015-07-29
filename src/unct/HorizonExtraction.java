/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

package unct;

import edu.mines.jtk.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.awt.ColorMap;
import static edu.mines.jtk.util.ArrayMath.*;

import ipfx.*;
import util.*;

/**
 * From a precomputed horizon volume x1(u1,x2,x3) or 
 * a relative geologic time volume u1(x1,x2,x3), seismic horizons can be easily 
 * extracted by simply horizontal slicing in the horizon volume at each relative 
 * geologic time u1
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.05.27
 */

public class HorizonExtraction {
  /**
   * Sets horizon volume and relative geologc time volume for horizon extraction.
   * @param x1 an array of horizon volume x1(u1,x2,x3).
   * @param u1 an array of relative geologic time volume u1(x1,x2,x3).
   */
  public HorizonExtraction(
    Sampling s1, Sampling s2, Sampling s3, float[][][] x1, float[][][] u1) 
  {
    _s1 = s1;
    _s2 = s2;
    _s3 = s3;
    _u1 = u1;
    if(x1!=null){_x1=x1;}
    else {_x1=horizonVolumeFromRgt(s1,u1);}
  }

  public float[][][] trigSurfaces(
    float color, float[][][] surfs, FaultSkin[] sks, float[][][] f) 
  {
    int ns = surfs.length;
    int n3 = _u1.length;
    int n2 = _u1[0].length;
    int n1 = _u1[0][0].length;
    short[][][] mk = zeroOnFaults(n1,n2,n3,sks);
    float[][][] sfs = new float[2][][];
    for (int is=0; is<ns; ++is) {
      float[][] fi = null;
      if(f!=null) {fi=f[is];}
      sfs[is] = buildTrigs(_s3,_s2,color,surfs[is],fi,mk);
    }
    return sfs;
  }

  public float[][][] ulOnSurface(float[][][] surfs, float[][][] f) {
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    int ns = surfs.length;
    float[][][] fsurf = new float[ns][n3][n2];
    SincInterpolator si = new SincInterpolator();
    for (int is=0; is<ns; ++is) {
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      float x3 = i3;
      float x2 = i2;
      float x1 = surfs[is][i3][i2];
      fsurf[is][i3][i2] = si.interpolate(s1,s2,s3,f,x1,x2,x3);
    }}}
    return fsurf;
  }

  public float[][][] horizonCurves(
    float[][][] uncs, int k2, int k3, FaultSkin[] sks) 
  {
    int n3 = _u1.length;
    int n2 = _u1[0].length;
    int n1 = _u1[0][0].length;
    int n1s = uncs.length;
    float[][][] xbs = new float[n1s*200][2][];
    float[][] surfs = new float[n1s*200][max(n2,n3)*3];
    short[][][] mk = zeroOnFaults(n1,n2,n3,sks);
    float rgtMin = min(_u1)+5;
    float rgtMax = max(_u1)-5;
    ColorMap cp = new ColorMap(rgtMin,rgtMax,ColorMap.JET);

    int k = 0;
    for (int i1s=0; i1s<n1s; ++i1s) {
      float u1i = i1s;
      int p = 0;
      for (int i3=0; i3<n3; ++i3) {
        float x3 = i3;
        float x2 = k2;
        float x1 = uncs[i1s][i3][k2];
        if(x1>=n1) {continue;}
        if(onFault(x3,x2,x1,mk)){
          float[] vs = fillfloat(u1i,p);
          xbs[k][0] = copy(p,surfs[k]);
          xbs[k][1] = cp.getRgbFloats(vs);
          if(p>0){k ++;}
          p = 0;
          continue;
        }
        surfs[k][p++] = x3;
        surfs[k][p++] = x2;
        surfs[k][p++] = x1;
      }
      float[] vs = fillfloat(u1i,p);
      xbs[k][0] = copy(p,surfs[k]);
      xbs[k][1] = cp.getRgbFloats(vs);
      if(p>0){k ++;}
      p = 0;
      for (int i2=0; i2<n2; ++i2) {
        float x3 = k3;
        float x2 = i2;
        float x1 = uncs[i1s][k3][i2];
        if(x1>=n1) {continue;}
        if(onFault(x3,x2,x1,mk)){
          vs = fillfloat(u1i,p);
          xbs[k][0] = copy(p,surfs[k]);
          xbs[k][1] = cp.getRgbFloats(vs);
          if(p>0){k ++;}
          p = 0;
          continue;
        }
        surfs[k][p++] = x3;
        surfs[k][p++] = x2;
        surfs[k][p++] = x1;
      }
      vs = fillfloat(u1i,p);
      xbs[k][0] = copy(p,surfs[k]);
      xbs[k][1] = cp.getRgbFloats(vs);
      if(p>0){k ++;}
      p = 0;
    }
    float[][][] sfs = new float[k][2][];
    for (int i=0; i<k; ++i) {
      int np = xbs[i][0].length;
      sfs[i][0] = copy(np,0,xbs[i][0]); 
      sfs[i][1] = copy(np,0,xbs[i][1]); 
      int j = 0;
      for (int ip=0; ip<np; ip+=3) {
        sfs[i][1][j++] = 1;
        sfs[i][1][j++] = 0;
        sfs[i][1][j++] = 1;
      }
    }
    return sfs;
  }


  public float[][][] horizonCurves(
    Sampling t1, int k1, int k2, int k3, FaultSkin[] sks, float[][][] uncs) 
  {
    int n3 = _u1.length;
    int n2 = _u1[0].length;
    int n1 = _u1[0][0].length;
    double[] u1s = t1.getValues();
    int n1s = u1s.length;
    float[][][] xbs = new float[n1s*200][2][];
    float[][] surfs = new float[n1s*200][max(n2,n3)*3];
    float d1 = (float)_s1.getDelta();
    float f1 = (float)_s1.getFirst();
    short[][][] uk = zeroUncs(n1,n2,n3,uncs);
    short[][][] mk = zeroOnFaults(n1,n2,n3,sks);
    float rgtMin = min(_u1)+5;
    float rgtMax = max(_u1)-5;
    ColorMap cp = new ColorMap(rgtMin,rgtMax,ColorMap.JET);

    int k = 0;
    for (int i1s=0; i1s<n1s; ++i1s) {
      float u1i = (float)u1s[i1s];
      int i1 = round((u1i-f1)/d1);
      int p = 0;
      for (int i3=0; i3<n3; ++i3) {
        float x3 = i3;
        float x2 = k2;
        float x1 = _x1[i3][k2][i1];
        if(x1>=n1) {continue;}
        if(onFault(x3,x2,x1,mk)||onUnc(x3,x2,x1,uk)){
          float[] vs = fillfloat(u1i,p);
          xbs[k][0] = copy(p,surfs[k]);
          xbs[k][1] = cp.getRgbFloats(vs);
          if(p>6){k ++;}
          p = 0;
          continue;
        }
        surfs[k][p++] = x3;
        surfs[k][p++] = x2;
        surfs[k][p++] = x1;
      }
      float[] vs = fillfloat(u1i,p);
      xbs[k][0] = copy(p,surfs[k]);
      xbs[k][1] = cp.getRgbFloats(vs);
      if(p>6){k ++;}
      p = 0;
      for (int i2=0; i2<n2; ++i2) {
        float x3 = k3;
        float x2 = i2;
        float x1 = _x1[k3][i2][i1];
        if(x1>=n1) {continue;}
        if(onFault(x3,x2,x1,mk)||onUnc(x3,x2,x1,uk)){
          vs = fillfloat(u1i,p);
          xbs[k][0] = copy(p,surfs[k]);
          xbs[k][1] = cp.getRgbFloats(vs);
          if(p>6){k ++;}
          p = 0;
          continue;
        }
        surfs[k][p++] = x3;
        surfs[k][p++] = x2;
        surfs[k][p++] = x1;
      }
      vs = fillfloat(u1i,p);
      xbs[k][0] = copy(p,surfs[k]);
      xbs[k][1] = cp.getRgbFloats(vs);
      if(p>6){k ++;}
      p = 0;
    }
    float[][][] sfs = new float[k][2][];
    for (int i=0; i<k; ++i) {
      int np = xbs[i][0].length;
      sfs[i][0] = copy(np,0,xbs[i][0]); 
      sfs[i][1] = copy(np,0,xbs[i][1]); 
    }
    return sfs;
  }

  public float[][][] zContours(int k1, FaultSkin[] sks) {
    int n3 = _u1.length;
    int n2 = _u1[0].length;
    int n1 = _u1[0][0].length;
    short[][][] mk = zeroOnFaults(n1,n2,n3,sks);
    float[][] r = new float[n3][n2];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      r[i3][i2] = _u1[i3][i2][k1];
    }}
    float dr = 2.0f;
    float rmin = min(r);
    float rmax = max(r);
    System.out.println("rmin="+rmin);
    System.out.println("rmax="+rmax);
    int nr = round((rmax-rmin)/dr)-1;
    Sampling sr = new Sampling(nr,dr,rmin);
    double[] rs = sr.getValues();
    ContourMaker cm = new ContourMaker(r);
    float fc = 0.5f*(rmin+rmax)+16;
    ContourMaker.Contour ct = cm.getContour(fc);
    int ns = ct.ns;
    float rgtMin = min(_u1)+5;
    float rgtMax = max(_u1)-5;
    ColorMap cp = new ColorMap(rgtMin,rgtMax,ColorMap.JET);
    float[][][] surfs = new float[ns][2][n2*n3*10];
    for (int is=0; is<ns; ++is) {
      int c1 = 0;
      int c2 = 0;
      float[] x1 = ct.x1.get(is);
      float[] x2 = ct.x2.get(is);
      int np = x1.length;
      for (int ip=0; ip<np; ++ip) {
        surfs[is][0][c1++] = x2[ip];
        surfs[is][0][c1++] = x1[ip];
        surfs[is][0][c1++] = k1;
        surfs[is][1][c2++] = fc;
        surfs[is][1][c2++] = fc;
        surfs[is][1][c2++] = fc;
      }
      surfs[is] = copy(c1,2,surfs[is]);
      surfs[is][1] = cp.getRgbFloats(surfs[is][1]);
    }
    return surfs;
  }

  /*
  public float[][][] horizonLines(
    Sampling t1, int k2, int k3, FaultSkin[] sks) 
  {
    double[] u1s = t1.getValues();
    int n1s = u1s.length;
    float[][][] surfs = new float[n1s*2][][];
    int k = 0;
    for (int i1s=0; i1s<n1s; i1s++) {
      float u1i = (float)u1s[i1s];
      surfs[k++] = horizonLine2(k2,u1i,sks);
      surfs[k++] = horizonLine3(k3,u1i,sks);
    }
    return surfs;
  }
  */

  public float[][] horizonLine2(int k2, float u1i, FaultSkin[] sks) {
    int n3 = _u1.length;
    int n2 = _u1[0].length;
    int n1 = _u1[0][0].length;
    float d1 = (float)_s1.getDelta();
    float f1 = (float)_s1.getFirst();
    int i1 = round((u1i-f1)/d1);
    float[] xyz = new float[n3*3];
    int k = 0;
    short[][][] mk = zeroOnFaults(n1,n2,n3,sks);
    for (int i3=0; i3<n3; ++i3) {
      float x3 = i3;
      float x2 = k2;
      float x1 = _x1[i3][k2][i1];
      if(x1>=n1) {continue;}
      if(onFault(x3,k2,x1,mk)){continue;}
      xyz[k++] = x3;
      xyz[k++] = x2;
      xyz[k++] = x1;
    }
    xyz = copy(k,0,xyz);
    int nt = xyz.length;
    float[] u1s = fillfloat(u1i,nt);
    float rgtMin = min(_u1)+5;
    float rgtMax = max(_u1)-5;
    ColorMap cp = new ColorMap(rgtMin,rgtMax,ColorMap.JET);
    float[] rgb = cp.getRgbFloats(u1s);
    return new float[][]{xyz,rgb};
  }

  public float[][] horizonLine3(int k3, float u1i, FaultSkin[] sks) {
    int n3 = _u1.length;
    int n2 = _u1[0].length;
    int n1 = _u1[0][0].length;
    float d1 = (float)_s1.getDelta();
    float f1 = (float)_s1.getFirst();
    int i1 = round((u1i-f1)/d1);
    float[] xyz = new float[n2*3];
    int k = 0;
    short[][][] mk = zeroOnFaults(n1,n2,n3,sks);
    for (int i2=0; i2<n2; ++i2) {
      float x3 = k3;
      float x2 = i2;
      float x1 = _x1[k3][i2][i1];
      if(x1>=n1) {continue;}
      if(onFault(x3,x2,x1,mk)){continue;}
      xyz[k++] = x3;
      xyz[k++] = x2;
      xyz[k++] = x1;
    }
    xyz = copy(k,0,xyz);
    int nt = xyz.length;
    float[] u1s = fillfloat(u1i,nt);
    float rgtMin = min(_u1)+5;
    float rgtMax = max(_u1)-5;
    ColorMap cp = new ColorMap(rgtMin,rgtMax,ColorMap.JET);
    float[] rgb = cp.getRgbFloats(u1s);
    return new float[][]{xyz,rgb};
  }

  public float[][][] multipleHorizons(Sampling t1, FaultSkin[] sks) {
    double[] u1s = t1.getValues();
    int n1s = u1s.length;
    float[][][] surfs = new float[n1s][][];
    for (int i1s=0; i1s<n1s; ++i1s) {
      float u1i = (float) u1s[i1s];
      surfs[i1s] = singleHorizon(u1i,sks);
    }
    return surfs;
  }


  public float[][] singleHorizon(float u1i, FaultSkin[] sks) {
    int n3 = _u1.length;
    int n2 = _u1[0].length;
    int n1 = _u1[0][0].length;
    float d1 = (float)_s1.getDelta();
    float f1 = (float)_s1.getFirst();
    int i1 = round((u1i-f1)/d1);
    float[][] sfi = new float[n3][n2];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      sfi[i3][i2] = _x1[i3][i2][i1];
    }}

    short[][][] mk = zeroOnFaults(n1,n2,n3,sks);
    float[] xyz = buildTrigs(_s3,_s2,sfi,mk);
    int nt = xyz.length;
    float[] u1s = fillfloat(u1i,nt);
    float rgtMin = min(_u1)+5;
    float rgtMax = max(_u1)-5;
    ColorMap cp = new ColorMap(rgtMin,rgtMax,ColorMap.JET);
    float[] rgb = cp.getRgbFloats(u1s);
    return new float[][]{xyz,rgb};
  }

  private float[] buildTrigs(
    Sampling sx, Sampling sy, float[][] z, short[][][] mk) 
  {
    int i = 0;
    int nz = mk[0][0].length;
    int nx = sx.getCount();
    int ny = sy.getCount();
    float[][] sc = fillfloat(1.0f,ny,nx);
    for (int ix=0;ix<nx; ++ix) {
    for (int iy=0;iy<ny; ++iy) {
      int iz = round(z[ix][iy]);
      if(iz<0){iz=0;} if(iz>=nz){iz=nz-1;}
      if(mk[ix][iy][iz]==1) {
        sc[ix][iy] = 0.0f;
      }
      
    }}
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(2.0f,sc,copy(z),z);
    float[] xyz = new float[nx*ny*6*3];
    for (int ix=4;ix<nx-4; ++ix) {
      float x0 = (float)sx.getValue(ix  );
      float x1 = (float)sx.getValue(ix+1);
      for (int iy=4; iy<ny-4; ++iy) {
        float y0 = (float)sy.getValue(iy  );
        float y1 = (float)sy.getValue(iy+1);
        float z1 = z[ix  ][iy  ];
        float z2 = z[ix  ][iy+1];
        float z3 = z[ix+1][iy  ];
        float z4 = z[ix+1][iy  ];
        float z5 = z[ix  ][iy+1];
        float z6 = z[ix+1][iy+1];
        if(abs(z1-z2)>1f){continue;}
        if(abs(z1-z3)>1f){continue;}
        if(abs(z2-z3)>1f){continue;}
        if(abs(z4-z5)>1f){continue;}
        if(abs(z4-z6)>1f){continue;}
        if(abs(z5-z6)>1f){continue;}
        if(onFault(x0,y0,z1,mk)){continue;}
        if(onFault(x0,y1,z2,mk)){continue;}
        if(onFault(x1,y0,z3,mk)){continue;}
        if(onFault(x1,y0,z4,mk)){continue;}
        if(onFault(x0,y1,z5,mk)){continue;}
        if(onFault(x1,y1,z6,mk)){continue;}
        if(z1<0||z2<0||z3<0){continue;}
        if(z4<0||z5<0||z6<0){continue;}
        if(z1>nz||z2>nz||z3>nz){continue;}
        if(z4>nz||z5>nz||z6>nz){continue;}

        xyz[i++] = x0;  xyz[i++] = y0;  xyz[i++] = z[ix  ][iy  ];
        xyz[i++] = x0;  xyz[i++] = y1;  xyz[i++] = z[ix  ][iy+1];
        xyz[i++] = x1;  xyz[i++] = y0;  xyz[i++] = z[ix+1][iy  ];
        xyz[i++] = x1;  xyz[i++] = y0;  xyz[i++] = z[ix+1][iy  ];
        xyz[i++] = x0;  xyz[i++] = y1;  xyz[i++] = z[ix  ][iy+1];
        xyz[i++] = x1;  xyz[i++] = y1;  xyz[i++] = z[ix+1][iy+1];
      }
    }
    return copy(i,0,xyz);
  }


  public float[][] buildTrigs(
    Sampling sx, Sampling sy, 
    float color, float[][] z, float[][] f, short[][][] mk) 
  {
    int i = 0;
    int k = 0;
    int c = 0;
    int nx = sx.getCount();
    int ny = sy.getCount();
    int nz = mk[0][0].length;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(1.0f);
    ref.apply(z,z);
    float[] zas = new float[nx*ny*6];
    float[] zfs = new float[nx*ny*6];
    float[] xyz = new float[nx*ny*6*3];
    for (int ix=0;ix<nx-1; ++ix) {
      float x0 = (float)sx.getValue(ix  );
      float x1 = (float)sx.getValue(ix+1);
      for (int iy=0; iy<ny-1; ++iy) {
        float y0 = (float)sy.getValue(iy  );
        float y1 = (float)sy.getValue(iy+1);
        float z1 = z[ix  ][iy  ];
        float z2 = z[ix  ][iy+1];
        float z3 = z[ix+1][iy  ];
        float z4 = z[ix+1][iy  ];
        float z5 = z[ix  ][iy+1];
        float z6 = z[ix+1][iy+1];
        if(abs(z1-z2)>1f){continue;}
        if(abs(z1-z3)>1f){continue;}
        if(abs(z2-z3)>1f){continue;}
        if(abs(z4-z5)>1f){continue;}
        if(abs(z4-z6)>1f){continue;}
        if(abs(z5-z6)>1f){continue;}
        if(Float.isNaN(z1)){continue;}
        if(Float.isNaN(z2)){continue;}
        if(Float.isNaN(z3)){continue;}
        if(Float.isNaN(z4)){continue;}
        if(Float.isNaN(z5)){continue;}
        if(Float.isNaN(z6)){continue;}
        if(onFault(x0,y0,z1,mk)){continue;}
        if(onFault(x0,y1,z2,mk)){continue;}
        if(onFault(x1,y0,z3,mk)){continue;}
        if(onFault(x1,y0,z4,mk)){continue;}
        if(onFault(x0,y1,z5,mk)){continue;}
        if(onFault(x1,y1,z6,mk)){continue;}
        if(z1<0||z2<0||z3<0){continue;}
        if(z4<0||z5<0||z6<0){continue;}
        if(z1>nz||z2>nz||z3>nz){continue;}
        if(z4>nz||z5>nz||z6>nz){continue;}
        //float v1 = sincInterp(z1,y0,x0,gx);
        //float v2 = sincInterp(z2,y1,x0,gx);
        //float v3 = sincInterp(z3,y0,x1,gx);
        //float v4 = sincInterp(z4,y0,x1,gx);
        //float v5 = sincInterp(z5,y1,x0,gx);
        //float v6 = sincInterp(z6,y1,x1,gx);
        //zas[k++] = v1;  zas[k++] = v2;  zas[k++] =v3;
        //zas[k++] = v4;  zas[k++] = v5;  zas[k++] =v6;
        zas[k++] = z1;  zas[k++] = z2;  zas[k++] =z3;
        zas[k++] = z4;  zas[k++] = z5;  zas[k++] =z6;
        if(f!=null) {
          zfs[c++] = f[ix  ][iy  ];
          zfs[c++] = f[ix  ][iy+1];
          zfs[c++] = f[ix+1][iy  ];
          zfs[c++] = f[ix+1][iy  ];
          zfs[c++] = f[ix  ][iy+1];
          zfs[c++] = f[ix+1][iy+1];
        }
        xyz[i++] = x0;  xyz[i++] = y0;  xyz[i++] = z[ix  ][iy  ];
        xyz[i++] = x0;  xyz[i++] = y1;  xyz[i++] = z[ix  ][iy+1];
        xyz[i++] = x1;  xyz[i++] = y0;  xyz[i++] = z[ix+1][iy  ];
        xyz[i++] = x1;  xyz[i++] = y0;  xyz[i++] = z[ix+1][iy  ];
        xyz[i++] = x0;  xyz[i++] = y1;  xyz[i++] = z[ix  ][iy+1];
        xyz[i++] = x1;  xyz[i++] = y1;  xyz[i++] = z[ix+1][iy+1];
      }
    }
    float[] rgb;
    zas = copy(k,0,zas);
    xyz = copy(i,0,xyz);
    float zmin = Float.MAX_VALUE;
    float zmax = -Float.MAX_VALUE;
    for (int ix=0; ix<nx; ++ix) {
    for (int iy=0; iy<ny; ++iy) {
      float zi = z[ix][iy];
      if (Float.isNaN(zi)) {continue;}
      if (zi<zmin) {zmin=zi;}
      if (zi>zmax) {zmax=zi;}
    }}
    if(color>0.0f) {
      zero(zas);
      add(zas,color,zas);
      ColorMap cp = new ColorMap(0.0f,1.0f,ColorMap.JET);
      rgb = cp.getRgbFloats(zas);
    } else if (f==null) {
      ColorMap cp = new ColorMap(-zmax,-zmin,ColorMap.JET);
      rgb = cp.getRgbFloats(mul(zas,-1f));
    } else {
      ColorMap cp = new ColorMap(0.3,1.0,ColorMap.JET);
      rgb = cp.getRgbFloats(zfs);
    }
    //ColorMap cp = new ColorMap(min(zas),max(zas),ColorMap.RED_WHITE_BLUE);
    //float[] rgb = cp.getRgbFloats(zas);
    return new float[][]{xyz,rgb};
  }


  private boolean onFault(float x3, float x2, float x1, short[][][] mk) {
    int n3 = mk.length;
    int n2 = mk[0].length;
    int n1 = mk[0][0].length;
    int i1 = round(x1);
    int i2 = round(x2);
    int i3 = round(x3);
    if(i1<0){i1=0;}if(i1>=n1){i1=n1-1;}
    if(i2<0){i2=0;}if(i2>=n2){i2=n2-1;}
    if(i3<0){i3=0;}if(i3>=n3){i3=n3-1;}
    if(mk[i3][i2][i1]==1){return true;}
    else                 {return false;}
  }

  private boolean onUnc(float x3, float x2, float x1, short[][][] mk) {
    int n3 = mk.length;
    int n2 = mk[0].length;
    int n1 = mk[0][0].length;
    int i1 = round(x1);
    int i2 = round(x2);
    int i3 = round(x3);
    int m1 = i1;
    int p1 = i1;
    if(i1>x1) {m1=i1-1;}
    if(i1<x1) {p1=i1+1;}
    if(i1<0){i1=0;}if(i1>=n1){i1=n1-1;}
    if(m1<0){m1=0;}if(m1>=n1){m1=n1-1;}
    if(p1<0){p1=0;}if(p1>=n1){p1=n1-1;}
    if(i2<0){i2=0;}if(i2>=n2){i2=n2-1;}
    if(i3<0){i3=0;}if(i3>=n3){i3=n3-1;}
    //if(mk[i3][i2][i1]==1){return true;}
    if(mk[i3][i2][m1]==1){return true;}
    if(mk[i3][i2][p1]==1){return true;}
    return false;
  }


  private short[][][] zeroOnFaults(
    int n1, int n2, int n3, FaultSkin[] sks) 
  {
    short[][][] mk = new short[n3][n2][n1];
    for (FaultSkin sk:sks) {
      for (FaultCell fc:sk) {
        float s1 = fc.getS1();
        if(s1<1.0f){continue;}
        int[] is = fc.getI();
        int i1i = is[0];
        int i2i = is[1];
        int i3i = is[2];
        mk[i3i][i2i][i1i] = 1;
      }
    }
    return mk;
  }

  private short[][][] zeroUncs(
    int n1, int n2, int n3, float[][][] uncs) 
  {
    int ns = uncs.length;
    short[][][] mk = new short[n3][n2][n1];
    for (int is=0; is<ns; ++is) {
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
        int i2i = i2;
        int i3i = i3;
        int i1i = round(uncs[is][i3][i2]);
        if(i1i<0)   {i1i=0;}
        if(i1i>=n1) {i1i=n1-1;}
        mk[i3i][i2i][i1i] = 1;
    }}}
    return mk;
  }


  private float[][][] horizonVolumeFromRgt(
    final Sampling s1, final float[][][] u1) 
  {
    cleanRGT(u1);
    final int n3 = u1.length;
    final int n2 = u1[0].length;
    final int n1 = u1[0][0].length;
    final float[][][] x1 = new float[n3][n2][n1];
    final InverseInterpolator ii = new InverseInterpolator(s1,s1);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) 
        ii.invert(u1[i3][i2],x1[i3][i2]);
    }});
    return x1;
  }

  private void cleanRGT(float[][][] u) {
    int n3 = u.length;
    int n2 = u[0].length;
    int n1 = u[0][0].length;
    float dMin = 0.0001f*(float)_s1.getDelta();
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          float ui = u[i3][i2][i1  ];
          float um = u[i3][i2][i1-1];
          float ut = um+dMin;
          if (ui<ut) {u[i3][i2][i1]=ut;}
        }
      }
    }
  }

  private Sampling _s1;
  private Sampling _s2;
  private Sampling _s3;
  private float[][][] _x1;
  private float[][][] _u1;
}
