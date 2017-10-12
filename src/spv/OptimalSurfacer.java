/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

package spv;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.awt.*;
import static edu.mines.jtk.util.ArrayMath.*;

import util.*;

/**
 * Enhance fault attributes and estimates fault strikes, and dips, 
 * by optimal surface voting.
 *
 * @author Xinming Wu, Univerity of Texas at Austin.
 * @version 2017.07.16
 */
public class OptimalSurfacer {

  public OptimalSurfacer(int ru, int rv, int rw) {
    _ru = ru;
    _rv = rv;
    _rw = rw;
    _lmin = -ru;
    _lmax =  ru;
    _nl = 1+_lmax-_lmin;
    updateShiftRanges();
  }
 

    /**
   * Sets bound on fault surface slopes in 1st and 2nd dimensions.
   * @param strainMax1 bound on strain in the 1st dimension.
   * @param strainMax2 bound on strain in the 2nd dimension.
   */
  public void setStrainMax(double strainMax1, double strainMax2) {
    Check.argument(strainMax1<=1.0,"strainMax1<=1.0");
    Check.argument(strainMax2<=1.0,"strainMax2<=1.0");
    Check.argument(strainMax1>0.0,"strainMax1>0.0");
    Check.argument(strainMax2>0.0,"strainMax2>0.0");
    _bstrain1 = (int)ceil(1.0/strainMax1);
    _bstrain2 = (int)ceil(1.0/strainMax2);
  }

  /**
   * Sets the number of nonlinear smoothings of fault attributes.
   * The default number of smoothings is one.
   * @param esmooth number of nonlinear smoothings.
   */
  public void setAttributeSmoothing(int esmooth) {
    _esmooth = esmooth;
  }

    /**
   * Sets extents of smoothing filters used to smooth an extracted fault surface.
   * Half-widths of smoothing filters are inversely proportional to
   * strain limits, and are scaled by the specified factors. Default 
   * factors are zero, for no smoothing.
   * @param usmooth1 extent of smoothing filter in 1st dimension.
   * @param usmooth2 extent of smoothing filter in 2nd dimension.
   */
  public void setSurfaceSmoothing(double usmooth1, double usmooth2) {
    _usmooth1 = usmooth1;
    _usmooth2 = usmooth2;
    updateSmoothingFilters();
  }

  public float[][][] transpose21(float[][][] fx) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[][][] gx = new float[n3][n1][n2];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      gx[i3][i1][i2] = fx[i3][i2][i1];
    }}}
    return gx;
  }


  public float[][][] shear1(float theta1, float theta2, float[][][] fx) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[][][] gx = new float[n3][n2][n1];
    float pi = (float)Math.PI;
    float ft1 = pi*theta1/180f;
    float ft2 = pi*theta2/180f;
    SincInterpolator si = new SincInterpolator();
    int r3 = round(n3*0.5f);
    for (int i3=r3; i3<n3; ++i3) {
      float[] x1 = rampfloat(0f,1f,n1);
      x1 = sub(x1,(i3-r3)*tan(ft1));
    for (int i2=0; i2<n2; ++i2) {
      si.interpolate(n1,1,0,fx[i3][i2],n1,x1,gx[i3][i2]);
    }}
    int d3 = 0;
    for (int i3=r3; i3>=0; --i3) {
      float[] x1 = rampfloat(0f,1f,n1);
      x1 = add(x1,d3*tan(ft2));
      d3++;
    for (int i2=0; i2<n2; ++i2) {
      si.interpolate(n1,1,0,fx[i3][i2],n1,x1,gx[i3][i2]);
    }}

    return gx;
  }

  public float[][][] getUvwBox(
    int c1, int c2, int c3, float[] u, float[] v, float[] w, float[][][] fx) {
    int nu = _nl;
    int nv = _rv*2+1;
    int nw = _rw*2+1;
    float[][] rws = new float[3][nw];
    float[][] rvs = new float[3][nv];
    float[][] rus = new float[3][nu];
    updateVectorMap(_ru,u,rus[0],rus[1],rus[2]);
    updateVectorMap(_rv,v,rvs[0],rvs[1],rvs[2]);
    updateVectorMap(_rw,w,rws[0],rws[1],rws[2]);
    float[][][] fb = new float[nw][nv][nu];
    samplesInUvwBox(c1,c2,c3,rus,rvs,rws,fb,fx);
    return fb;
  }

  private void updateVectorMap(
    int r, float[] u, float[] ru1, float[] ru2, float[] ru3) {
    for (int i=1; i<=r; ++i) {
      int kp = r+i;
      int km = r-i;
      float iu1 = i*u[0];
      float iu2 = i*u[1];
      float iu3 = i*u[2];
      ru1[kp] =  iu1;
      ru2[kp] =  iu2;
      ru3[kp] =  iu3;
      ru1[km] = -iu1;
      ru2[km] = -iu2;
      ru3[km] = -iu3;
    }
  }



  /**
   * Extract optimal fault surface from an input fault attribute image.
   * @param fx input array for the fault attribute image.
   */
  public float[][] findSurface(int c1, int c2, int c3, float[][][] fx) {
    final int n2 = fx.length;
    final int n1 = fx[0].length;
    final int nl = fx[0][0].length;
    final float[][] u = new float[n2][n1];
    final float[][] uf = u;
    for (int is=0; is<_esmooth; ++is) {
      smoothFaultAttributes(fx,fx);
      //setControlPoint(c1,c2,c3,fx);
    }
    float[][] d = new float[n1][nl];
    for (int i2=0; i2<n2; ++i2) {
      accumulateForward(fx[i2],d);
      backtrackReverse(d,fx[i2],uf[i2]);
    }
    sub(u,_lmin,u);
    smoothSurface(u,u);
    return u;
  }

  public FaultCell[] pickSeeds(
    int c1,
    int d2, int d3, float fm, float[][][] ft, float[][][] pt, float[][][] tt) {
    final int n3 = ft.length;
    final int n2 = ft[0].length;
    final int n1 = ft[0][0].length;
    final ArrayList<FaultCell> cs = new ArrayList<FaultCell>();
    for (int i3=_rw; i3<n3-_rw; i3++) {
    for (int i2=_rv; i2<n2-_rv; i2++) {
      float fti = ft[i3][i2][c1];
      float pti = pt[i3][i2][c1];
      float tti = tt[i3][i2][c1];
      if(fti>fm) {
        FaultCell cell = new FaultCell(c1,i2,i3,fti,pti,tti);
        cs.add(cell);
      }
    }}
    int np = cs.size();
    int[] is = new int[np];
    float[] fs = new float[np];
    for (int ip=0; ip<np; ++ip) {
      is[ip] = ip;
      fs[ip] = cs.get(ip).getFl();
    }
    quickIndexSort(fs,is);
    int[][][] mark = new int[n3][n2][n1];
    ArrayList<FaultCell> seeds = new ArrayList<FaultCell>();
    for (int ip=np-1; ip>=0; --ip) {
      FaultCell cell = cs.get(is[ip]);
      int i1 = cell.getI1();
      int i2 = cell.getI2();
      int i3 = cell.getI3();
      int b2 = i2-d2; b2=max(b2,0);
      int b3 = i3-d3; b3=max(b3,0);
      int e2 = i2+d2; e2=min(e2,n2-1);
      int e3 = i3+d3; e3=min(e3,n3-1);
      boolean ok = true;
      for (int k3=b3;k3<=e3;k3++) {
      for (int k2=b2;k2<=e2;k2++) {
        if(mark[k3][k2][i1]==1) {
          ok=false;
          break;
        }
      }}
      if(ok) {
        seeds.add(cell);
        mark[i3][i2][i1] = 1;
      }
    }
    return seeds.toArray(new FaultCell[0]);
  }


  public void setControlPoint(int c1, int c2, int c3, float[][][] fb) {
    int nw = fb.length;
    int nv = fb[0].length;
    int nu = fb[0][0].length;
    for (int iu=0; iu<c1; ++iu) {
      int dm = round((c1-iu)*0.5f);
      int vb = max(c2-dm,0);
      int ve = min(c2+dm,nv-1);
      int wb = max(c3-dm,0);
      int we = min(c3+dm,nw-1);
      for (int iw=wb; iw<=we; ++iw) {
      for (int iv=vb; iv<=ve; ++iv) {
        int dw = c3-iw;
        int dv = c2-iv;
        float ri = sqrt(dw*dw+dv*dv);
        if(ri<dm) fb[iw][iv][iu] = 1f;
      }}
    }
    for (int iu=nu-1; iu>c1; --iu) {
      int dm = round((iu-c1)*0.5f);
      int vb = max(c2-dm,0);
      int ve = min(c2+dm,nv-1);
      int wb = max(c3-dm,0);
      int we = min(c3+dm,nw-1);
      for (int iw=wb; iw<=we; ++iw) {
      for (int iv=vb; iv<=ve; ++iv) {
        int dw = c3-iw;
        int dv = c2-iv;
        float ri = sqrt(dw*dw+dv*dv);
        if(ri<dm) fb[iw][iv][iu] = 1f;
      }}
    }

  }

  public float[] getSurfaceTriangles(
    float c1, float c2, float c3, 
    float[] u, float[] v, float[] w,
    float[][] z) 
  {
    int i = 0;
    int nx = z.length;
    int ny = z[0].length;
    float[] xyz = new float[nx*ny*6*3];
    for (int ix=1;ix<nx-1; ++ix) {
      float w0 = ix-_rw;
      float w1 = ix-_rw+1;
      for (int iy=1; iy<ny-1; ++iy) {
        float v0 = iy-_rv;
        float v1 = iy-_rv+1;
        float u00 = z[ix  ][iy  ]-_ru;
        float u01 = z[ix  ][iy+1]-_ru;
        float u10 = z[ix+1][iy  ]-_ru;
        float u11 = z[ix+1][iy+1]-_ru;
        float x1 = w0*w[2]+v0*v[2]+u00*u[2]+c3;
        float y1 = w0*w[1]+v0*v[1]+u00*u[1]+c2;
        float z1 = w0*w[0]+v0*v[0]+u00*u[0]+c1;

        float x2 = w0*w[2]+v1*v[2]+u01*u[2]+c3;
        float y2 = w0*w[1]+v1*v[1]+u01*u[1]+c2;
        float z2 = w0*w[0]+v1*v[0]+u01*u[0]+c1;

        float x3 = w1*w[2]+v0*v[2]+u10*u[2]+c3;
        float y3 = w1*w[1]+v0*v[1]+u10*u[1]+c2;
        float z3 = w1*w[0]+v0*v[0]+u10*u[0]+c1;

        float x4 = w1*w[2]+v1*v[2]+u11*u[2]+c3;
        float y4 = w1*w[1]+v1*v[1]+u11*u[1]+c2;
        float z4 = w1*w[0]+v1*v[0]+u11*u[0]+c1;


        xyz[i++] = x1;  xyz[i++] = y1;  xyz[i++] = z1;
        xyz[i++] = x2;  xyz[i++] = y2;  xyz[i++] = z2;
        xyz[i++] = x3;  xyz[i++] = y3;  xyz[i++] = z3;
        xyz[i++] = x3;  xyz[i++] = y3;  xyz[i++] = z3;
        xyz[i++] = x2;  xyz[i++] = y2;  xyz[i++] = z2;
        xyz[i++] = x4;  xyz[i++] = y4;  xyz[i++] = z4;
      }
    }
    return copy(i,0,xyz);
  }

  public float[][] getSurfaceTriangles(
    boolean smooth,
    float c1, float c2, float c3, 
    float fmin, float fmax,
    float[] u, float[] v, float[] w,
    float[][] z, float[][][] fx) 
  {
    int i = 0;
    int nx = z.length;
    int ny = z[0].length;
    float[] xyz = new float[nx*ny*6*3];
    float[][] zf = new float[nx][ny];
    float[] zfs = new float[nx*ny*6];
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    SincInterpolator si = new SincInterpolator();
    for (int ix=0;ix<nx; ++ix) {
      float wi = ix-_rw;
    for (int iy=0; iy<ny; ++iy) {
        float vi = iy-_rv;
        float ui = z[ix][iy]-_ru;
        float xi = wi*w[2]+vi*v[2]+ui*u[2]+c3;
        float yi = wi*w[1]+vi*v[1]+ui*u[1]+c2;
        float zi = wi*w[0]+vi*v[0]+ui*u[0]+c1;
        zf[ix][iy] = si.interpolate(s1,s2,s3,fx,zi,yi,xi);
    }}
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(6);
    if(smooth) {
      rgf.apply0X(zf,zf);
      rgf.applyX0(zf,zf);
    }

    int c = 0;
    for (int ix=1;ix<nx-1; ++ix) {
      float w0 = ix-_rw;
      float w1 = ix-_rw+1;
      for (int iy=1; iy<ny-1; ++iy) {
        float v0 = iy-_rv;
        float v1 = iy-_rv+1;
        float u00 = z[ix  ][iy  ]-_ru;
        float u01 = z[ix  ][iy+1]-_ru;
        float u10 = z[ix+1][iy  ]-_ru;
        float u11 = z[ix+1][iy+1]-_ru;
        float x1 = w0*w[2]+v0*v[2]+u00*u[2]+c3;
        float y1 = w0*w[1]+v0*v[1]+u00*u[1]+c2;
        float z1 = w0*w[0]+v0*v[0]+u00*u[0]+c1;

        float x2 = w0*w[2]+v1*v[2]+u01*u[2]+c3;
        float y2 = w0*w[1]+v1*v[1]+u01*u[1]+c2;
        float z2 = w0*w[0]+v1*v[0]+u01*u[0]+c1;

        float x3 = w1*w[2]+v0*v[2]+u10*u[2]+c3;
        float y3 = w1*w[1]+v0*v[1]+u10*u[1]+c2;
        float z3 = w1*w[0]+v0*v[0]+u10*u[0]+c1;

        float x4 = w1*w[2]+v1*v[2]+u11*u[2]+c3;
        float y4 = w1*w[1]+v1*v[1]+u11*u[1]+c2;
        float z4 = w1*w[0]+v1*v[0]+u11*u[0]+c1;


        xyz[i++] = x1;  xyz[i++] = y1;  xyz[i++] = z1;
        xyz[i++] = x2;  xyz[i++] = y2;  xyz[i++] = z2;
        xyz[i++] = x3;  xyz[i++] = y3;  xyz[i++] = z3;
        xyz[i++] = x3;  xyz[i++] = y3;  xyz[i++] = z3;
        xyz[i++] = x2;  xyz[i++] = y2;  xyz[i++] = z2;
        xyz[i++] = x4;  xyz[i++] = y4;  xyz[i++] = z4;

        zfs[c++] = zf[ix  ][iy  ];
        zfs[c++] = zf[ix  ][iy+1];
        zfs[c++] = zf[ix+1][iy  ];
        zfs[c++] = zf[ix+1][iy  ];
        zfs[c++] = zf[ix  ][iy+1];
        zfs[c++] = zf[ix+1][iy+1];

      }
    }
    xyz =  copy(i,0,xyz);
    zfs =  copy(c,0,zfs);
    ColorMap cp = new ColorMap(fmin,fmax,ColorMap.JET);
    float[] rgb = cp.getRgbFloats(zfs);
    return new float[][]{xyz,rgb};
  }



  public float[] getSurfaceTriangles(
    Sampling sx, Sampling sy, float[][] z) 
  {
    int i = 0;
    int nx = sx.getCount();
    int ny = sy.getCount();
    float[] xyz = new float[nx*ny*6*3];
    for (int ix=1;ix<nx-1; ++ix) {
      float x0 = (float)sx.getValue(ix  );
      float x1 = (float)sx.getValue(ix+1);
      for (int iy=1; iy<ny-1; ++iy) {
        float y0 = (float)sy.getValue(iy  );
        float y1 = (float)sy.getValue(iy+1);
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

  public float[] getSurfaceTrianglesT(
    Sampling sx, Sampling sy, float[][] z) 
  {
    int i = 0;
    int nx = sx.getCount();
    int ny = sy.getCount();
    float[] xyz = new float[nx*ny*6*3];
    for (int ix=1;ix<nx-1; ++ix) {
      float x0 = (float)sx.getValue(ix  );
      float x1 = (float)sx.getValue(ix+1);
      for (int iy=1; iy<ny-1; ++iy) {
        float y0 = (float)sy.getValue(iy  );
        float y1 = (float)sy.getValue(iy+1);
        xyz[i++] = x0;  
        xyz[i++] = z[ix  ][iy  ];
        xyz[i++] = y0;  

        xyz[i++] = x0;  
        xyz[i++] = z[ix  ][iy+1];
        xyz[i++] = y1;  

        xyz[i++] = x1;  
        xyz[i++] = z[ix+1][iy  ];
        xyz[i++] = y0;  

        xyz[i++] = x1;  
        xyz[i++] = z[ix+1][iy  ];
        xyz[i++] = y0;  

        xyz[i++] = x0;  
        xyz[i++] = z[ix  ][iy+1];
        xyz[i++] = y1;  

        xyz[i++] = x1;  
        xyz[i++] = z[ix+1][iy+1];
        xyz[i++] = y1;  
      }
    }
    return copy(i,0,xyz);
  }

  public float[][] planarityOnFault(float[][][] fx, float[][] sf) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[][] sa = new float[n3][n2];
    SincInterpolator si = new SincInterpolator();
    for (int i3=0; i3<n3; ++i3) { 
    for (int i2=0; i2<n2; ++i2) {
      sa[i3][i2] = si.interpolate(n1,1,0,fx[i3][i2],sf[i3][i2]);
    }}
    return sa;
  }

  public void surfaceToImage(float[][] sf, float[][] sx, float[][][] fx) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    for(int i3=0; i3<n3; ++i3) {
    for(int i2=0; i2<n2; ++i2) {
      int i1 = round(sf[i3][i2]);
      float sxi = sx[i3][i2];
      fx[i3][i2][i1] = sxi;
      if(i3==72) fx[i3][i2][i1-1] = sxi;
    }}

  }


  public float[][] getSurfaceTrianglesT(
    float fmin, float fmax,
    Sampling sx, Sampling sy, float[][] z, float[][] f) 
  {
    int i = 0;
    int c = 0;
    int nx = sx.getCount();
    int ny = sy.getCount();
    float[] xyz = new float[nx*ny*6*3];
    float[] zfs = new float[nx*ny*6];
    for (int ix=1;ix<nx-1; ++ix) {
      float x0 = (float)sx.getValue(ix  );
      float x1 = (float)sx.getValue(ix+1);
      for (int iy=1; iy<ny-1; ++iy) {
        float y0 = (float)sy.getValue(iy  );
        float y1 = (float)sy.getValue(iy+1);
        zfs[c++] = f[ix  ][iy  ];
        zfs[c++] = f[ix  ][iy+1];
        zfs[c++] = f[ix+1][iy  ];
        zfs[c++] = f[ix+1][iy  ];
        zfs[c++] = f[ix  ][iy+1];
        zfs[c++] = f[ix+1][iy+1];

        xyz[i++] = x0;  
        xyz[i++] = z[ix  ][iy  ];
        xyz[i++] = y0;  

        xyz[i++] = x0;  
        xyz[i++] = z[ix  ][iy+1];
        xyz[i++] = y1;  

        xyz[i++] = x1;  
        xyz[i++] = z[ix+1][iy  ];
        xyz[i++] = y0;  

        xyz[i++] = x1;  
        xyz[i++] = z[ix+1][iy  ];
        xyz[i++] = y0;  

        xyz[i++] = x0;  
        xyz[i++] = z[ix  ][iy+1];
        xyz[i++] = y1;  

        xyz[i++] = x1;  
        xyz[i++] = z[ix+1][iy+1];
        xyz[i++] = y1;  
      }
    }
    xyz =  copy(i,0,xyz);
    zfs =  copy(c,0,zfs);
    ColorMap cp = new ColorMap(fmin,fmax,ColorMap.JET);
    float[] rgb = cp.getRgbFloats(zfs);
    return new float[][]{xyz,rgb};
  }




 /**
   * Smooths the specified shifts. Smoothing can be performed 
   * in place; input and output arrays can be the same array.
   * @param u input array of shifts to be smoothed.
   * @param us output array of smoothed shifts.
   */
  public void smoothSurface(float[][] u, float[][] us) {
    if (_ref1!=null) {
      _ref1.apply1(u,us);
    } else {
      copy(u,us);
    }
    if (_ref2!=null)
      _ref2.apply2(us,us);
  }

    /**
   * Smooths (and normalizes) alignment errors.
   * Input and output arrays can be the same array.
   * @param e input array[n2][n1][nl] of alignment errors.
   * @param es output array[n2][n1][nl] of smoothed errors.
   */
  public void smoothFaultAttributes(float[][][] fx, float[][][] fs) {
    smoothFaultAttributes1(_bstrain1,fx,fs);
    smoothFaultAttributes2(_bstrain2,fs,fs);
    normalizeErrors(fs);
  }

  public void smoothAttribute1(float[][][] fx, float[][][] fs) {
    smoothFaultAttributes1(_bstrain1,fx,fs);
    normalizeErrors(fs);
  }

  public void smoothAttribute2(float[][][] fx, float[][][] fs) {
    smoothFaultAttributes2(_bstrain2,fx,fs);
    normalizeErrors(fs);
  }



  /**
   * Accumulates alignment errors in forward direction.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   */
  public void accumulateForward(float[][] e, float[][] d) {
    accumulate( 1,_bstrain1,e,d);
  }

  /**
   * Returns shifts found by backtracking in reverse.
   * @param d array of accumulated errors.
   * @param e array of alignment errors.
   */
  public float[] backtrackReverse(float[][] d, float[][] e) {
    float[] u = new float[d.length];
    backtrackReverse(d,e,u);
    return u;
  }

  private void samplesInUvwBox(
    int c1, int c2, int c3,
    float[][] dus, float[][] dvs, float[][] dws, 
    float[][][] fb, float[][][] fx) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    int nw = fb.length;
    int nv = fb[0].length;
    int nu = fb[0][0].length;
    SincInterpolator si = new SincInterpolator();
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    for (int kw=0; kw<nw; kw++) {
      float dw1 = dws[0][kw]+c1;
      float dw2 = dws[1][kw]+c2;
      float dw3 = dws[2][kw]+c3;
      for (int kv=0; kv<nv; kv++) {
        float dv1 = dw1+dvs[0][kv];
        float dv2 = dw2+dvs[1][kv];
        float dv3 = dw3+dvs[2][kv];
        for (int ku=0; ku<nu; ku++) {
          float x1 = dv1+dus[0][ku];
          float x2 = dv2+dus[1][ku];
          float x3 = dv3+dus[2][ku];
          fb[kw][kv][ku] = si.interpolate(s1,s2,s3,fx,x1,x2,x3);
        }
      }
    }
  }


    /**
   * Computes shifts by backtracking in reverse direction.
   * @param d input array of accumulated errors.
   * @param e input array of alignment errors.
   * @param u output array of shifts.
   */
  public void backtrackReverse(float[][] d, float[][] e, float[] u) {
    backtrack(-1,_bstrain1,_lmin,d,e,u);
  }

  /**
   * Finds shifts by backtracking in accumulated alignment errors.
   * Backtracking must be performed in the direction opposite to
   * that for which accumulation was performed.
   * @param dir backtrack direction, positive or negative.
   * @param b sample offset used to constrain changes in lag.
   * @param lmin minimum lag corresponding to lag index zero.
   * @param d input array[ni][nl] of accumulated errors.
   * @param e input array[ni][nl] of alignment errors.
   * @param u output array[ni] of computed shifts.
   */
  private static void backtrack(
    int dir, int b, int lmin, float[][] d, float[][] e, float[] u) 
  {
    float ob = 1.0f/b;
    int nl = d[0].length;
    int ni = d.length;
    int nlm1 = nl-1;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1;
    int ie = (dir>0)?nim1:0;
    int is = (dir>0)?1:-1;
    int ii = ib;
    int il = max(0,min(nlm1,-lmin));
    float dl = d[ii][il];
    for (int jl=1; jl<nl; ++jl) {
      if (d[ii][jl]<dl) {
        il = jl;
        dl = d[ii][jl];
      }
    }
    u[ii] = il+lmin;
    while (ii!=ie) {
      int ji = max(0,min(nim1,ii+is));
      int jb = max(0,min(nim1,ii+is*b));
      int ilm1 = il-1; if (ilm1==-1) ilm1 = 0;
      int ilp1 = il+1; if (ilp1==nl) ilp1 = nlm1;
      float dm = d[jb][ilm1];
      float di = d[ji][il  ];
      float dp = d[jb][ilp1];
      for (int kb=ji; kb!=jb; kb+=is) {
        dm += e[kb][ilm1];
        dp += e[kb][ilp1];
      }
      dl = min3(dm,di,dp);
      if (dl!=di) {
        if (dl==dm) {
          il = ilm1;
        } else {
          il = ilp1;
        }
      }
      ii += is;
      u[ii] = il+lmin;
      if (il==ilm1 || il==ilp1) {
        float du = (u[ii]-u[ii-is])*ob;
        u[ii] = u[ii-is]+du;
        for (int kb=ji; kb!=jb; kb+=is) {
          ii += is;
          u[ii] = u[ii-is]+du;
        }
      }
    }
  }

  /**
   * Smooths fault attributes in 1st dimension.
   * @param b strain parameter in 1st dimension.
   * @param e input array of alignment errors to be smooothed.
   * @param es output array of smoothed alignment errors.
   */
  private static void smoothFaultAttributes1(
    int b, float[][][] e, float[][][] es) {
    final int n2 = e.length;
    final int bf = b;
    final float[][][] ef = e;
    final float[][][] esf = es;
    for (int i2=0; i2<n2; ++i2)
      smoothFaultAttributes1(bf,ef[i2],esf[i2]);
  }

    /**
   * Smooths fault attributes in 1st dimension.
   * Does not normalize errors after smoothing.
   * @param b strain parameter in 1st dimension.
   * @param e input array of alignment errors to be smooothed.
   * @param es output array of smoothed alignment errors.
   */
  private static void smoothFaultAttributes1(
    int b, float[][] e, float[][] es) {
    int nl = e[0].length;
    int n1 = e.length;
    float[][] ef = new float[n1][nl];
    float[][] er = new float[n1][nl];
    accumulate( 1,b,e,ef);
    accumulate(-1,b,e,er);
    for (int i1=0; i1<n1; ++i1)
      for (int il=0; il<nl; ++il)
        es[i1][il] = ef[i1][il]+er[i1][il]-e[i1][il];
  }

 /**
   * Smooths fault attributes in 2nd dimension.
   * @param b strain parameter in 2nd dimension.
   * @param e input array of alignment errors to be smooothed.
   * @param es output array of smoothed alignment errors.
   */
  private static void smoothFaultAttributes2(
    int b, float[][][] e, float[][][] es) {
    int bf = b;
    int n2 = e.length;
    int n1 = e[0].length;
    int nl = e[0][0].length;
    float[][] e1  = new float[n2][nl];
    float[][] es1 = new float[n2][nl];
    float[][] ef1 = new float[n2][nl];
    float[][] er1 = new float[n2][nl];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
         e1[i2] =  e[i2][i1];
        es1[i2] = es[i2][i1];
      }
      accumulate( 1,bf,e1,ef1);
      accumulate(-1,bf,e1,er1);
      for (int i2=0; i2<n2; ++i2)
      for (int il=0; il<nl; ++il)
        es1[i2][il] = ef1[i2][il]+er1[i2][il]-e1[i2][il];
    }
  }

  /**
   * Non-linear accumulation of alignment errors.
   * @param dir accumulation direction, positive or negative.
   * @param b sample offset used to constrain changes in lag.
   * @param e input array[ni][nl] of alignment errors.
   * @param d output array[ni][nl] of accumulated errors.
   */
  private static void accumulate(int dir, int b, float[][] e, float[][] d) {
    int nl = e[0].length;
    int ni = e.length;
    int nlm1 = nl-1;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1;
    int ie = (dir>0)?ni:-1;
    int is = (dir>0)?1:-1;
    for (int il=0; il<nl; ++il)
      d[ib][il] = 0.0f;
    for (int ii=ib; ii!=ie; ii+=is) {
      int ji = max(0,min(nim1,ii-is));
      int jb = max(0,min(nim1,ii-is*b));
      for (int il=0; il<nl; ++il) {
        int ilm1 = il-1; if (ilm1==-1) ilm1 = 0;
        int ilp1 = il+1; if (ilp1==nl) ilp1 = nlm1;
        float dm = d[jb][ilm1];
        float di = d[ji][il  ];
        float dp = d[jb][ilp1];
        for (int kb=ji; kb!=jb; kb-=is) {
          dm += e[kb][ilm1];
          dp += e[kb][ilp1];
        }
        d[ii][il] = min3(dm,di,dp)+e[ii][il];
      }
    }
  }

  private static float min3(float a, float b, float c) {
    return b<=a?(b<=c?b:c):(a<=c?a:c); // if equal, choose b
  }

  private void updateSmoothingFilters() {
    _ref1 = (_usmooth1<=0.0) ? null :
      new RecursiveExponentialFilter(_usmooth1*_bstrain1);
    _ref2 = (_usmooth2<=0.0) ? null :
      new RecursiveExponentialFilter(_usmooth2*_bstrain2);
  }

  private void updateShiftRanges() {
    int nw = _rw*2+1;
    int nv = _rv*2+1;
    _lmins = new int[nw][nv];
    _lmaxs = new int[nw][nv];
    for (int iw=-_rw; iw<=_rw; ++iw) {
    for (int iv=-_rv; iv<=_rv; ++iv) {
      float wv = sqrt(iw*iw+iv*iv);
      if(wv>2) {
        _lmins[iw+_rw][iv+_rv] = max(-round(wv),_lmin);
        _lmaxs[iw+_rw][iv+_rv] = min( round(wv),_lmax);
      }
    }}
  }

    /**
   * Normalizes alignment errors to be in range [0,1].
   * @param e input/output array of alignment errors.
   */
  public static void normalizeErrors(float[][][] e) {
    final float[][][] ef = e;
    int n2 = e.length;
    MinMax mm = Parallel.reduce(n2,new Parallel.ReduceInt<MinMax>() {
    public MinMax compute(int i2) {
      int nl = ef[i2][0].length;
      int n1 = ef[i2].length;
      float emin =  Float.MAX_VALUE;
      float emax = -Float.MAX_VALUE;
      for (int i1=0; i1<n1; ++i1) {
        for (int il=0; il<nl; ++il) {
          float ei = ef[i2][i1][il];
          if (ei<emin) emin = ei;
          if (ei>emax) emax = ei;
        }
      }
      return new MinMax(emin,emax);
    }
    public MinMax combine(MinMax mm1, MinMax mm2) {
      return new MinMax(min(mm1.emin,mm2.emin),max(mm1.emax,mm2.emax));
    }});
    shiftAndScale(mm.emin,mm.emax,e);
  }

  private static class MinMax {
    float emin,emax;
    MinMax(float emin, float emax) {
      this.emin = emin;
      this.emax = emax;
    }
  }


    /**
   * Shifts and scales alignment errors to be in range [0,1].
   * @param emin minimum alignment error before normalizing.
   * @param emax maximum alignment error before normalizing.
   * @param e input/output array of alignment errors.
   */
  private static void shiftAndScale(float emin, float emax, float[][][] e) {
    final int n2 = e.length;
    final float eshift = emin;
    final float escale = (emax>emin)?1.0f/(emax-emin):1.0f;
    final float[][][] ef = e;
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      int nl = ef[i2][0].length;
      int n1 = ef[i2].length;
      for (int i1=0; i1<n1; ++i1) {
        for (int il=0; il<nl; ++il) {
          ef[i2][i1][il] = (ef[i2][i1][il]-eshift)*escale;
        }
      }
    }});
  }




  ///////////////////////////////////////////////////////////////////////////
  // private
  private int _nl; // number of lags
  private int _lmin,_lmax; // min,max lags
  private int[][] _lmins,_lmaxs;
  private int _ru,_rv,_rw;
  private int _esmooth = 1; // number of nonlinear smoothings of attributes
  private int _bstrain1 = 4; // inverse of bound on slope in 1st dimension
  private int _bstrain2 = 4; // inverse of bound on slope in 2nd dimension
  private double _usmooth1 = 2.0; // extent of smoothing shifts in 1st dim
  private double _usmooth2 = 2.0; // extent of smoothing shifts in 2nd dim
  private RecursiveExponentialFilter _ref1; // for smoothing shifts
  private RecursiveExponentialFilter _ref2; // for smoothing shifts

}
