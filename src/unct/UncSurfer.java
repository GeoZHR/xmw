/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package unct;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.util.ArrayMath;

import java.util.*;

/**
  *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.04.08
 */
public class UncSurfer {

  public void setForLof(double sigma1, double sigma2) {
    _lof = new LocalOrientFilter(sigma1,sigma2);
  }

  // set horizontal and vertical intervals for sparse sampling
  public void setSampling(int dh, int dv) {
    _dh = dh;
    _dv = dv;
  }

  public void setForLofu(double sigma1, double sigma2, int niter) {
    _niter = niter;
    _lofuSigma1  = sigma1;
    _lofuSigma2  = sigma2;
    //_lofu = new LocalOrientFilterUnc(_lofuSigma1,_lofuSigma2,_niter);
    _lofu = new LocalOrientFilterUM(_lofuSigma1,_lofuSigma2,_niter);
  }

  public float[][][] likelihood(float[][][] f) {
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    int n3s = reSampling(_dh,n3);
    int n2s = reSampling(_dh,n2);
    int nms = n3s;if(nms<n2s){nms=n2s;}
    int nit = nms*2+200;
    double sig1u = 0.75;
    double sig2u = (double)nms*2.0;
    float[][][] uls = likelihood(nit,sig1u,sig2u,f);
    float[][][] uli = new float[n3][n2][n1];
    interp(_dh,_dv,copy(uls),uli);
    return uli;
  }

  // compute unconformity likelihoods
  public float[][][] likelihood(
    int nit, double sig1u, double sig2u, float[][][] f) {
    int n3 = f.length; 
    int n2 = f[0].length; 
    int n1 = f[0][0].length; 
    int n3s = reSampling(_dh,n3);
    int n2s = reSampling(_dh,n2);
    int n1s = reSampling(_dv,n1);
    float[][][] u   = new float[n3s][n2s][n1s];
    float[][][] uc1 = new float[n3s][n2s][n1s];
    float[][][] uc2 = new float[n3s][n2s][n1s];
    float[][][] uc3 = new float[n3s][n2s][n1s];
    float[][][] ua1 = new float[n3s][n2s][n1s];
    float[][][] ua2 = new float[n3s][n2s][n1s];
    float[][][] ua3 = new float[n3s][n2s][n1s];
    EigenTensors3 et = _lof.applyForTensors(f,true);
    _lofu = new LocalOrientFilterUM(sig1u,sig2u,nit);
    _lofu.setSampling(_dh,_dv);
    _lofu.applyForNormal(et,f,uc1,uc2,uc3,ua1,ua2,ua3); 
    System.out.println("comupte normals finish~");
    mul(uc1,ua1,uc1);
    mul(uc2,ua2,uc2);
    mul(uc3,ua3,uc3);
    add(uc1,uc2,u);
    add(uc3,u  ,u);
    sub(u,min(u),u);
    div(u,max(u),u);
    return u;
  }

  // find all ridges of unconformity likelihoods 
  public void thin(float th,float[][] u, float[][] ut) {
    thinParallel(th,u,ut);
  }
  public void thin(float th, float[][][] u, float[][][] ut) {
    thinParallel(th,u,ut);
  }
  public void thin(float th, float[] u, float[] ut) {
    int n1 = u.length;
    for (int i1=1; i1<n1-1; ++i1) {
      float ui = u[i1];
      float dum = ui-u[i1-1]; 
      float dup = ui-u[i1+1]; 
      float sum1 = dum+dup;
      float sum2 = abs(dum)+abs(dup);
      if (sum1==sum2 && ui>th) {
        ut[i1] = ui;
      }
    }
  }

  // connect ajacent points to form unconformity surfaces for 3D images
  // x thinned unconformity likelihood with only 1 sample vertically on the ridges
  // y unconformity likelihoods
  // output: unconformity surfaces saved in surf[ns][n3][n2]
  // ns: number of surfaces 
  // surf[is]: depth of samples on a surface
  public float[][][] surfer(
    int n2d, int n3d,float th, int pMin, float[][][] x, float[][][] y) 
  {
    float[][][] u = copy(x);
    int n3 = u.length;
    int n2 = u[0].length;
    int n1 = u[0][0].length;
    int is = 0;
    int[][] mark = new int[n3][n2];
    float[][][] surf = new float[100][n3][n2];
    ArrayMath.fill(Float.NaN,surf);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=1; i1<n1-1; ++i1) {
          float ui = u[i3][i2][i1];  
          int ip = 0;
          if(ui>0.0f && ui<th) {
            u[i3][i2][i1] = 0.0f;
            mark[i3][i2] = i1;
            double ym = y[i3][i2][i1-1];
            double yi = y[i3][i2][i1  ];
            double yp = y[i3][i2][i1+1];
            surf[is][i3][i2] = findPeak2(i1,_dv,ym,yi,yp);
            ip ++;
            while (sum(mark)>0) {
              int[] ind = findMark(mark);
              int i3t = ind[0]; 
              int i2t = ind[1]; 
              int i1t = ind[2]; 
              mark[i3t][i2t] = 0;
              ip = floodFill(ip,i1t,i2t,i3t,2,2,2,th,u,y,mark,surf[is]);
            }
            if(ip>=pMin) { 
              is = is + 1;
              System.out.println("ip="+ip);
              System.out.println("is="+is);
            }else{
              ArrayMath.fill(Float.NaN,surf[is]);
            }
          }
        }
      }
    }
    float[][][] sfs = copy(n2,n3,is,0,0,0,surf);
    int ns = sfs.length;
    float[][][] rsfs = fillfloat(Float.NaN,n2d,n3d,ns);
    refineSample(_dh,sfs,rsfs);
    return rsfs;
  }

  public void surfaceUpdate(double sigma1, double sigma2, 
    float[][][] f, float[][][] sf) 
  {
    int ns = sf.length;
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    float[][][] p2 = new float[n3][n2][n1];
    float[][][] p3 = new float[n3][n2][n1];
    float[][][] ep = new float[n3][n2][n1];
    System.out.println("sigma2="+sigma2);
    LocalSlopeFinder lsf = new LocalSlopeFinder(sigma1,sigma2,20.0);
    lsf.findSlopes(f,p2,p3,ep);
    ep = pow(ep,12.0f);
    replaceBounds(ep);
    replaceBounds(p2);
    replaceBounds(p3);
    SurfaceExtractor se = new SurfaceExtractor();
    se.setCG(0.01f,200);
    se.setWeights(1.0f);
    se.setSmoothings(4.0f,4.0f);
    for (int is=0; is<ns; ++is) {
      boolean NaNs = checkNaNs(sf[is]);
      for (int iter=0; iter<10; iter++) {
        float ad = se.surfaceUpdateFromSlopes(ep,p2,p3,sf[is],NaNs);
        System.out.println(" Average adjustments per sample = "+ad);
        if(ad<0.02f) {break;}
      }
    }
  }

  private boolean checkNaNs(float[][] x) {
    int n2 = x.length;
    int n1 = x[0].length;
    for(int i2=0; i2<n2; ++i2) {
      for(int i1=0; i1<n1; ++i1) {
        float xi = x[i2][i1];
        if (Float.isNaN(xi)) {return true;}
      }
    }
    return false;
  }

  private void replaceBounds(float[][][] x) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    for(int i3=0; i3<n3; ++i3) {
      for(int i1=0; i1<n1; ++i1) {
        x[i3][0   ][i1] = x[i3][2   ][i1];
        x[i3][1   ][i1] = x[i3][2   ][i1];
        x[i3][n2-1][i1] = x[i3][n2-3][i1];
        x[i3][n2-2][i1] = x[i3][n2-3][i1];
      }
    }
    for(int i2=0; i2<n2; ++i2) {
      for(int i1=0; i1<n1; ++i1) {
        x[0   ][i2][i1] = x[2   ][i2][i1];
        x[1   ][i2][i1] = x[2   ][i2][i1];
        x[n3-1][i2][i1] = x[n3-3][i2][i1];
        x[n3-2][i2][i1] = x[n3-3][i2][i1];
      }
    }
  }

  private void refineSample(int d, float[][][] sfs, float[][][] rsfs) {
    int ns  = sfs.length;
    int n3  =  sfs[0].length;
    int n2  =  sfs[0][0].length;
    int n3r =  rsfs[0].length;
    int n2r =  rsfs[0][0].length;
//    Sampling s2 = new Sampling(n2r,1.0,0.0);
//    Sampling s3 = new Sampling(n3r,1.0,0.0);
    for (int is=0; is<ns; ++is) {
      boolean[][] mkp = boolArray(n2r,n3r,false);
      boolean[][] mki = boolArray(n2r,n3r,false);
      ArrayList<Float> x2l = new ArrayList<Float>();
      ArrayList<Float> x3l = new ArrayList<Float>();
      ArrayList<Float> fvl = new ArrayList<Float>();
      float[][] sfi = copy(sfs[is]);
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          int i2d = i2*d;
          int i3d = i3*d;
          float fvi = sfi[i3][i2];
          if(Float.isNaN(fvi)) {
            continue;
          }else{
            fvl.add(fvi);
            x2l.add((float)i2);
            x3l.add((float)i3);
            mkp[i3d][i2d] = true;
            rsfs[is][i3d][i2d] = fvi;
          }
        }
      }
      int n = x2l.size();
      float[] x2 = new float[n];getArray(x2l,x2);
      float[] x3 = new float[n];getArray(x3l,x3);
      float[] fv = new float[n];getArray(fvl,fv);
      if(n<n2*n3){fillHoles(x2,x3,fv,sfi);}
      for (int i3r=0; i3r<n3r; ++i3r) {
        for (int i2r=0; i2r<n2r; ++i2r) {
          if(mkp[i3r][i2r]) {
            continue;
          }else{
            mki[i3r][i2r]=markToInterp(d,i2r,i3r,n2r,n3r,mkp);
          }
        }
      }
      SincInterp sfsi =  new SincInterp();
      sfsi.setExtrapolation(SincInterp.Extrapolation.CONSTANT);
      for (int i3r=0; i3r<n3r; i3r++) {
        for (int i2r=0; i2r<n2r; i2r++) {
          if(mki[i3r][i2r]) {
            float x2i = (float)i2r;
            float x3i = (float)i3r;
            rsfs[is][i3r][i2r] = sfsi.interpolate(n2,d,0.0,n3,d,0.0,sfi,x2i,x3i);
          }
        }
      }
    }
  }

  private void fillHoles(float[] x2, float[] x3, float[] fv, float[][] sf) {
    int n3 = sf.length;
    int n2 = sf[0].length;
    Sampling s2 = new Sampling(n2,1.0f,0.0f);
    Sampling s3 = new Sampling(n3,1.0f,0.0f);
    NearestGridder2 ng = new NearestGridder2(fv,x2,x3);
    float[][] surf = ng.grid(s2,s3);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        sf[i3][i2] = surf[i3][i2];
      }
    }
  }

  private boolean[][] boolArray(int n2, int n3, boolean bool) {
    boolean[][] ba = new boolean[n3][n2];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        ba[i3][i2] = bool;
      }
    }
    return ba;
  }

  private void getArray(ArrayList<Float> list, float[] a) {
    int n = list.size();
    for (int i=0; i<n; ++i) {
      a[i] = list.get(i);
    }
  }
  private boolean markToInterp(int d, int i2, int i3, int n2, int n3, boolean[][] mk) {
    int i3b = i3-d+1;if(i3b<0)   {i3b=0;}
    int i3e = i3+d-1;if(i3e>n3-1){i3e=n3-1;}
    int i2b = i2-d+1;if(i2b<0)   {i2b=0;}
    int i2e = i2+d-1;if(i2e>n2-1){i2e=n2-1;}
    int sum = 0;
    int cot1 = 0;
    for (int i3i=i3b; i3i<=i3e; ++i3i) {
      if(mk[i3i][i2]) {cot1++;}
      if(i2>0 && i2<n2-1) {
        if(i3i==0||i3i==n3-1)
          cot1++;
      }
    }
    int cot2 = 0;
    for (int i2i=i2b; i2i<=i2e; ++i2i) {
      if(mk[i3][i2i]) {cot2++;}
      if(i3>0 && i3<n3-1) {
        if(i2i==0||i2i==n2-1)
          cot2++;
      }
    }
    int cot3 = 0;
    for (int i3i=i3b; i3i<=i3e; ++i3i) {
      for (int i2i=i2b; i2i<=i2e; ++i2i) {
        if(mk[i3i][i2i])     {cot3++;}
        if(i2i==0||i2i==n2-1){cot3++;}
        if(i3i==0||i3i==n3-1){cot3++;}
      }
    }
    if(cot1>=2 && i2>0 && i2<n2-1) {sum ++;}
    if(cot2>=2 && i3>0 && i3<n3-1) {sum ++;}
    if(cot3>=4                   ) {sum ++;}
    if(sum>0) {return true;}
    else      {return false;}
  }

  private void thinParallel(
     final float th, final float[][] u, final float[][] ut)
  {
    final int n2 = u.length;
    Parallel.loop(n2,new Parallel.LoopInt() {
      public void compute(int i2) {
        thin(th,u[i2],ut[i2]);
      }
    });
  }

  private void thinParallel(
    final float th,final float[][][] u,final float[][][] ut)
  {
    final int n3 = u.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        thin(th,u[i3],ut[i3]);
      }
    });
  }
    
  private int floodFill(
    int ip, int i1t, int i2t, int i3t, 
    int d1, int d2, int d3, float th, 
    float[][][] u, float[][][] y, int[][] mark, float[][] surf) 
  {
    int n3 = u.length;
    int n2 = u[0].length;
    int n1 = u[0][0].length;
    int i1b = i1t-d1; if(i1b<0) i1b=0;
    int i2b = i2t-d2; if(i2b<0) i2b=0;
    int i3b = i3t-d3; if(i3b<0) i3b=0;
    int i1e = i1t+d1; if(i1e>n1-1) i1e=n1-1;
    int i2e = i2t+d2; if(i2e>n2-1) i2e=n2-1;
    int i3e = i3t+d3; if(i3e>n3-1) i3e=n3-1;
    float d = 4.f;
    for (int i3=i3b; i3<=i3e; ++i3) { 
      for (int i2=i2b; i2<=i2e; ++i2) { 
        for (int i1=i1b; i1<=i1e; ++i1) { 
          float ui = u[i3][i2][i1];
          if(ui>0.0f && ui<th) {
            float d1i = i1-i1t;
            float d2i = i2-i2t;
            float d3i = i3-i3t;
            float di  = d1i*d1i+d2i*d2i+d3i*d3i;
            if(di<=d) {
              double ym = y[i3][i2][i1-1];
              double yi = y[i3][i2][i1  ];
              double yp = y[i3][i2][i1+1];
              surf[i3][i2] = findPeak2(i1,_dv,ym,yi,yp);
              mark[i3][i2] = i1;
              u[i3][i2][i1] = 0.0f;
              ip++;
            }
          }
        }
      }
    }
    return ip;
  }
  private int[] findMark(int[][]mark) {
    int n3 = mark.length;
    int n2 = mark[0].length;
    int[] ind = new int[3];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        int i1 = mark[i3][i2];
        if (i1>0) {
          ind[0] = i3;
          ind[1] = i2;
          ind[2] = i1;
          break;
        }
      }
    }
    return ind;
  }
  /*
  private float findPeak(int i1, double u1, double u2, double u3) {
    System.out.println("test");
    float z = (float) i1;
    double a = u1-u3;
    double b = 2.0*(u1+u3)-4.0*u2;
    double d = a/b;
    return (z+(float)d);
  }
  */

  private float findPeak2(int i1, double di, double xm, double xi, double xp) {
    double z = (double) i1*di;
    double a = (xm-xp)*di;
    double b = 2.0*(xp+xm)-4.0*xi;
    double d = a/b;
    return (float)(z+d);
  }

  private int reSampling(int d, int n) {
    int ns = 0; 
    for (int i=0; i<n; i+=d) {
      ns += 1;
    }
    return ns;
  }

  private static void interp(int dh, int dv, float[][][] x, float[][][] y) {
    int n3x = x.length;
    int n3y = y.length;
    int n2x = x[0].length;
    int n2y = y[0].length;
    int n1x = x[0][0].length;
    int n1y = y[0][0].length;
    SincInterp xi = new SincInterp();
    //xi.setUniform(n1x,dv,0.0,n2x,dh,0.0,n3x,dh,0.0,x);
    xi.setExtrapolation(SincInterp.Extrapolation.CONSTANT);
    for (int i3=0; i3<n3y; ++i3) {
      for (int i2=0; i2<n2y; ++i2) {
        for (int i1=0; i1<n1y; ++i1) {
          float x1 = (float) i1;
          float x2 = (float) i2;
          float x3 = (float) i3;
          y[i3][i2][i1] = 
            1.0f-xi.interpolate(n1x,dv,0.0,n2x,dh,0.0,n3x,dh,0.0,x,x1,x2,x3);
        }
      }
    }
  }



  private int _dh = 1;
  private int _dv = 1;
  private double _lofuSigma1 = 0.95;
  private double _lofuSigma2 = 128.0;
  private int _niter = 400;
  private LocalOrientFilterUM _lofu = 
    new LocalOrientFilterUM(_lofuSigma1,_lofuSigma2,_niter);
  private LocalOrientFilter _lof = new LocalOrientFilter(1.0,4.0);

}
