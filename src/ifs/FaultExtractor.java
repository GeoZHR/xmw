/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ifs;

import vec.*;
import util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Regrid fault cells from known cells.
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.01.21
 */

public class FaultExtractor {

  public FaultExtractor(int n1, int n2, int n3, FaultCell[] fcs) {
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
    _fcs = fcs;
    _kt = setKdTree();
  }

  public float[][] faultExtract(FaultCell fc) {
    float i1 = fc.i1;
    float i2 = fc.i2;
    float i3 = fc.i3;
    float w1 = fc.w1;
    float w2 = fc.w2;
    float w3 = fc.w3;
    float fp = fc.fp;
    if(abs(w2)>abs(w3)) {
      float[] ws = new float[]{w2,w1,w3};
      float[] xs = new float[]{i2,i1,i3};
      float[][] sf = new float[_n3][_n1];
      faultInitial(xs,ws,sf);
      float[] k2 = new float[]{i1};
      float[] k3 = new float[]{i3};
      surfaceUpdateFromSlopes(_n2,fp,k2,k3,sf);
      return sf;
    } else {
      float[] ws = new float[]{w3,w1,w2};
      float[] xs = new float[]{i3,i1,i2};
      float[][] sf = new float[_n2][_n1];
      faultInitial(xs,ws,sf);
      float[] k2 = new float[]{i1};
      float[] k3 = new float[]{i2};
      surfaceUpdateFromSlopes(_n3,fp,k2,k3,sf);
      return sf;
    }
  }

  public float[][] faultInitial(float[] xs, float[] ws, float[][] sf) {
    float x1 = xs[0];
    float x2 = xs[1];
    float x3 = xs[2];
    float w1 = ws[0];
    float w2 = ws[1];
    float w3 = ws[2];
    if(w1<0f){w1=-w1;w2=-w2;w3=-w3;}
    int n3 = sf.length;
    int n2 = sf[0].length;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        float x2s = w2*(i2-x2);
        float x3s = w3*(i3-x3);
        float x1s = x1-(x2s+x3s)/w1; 
        if(x1s>=_n1-1) {sf[i3][i2]=_n1-1;}
        else if(x1s<0f){sf[i3][i2]=0f;}
        else           {sf[i3][i2]=x1s;}
      }
    }
    return sf;
  }

  public float[][][][] faultSlopeField() {
    float[][][] p2 = new float[_n3][_n2][_n1];
    float[][][] p3 = new float[_n3][_n2][_n1];
    float[][][] ep = new float[_n3][_n2][_n1];
    float[][][] et = new float[_n3][_n2][_n1];
    for (FaultCell fc:_fcs) {
      int i1 = fc.i1;
      int i2 = fc.i2;
      int i3 = fc.i3;
      ep[i3][i2][i1] = fc.fl;
    }
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(3.0);
    ref.apply(ep,ep);
    div(ep,max(ep),ep);
    LocalSlopeFinder lsf = new LocalSlopeFinder(8.0f,2.0f,10f);
    lsf.findSlopes(ep,p2,p3,et);
    return new float[][][][]{p2,p3,ep};
  }

  public float[][][][] faultSlopeField2(float sig1,float sig2,float sig3) {
    float[][][] p2 = new float[_n3][_n1][_n2];
    float[][][] p3 = new float[_n3][_n1][_n2];
    float[][][] ep = new float[_n3][_n1][_n2];
    float[][][] et = new float[_n3][_n1][_n2];
    for (FaultCell fc:_fcs) {
      int i1 = fc.i1;
      int i2 = fc.i2;
      int i3 = fc.i3;
      float w2 = fc.w2;
      float w3 = fc.w3;
      if(abs(w2)<abs(w3)){continue;}
      ep[i3][i1][i2] = fc.fl;
    }
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sig1);
    ref.apply(ep,ep);
    div(ep,max(ep),ep);
    LocalSlopeFinder lsf = new LocalSlopeFinder(2.0,8.0,2.0,1f);
    lsf.findSlopes(ep,p2,p3,et);
    return new float[][][][]{ep,p2,p3};
  }

  public float[][][][] faultSlopeField3(float sig1,float sig2,float sig3) {
    float[][][] p2 = new float[_n2][_n1][_n3];
    float[][][] p3 = new float[_n2][_n1][_n3];
    float[][][] ep = new float[_n2][_n1][_n3];
    float[][][] et = new float[_n2][_n1][_n3];
    for (FaultCell fc:_fcs) {
      int i1 = fc.i1;
      int i2 = fc.i2;
      int i3 = fc.i3;
      float w2 = fc.w2;
      float w3 = fc.w3;
      if(abs(w3)<abs(w2)){continue;}
      ep[i2][i1][i3] = fc.fl;
    }
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sig1);
    ref.apply(ep,ep);
    div(ep,max(ep),ep);
    LocalSlopeFinder lspf = new LocalSlopeFinder(2.0,8.0,2.0,1f);
    lspf.findSlopes(ep,p2,p3,et);
    return new float[][][][]{ep,p2,p3};
  }


  public float[] makeVertices2(Sampling sx, Sampling sy, float[][] z) {
    int nx = sx.getCount()-1;
    int ny = sy.getCount()-1;
    float[] xyz = new float[3*6*nx*ny];
    for (int ix=0,i=0; ix<nx; ++ix) {
      float x0 = (float)sx.getValue(ix  );
      float x1 = (float)sx.getValue(ix+1);
      for (int iy=0; iy<ny; ++iy) {
        float y0 = (float)sy.getValue(iy  );
        float y1 = (float)sy.getValue(iy+1);
        xyz[i++] = x0;  xyz[i++] = z[ix  ][iy  ]; xyz[i++] = y0;
        xyz[i++] = x0;  xyz[i++] = z[ix  ][iy+1]; xyz[i++] = y1;  
        xyz[i++] = x1;  xyz[i++] = z[ix+1][iy  ]; xyz[i++] = y0; 
        xyz[i++] = x1;  xyz[i++] = z[ix+1][iy  ]; xyz[i++] = y0; 
        xyz[i++] = x0;  xyz[i++] = z[ix  ][iy+1]; xyz[i++] = y1; 
        xyz[i++] = x1;  xyz[i++] = z[ix+1][iy+1]; xyz[i++] = y1; 
      }
    }
    return xyz;
  }

  public float[] makeVertices3(Sampling sx, Sampling sy, float[][] z) {
    int nx = sx.getCount()-1;
    int ny = sy.getCount()-1;
    float[] xyz = new float[3*6*nx*ny];
    for (int ix=0,i=0; ix<nx; ++ix) {
      float x0 = (float)sx.getValue(ix  );
      float x1 = (float)sx.getValue(ix+1);
      for (int iy=0; iy<ny; ++iy) {
        float y0 = (float)sy.getValue(iy  );
        float y1 = (float)sy.getValue(iy+1);
        xyz[i++] = z[ix  ][iy  ]; xyz[i++] = x0;  xyz[i++] = y0;
        xyz[i++] = z[ix  ][iy+1]; xyz[i++] = x0;  xyz[i++] = y1;  
        xyz[i++] = z[ix+1][iy  ]; xyz[i++] = x1;  xyz[i++] = y0; 
        xyz[i++] = z[ix+1][iy  ]; xyz[i++] = x1;  xyz[i++] = y0; 
        xyz[i++] = z[ix  ][iy+1]; xyz[i++] = x0;  xyz[i++] = y1; 
        xyz[i++] = z[ix+1][iy+1]; xyz[i++] = x1;  xyz[i++] = y1; 
      }
    }
    return xyz;
  }

  public float[] makeColors(float[][] r, float[][] g, float[][] b) {

    int nx = r.length-1;
    int ny = r[0].length-1;
    float[] rgb = new float[3*6*nx*ny];
    for (int ix=0,i=0; ix<nx; ++ix) {
      for (int iy=0; iy<ny; ++iy) {
        rgb[i++] = r[ix  ][iy  ];
        rgb[i++] = g[ix  ][iy  ];
        rgb[i++] = b[ix  ][iy  ];
        rgb[i++] = r[ix  ][iy+1];
        rgb[i++] = g[ix  ][iy+1];
        rgb[i++] = b[ix  ][iy+1];
        rgb[i++] = r[ix+1][iy  ];
        rgb[i++] = g[ix+1][iy  ];
        rgb[i++] = b[ix+1][iy  ];
        rgb[i++] = r[ix+1][iy  ];
        rgb[i++] = g[ix+1][iy  ];
        rgb[i++] = b[ix+1][iy  ];
        rgb[i++] = r[ix  ][iy+1];
        rgb[i++] = g[ix  ][iy+1];
        rgb[i++] = b[ix  ][iy+1];
        rgb[i++] = r[ix+1][iy+1];
        rgb[i++] = g[ix+1][iy+1];
        rgb[i++] = b[ix+1][iy+1];
      }
    }
    return rgb;
  }

  public void surfaceUpdateFromSlopes
    (int n1, float fp, float[] k2, float[] k3,float[][] surf)
  {	
    int n3 = surf.length; 
    int n2 = surf[0].length; 
    float lmt = (float)n1-1.f;
    float[][] surft = copy(surf);
    float[][] b   = new float[n3][n2]; 
    float[][] pi1 = new float[n3][n2]; 
    float[][] qi1 = new float[n3][n2]; 
    float[][] wi1 = new float[n3][n2]; 
    checkControlPoints(k2, k3, surf); 
    int niter = 50;
    for (int n=1; n<=_exniter; n++){
      System.out.println(" Iteration "+n+"......");
      VecArrayFloat2 vb    = new VecArrayFloat2(b);
      VecArrayFloat2 vsurf = new VecArrayFloat2(surf);
      updateSlopesAndWeights(fp,surf,pi1,qi1,wi1);
      A2 a2 = new A2(wi1,_weight);
      M2 m2 = new M2(_sigma1,_sigma2,wi1,k2,k3);
      if(n>5) {niter=_niter;}
      CgSolver cs = new CgSolver(_small,niter);
      vb.zero();
      makeRhs(wi1,pi1,qi1,b);
      cs.solve(a2,m2,vb,vsurf);
      surf = vsurf.getArray();
      surfaceCorrect(surf,lmt);
      float ad = sum(abs(sub(surft,surf)))/(n3*n2); 
      System.out.println(" Average adjustments per sample = "+ad);
      if (ad<0.02f) break;
      surft = copy(surf);
    }
  
  }

  private void updateSlopesAndWeights(final float fp,
    final float[][] sf, final float[][] p2, 
    final float[][] p3, final float[][] wp) 
  {
    final int dx = 10;
    final int n3 = sf.length;
    final int n2 = sf[0].length;
    Parallel.loop(n3,new Parallel.LoopInt(){
    public void compute(int i3) {
      float[] xmin = new float[3];
      float[] xmax = new float[3];
      for (int i2=0; i2<n2; ++i2) {
        float x1i,x2i,x3i;
        if(_n3==n3) {
          x1i = i2;
          x3i = i3;
          x2i = sf[i3][i2];
        } else {
          x1i = i2;
          x2i = i3;
          x3i = sf[i3][i2];
        }
        //int ip = _kt.findNearest(new float[]{x1i,x2i,x3i});
        //float fpt = _fcs[ip].fp;
        xmin[0] = x1i-dx; xmax[0] = x1i+dx;
        xmin[1] = x2i-dx; xmax[1] = x2i+dx;
        xmin[2] = x3i-dx; xmax[2] = x3i+dx;
        int[] id = _kt.findInRange(xmin,xmax);
        int nd = id.length;
        if(nd<1) {
          p2[i3][i2] = 0f;
          p3[i3][i2] = 0f;
          wp[i3][i2] = 0f;
        } else {
          float g11 = 0;
          float g12 = 0;
          float g13 = 0;
          float g22 = 0;
          float g23 = 0;
          float g33 = 0;
          float wss = 0;
          float fls = 0;
          for (int ik=0; ik<nd; ++ik) {
            int ic = id[ik];
            FaultCell fc = _fcs[ic];
            float dp = fc.fp-fp;
            dp = min(abs(dp),abs(dp+360),abs(dp-360));
            if(dp>15f){continue;}
            float fli = fc.fl;
            float w11 = fc.w11;
            float w22 = fc.w22;
            float w33 = fc.w33;
            float w12 = fc.w12;
            float w23 = fc.w23;
            float w13 = fc.w13;
            float dx1 = x1i-fc.x1;
            float dx2 = x2i-fc.x2;
            float dx3 = x3i-fc.x3;
            float dxs = dx1*dx1+dx2*dx2+dx3*dx3;
            float wpi = gaussian(dxs,dx/4);
            wss += wpi;
            fls += wpi*fli;
            g11 += wpi*w11;
            g22 += wpi*w22;
            g33 += wpi*w33;
            g12 += wpi*w12;
            g13 += wpi*w13;
            g23 += wpi*w23;
          }
          if(wss==0){
            p2[i3][i2] = 0f;
            p3[i3][i2] = 0f;
            wp[i3][i2] = 0f;
          } else {
            float[] us = solveEigenproblem(g11,g12,g13,g22,g23,g33);
            float u1 = us[0];
            float u2 = us[1];
            float u3 = us[2];
            wp[i3][i2] = fls/wss;
            if(n3==_n3) {
              if(u2<0f){u1=-u1;u2=-u2;u3=-u3;}
              p2[i3][i2] = -u1/u2;
              p3[i3][i2] = -u3/u2;
            } else {
              if(u3<0f){u1=-u1;u2=-u2;u3=-u3;}
              p2[i3][i2] = -u1/u3;
              p3[i3][i2] = -u2/u3;
            }
          }
        }
      }
    }});
  }

  private float[] solveEigenproblem(
    float g11, float g12, float g13,
    float g22, float g23, float g33)
  {
    double[] e = new double[3];
    double[][] a = new double[3][3];
    double[][] z = new double[3][3];
    a[0][0] = g11;
    a[0][1] = g12;
    a[0][2] = g13;
    a[1][0] = g12;
    a[1][1] = g22;
    a[1][2] = g23;
    a[2][0] = g13;
    a[2][1] = g23;
    a[2][2] = g33;
    Eigen.solveSymmetric33(a,z,e);
    float u1i = (float)z[0][0];
    float u2i = (float)z[0][1];
    float u3i = (float)z[0][2];
    return new float[]{u1i,u2i,u3i};
  }


  private float gaussian(float xs, float sig) {
    return exp(-xs/(sig*sig));
  }

  private KdTree setKdTree() {
    int nc = _fcs.length;
    float[][] xc = new float[3][nc];
    for (int ic=0; ic<nc; ++ic) {
      FaultCell fc = _fcs[ic];
      xc[0][ic] = fc.x1;
      xc[1][ic] = fc.x2;
      xc[2][ic] = fc.x3;
    }
    return new KdTree(xc);
  }


  private void surfaceCorrect(float[][] surf, float lmt) {
    int n1 = surf[0].length;
    int n2 = surf.length;
    for(int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        if (surf[i2][i1]<0.f) surf[i2][i1]=0.f;
        if (surf[i2][i1]>lmt) surf[i2][i1]=lmt;
      }
    }
  }

  private static class A2 implements CgSolver.A{
    A2(float[][] wp, float w){
      _w  = w;
      _wp = wp;
    }  
    public void apply(Vec vx, Vec vy){
      VecArrayFloat2 v2x = (VecArrayFloat2) vx;
      VecArrayFloat2 v2y = (VecArrayFloat2) vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      int n1 = y[0].length; int n2 = y.length;
      float[][] yy = new float[n2][n1];
      float[][] yt = new float[n2][n1];
      VecArrayFloat2 v2yy = new VecArrayFloat2(yy);
      v2y.zero();
      v2yy.zero();
      applyLhs(_wp,x,y);
      if (_w>0.0f) {
        applyLhs(_wp,x,yt);
        applyLhs(_wp,yt,yy);
        v2y.add(1.f,v2yy,_w);
      }
    }
    private float[][] _wp;
    private float _w;
  }

   // Preconditioner; includes smoothers and (optional) constraints.
  private static class M2 implements CgSolver.A {
    M2(float sigma1, float sigma2, float[][] wp, float[] k2, float[] k3) {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _wp = wp;
      if (k2!=null && k3!=null) {
        _k2 = copy(k2);
        _k3 = copy(k3);
      }
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      copy(x,y);
      constrain(_k2,_k3,y);
      smooth2(_sigma2,_wp,y);
      smooth1(2.f*_sigma1,_wp,y);
      smooth2(_sigma2,_wp,y);
      constrain(_k2,_k3,y);
    }
    private float _sigma1,_sigma2;
    private float[][] _wp;
    private float[] _k2,_k3;
  }

  private static void checkControlPoints(float[] k2, float[] k3, float[][] f) {
    if (k2!=null && k3!=null) {
      int np = k2.length;
      for (int ip=0; ip<np; ++ip) {
        int i2 = (int)k2[ip];
        int i3 = (int)k3[ip];
        System.out.println(" i2="+i2+" i3="+i3+" f1="+f[i3][i2]);
      }
    }
  }

  private static void constrain(float[] k2, float[] k3, float[][] x) {
    if (k2!=null && k3!=null) {
      int np = k2.length;
      for (int ip=0; ip<np; ++ip) {
        int i2 = (int)k2[ip]; 
        int i3 = (int)k3[ip]; 
        x[i3][i2] = 0.f;
      }
    }
  }

  // Smoothing for dimension 1
  private static void smooth1(float sigma, float[][] s, float[][] x){
    if (sigma<1.0f)
      return;
    float c = 0.5f*sigma*sigma;
    int n1 = x[0].length;
    int n2 = x.length;
    float[] xt = zerofloat(n1);
    float[] yt = zerofloat(n1);
    float[] st = zerofloat(n1);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        xt[i1] = x[i2][i1];
        st[i1] = s[i2][i1];
      }
      lsf.apply(c,st,xt,yt);
      for (int i1=0; i1<n1; ++i1)
        x[i2][i1] = yt[i1];
    }
  }

  // Smoothing for dimension 2
  private static void smooth2(float sigma, float[][] s, float[][] x){
    if (sigma<1.0f)
      return;
    float c = 0.5f*sigma*sigma;
    int n1 = x[0].length;
    int n2 = x.length;
    float[] xt = zerofloat(n2);
    float[] yt = zerofloat(n2);
    float[] st = zerofloat(n2);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        xt[i2] = x[i2][i1];
        st[i2] = s[i2][i1];
      }
      lsf.apply(c,st,xt,yt);
      for (int i2=0; i2<n2; ++i2)
        x[i2][i1] = yt[i2];
    }
  }

  private static void applyLhs( float[][] wp, float[][] x, float[][] y) {
    zero(y);
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float wpi = (wp!=null)?wp[i2][i1]:1.000f;
        if(wpi<0.05f) {wpi=0.05f;}
        float wps = wpi*wpi;
        float d11 = wps;
        float d22 = wps;
        float xa = 0.0f;
        float xb = 0.0f;
        xa += x[i2  ][i1  ];
        xb -= x[i2  ][i1-1];
        xb += x[i2-1][i1  ];
        xa -= x[i2-1][i1-1];
        float x1 = 0.5f*(xa+xb);
        float x2 = 0.5f*(xa-xb);
        float y1 = d11*x1;
        float y2 = d22*x2;
        float ya = 0.5f*(y1+y2);
        float yb = 0.5f*(y1-y2);
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }
  
  //private static void makeRhs
  private static void makeRhs(float[][] wp, float[][] p2, float[][] p3, float[][] y) {
    zero(y);
    int n1 = y[0].length;
    int n2 = y.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float wpi = (wp!=null)?wp[i2][i1]:1.000f;
        if(wpi<0.05f) {wpi=0.05f;}
        float p2i = p2[i2][i1];
        float p3i = p3[i2][i1];
        float b11 = wpi;
        float b22 = wpi;
        float x1 = wpi*p2i;
        float x2 = wpi*p3i;
        float y1 = b11*x1;
        float y2 = b22*x2;
        float ya = 0.5f*(y1+y2);
        float yb = 0.5f*(y1-y2);
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  private float _weight = 0.5f; // -(fxx+fyy)+weight*(fxxxx+2fxxyy+fyyyy)=px+qy
  private float _sigma1 = 6.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 6.0f; // precon smoothing extent for 2nd dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 200; // maximum number of CG iterations
  private int _exniter = 10; // external iterations of surface updating

  private KdTree _kt;
  private int _n1,_n2,_n3;
  private FaultCell[] _fcs;

}
