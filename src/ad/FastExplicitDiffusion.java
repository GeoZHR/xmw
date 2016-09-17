/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ad;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import util.*;

/**
 * Fast explicit diffusion filter. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.12.05
 */

public class FastExplicitDiffusion {

  /**
   * Construction.
   * @param m number of cycles
   * @param d stability limit
   */
  public void setCycles(int m, float d) {
    _m = m;
    _d = d;
  }

  public float[][] apply(
    float sigma, EigenTensors2 et, float[][] fx) 
  {
    float t = sigma*sigma*0.5f;
    FedStep fs = new FedStep(t,_m,_d);
    float[] ts = fs.getSteps(true);
    int nc = ts.length;
    //et.setEigenvalues(0.001f,1.0f);
    float[][] gx = copy(fx);
    for (int m=0; m<_m; ++m) {
    for (int ic=0; ic<nc; ++ic) {
      applyLaplacian(et,-ts[ic],copy(gx),gx);
    }}
    return gx;
  }

  public float[][] apply(
    float sigma, float lambda, EigenTensors2 et, float[][] fx) 
  {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] p2 = new float[n2][n1];
    slopesFromTensors(et,p2);
    float t = sigma*sigma*0.5f;
    FedStep fs = new FedStep(t,_m,_d);
    float[] ts = fs.getSteps(true);
    int nc = ts.length;
    float[][] gx = copy(fx);
    et.setEigenvalues(0.001f,1.0f);
    for (int m=0; m<_m; ++m) {
      setNonlinearDiffusion(lambda,p2,gx,et);
      for (int ic=0; ic<nc; ++ic) {
        applyLaplacian(et,-ts[ic],copy(gx),gx);
      }
    }
    return gx;
  }

    /**
   * Applies a simple 3x3x3 weighted-average smoothing filter S.
   * Input and output arrays x and y may be the same array.
   * @param x input array.
   * @param y output array.
   */
  public void applySmoothS(float[][][] x, float[][][] y) {
    smoothS(x,y);
  }


  public float[][][] apply(
    float sigma, EigenTensors3 et,float[][][] fx) 
  {
    float[][][] gx = copy(fx);
    float t = sigma*sigma*0.5f;
    FedStep fs = new FedStep(t,_m,_d);
    float[] ts = fs.getSteps(true);
    int nc = ts.length;
    Stopwatch sw = new Stopwatch();
    sw.start();
    int ik = 0;
    int nk = _m*nc;
    for (int m=0; m<_m; ++m) {
    for (int ic=0; ic<nc; ++ic) {
      if (ik>0) {
        double timeUsed = sw.time();
        double timeLeft = ((double)nk/(double)ik-1.0)*timeUsed;
        int timeLeftSec = 1+(int)timeLeft;
        trace("Linear diffusion: done in "+timeLeftSec+" seconds");
      }
      applyLaplacian(et,-ts[ic],copy(gx),gx);
      ik++;
    }}
    sw.stop();
    trace("Linear diffusion: done");
    return gx;
  }

  public float[][][] apply(
    float sigma, EigenTensors3 et,float[][][] wp, float[][][] fx) 
  {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float wpi = wp[i3][i2][i1];
      et.setEigenvalues(i1,i2,i3,0.0001f,wpi,1f);
    }}}
    float[][][] gx = copy(fx);
    float t = sigma*sigma*0.5f;
    FedStep fs = new FedStep(t,_m,_d);
    float[] ts = fs.getSteps(true);
    int nc = ts.length;
    Stopwatch sw = new Stopwatch();
    sw.start();
    int ik = 0;
    int nk = _m*nc;
    for (int m=0; m<_m; ++m) {
    for (int ic=0; ic<nc; ++ic) {
      if (ik>0) {
        double timeUsed = sw.time();
        double timeLeft = ((double)nk/(double)ik-1.0)*timeUsed;
        int timeLeftSec = 1+(int)timeLeft;
        trace("Linear diffusion: done in "+timeLeftSec+" seconds");
      }
      applyLaplacian(et,-ts[ic],copy(gx),gx);
      ik++;
    }}
    sw.stop();
    trace("Linear diffusion: done");
    return gx;
  }

  public float[][][][] apply(
    float sigma, float lambda, EigenTensors3 et, float[][][] fx)
  {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[][][] gx = copy(fx);
    float[][][] p2 = new float[n3][n2][n1];
    float[][][] p3 = new float[n3][n2][n1];
    float[][][] sc = new float[n3][n2][n1];
    slopesFromTensors(et,p2,p3);
    float t = sigma*sigma*0.5f;
    FedStep fs = new FedStep(t,_m,_d);
    float[] ts = fs.getSteps(true);
    int nc = ts.length;
    Stopwatch sw = new Stopwatch();
    sw.start();
    int ik = 0;
    int nk = _m*nc;
    for (int m=0; m<_m; ++m) {
    setNonlinearDiffusion(lambda,p2,p3,gx,sc,et);
    for (int ic=0; ic<nc; ++ic) {
      if (ik>0) {
        double timeUsed = sw.time();
        double timeLeft = ((double)nk/(double)ik-1.0)*timeUsed;
        int timeLeftSec = 1+(int)timeLeft;
        trace("Nonlinear diffusion: done in "+timeLeftSec+" seconds");
      }
      applyLaplacian(et,-ts[ic],copy(gx),gx);
      ik++;
   }}
   sw.stop();
   trace("Nonlinear diffusion: done");
   return new float[][][][]{sc,gx};
  }

  public float[][][][] apply(
    float sigma, float lambda, float av, EigenTensors3 et, float[][][] fx)
  {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[][][] gx = copy(fx);
    float[][][] p2 = new float[n3][n2][n1];
    float[][][] p3 = new float[n3][n2][n1];
    float[][][] sc = new float[n3][n2][n1];
    slopesFromTensors(et,p2,p3);
    float t = sigma*sigma*0.5f;
    FedStep fs = new FedStep(t,_m,_d);
    float[] ts = fs.getSteps(true);
    int nc = ts.length;
    Stopwatch sw = new Stopwatch();
    sw.start();
    int ik = 0;
    int nk = _m*nc;
    for (int m=0; m<_m; ++m) {
    setNonlinearDiffusion(lambda,av,p2,p3,gx,sc,et);
    for (int ic=0; ic<nc; ++ic) {
      if (ik>0) {
        double timeUsed = sw.time();
        double timeLeft = ((double)nk/(double)ik-1.0)*timeUsed;
        int timeLeftSec = 1+(int)timeLeft;
        trace("Nonlinear diffusion: done in "+timeLeftSec+" seconds");
      }
      applyLaplacian(et,-ts[ic],copy(gx),gx);
      ik++;
   }}
   sw.stop();
   trace("Nonlinear diffusion: done");
   return new float[][][][]{sc,gx};
  }


  public float[][][] getWeights(float lambda, EigenTensors3 et, 
    float[][][] fx) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[][][] p2 = new float[n3][n2][n1];
    float[][][] p3 = new float[n3][n2][n1];
    float[][][] ws = new float[n3][n2][n1];
    slopesFromTensors(et,p2,p3);
    final float[][][] g2 = new float[n3][n2][n1];
    final float[][][] g3 = new float[n3][n2][n1];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      applyForDirectionalDerivative(p2[i3],fx[i3],g2[i3]);
    }});
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] p32 = new float[n3][n1];
      float[][] fx2 = new float[n3][n1];
      float[][] g32 = new float[n3][n1];
      for (int i3=0; i3<n3; ++i3) {
        p32[i3] = p3[i3][i2];
        fx2[i3] = fx[i3][i2];
      }
      applyForDirectionalDerivative(p32,fx2,g32);
      for (int i3=0; i3<n3; ++i3)
        g3[i3][i2] = g32[i3];
    }});
    final float[][][] gs = new float[n3][n2][n1];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float g2i = g2[i3][i2][i1];
        float g3i = g3[i3][i2][i1];
        gs[i3][i2][i1] = g2i*g2i+g3i*g3i;
      }}
    }});
    final float ls = lambda*lambda;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float gsi = gs[i3][i2][i1];
        float gli = gsi/ls;
        gli *= gli;
        gli *= gli;
        float av = 1.0f;
        if(gsi>0f){av=1-exp(-3.315f/gli);}
        ws[i3][i2][i1] = av;
      }}
    }});
    return ws;
  }

  private void slopesFromTensors(EigenTensors2 et, float[][] p2) {
    int n2 = p2.length;
    int n1 = p2[0].length;
    float p2min = -5f;
    float p2max =  5f;
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float[] u = et.getEigenvectorU(i1,i2);
      float u1i = u[0];
      float u2i = u[1];
      if (u1i<0f) {u1i=-u1i; u2i=-u2i;}
      if (-u2i<p2min*u1i) u2i = -p2min*u1i;
      if (-u2i>p2max*u1i) u2i = -p2max*u1i;
      if (u1i==0.0f) {
        p2[i2][i1] = (u2i<0.0f)?p2max:p2min;
      } else {
        p2[i2][i1] = -u2i/u1i;
      }
    }}
  }

  private void slopesFromTensors(
    EigenTensors3 et, float[][][] p2, float[][][] p3) 
  {
    int n3 = p3.length;
    int n2 = p2[0].length;
    int n1 = p2[0][0].length;
    float p2min = -5f;
    float p3min = -5f;
    float p2max =  5f;
    float p3max =  5f;
    // Compute slopes from normal vectors.
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float[] u = et.getEigenvectorU(i1,i2,i3);
      float u1i = u[0];
      float u2i = u[1];
      float u3i = u[2];
      if (u1i<0f){u1i=-u1i;u2i=-u2i;u3i=-u3i;}
      if (-u2i<p2min*u1i) u2i = -p2min*u1i;
      if (-u2i>p2max*u1i) u2i = -p2max*u1i;
      if (-u3i<p3min*u1i) u3i = -p3min*u1i;
      if (-u3i>p3max*u1i) u3i = -p3max*u1i;
      if (u1i==0.0f) {
        p2[i3][i2][i1] = (u2i<0.0f)?p2max:p2min;
        p3[i3][i2][i1] = (u3i<0.0f)?p3max:p3min;
      } else {
        p2[i3][i2][i1] = -u2i/u1i;
        p3[i3][i2][i1] = -u3i/u1i;
      }
    }}}
  }


  private void setNonlinearDiffusion(
    float lambda, float[][] p2, float[][] fx, EigenTensors2 et) 
  {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] gx = new float[n2][n1];
    float[][] sx = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
    }}
    applyForDirectionalDerivative(p2,fx,gx);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float gi = gx[i2][i1];
      float gs = gi*gi;
      float ls = lambda*lambda;
      float gl = gs/ls;
      gl *= gl;
      gl *= gl;
      if(gi==0f) { sx[i2][i1] = 1f;
      } else {sx[i2][i1] = 1f-exp(-3.315f/gl);}
    }}
    et.setEigenvalues(1.0f,0.05f); //smooth in seismic normal direction
    FedStep fs = new FedStep(4,_m,_d);
    float[] ts = fs.getSteps(true);
    int nc = ts.length;
    for (int m=0; m<_m; ++m) {
    for (int ic=0; ic<nc; ++ic) {
      applyLaplacian(et,-ts[ic],copy(sx),sx);
    }}
    RecursiveGaussianFilterP rgf =  new RecursiveGaussianFilterP(1.0);
    rgf.applyX0(sx,sx);
    float[] sm = new float[n1];
    float[] tm = new float[n1];
    float[] sp = new float[n1];
    float[] tp = new float[n1];
    for (int i2=0; i2<n2; ++i2) {
      int i2m = max(i2-1,0);
      int i2p = min(i2+1,n2-1);
      float dx2 = 0.5f*(i2p-i2m);
      for (int i1=0; i1<n1; ++i1) {
        tm[i1] = i1-p2[i2][i1]*dx2;
        tp[i1] = i1+p2[i2][i1]*dx2;
      }
      _si.interpolate(n1,1.0,0.0,sx[i2m],n1,tm,sm);
      _si.interpolate(n1,1.0,0.0,sx[i2p],n1,tp,sp);
      for (int i1=0; i1<n1; ++i1) {
        float mi = sm[i1];
        float pi = sp[i1];
        float si = sx[i2][i1];
        if(si<=mi&&si<=pi&&si<1f) {
          float s2 = si*si;
          float s4 = s2*s2;
          si = s2*s4;
          mi *= mi;
          mi *= mi;
          int k2  = min(i2+1,n2-1);
          int k2m = min(i2m+1,n2-1);
          int k2p = min(i2p+1,n2-1);
          et.setEigenvalues(i1,k2, 0.001f,si);
          et.setEigenvalues(i1,k2m,0.001f,mi);
          et.setEigenvalues(i1,k2p,0.001f,pi);
        } else {
          int k2 = min(i2+1,n2-1);
          et.setEigenvalues(i1,k2, 0.001f,1.0f);
        }
      }
    }
  }

  // for plots only
  public void applyForWeightsP0(
    float lambda, EigenTensors2 et, float[][] f, float[][] w) 
  {
    int n2 = f.length;
    int n1 = f[0].length;
    float[][] g = new float[n2][n1];
    float[][] p = new float[n2][n1];
    slopesFromTensors(et,p);
    applyForDirectionalDerivative(p,f,g);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float gi = g[i2][i1];
      float gs = gi*gi;
      float ls = lambda*lambda;
      float gl = gs/ls;
      gl *= gl;
      gl *= gl;
      if(gi==0f) { w[i2][i1] = 1f;
      } else {w[i2][i1] = 1f-exp(-3.315f/gl);}
    }}
  }

  // for plots only
  public void applyForWeightsP1(
    float lambda, EigenTensors2 et, float[][] f, float[][] w) 
  {
    int n2 = f.length;
    int n1 = f[0].length;
    float[][] g = new float[n2][n1];
    float[][] p = new float[n2][n1];
    slopesFromTensors(et,p);
    applyForDirectionalDerivative(p,f,g);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float gi = g[i2][i1];
      float gs = gi*gi;
      float ls = lambda*lambda;
      float gl = gs/ls;
      gl *= gl;
      gl *= gl;
      if(gi==0f) {w[i2][i1] = 1f;} 
      else {w[i2][i1] = 1f-exp(-3.315f/gl);}
    }}
    et.setEigenvalues(1.0f,0.05f);
    FedStep fs = new FedStep(4,_m,_d);
    float[] ts = fs.getSteps(true);
    int nc = ts.length;
    for (int m=0; m<_m; ++m) {
    for (int ic=0; ic<nc; ++ic) {
      applyLaplacian(et,-ts[ic],copy(w),w);
    }}
    RecursiveGaussianFilterP rgf =  new RecursiveGaussianFilterP(1.0);
    rgf.applyX0(w,w);
  }

  // for plots only
  public void applyForWeightsP(
    float lambda, EigenTensors2 et, float[][] f, float[][] w) 
  {
    int n2 = f.length;
    int n1 = f[0].length;
    float[][] g = new float[n2][n1];
    float[][] s = new float[n2][n1];
    float[][] p = new float[n2][n1];
    slopesFromTensors(et,p);
    applyForDirectionalDerivative(p,f,g);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float gi = g[i2][i1];
      float gs = gi*gi;
      float ls = lambda*lambda;
      float gl = gs/ls;
      gl *= gl;
      gl *= gl;
      if(gi==0f) { s[i2][i1] = 1f;
      } else {s[i2][i1] = 1f-exp(-3.315f/gl);}
      w[i2][i1] = 1f; // reset w
    }}
    et.setEigenvalues(1.0f,0.05f);
    FedStep fs = new FedStep(4,_m,_d);
    float[] ts = fs.getSteps(true);
    int nc = ts.length;
    for (int m=0; m<_m; ++m) {
    for (int ic=0; ic<nc; ++ic) {
      applyLaplacian(et,-ts[ic],copy(s),s);
    }}
    et.setEigenvalues(0.001f,1.0f);
    RecursiveGaussianFilterP rgf =  new RecursiveGaussianFilterP(1.0);
    rgf.applyX0(s,s);
    float[] sm = new float[n1];
    float[] tm = new float[n1];
    float[] sp = new float[n1];
    float[] tp = new float[n1];
    for (int i2=0; i2<n2; ++i2) {
      int i2m = max(i2-1,0);
      int i2p = min(i2+1,n2-1);
      float dx2 = 0.5f*(i2p-i2m);
      for (int i1=0; i1<n1; ++i1) {
        tm[i1] = i1-p[i2][i1]*dx2;
        tp[i1] = i1+p[i2][i1]*dx2;
      }
      _si.interpolate(n1,1.0,0.0,s[i2m],n1,tm,sm);
      _si.interpolate(n1,1.0,0.0,s[i2p],n1,tp,sp);
      for (int i1=0; i1<n1; ++i1) {
        float si = s[i2][i1];
        float mi = sm[i1];
        float pi = sp[i1];
        if(si<=mi&&si<=pi&&si<1f) {
          float s2 = si*si;
          float s4 = s2*s2;
          si = s2*s4;
          int p2 = min(i2+1,n2-1);
          w[p2][i1]  = si;
        }
      }
    }
  }

  public float[][][] thin(float sig1, float sig2, 
      final float flmax, final float[][][] fl) {
    int n3 = fl.length;
    int n2 = fl[0].length;
    int n1 = fl[0][0].length;
    LocalOrientFilter lof = new LocalOrientFilter(sig1,sig2);
    final float[][][] u1 = new float[n3][n2][n1];
    final float[][][] u2 = new float[n3][n2][n1];
    final float[][][] u3 = new float[n3][n2][n1];
    final float[][][] ft = fillfloat(1f,n1,n2,n3);
    lof.applyForNormal(fl,u1,u2,u3);
    final Sampling s1 = new Sampling(n1);
    final Sampling s2 = new Sampling(n2);
    final Sampling s3 = new Sampling(n3);
    final SincInterpolator si = new SincInterpolator();
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float u1i = u1[i3][i2][i1];
        float u2i = u2[i3][i2][i1];
        float u3i = u3[i3][i2][i1];
        float x1p = i1+u1i;
        float x2p = i2+u2i;
        float x3p = i3+u3i;
        float x1m = i1-u1i;
        float x2m = i2-u2i;
        float x3m = i3-u3i;
        float fli = fl[i3][i2][i1];
        float flp = si.interpolate(s1,s2,s3,fl,x1p,x2p,x3p);
        float flm = si.interpolate(s1,s2,s3,fl,x1m,x2m,x3m);
        if (flp>fli && flm>fli && fli<flmax) {
          ft[i3][i2][i1] = fl[i3][i2][i1];
        }
      }}
    }});
    return ft;
  }

  private void setNonlinearDiffusion(
    final float lambda, final float av, final float[][][] p2, final float[][][] p3, 
    final float[][][] fx, final float[][][] sc, final EigenTensors3 et) 
  {
    final int n3 = fx.length;
    final int n2 = fx[0].length;
    final int n1 = fx[0][0].length;
    final float[][][] g2 = new float[n3][n2][n1];
    final float[][][] g3 = new float[n3][n2][n1];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      applyForDirectionalDerivative(p2[i3],fx[i3],g2[i3]);
    }});
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] p32 = new float[n3][n1];
      float[][] fx2 = new float[n3][n1];
      float[][] g32 = new float[n3][n1];
      for (int i3=0; i3<n3; ++i3) {
        p32[i3] = p3[i3][i2];
        fx2[i3] = fx[i3][i2];
      }
      applyForDirectionalDerivative(p32,fx2,g32);
      for (int i3=0; i3<n3; ++i3)
        g3[i3][i2] = g32[i3];
    }});
    final float[][][] gs = new float[n3][n2][n1];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float g2i = g2[i3][i2][i1];
        float g3i = g3[i3][i2][i1];
        gs[i3][i2][i1] = g2i*g2i+g3i*g3i;
      }}
    }});
    final float ls = lambda*lambda;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float gsi = gs[i3][i2][i1];
        float gli = gsi/ls;
        gli *= gli;
        gli *= gli;
        float avi = 1.0f;
        if(gsi>0f){avi=1-exp(-3.315f/gli);}
        sc[i3][i2][i1] = avi;
        et.setEigenvalues(i1,i2,i3,0.0001f,avi*avi*av,1f);
      }}
    }});
  }


  private void setNonlinearDiffusion(
    final float lambda, final float[][][] p2, final float[][][] p3, 
    final float[][][] fx, final float[][][] sc, final EigenTensors3 et) 
  {
    final int n3 = fx.length;
    final int n2 = fx[0].length;
    final int n1 = fx[0][0].length;
    final float[][][] g2 = new float[n3][n2][n1];
    final float[][][] g3 = new float[n3][n2][n1];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      applyForDirectionalDerivative(p2[i3],fx[i3],g2[i3]);
    }});
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] p32 = new float[n3][n1];
      float[][] fx2 = new float[n3][n1];
      float[][] g32 = new float[n3][n1];
      for (int i3=0; i3<n3; ++i3) {
        p32[i3] = p3[i3][i2];
        fx2[i3] = fx[i3][i2];
      }
      applyForDirectionalDerivative(p32,fx2,g32);
      for (int i3=0; i3<n3; ++i3)
        g3[i3][i2] = g32[i3];
    }});
    final float[][][] gs = new float[n3][n2][n1];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float g2i = g2[i3][i2][i1];
        float g3i = g3[i3][i2][i1];
        gs[i3][i2][i1] = g2i*g2i+g3i*g3i;
      }}
    }});
    final float ls = lambda*lambda;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float gsi = gs[i3][i2][i1];
        float gli = gsi/ls;
        gli *= gli;
        gli *= gli;
        float av = 1.0f;
        if(gsi>0f){av=1-exp(-3.315f/gli);}
        sc[i3][i2][i1] = av;
        float[] eu = et.getEigenvalues(i1,i2,i3);
        et.setEigenvalues(i1,i2,i3,0.0001f,av*eu[1],eu[2]);
      }}
    }});
  }

  private void applyForDirectionalDerivative(
    float[][] p, float[][] f, float[][] g)
  {
    int n2 = f.length;
    int n1 = f[0].length;
    float[] fm = new float[n1];
    float[] tm = new float[n1];
    float[] fp = new float[n1];
    float[] tp = new float[n1];
    for (int i2=0; i2<n2; ++i2) {
      int i2m = max(i2-1,0);
      int i2p = min(i2+1,n2-1);
      float dx2 = 0.5f*(i2p-i2m);
      for (int i1=0; i1<n1; ++i1) {
        tm[i1] = i1-p[i2][i1]*dx2;
        tp[i1] = i1+p[i2][i1]*dx2;
      }
      _si.interpolate(n1,1.0,0.0,f[i2m],n1,tm,fm);
      _si.interpolate(n1,1.0,0.0,f[i2p],n1,tp,fp);
      float scale = 0.5f/dx2;
      for (int i1=0; i1<n1; ++i1)
        g[i2][i1] = scale*(fp[i1]-fm[i1]);
    }
  }

  private void applyLaplacian(
    EigenTensors2 et, float s, float[][] x, float[][] y){
    s *= 0.25f;
    float[] ds = fillfloat(1.0f,3);
    ds[1] = 0.0f;
    int n2 = x.length;
    int n1 = x[0].length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        if(et!=null){et.getTensor(i1,i2,ds);}
        float d11 = ds[0];
        float d12 = ds[1];
        float d22 = ds[2];
        float xa = 0.0f;
        float xb = 0.0f;
        xa += x[i2  ][i1  ];
        xb -= x[i2  ][i1-1];
        xb += x[i2-1][i1  ];
        xa -= x[i2-1][i1-1];
        float x1 = (xa+xb);
        float x2 = (xa-xb);
        float y1 = (d11*x1+d12*x2)*s;
        float y2 = (d12*x1+d22*x2)*s;
        float ya = (y1+y2);
        float yb = (y1-y2);
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }


  private void applyLaplacian(final Tensors3 d, final float s, 
    final float[][][] x, final float[][][] y) 
  {
    int i3start = 1; 
    final int i3step = 2; 
    final int i3stop = x.length;
    for (int i3pass=0; i3pass<i3step; ++i3pass,++i3start) {
      Parallel.loop(i3start,i3stop,i3step,new Parallel.LoopInt() {
        public void compute(int i3) {apply22(i3,d,s,x,y);}});
    }
  }

  private void apply22(int i3, Tensors3 d, float s, 
    float[][][] x, float[][][] y) 
  {
    s *= 0.0625f;
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    float[] di = new float[6];
    for (int i2=1; i2<n2; ++i2) {
      float[] x00 = x[i3  ][i2  ];
      float[] x0m = x[i3  ][i2-1];
      float[] xm0 = x[i3-1][i2  ];
      float[] xmm = x[i3-1][i2-1];
      float[] y00 = y[i3  ][i2  ];
      float[] y0m = y[i3  ][i2-1];
      float[] ym0 = y[i3-1][i2  ];
      float[] ymm = y[i3-1][i2-1];
      for (int i1=1,m1=0; i1<n1; ++i1,++m1) {
        d.getTensor(i1,i2,i3,di);
        float d11 = di[0];
        float d12 = di[1];
        float d13 = di[2];
        float d22 = di[3];
        float d23 = di[4];
        float d33 = di[5];
        float xa = x00[i1]-xmm[m1];
        float xb = x00[m1]-xmm[i1];
        float xc = x0m[i1]-xm0[m1];
        float xd = xm0[i1]-x0m[m1];
        float x1 = xa-xb+xc+xd;
        float x2 = xa+xb-xc+xd;
        float x3 = xa+xb+xc-xd;
        float y1 = (d11*x1+d12*x2+d13*x3)*s;
        float y2 = (d12*x1+d22*x2+d23*x3)*s;
        float y3 = (d13*x1+d23*x2+d33*x3)*s;
        float ya = y1+y2+y3; y00[i1] += ya; ymm[m1] -= ya;
        float yb = y1-y2+y3; y0m[i1] += yb; ym0[m1] -= yb;
        float yc = y1+y2-y3; ym0[i1] += yc; y0m[m1] -= yc;
        float yd = y1-y2-y3; ymm[i1] += yd; y00[m1] -= yd;
      }
    }
  }

    /*
   * Computes y = S'Sx. Arrays x and y may be the same array.
   */
  private static void smoothS(float[][][] x, float[][][] y) {
    int n1 = x[0][0].length;
    int n2 = x[0].length;
    int n3 = x.length;
    int n1m = n1-1;
    int n2m = n2-1;
    int n3m = n3-1;
    float[][][] t = new float[3][n2][n1];
    scopy(x[0],t[0]);
    scopy(x[0],t[1]);
    for (int i3=0; i3<n3; ++i3) {
      int i3m = (i3>0)?i3-1:0;
      int i3p = (i3<n3m)?i3+1:n3m;
      int j3m = i3m%3;
      int j3  = i3%3;
      int j3p = i3p%3;
      scopy(x[i3p],t[j3p]);
      float[][] x3m = t[j3m];
      float[][] x3p = t[j3p];
      float[][] x30 = t[j3];
      float[][] y30 = y[i3];
      for (int i2=0; i2<n2; ++i2) {
        int i2m = (i2>0)?i2-1:0;
        int i2p = (i2<n2m)?i2+1:n2m;
        float[] x3m2m = x3m[i2m];
        float[] x3m20 = x3m[i2 ];
        float[] x3m2p = x3m[i2p];
        float[] x302m = x30[i2m];
        float[] x3020 = x30[i2 ];
        float[] x302p = x30[i2p];
        float[] x3p2m = x3p[i2m];
        float[] x3p20 = x3p[i2 ];
        float[] x3p2p = x3p[i2p];
        float[] y3020 = y30[i2 ];
        for (int i1=0; i1<n1; ++i1) {
          int i1m = (i1>0)?i1-1:0;
          int i1p = (i1<n1m)?i1+1:n1m;
          y3020[i1] = 0.125000f*(x3020[i1 ]) +
                      0.062500f*(x3020[i1m]+x3020[i1p]+
                                 x302m[i1 ]+x302p[i1 ]+
                                 x3m20[i1 ]+x3p20[i1 ]) +
                      0.031250f*(x3m20[i1m]+x3m20[i1p]+
                                 x3m2m[i1 ]+x3m2p[i1 ]+
                                 x302m[i1m]+x302m[i1p]+
                                 x302p[i1m]+x302p[i1p]+
                                 x3p20[i1m]+x3p20[i1p]+
                                 x3p2m[i1 ]+x3p2p[i1 ]) +
                      0.015625f*(x3m2m[i1m]+x3m2m[i1p]+
                                 x3m2p[i1m]+x3m2p[i1p]+
                                 x3p2m[i1m]+x3p2m[i1p]+
                                 x3p2p[i1m]+x3p2p[i1p]);
        }
      }
    }
  }

    // Copys array x to array y.
  private static void scopy(float[] x, float[] y) {
    copy(x,y);
  }
  private static void scopy(float[][] x, float[][] y) {
    copy(x,y);
  }
  private static void scopy(final float[][][] x, final float[][][] y) {
    final int n3 = x.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        scopy(x[i3],y[i3]);
      }
    });
  }



  private static void trace(String s) {
    System.out.println(s);
  }

  private int _m = 5; //number of cycles
  private float _d = 0.5f; //stability limit
  private SincInterpolator _si = new SincInterpolator();
}
