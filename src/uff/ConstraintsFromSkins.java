/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package uff;

import ipf.*;
import java.util.*;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.SincInterpolator;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Set weights and constraints for flattening: 
 * set zero weights at fault, 
 * set hard constraints using fault slip vectors.
 * @author Xinming Wu
 * @version 2014.09.18
 */

public class ConstraintsFromSkins {

  public ConstraintsFromSkins(FaultSkin[] fss, int wse, int cse, float[][][] p, float[][][] q, float[][][] w) {
    _p = p;
    _q = q;
    _w = w;
    _fss = fss;
    _wse = wse;
    _cse = cse;
    _n3 = p.length;
    _n2 = p[0].length;
    _n1 = p[0][0].length;
  }

  public void setForHorizonExtraction(float sigma1, float sigma2, float small, int niter) {
    _small = small;
    _niter = niter;
    _sigma1 = sigma1;
    _sigma2 = sigma2;
    _se.setCG(_small,_niter);
    _se.setSmoothings(_sigma1,_sigma2);
  }

  public float[][][] getWeightsAndConstraints(float[][][] ws) {
    int n3 = ws.length;
    int n2 = ws[0].length;
    ArrayList<float[][]> cl = new ArrayList<float[][]>();
    for (FaultSkin fs:_fss) {
      for (FaultCell fc:fs) {
        float[] pi = fc.getX();
        float[] wi = fc.getW();
        float[] si = new float[3];
        si[0] = fc.getS1();
        si[1] = fc.getS2();
        si[2] = fc.getS3();
        float[] pp = add(pi,si);
        int pi1 = round(pi[0]);
        int pi2 = round(pi[1]);
        int pi3 = round(pi[2]);
        int p2m = pi2-1;if(p2m<0){p2m=0;}
        int p3m = pi3-1;if(p3m<0){p3m=0;}
        int p2p = pi2+1;if(p2p>n2-1){p2p=n2-1;}
        int p3p = pi3+1;if(p3p>n3-1){p3p=n3-1;}
        ws[pi3][pi2][pi1] = 0.0f;
        ws[pi3][p2m][pi1] = 0.0f;
        ws[p3m][pi2][pi1] = 0.0f;
        /*
        ws[pi3][p2p][pi1] = 0.0f;
        ws[p3m][pi2][pi1] = 0.0f;
        ws[p3p][pi2][pi1] = 0.0f;
        */
        /*
        boolean bd = boundCheck(pi,pp);
        if(bd){continue;}
        int[] bi = new int[3];
        int[] bp = new int[3];
        float[][] sfi = subSurfer(pi,bi);
        float[][] sfp = subSurfer(pp,bp);
        float[][] ps = findControlPoint(pi,pp,bi,bp,sfi,sfp);
        */
        //float[][] ps = getControlPoints(pi,pp);
        float[][] ps = getControlPoints(pi,pp,wi);
        if (ps==null){continue;}
        cl.add(ps);
      }
    }
    int np = cl.size();
    System.out.println("sets of constraints before:"+np);
    int nn = checkConstraints(cl);
    np = cl.size();
    int nc = np-nn;
    System.out.println("sets of constraints after:"+np);
    System.out.println("sets of constraints after:"+nc);
    float[][][] cs = new float[3][nc][];
    int ic = 0;
    for (int ip=0; ip<np; ++ip) {
      float[][] ps = cl.get(ip);
      if(ps!=null) {
        int ns = ps.length;
        cs[0][ic] = new float[ns];
        cs[1][ic] = new float[ns];
        cs[2][ic] = new float[ns];
        for (int is=0; is<ns; ++is) {
          cs[0][ic][is] = ps[is][0];
          cs[1][ic][is] = ps[is][1];
          cs[2][ic][is] = ps[is][2];
        }
        ic++;
      }
    }
    return cs;
  }

  private boolean boundCheck(float[] pi, float[] pp) {
    if(pi[0]>=_n1)  {return true;}
    if(pp[0]>=_n1)  {return true;}
    if(pi[1]<=1)    {return true;}
    if(pp[1]<=1)    {return true;}
    if(pi[2]<=1)    {return true;}
    if(pp[2]<=1)    {return true;}
    if(pi[1]+2>=_n2){return true;}
    if(pp[1]+2>=_n2){return true;}
    if(pi[2]+2>=_n3){return true;}
    if(pp[2]+2>=_n3){return true;}
    return false;
  }

  private float[][] getControlPoints(float[] pi, float[] pp) {
    int pi1 = round(pi[0]);
    int pi2 = round(pi[1]);
    int pi3 = round(pi[2]);
    int pp1 = round(pp[0]);
    int pp2 = round(pp[1]);
    int pp3 = round(pp[2]);
    if (pp2>pi2) {pp2 +=_cse;pi2 -=_cse;} 
    else         {pp2 -=_cse;pi2 +=_cse;}
    if (pp3>pi3) {pp3 +=_cse;pi3 -=_cse;} 
    else         {pp3 -=_cse;pi3 +=_cse;}
    if (pi1<0||pi1>=_n1) {return null;}
    if (pp1<0||pp1>=_n1) {return null;}
    if (pi2<0||pi2>=_n2) {return null;}
    if (pp2<0||pp2>=_n2) {return null;}
    if (pi3<0||pi3>=_n3) {return null;}
    if (pp3<0||pp3>=_n3) {return null;}
    if (pi1!=pp1 || pi2!=pp2 || pi3!=pp3) { 
      float[] p1 = new float[]{pi1,pi2,pi3};
      float[] p2 = new float[]{pp1,pp2,pp3};
      return new float[][]{p1,p2};
    } else {return null;}
  }

  private int checkConstraints(ArrayList<float[][]> cl) {
    int n12 = _n1*_n2;
    int np  = n12*_n3;
    int nc = cl.size();
    int[][][] mp = new int[np][2][10];
    for (int ic=0; ic<nc; ic++) {
      float[][] ps = cl.get(ic);
      for (int is=0; is<=1; ++is) {
        int i1 = (int)ps[is][0];
        int i2 = (int)ps[is][1];
        int i3 = (int)ps[is][2];
        int n2 = i2*_n1;
        int n3 = i3*n12;
        int pi = i1+n2+n3;
        mp[pi][0][0] += 1;
        int id = mp[pi][0][0];
        mp[pi][0][id] = ic;
        mp[pi][1][id] = is;
      }
    }
    for (int ip=0; ip<np; ++ip) {
      int mpi = mp[ip][0][0]+1;
      float[][] cs = new float[mpi][3];
      if (mpi>2) {
        for (int i=1; i<mpi; i++) {
          int ic = mp[ip][0][i];
          int is = mp[ip][1][i];
          float[][] ps = cl.get(ic);
          if (is==0) {cs[0]=ps[0];cs[i] = ps[1];} 
          else       {cs[0]=ps[1];cs[i] = ps[0];}
        }
        cl.add(cs);
      }
    }
    int nn = 0;
    for (int ip=0; ip<np; ++ip) {
      int mpi = mp[ip][0][0]+1;
      if (mpi>2) {
        for (int i=1; i<mpi; i++) {
          int ic = mp[ip][0][i];
          float[][] ps = cl.get(ic);
          if(ps!=null){
            cl.set(ic,null);
            nn++;
          }
        }
      }
    }
    System.out.println("nn="+nn);
    System.out.println("listSize="+cl.size());
    return nn;
  }

  private float[][] getControlPoints(float[] pi, float[] pp, float[] wi) {
    float ds = 1.5f;
    float w2 = wi[1];
    float w3 = wi[2];
    float sc = 1.0f/sqrt(w2*w2+w3*w3);
    float dw2 = ds*w2*sc;
    float dw3 = ds*w3*sc;
    int pi1 = round(pi[0]);
    int pp1 = round(pp[0]);
    int pi2 = round(pi[1]-dw2);
    int pi3 = round(pi[2]-dw3);
    int pp2 = round(pp[1]+dw2);
    int pp3 = round(pp[2]+dw3);
    if(pi1<0){pi1=0;} if(pi1>=_n1){pi1=_n1-1;}
    if(pi2<0){pi2=0;} if(pi2>=_n2){pi2=_n2-1;}
    if(pi3<0){pi3=0;} if(pi3>=_n3){pi3=_n3-1;}
    if(pp1<0){pp1=0;} if(pp1>=_n1){pp1=_n1-1;}
    if(pp2<0){pp2=0;} if(pp2>=_n2){pp2=_n2-1;}
    if(pp3<0){pp3=0;} if(pp3>=_n3){pp3=_n3-1;}
    float p2i = _p[pi3][pi2][pi1]; 
    float p3i = _q[pi3][pi2][pi1]; 
    float p2p = _p[pp3][pp2][pp1]; 
    float p3p = _q[pp3][pp2][pp1]; 
    float di1 = dw2*p2i+dw3*p3i;
    float dp1 = dw2*p2p+dw3*p3p;
    pi1 = round(pi[0]-di1);
    pp1 = round(pp[0]+dp1);
    if (pi1<0||pi1>=_n1) {return null;}
    if (pp1<0||pp1>=_n1) {return null;}
    if (pi2<0||pi2>=_n2) {return null;}
    if (pp2<0||pp2>=_n2) {return null;}
    if (pi3<0||pi3>=_n3) {return null;}
    if (pp3<0||pp3>=_n3) {return null;}
    if (pi1==pp1 && pi2==pp2 && pi3==pp3) { 
      return null;
    } else {
      float[] p1 = new float[]{pi1,pi2,pi3};
      float[] p2 = new float[]{pp1,pp2,pp3};
      return new float[][]{p1,p2};
    }
  }


  private float[][] subSurfer(float[] cp, int[] bs) {
    int k1 = round(cp[0]);
    int k2 = round(cp[1]);
    int k3 = round(cp[2]);
    int[][] ind = getBox(k1,k2,k3,_n1,_n2,_n3,_cse);
    int n1s = ind[0][1] - ind[0][0]+1;
    int n2s = ind[1][1] - ind[1][0]+1;
    int n3s = ind[2][1] - ind[2][0]+1;
    float[][][] ws = new float[n3s][n2s][n1s];
    float[][][] ps = new float[n3s][n2s][n1s];
    float[][][] qs = new float[n3s][n2s][n1s];
    int j1x = bs[0]= ind[0][0];
    int j2x = bs[1]= ind[1][0];
    int j3x = bs[2]= ind[2][0];
    copy(n1s,n2s,n3s,j1x,j2x,j3x,_w,0,0,0,ws);
    copy(n1s,n2s,n3s,j1x,j2x,j3x,_p,0,0,0,ps);
    copy(n1s,n2s,n3s,j1x,j2x,j3x,_q,0,0,0,qs);
    float[][] cps = new float[3][1];
    cps[0][0] = cp[0] - (float)j1x;
    cps[1][0] = cp[1] - (float)j2x;
    cps[2][0] = cp[2] - (float)j3x;
    return horizonSurfer(12,cps,ps,qs,ws,false);
  }

  private float[][] findControlPoint(
    float[] pi, float[] pp, int[] bi, int[] bp, float[][] sfi, float[][] sfp) 
  {
    float pi2 = pi[1];
    float pi3 = pi[2];
    float pp2 = pp[1];
    float pp3 = pp[2];
    int ni3 = sfi.length;
    int np3 = sfp.length;
    int ni2 = sfi[0].length;
    int np2 = sfp[0].length;
    int fi2,li2,fi3,li3;
    int fp2,lp2,fp3,lp3;
    if (pi2==pp2||pi3==pp3){return null;}
    if (pi2<pp2) {
      fi2 = 0;
      li2 = round(pi2)-bi[1]-1;
      fp2 = round(pp2)-bp[1]+1;
      lp2 = np2-1;
    } else {
      fp2 = 0;
      lp2 = round(pp2)-bp[1]-1;
      fi2 = round(pi2)-bi[1]+1;
      li2 = ni2-1;
    }
    if (pi3<pp3) {
      fi3 = 0;
      li3 = round(pi3)-bi[2]-1;
      fp3 = round(pp3)-bp[2]+1;
      lp3 = np3-1;
    } else {
      fp3 = 0;
      lp3 = round(pp3)-bp[2]-1;
      fi3 = round(pi3)-bi[2]+1;
      li3 = ni3-1;
    }
    float bg = 1000.0f;
    float[][] ps = new float[2][3];
    for (int i3=fi3; i3<=li3; ++i3) {
      for (int i2=fi2; i2<=li2; ++i2) {
        float hi = sfi[i3][i2];
        int   hr = round(hi);
        float df = abs(hi-hr);
        int i1p = hr+(int)bi[0];
        int i2p = i2+(int)bi[1];
        int i3p = i3+(int)bi[2];
        if(df<bg) {
          bg = df;
          ps[0][0] = i1p;
          ps[0][1] = i2p;
          ps[0][2] = i3p;
        }
      }
    }
    bg = 1000.0f;
    for (int i3=fp3; i3<=lp3; ++i3) {
      for (int i2=fp2; i2<=lp2; ++i2) {
        float hi = sfp[i3][i2];
        int   hr = round(hi);
        float df = abs(hi-hr);
        int i1p = hr+(int)bp[0];
        int i2p = i2+(int)bp[1];
        int i3p = i3+(int)bp[2];
        if(df<bg) {
          bg = df;
          ps[1][0] = i1p;
          ps[1][1] = i2p;
          ps[1][2] = i3p;
        }
      }
    }
    float pi1 = ps[0][0];
    float pp1 = ps[1][0];
    if(pi1<=1||pi1>=_n1-2) {return null;}
    if(pp1<=1||pp1>=_n1-2) {return null;}
    if(abs(pi1-pp1)>10.0f) {return null;} 
    return ps;
  }


  public float[][] horizonSurfer(int niter, float[][] ks, 
    float[][][] p, float[][][] q, float[][][] w, boolean NaNs) 
  {
    int n3 = p.length;
    int n2 = p[0].length;
    _se.setConstraints(ks[0],ks[1],ks[2]);
    float[][] surf = _se.surfaceInitialization(n2,n3,p,NaNs);
    float adp = 10000000000000.f;
    for (int i=0; i<niter; i++) {
      float ad = _se.surfaceUpdateFromSlopes(w,p,q,surf,NaNs);
      //System.out.println("Average adjustment="+ad);
      if(ad<0.01f || ad>=adp) {break;}
      adp = ad;
    }
    return surf;
  }


  private int[][] getBox(int k1, int k2, int k3, int n1, int n2, int n3,int dn) {
    int j3x = k3-dn;if(j3x<0){j3x=0;}
    int j2x = k2-dn;if(j2x<0){j2x=0;}
    int j1x = k1-dn;if(j1x<0){j1x=0;}
    int k3x = k3+dn;if(k3x>n3-1){k3x=n3-1;}
    int k2x = k2+dn;if(k2x>n2-1){k2x=n2-1;}
    int k1x = k1+dn;if(k1x>n1-1){k1x=n1-1;}
    int[][] ind = new int[3][2];
    ind[0][0] = j1x;ind[0][1] = k1x;
    ind[1][0] = j2x;ind[1][1] = k2x;
    ind[2][0] = j3x;ind[2][1] = k3x;
    return ind;
  }


  private int _wse = 0;
  private int _cse = 0;
  private FaultSkin[] _fss;
  private int _n1,_n2,_n3;
  private float[][][] _p = null;
  private float[][][] _q = null;
  private float[][][] _w = null;
  private float _small = 0.01f;
  private int _niter = 200;
  private float _sigma1 = 6.0f;
  private float _sigma2 = 6.0f;
  private SurfaceExtractor _se = new SurfaceExtractor();
}

