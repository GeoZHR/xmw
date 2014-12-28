/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ifs;

import java.util.*;
import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;


/**
 * Remove noisy fault cells. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2014.12.16
 */
public class RemoveOutlierCells {

  public RemoveOutlierCells (int n1, int n2, int n3, FaultCell[] fc) {
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
    _fc = fc;
    _nc = fc.length;
    initialTensor();
  }

  public FaultCell[] apply(final int d, final float ept) {
    final HashSet<Integer> hsc = new HashSet<Integer>();
    for (int ic=0; ic<_nc; ++ic) hsc.add(ic);
    final KdTree kd = setKDTree();
    final int[] ns = new int[]{_n1,_n2,_n3};
    for (int ic=0; ic<_nc; ++ic) {
      float x1 = _fc[ic].x1;
      float x2 = _fc[ic].x2;
      float x3 = _fc[ic].x3;
      float w1 = _fc[ic].w1;
      float w2 = _fc[ic].w2;
      float w3 = _fc[ic].w3;
      float d1 = d*(1f-w1); if(d1<2f) d1=2f;
      float d2 = d*(1f-w2); if(d2<2f) d2=2f;
      float d3 = d*(1f-w3); if(d3<2f) d3=2f;

      float[] xmin = new float[3];
      float[] xmax = new float[3];
      float[] ds = new float[]{d1,d2,d3};
      float[] xs = new float[]{x1,x2,x3};
      getRange(ds,xs,ns,xmin,xmax);
      int[] id = kd.findInRange(xmin,xmax);
      int nd = id.length;
      if(nd<1){continue;}
      else {
        float[] u = new float[3];
        float ep = tensorVoteS(d,xs,id,u);
        if(ep<ept) hsc.remove(ic);
      }
    }
    int ik = -1;
    int nc = hsc.size();
    FaultCell[] fc = new FaultCell[nc];
    for (int ic:hsc) {
      ik++;
      fc[ik] = _fc[ic];
    }
    return fc;
  }

  private KdTree setKDTree() {
    float[][] xc = new float[3][_nc];
    for (int ic=0; ic<_nc; ++ic) {
      xc[0][ic] = _fc[ic].x1;
      xc[1][ic] = _fc[ic].x2;
      xc[2][ic] = _fc[ic].x3;
    }
    return new KdTree(xc);
  }

  private float tensorVoteS(float d, float[] x, int[] id, float[] u) {
    float ck = 3.5f;
    int nd = id.length;
    float ss = 0.5f/(d*d);
    float tm = (float)Math.PI/6.0f;
    double[][] a = new double[3][3];
    for (int ik=0; ik<nd; ++ik) {
      int ic = id[ik];
      float d1 = _fc[ic].x1-x[0];
      float d2 = _fc[ic].x2-x[1];
      float d3 = _fc[ic].x3-x[2];
      float ds = d1*d1+d2*d2+d3*d3;
      float dl = sqrt(ds);
      float fl = pow(_fc[ic].fl,4.0f);
      if(ds==0.0f) {
        a[0][0] += fl*_st[ic][0]; //s11
        a[0][1] += fl*_st[ic][1]; //s12
        a[0][2] += fl*_st[ic][2]; //s13
        a[1][0] += fl*_st[ic][1]; //s12
        a[1][1] += fl*_st[ic][3]; //s22
        a[1][2] += fl*_st[ic][4]; //s23
        a[2][0] += fl*_st[ic][2]; //s13
        a[2][1] += fl*_st[ic][4]; //s23
        a[2][2] += fl*_st[ic][5]; //s33
      } else {
        float w1 = _fc[ic].w1;
        float w2 = _fc[ic].w2;
        float w3 = _fc[ic].w3;
        float wd = w1*d1+w2*d2+w3*d3;
        float ci = -2.0f*wd/ds;
        float v1 = w1+ci*d1;
        float v2 = w2+ci*d2;
        float v3 = w3+ci*d3;
        float vs = sqrt(v1*v1+v2*v2+v3*v3);
        if(vs==0.0f){continue;}
        vs = 1.0f/vs;
        float st = wd/dl;
        if(st>sin(tm)){continue;}
        v1 *= vs;v2 *= vs;v3 *= vs;
        float ki = 2.0f*st/dl;
        float si = dl*asin(st)/st;
        float sc = fl*exp(-(si*si+ck*ki*ki)*ss);
        float s11 = sc*v1*v1;
        float s12 = sc*v1*v2;
        float s13 = sc*v1*v3;
        float s22 = sc*v2*v2;
        float s23 = sc*v2*v3;
        float s33 = sc*v3*v3;
        a[0][0] += s11; //s11
        a[0][1] += s12; //s12
        a[0][2] += s13; //s13
        a[1][0] += s12; //s12
        a[1][1] += s22; //s22
        a[1][2] += s23; //s23
        a[2][0] += s13; //s13
        a[2][1] += s23; //s23
        a[2][2] += s33; //s33
      }
    }
    float[][] ue = solveEigenproblems(a);
    u[0] = ue[0][0];
    u[1] = ue[0][1];
    u[2] = ue[0][2];
    float eu = ue[1][0];
    float ev = ue[1][1];
    float es = (eu>0.0f)?1.0f/eu:1.0f;
    float ep = (eu-ev)*es;
    return ep;
  }

  private float tensorVote(float[] x, int[] id, float[] u) {
    int nd = id.length;
    double[][] a = new double[3][3];
    double[][] r = new double[3][3];
    for (int ik=0; ik<nd; ++ik) {
       int ic = id[ik];
       float v1 = (float)x[0]-_fc[ic].x1;
       float v2 = (float)x[1]-_fc[ic].x2;
       float v3 = (float)x[2]-_fc[ic].x3;
       float[] v = new float[]{v1,v2,v3};
      if(sum(abs(v))==0.0f) {
        a[0][0] += _st[ic][0]; //s11
        a[0][1] += _st[ic][1]; //s12
        a[0][2] += _st[ic][2]; //s13
        a[1][0] += _st[ic][1]; //s12
        a[1][1] += _st[ic][3]; //s22
        a[1][2] += _st[ic][4]; //s23
        a[2][0] += _st[ic][2]; //s13
        a[2][1] += _st[ic][4]; //s23
        a[2][2] += _st[ic][5]; //s33
      } else {
        float[] t = new float[1];
        float w1 = _fc[ic].w1;
        float w2 = _fc[ic].w2;
        float w3 = _fc[ic].w3;
        float[] w = new float[]{w1,w2,w3};
        float sc = computeScale(t,v,_st[ic],w,r);
        if(sc==0.0f){continue;}
        if(Float.isNaN(sc)){continue;}
        double[][] sr = computeRSR(_st[ic],r);
        a[0][0] += sc*sr[0][0]; 
        a[0][1] += sc*sr[0][1];
        a[0][2] += sc*sr[0][2];
        a[1][0] += sc*sr[1][0]; 
        a[1][1] += sc*sr[1][1];
        a[1][2] += sc*sr[1][2];
        a[2][0] += sc*sr[2][0]; 
        a[2][1] += sc*sr[2][1];
        a[2][2] += sc*sr[2][2];
      }
    }
    float[][] ue = solveEigenproblems(div(a,nd));
    u[0] = ue[0][0];
    u[1] = ue[0][1];
    u[2] = ue[0][2];
    float eu = ue[1][0];
    float ev = ue[1][1];
    float es = (eu>0.0f)?1.0f/eu:1.0f;
    return (eu-ev)*es;
  }


  private float computeScale(
    float[] k, float[] v, float[] s, float[] w, double[][] r) 
  {
    float b = 0.0f;
    float sigma = 20.0f;
    float v1 = v[0];
    float v2 = v[1];
    float v3 = v[2];
    float s11 = s[0];
    float s12 = s[1];
    float s13 = s[2];
    float s22 = s[3];
    float s23 = s[4];
    float s33 = s[5];
    float sds = v1*v1+v2*v2+v3*v3;

    float sv1 = s11*v1+s12*v2+s13*v3;
    float sv2 = s12*v1+s22*v2+s23*v3;
    float sv3 = s13*v1+s23*v2+s33*v3;
    float sns = sv1*v1+sv2*v2+sv3*v3;
    float vsr = sqrt(sds); 
    float sin = sqrt(sns/sds);
    float tha = asin(sin);
    float arc = vsr*tha/sin;
    if(sin==0.0f){arc=vsr;}
    float cur = 2.0f*sin/vsr;
    k[0] = cur;
    if(tha>Math.PI/4.0f) {
      return 0.0f;
    } else {
      rotateMatrix(tha,v,w,r);
      float ep = (arc*arc+b*cur*cur)/(sigma*sigma);
      return exp(-ep);
    }
  }

  private double[][] computeRSR(float[] s, double[][] r) {
    double[][] sm = new double[3][3];
    double[][] rt = new double[3][3];
    sm[0][0] = s[0];
    sm[1][1] = s[3];
    sm[2][2] = s[5];
    sm[0][1] = s[1];
    sm[1][0] = s[1];
    sm[0][2] = s[2];
    sm[2][0] = s[2];
    sm[1][2] = s[4];
    sm[2][1] = s[4];

    rt[0][0] = r[0][0];
    rt[1][1] = r[1][1];
    rt[2][2] = r[2][2];
    
    rt[0][1] = r[1][0];
    rt[1][0] = r[0][1];

    rt[0][2] = r[2][0];
    rt[2][0] = r[0][2];

    rt[1][2] = r[2][1];
    rt[2][1] = r[1][2];

    return mul(r,mul(sm,rt));
  }

  private void rotateMatrix(float tha, float[] v, float[] s, double[][] r) {
    float cti = cos(2.0f*tha);
    float sti = sin(2.0f*tha);
    float ctm = 1.0f-cti;
    float[] p = cross(v,s);

    float p1 = p[0];
    float p2 = p[1];
    float p3 = p[2];

    float p11 = p1*p1; 
    float p12 = p1*p2; 
    float p13 = p1*p3; 
    float p22 = p2*p2; 
    float p23 = p2*p3; 
    float p33 = p3*p3; 

    r[0][0] = cti+ctm*p11; 
    r[1][1] = cti+ctm*p22;
    r[2][2] = cti+ctm*p33;

    r[0][1] = ctm*p12-sti*p3; 
    r[1][0] = ctm*p12+sti*p3;

    r[0][2] = ctm*p13+sti*p2; 
    r[2][0] = ctm*p13-sti*p2; 

    r[1][2] = ctm*p23-sti*p1;
    r[2][1] = ctm*p23+sti*p1;
  }

  private float[] cross(float[] u, float[] v) {
    float u1 = u[0];
    float u2 = u[1];
    float u3 = u[2];
    float v1 = v[0];
    float v2 = v[1];
    float v3 = v[2];

    float uv1 = u2*v3-u3*v2;
    float uv2 = u3*v1-u1*v3;
    float uv3 = u1*v2-u2*v1;
    float uvs = 1.0f/sqrt(uv1*uv1+uv2*uv2+uv3*uv3);
    return new float[]{uv1*uvs,uv2*uvs,uv3*uvs};
  }

  
  private float[][] solveEigenproblems(double[][] a) {
    double[] e = new double[3];
    double[][] z = new double[3][3];
    Eigen.solveSymmetric33(a,z,e);
    float eui = (float)e[0];
    float evi = (float)e[1];
    float ewi = (float)e[2];
    if (ewi<0.0f) ewi = 0.0f;
    if (evi<ewi) evi = ewi;
    if (eui<evi) eui = evi;
    float u1i = (float)z[0][0];
    float u2i = (float)z[0][1];
    float u3i = (float)z[0][2];
    float[] es = new float[]{eui,evi,ewi};
    float[] us = new float[]{u1i,u2i,u3i};
    return new float[][]{us,es};
  }

  private static float distance(float[] x, float[] y) {
    float d1 = y[0]-x[0];
    float d2 = y[1]-x[1];
    float d3 = y[2]-x[2];
    return sqrt(d1*d1+d2*d2+d3*d3);
  }


  private static void getRange(float[] ds, float[] xs, int[] ns, 
    float[] xmin, float[] xmax) {
    int i1m = round(xs[0]-ds[0]); if(i1m<0){i1m=0;}
    int i2m = round(xs[1]-ds[1]); if(i2m<0){i2m=0;}
    int i3m = round(xs[2]-ds[2]); if(i3m<0){i3m=0;}
    int i1p = round(xs[0]+ds[0]); if(i1p>=ns[0]){i1p=ns[0]-1;}
    int i2p = round(xs[1]+ds[1]); if(i2p>=ns[1]){i2p=ns[1]-1;}
    int i3p = round(xs[2]+ds[2]); if(i3p>=ns[2]){i3p=ns[2]-1;}
    xmin[0] = i1m;
    xmin[1] = i2m;
    xmin[2] = i3m;
    xmax[0] = i1p;
    xmax[1] = i2p;
    xmax[2] = i3p;
  }


  private float[][][] setKdTreePoints() {
    int ik = 0;
    float[][][] xu = new float[2][3][_nc];
    for (int ic=0; ic<_nc; ic+=1) {
      xu[0][0][ik] = _fc[ic].i1;
      xu[0][1][ik] = _fc[ic].i2;
      xu[0][2][ik] = _fc[ic].i3;
      xu[1][0][ik] = _fc[ic].w1;
      xu[1][1][ik] = _fc[ic].w2;
      xu[1][2][ik] = _fc[ic].w3;
      ik++;
    }
    return xu;
  }


  private void initialTensor() {
    _st = new float[_nc][6];
    for (int ic=0; ic<_nc; ++ic) {
      float u1 = _fc[ic].w1;
      float u2 = _fc[ic].w2;
      float u3 = _fc[ic].w3;
      float u11 = u1*u1;
      float u12 = u1*u2;
      float u13 = u1*u3;
      float u22 = u2*u2;
      float u23 = u2*u3;
      float u33 = u3*u3;
      _st[ic][0] = u11;
      _st[ic][1] = u12;
      _st[ic][2] = u13;
      _st[ic][3] = u22;
      _st[ic][4] = u23;
      _st[ic][5] = u33;
    }
  }

  private int[] getBounds(int n, float[] x) {
    int[] bs = new int[2];
    int n1m = (int)min(x)-5; 
    int n1p = (int)max(x)+5; 
    if(n1m<0){n1m=0;}
    if(n1p>n){n1p=n;}
    bs[0] = n1m;
    bs[1] = n1p;
    return bs;
  }


  private int _n1;
  private int _n2;
  private int _n3;
  private int _nc;
  private float[][] _st;
  private FaultCell[] _fc;

}
