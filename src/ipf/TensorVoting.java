/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ipf;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.sgl.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Computes fault blocks. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2014.12.16
 */
public class TensorVoting {

  public TensorVoting (int n1, int n2, int n3, FaultCell[] fc) {
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
    _fc = fc;
    _nc = fc.length;
  }

  public float[][][][] tensorVoteS(
    final float[][][] fl, final float[][][] fp, final float[][][] ft) {
    final float[][][] fpp = copy(fp);
    final float[][][] ftp = copy(ft);
    final FaultGeometry fg = new FaultGeometry();
    final float[][] st = initialTensor();
    final float[][][] xu = setKdTreePoints();
    final KdTree kt = new KdTree(xu[0]);
    final int[] ds = new int[]{15,10,10};
    final int[] ns = new int[]{_n1,_n2,_n3};
    final int[] bs1 = getBounds(_n1,xu[0][0]);
    final int[] bs2 = getBounds(_n2,xu[0][1]);
    final int[] bs3 = getBounds(_n3,xu[0][2]);
    final float[][][] sm =  new float[_n3][_n2][_n1];
    final float[][][] cm =  new float[_n3][_n2][_n1];
    final float[][][] jm =  new float[_n3][_n2][_n1];
    Parallel.loop(bs3[0],bs3[1],1,new Parallel.LoopInt() {
    public void compute(int i3) {
    //for (int i3=bs3[0]; i3<bs3[1]; ++i3) {
      System.out.println("i3="+i3);
      float[] xmin = new float[3];
      float[] xmax = new float[3];
      for (int i2=bs2[0]; i2<bs2[1]; ++i2) {
        for (int i1=bs1[0]; i1<bs1[1]; ++i1) {
          int[] is = new int[]{i1,i2,i3};
          float[] y = new float[]{i1,i2,i3};
          int ne = kt.findNearest(y);
          float x1 = xu[0][0][ne];
          float x2 = xu[0][1][ne];
          float x3 = xu[0][2][ne];
          float dd = distance(new float[]{x1,x2,x3},y);
          if(dd>15.0f){continue;}
          getRange(ds,is,ns,xmin,xmax);
          int[] id = kt.findInRange(xmin,xmax);
          int nd = id.length;
          if(nd<1) {continue;}
          double[][] a = new double[3][3];
          float fpa = 0.0f;
          float fta = 0.0f;
          float pts = 0.0f;

          for (int ik=0; ik<nd; ++ik) {
            int ic = id[ik];
            int j1 = (int)xu[0][0][ic];
            int j2 = (int)xu[0][1][ic];
            int j3 = (int)xu[0][2][ic];
            float[] v = new float[]{i1-j1,i2-j2,i3-j3};
            float gn = fl[j3][j2][j1];
            if(sum(abs(v))==0.0f) {
              a[0][0] += gn*st[ic][0];
              a[1][1] += gn*st[ic][3];
              a[2][2] += gn*st[ic][5];
              a[0][1] += gn*st[ic][1];
              a[1][0] += gn*st[ic][1];
              a[0][2] += gn*st[ic][2];
              a[2][0] += gn*st[ic][2];
              a[1][2] += gn*st[ic][4];
              a[2][1] += gn*st[ic][4];
              fpa += fpp[j3][j2][j1];
              fta += ftp[j3][j2][j1];
              pts += 1.0f;
            } else {
              float w1 = xu[1][0][ic];
              float w2 = xu[1][1][ic];
              float w3 = xu[1][2][ic];
              float[] t = new float[1];
              float[] w = new float[]{w1,w2,w3};
              float sc = computeScale(t,v,w);
              if(sc==0.0f){continue;}
              if(Float.isNaN(sc)){continue;}
              sc *= gn;
              float u1 = w1+v[0]*t[0];
              float u2 = w2+v[1]*t[0];
              float u3 = w3+v[2]*t[0];
              float u11 = u1*u1;
              float u12 = u1*u2;
              float u13 = u1*u3;
              float u22 = u2*u2;
              float u23 = u2*u3;
              float u33 = u3*u3;
              float uss = 1.0f/(u11+u22+u33);
              a[0][0] += sc*u11*uss;
              a[0][1] += sc*u12*uss;
              a[0][2] += sc*u13*uss;

              a[1][0] += sc*u12*uss;
              a[1][1] += sc*u22*uss;
              a[1][2] += sc*u23*uss;

              a[2][0] += sc*u13*uss;
              a[2][1] += sc*u23*uss;
              a[2][2] += sc*u33*uss;
              int k1 = _fc[ne].i1;
              int k2 = _fc[ne].i2;
              int k3 = _fc[ne].i3;
              float fpr = fpp[k3][k2][k1];
              float fpi = fpp[j3][j2][j1];
              if(abs(fpr-fpi)<40) {
                fpa += fpp[j3][j2][j1];
                fta += ftp[j3][j2][j1];
                pts += 1.0f;
              }
            }
          }
          float[][] ue = solveEigenproblems(a);
          float u1 = ue[0][0];
          float u2 = ue[0][1];
          float u3 = ue[0][2];
          if(u2==0.0f&&u3==0.0f){continue;}
          sm[i3][i2][i1] = (ue[1][0]-ue[1][1]);
          cm[i3][i2][i1] = (ue[1][1]-ue[1][2]);
          jm[i3][i2][i1] = (ue[1][2]);
          /*
          double fpo = fpp[i3][i2][i1];
          double fto = ftp[i3][i2][i1];
          float[] uo = fg.faultNormalVectorFromStrikeAndDip(fpo,fto);
          if((u1*uo[0]+u2*uo[1]+u3*uo[2])<0.0f){
            u1=-u1;
            u2=-u2;
            u3=-u3;
          }
          ft[i3][i2][i1] = fg.faultDipFromNormalVector(u1,u2,u3);
          fp[i3][i2][i1] = fg.faultStrikeFromNormalVector(u1,u2,u3);
          */
          ft[i3][i2][i1] = fta/pts;
          fp[i3][i2][i1] = fpa/pts;
        }
      }
      //}
    }});
    return new float[][][][]{sm,cm,jm};
  }


  public float[][][][] tensorVote(
    final float[][][] fl, final float[][][] fp, final float[][][] ft) {
    final float[][][] fpp = copy(fp);
    final float[][][] ftp = copy(ft);
    final FaultGeometry fg = new FaultGeometry();
    final float[][] st = initialTensor();
    final float[][][] xu = setKdTreePoints();
    final KdTree kt = new KdTree(xu[0]);
    final int[] ds = new int[]{10,10,10};
    final int[] ns = new int[]{_n1,_n2,_n3};
    final int[] bs1 = getBounds(_n1,xu[0][0]);
    final int[] bs2 = getBounds(_n2,xu[0][1]);
    final int[] bs3 = getBounds(_n3,xu[0][2]);
    final float[][][] sm =  new float[_n3][_n2][_n1];
    final float[][][] cm =  new float[_n3][_n2][_n1];
    final float[][][] jm =  new float[_n3][_n2][_n1];
    Parallel.loop(bs3[0],bs3[1],1,new Parallel.LoopInt() {
    public void compute(int i3) {
    //for (int i3=bs3[0]; i3<bs3[1]; ++i3) {
      System.out.println("i3="+i3);
      float[] xmin = new float[3];
      float[] xmax = new float[3];
      for (int i2=bs2[0]; i2<bs2[1]; ++i2) {
        for (int i1=bs1[0]; i1<bs1[1]; ++i1) {
          int[] is = new int[]{i1,i2,i3};
          float[] y = new float[]{i1,i2,i3};
          int ne = kt.findNearest(y);
          float x1 = xu[0][0][ne];
          float x2 = xu[0][1][ne];
          float x3 = xu[0][2][ne];
          float dd = distance(new float[]{x1,x2,x3},y);
          if(dd>20.0f){continue;}
          getRange(ds,is,ns,xmin,xmax);
          int[] id = kt.findInRange(xmin,xmax);
          int nd = id.length;
          if(nd<1) {continue;}
          double[][] a = new double[3][3];
          double[][] r = new double[3][3];
          float fpa = 0.0f;
          float fta = 0.0f;
          float pts = 0.0f;
          for (int ik=0; ik<nd; ++ik) {
            int ic = id[ik];
            int j1 = (int)xu[0][0][ic];
            int j2 = (int)xu[0][1][ic];
            int j3 = (int)xu[0][2][ic];
            float[] v = new float[]{i1-j1,i2-j2,i3-j3};
            float gn = fl[j3][j2][j1];
            if(sum(abs(v))==0.0f) {
              a[0][0] += gn*st[ic][0]; //s11
              a[0][1] += gn*st[ic][1]; //s12
              a[0][2] += gn*st[ic][2]; //s13
              a[1][0] += gn*st[ic][1]; //s12
              a[1][1] += gn*st[ic][3]; //s22
              a[1][2] += gn*st[ic][4]; //s23
              a[2][0] += gn*st[ic][2]; //s13
              a[2][1] += gn*st[ic][4]; //s23
              a[2][2] += gn*st[ic][5]; //s33
            } else {
              float[] t = new float[1];
              float w1 = _fc[ic].w1;
              float w2 = _fc[ic].w2;
              float w3 = _fc[ic].w3;
              float[] w = new float[]{w1,w2,w3};
              float sc = computeScale(t,v,st[ic],w,r);
              if(sc==0.0f){continue;}
              if(Float.isNaN(sc)){continue;}
              sc *= gn;
              double[][] sr = computeRSR(st[ic],r);
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
            int k1 = _fc[ne].i1;
            int k2 = _fc[ne].i2;
            int k3 = _fc[ne].i3;
            float fpr = fpp[k3][k2][k1];
            float fpi = fpp[j3][j2][j1];
            if(abs(fpr-fpi)<40) {
              fpa += fpp[j3][j2][j1];
              fta += ftp[j3][j2][j1];
              pts += 1.0f;
            }

          }
          //float[][] ue = solveEigenproblems(mul(a,1.0/cot));
          float[][] ue = solveEigenproblems(a);
          float u1 = ue[0][0];
          float u2 = ue[0][1];
          float u3 = ue[0][2];
          if(u2==0.0f&&u3==0.0f){continue;}
          sm[i3][i2][i1] = (ue[1][0]-ue[1][1]);
          cm[i3][i2][i1] = (ue[1][1]-ue[1][2]);
          jm[i3][i2][i1] = (ue[1][2]);
          ft[i3][i2][i1] = fta/pts;
          fp[i3][i2][i1] = fpa/pts;
          /*
          double fpo = fpp[i3][i2][i1];
          double fto = ftp[i3][i2][i1];
          float[] uo = fg.faultNormalVectorFromStrikeAndDip(fpo,fto);
          if(u1*uo[0]<0.0f){u1=-u1;}
          if(u2*uo[1]<0.0f){u2=-u2;}
          if(u3*uo[2]<0.0f){u3=-u3;}

          ft[i3][i2][i1] = fg.faultDipFromNormalVector(u1,u2,u3);
          fp[i3][i2][i1] = fg.faultStrikeFromNormalVector(u1,u2,u3);
          */
        }
      }
      //}
    }});
    return new float[][][][]{sm,cm,jm};
  }

  private float computeScale(float[] k, float[] v, float[] w) {
    float b = 0.0f;
    float sigma = 20.0f;
    float v1 = v[0];
    float v2 = v[1];
    float v3 = v[2];
    float w1 = w[0];
    float w2 = w[1];
    float w3 = w[2];
    float wv = abs(w1*v1+w2*v2+w3*v3);
    float vs = sqrt(v1*v1+v2*v2+v3*v3);
    float sin = sqrt(wv/vs);
    float tha = asin(sin);
    float arc = vs*tha/sin;
    if(sin==0.0f){arc=vs;}
    float cur = 2.0f*sin/vs;
    k[0] = cur;
    if(tha>Math.PI/4.0f) {
      return 0.0f;
    } else {
      float ep = (arc*arc+b*cur*cur)/(sigma*sigma);
      return exp(-ep);
    }
  }


  private float computeScale(float[] k, float[] v, float[] s, float[] w, double[][] r) {
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


  private static void getRange(int[] ds, int[] is, int[] ns, 
    float[] xmin, float[] xmax) {
    int i1m = is[0]-ds[0]; if(i1m<0){i1m=0;}
    int i2m = is[1]-ds[1]; if(i2m<0){i2m=0;}
    int i3m = is[2]-ds[2]; if(i3m<0){i3m=0;}
    int i1p = is[0]+ds[0]; if(i1p>=ns[0]){i1p=ns[0]-1;}
    int i2p = is[1]+ds[1]; if(i2p>=ns[1]){i2p=ns[1]-1;}
    int i3p = is[2]+ds[2]; if(i3p>=ns[2]){i3p=ns[2]-1;}
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


  private float[][] initialTensor() {
    float[][] st = new float[_nc][6];
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
      st[ic][0] = u11;
      st[ic][1] = u12;
      st[ic][2] = u13;
      st[ic][3] = u22;
      st[ic][4] = u23;
      st[ic][5] = u33;
    }
    return st;
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
  private FaultCell[] _fc;

}
