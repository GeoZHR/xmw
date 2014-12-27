/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package ifs;

import java.util.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Computes fault blocks. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2014.09.16
 */
public class PointSetSurface {

  public float[][][] findScalarField(int n1, int n2, int n3, FaultCell[] fc) {
    return scalarField(n1,n2,n3,fc);
  }

  public FaultCell[] getCells(FaultSkin[] fs) {
    FaultSkin[] fsi = new FaultSkin[2];
    fsi[0] = fs[2];
    fsi[1] = fs[3];
    return FaultSkin.getCells(fs);
  }

  private float[][][] setKdTreePoints(FaultCell[] fc) {
    int nc = fc.length;
    float[][][] xu = new float[2][3][nc];
    int ik = 0;
    for (int ic=0; ic<nc; ic+=1) {
      xu[0][0][ik] = fc[ic].i1;
      xu[0][1][ik] = fc[ic].i2;
      xu[0][2][ik] = fc[ic].i3;
      xu[1][0][ik] = fc[ic].w1;
      xu[1][1][ik] = fc[ic].w2;
      xu[1][2][ik] = fc[ic].w3;
      ik++;
    }
    return xu;
  }

  public void alignNormals(int n1, int n2, int n3, FaultCell[] fc) {
    int nc = fc.length;
    int ib = startCell(fc);
    int[] mk = new int[nc];
    float[][][] xu = setKdTreePoints(fc);
    HashSet<Integer> idh = new HashSet<Integer>();
    idh.add(ib);
    float[] xmin = new float[3];
    float[] xmax = new float[3];
    int[] ds = new int[]{5,5,5};
    int[] ns = new int[]{n1,n2,n3};
    KdTree kt = new KdTree(xu[0]);
    while (idh.size()>0) {
      for (int ic:idh) {
        mk[ic] = 1;
        idh.remove(ic);
        int[] is = fc[ic].getI();
        float[] us = fc[ic].getW();
        getRange(ds,is,ns,xmin,xmax);
        int[] ids = kt.findInRange(xmin,xmax);
        int nd = ids.length;
        if (nd<0) {continue;}
        for (int ik=0; ik<nd; ++ik) {
          int id = ids[ik];
          if(mk[id]==0) {
            float x1 = fc[id].i1-is[0];
            float x2 = fc[id].i2-is[1];
            float x3 = fc[id].i3-is[2];
            float sx = x1*x1+x2*x2+x3*x3;
            if(sx==0.0f){continue;}
            idh.add(id);
            float[] vs = fc[id].getW();
            float[] xs = new float[]{x1/sx,x2/sx,x3/sx};
            float[] ur = reflectionVector(xs,us);
            if(sum(mul(vs,ur))<0.0f) {
              fc[id].setNormal(-vs[0],-vs[1],-vs[2]);
            }
          }
        }
      }
    }
  }

  private float[] reflectionVector(float[] xi, float[] ui) {
    float u1 = ui[0];
    float u2 = ui[1];
    float u3 = ui[2];
    float x1 = xi[0];
    float x2 = xi[1];
    float x3 = xi[2];

    float ux = u1*x1+u2*x2+u3*x3;
    float w1 = u1-x1;
    float w2 = u2-x2;
    float w3 = u3-x3;
    float ws = 1.0f/(w1*w1+w2*w2+w3*w3);
    w1 *= ws; w2 *= ws; w3 *= ws;
    float wu = 2.0f*(w1*u1+w2*u2+w3*u3);
    float v1 = wu*w1-u1;
    float v2 = wu*w2-u2;
    float v3 = wu*w3-u3;
    return new float[]{v1,v2,v3};
  }

  private int startCell(FaultCell[] fc) {
    int id = 0;
    float flMax = 0.0f;
    int nc = fc.length;
    for (int ic=0; ic<nc; ++ic) {
      float fl = fc[ic].fl;
      if(fl>flMax) {id = ic; flMax = fl;}
    }
    return id;
  }


  public float[][][] scalarField(
    final int n1, final int n2, final int n3, FaultCell[] fc) 
  {
    final int d2=10;
    final int d3=10;
    final int d1=15;
    final float v = -30.f;
    final float[][][] xu = setKdTreePoints(fc);
    final int[] bs1 = setBounds(n1,xu[0][0]);
    final int[] bs2 = setBounds(n2,xu[0][1]);
    final int[] bs3 = setBounds(n3,xu[0][2]);
    final KdTree kt = new KdTree(xu[0]);
    final float[][][] sf = fillfloat(v,n1,n2,n3);
    Parallel.loop(bs3[0],bs3[1],1,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[] xmin = new float[3];
      float[] xmax = new float[3];
      System.out.println("i3="+i3);
      for (int i2=bs2[0]; i2<bs2[1]; ++i2) {
        for (int i1=bs1[0]; i1<bs1[1]; ++i1) {
          float[] y = new float[]{i1,i2,i3};
          int ne = kt.findNearest(y);
          float x1 = xu[0][0][ne];
          float x2 = xu[0][1][ne];
          float x3 = xu[0][2][ne];
          float dd = distance(new float[]{x1,x2,x3},y);
          if(dd>10.0f){continue;}
          getRange(d1,d2,d3,i1,i2,i3,n1,n2,n3,xmin,xmax);
          int[] id = kt.findInRange(xmin,xmax);
          int nd = id.length;
          if(nd<d2*d1/2){continue;}
          float[][] xf = new float[3][nd];
          float[][] uf = new float[3][nd];
          for (int ik=0; ik<nd; ++ik) {
            int ip = id[ik];
            xf[0][ik] = xu[0][0][ip];
            xf[1][ik] = xu[0][1][ip];
            xf[2][ik] = xu[0][2][ip];
            uf[0][ik] = xu[1][0][ip];
            uf[1][ik] = xu[1][1][ip];
            uf[2][ik] = xu[1][2][ip];
          }
          float[] r = scalarField(1.0f,y,xf,uf);
          sf[i3][i2][i1] = r[0];
        }
      }
    }});
    //cleanScalarField(sf,mk);
    return sf;
  }

  private int[] setBounds(int n, float[] x) {
    int[] bs = new int[2];
    int n1m = (int)min(x)-5; 
    int n1p = (int)max(x)+5; 
    if(n1m<0){n1m=0;}
    if(n1p>n){n1p=n;}
    bs[0] = n1m;
    bs[1] = n1p;
    return bs;
  }


  private void cleanScalarFieldM(float[][][] sf, float[][][] mk) {
    int n3 = sf.length;
    int n2 = sf[0].length;
    int n1 = sf[0][0].length;
    int nk = round(sum(mk));
    float[][] xf = new float[3][n3*n2*n1-nk];
    int ik = 0;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float mki = mk[i3][i2][i1];
          if(mki==0.0f) {
            xf[0][ik] = i1;
            xf[1][ik] = i2;
            xf[2][ik] = i3;
            ik++;
          }
        }
      }
    }
    KdTree kt = new KdTree(xf);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float mki = mk[i3][i2][i1];
          if(mki==1.0f) {
            int id = kt.findNearest(new float[]{i1,i2,i3});
            int x1 = round(xf[0][id]);
            int x2 = round(xf[1][id]);
            int x3 = round(xf[2][id]);
            sf[i3][i2][i1] = sf[x3][x2][x1];
          }
        }
      }
    }

  }

  private float scalarFieldS(float h, float[] y,float[][] xp, float[][] up) {
    float c = 0.0f;
    int np = xp[0].length;
    //float[] wp = computeWeights(h,dx,y,xp);
    float[] wp = fillfloat(1.0f,np);
    for (int ip=0; ip<np; ++ip) {
      float y1 = y[0];
      float y2 = y[1];
      float y3 = y[2];
      float x1 = xp[0][ip];
      float x2 = xp[1][ip];
      float x3 = xp[2][ip];
      float u1 = up[0][ip];
      float u2 = up[1][ip];
      float u3 = up[2][ip];
      float wi = wp[ip];
      c += wi*((y1-x1)*u1+(y2-x2)*u2+(y3-x3)*u3);
      c += ((y1-x1)*u1+(y2-x2)*u2+(y3-x3)*u3);
    }
    return c/sum(wp);
  }

  private void cleanScalarField(float[][][] sf, float[][][] mk) {
    int n3 = sf.length;
    int n2 = sf[0].length;
    int n1 = sf[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          if(mk[i3][i2][i1]==0.0f)
            continue;
          int cot = 0;
          float sfv = 0.0f;
          for(int d=1; d<5; ++d) {
            for(int d3=0; d3<d; ++d3) {
              for(int d2=0; d2<d; ++d2) {
                for(int d1=0; d1<d; ++d1) {
                  int i1p = i1+d1; if(i1p>=n1){i1p=n1-1;}
                  int i2p = i2+d2; if(i2p>=n2){i2p=n2-1;}
                  int i3p = i3+d3; if(i3p>=n3){i3p=n3-1;}
                  int i1m = i1-d1; if(i1m<0){i1m=0;}
                  int i2m = i2-d2; if(i2m<0){i2m=0;}
                  int i3m = i3-d3; if(i3m<0){i3m=0;}
                  if(mk[i3m][i2m][i1m]==0.0f) {
                    cot++;
                    sfv += sf[i3m][i2m][i1m];
                  }
                  if(mk[i3p][i2p][i1p]==0.0f) {
                    cot++;
                    sfv += sf[i3p][i2p][i1p];
                  }
                }
                if(cot>2) {
                  sf[i3][i2][i1] = sfv/cot;
                  break;
                }
              }
            }
          }
        }
      }
    }
  }

  private static float[] scalarField(float h, float[] xp, float[][] xf, float[][] uf) {
    int np = xf[0].length;
    float[] r = new float[2];
    float[] wp = fillfloat(1.0f,np);
    //float[] wp = computeWeights(h,dx,xp,xf);
    float[] wn = normalizeWeights(wp);
    float x0 = 1.0f;
    float x1 = xp[0];
    float x2 = xp[1];
    float x3 = xp[2];
    float x4 = x1*x1+x2*x2+x3*x3;
    float beta = 1.0f;
    float c1 = 0.0f;
    float d1 = 0.0f;
    float b1 = 0.0f;
    float[] wpu = new float[3];
    float[] wnu = new float[3];
    float[] wpx = new float[3];
    float[] wnx = new float[3];
    for (int ip =0; ip<np; ++ip) {
      float wpi = wp[ip];
      float wni = wn[ip];
      float xi1 = xf[0][ip];
      float xi2 = xf[1][ip];
      float xi3 = xf[2][ip];
      float ui1 = uf[0][ip];
      float ui2 = uf[1][ip];
      float ui3 = uf[2][ip];
      float xsi = xi1*xi1+xi2*xi2+xi3*xi3;
      d1 += wpi*xsi;
      b1 += wni*xsi;
      c1 += wpi*(xi1*ui1+xi2*ui2+xi3*ui3); 
      wnx[0] += wni*xi1;
      wnx[1] += wni*xi2;
      wnx[2] += wni*xi3;

      wpx[0] += wpi*xi1;
      wpx[1] += wpi*xi2;
      wpx[2] += wpi*xi3;

      wnu[0] += wni*ui1;
      wnu[1] += wni*ui2;
      wnu[2] += wni*ui3;

      wpu[0] += wpi*ui1;
      wpu[1] += wpi*ui2;
      wpu[2] += wpi*ui3;
    }
    float c2 = wnx[0]*wpu[0]+wnx[1]*wpu[1]+wnx[2]*wpu[2];
    float d2 = wnx[0]*wpx[0]+wnx[1]*wpx[1]+wnx[2]*wpx[2];
    float u4 = beta*(c1-c2)/(d1-d2);
    if(np==1){u4=0.0f;}
    if(d1-d2==0.0f||abs(u4)>1.0f){
      r[1] = 1.0f;
    }
    float u1 = wnu[0]-u4*wnx[0];
    float u2 = wnu[1]-u4*wnx[1];
    float u3 = wnu[2]-u4*wnx[2];
    u4 *= 0.5f;
    float u0 = -u1*wnx[0]-u2*wnx[1]-u3*wnx[2]-u4*b1;
    r[0] = x0*u0+x1*u1+x2*u2+x3*u3+x4*u4;
    return r;
  }


  private static float[] normalizeWeights(float[] wp) {
    int np = wp.length;
    float sum = 1.0f/sum(wp);
    float[] wn = new float[np];
    for (int ip=0; ip<np; ++ip) {
      wn[ip] = wp[ip]*sum;
    }
    return wn;
  }

  private static float[] computeWeights(float h, float[] dx, float[] xp, float[][] xf) {
    int np = xf[0].length;
    float[] wp = new float[np];
    MedianFinder mf = new MedianFinder(np);
    float md = mf.findMedian(abs(dx));
    md *= h;
    for (int ip=0; ip<np; ++ip) {
      float x = abs(dx[ip])/md;
      if(x<1.0f) {
        float wi = pow(x,2.f);
        wp[ip] = pow((1.0f-wi),4.0f);
      } else {
        wp[ip] = 0.0f;
      }
    }
    return wp;
  }

  private static float distance(float[] x, float[] y) {
    float d1 = y[0]-x[0];
    float d2 = y[1]-x[1];
    float d3 = y[2]-x[2];
    return sqrt(d1*d1+d2*d2+d3*d3);
  }

  private static float distance(float[] x, float[] u, float[] y) {
    float d1 = y[0]-x[0];
    float d2 = y[1]-x[1];
    float d3 = y[2]-x[2];
    float du = d1*u[0]+d2*u[1]+d3*u[2];
    return du;
  }


  private static void getRange(int d1, int d2, int d3, int i1, int i2, int i3, 
    int n1, int n2, int n3, float[] xmin, float[] xmax) 
  {
    int i1m = i1-d1; if(i1m<0){i1m=0;}
    int i2m = i2-d2; if(i2m<0){i2m=0;}
    int i3m = i3-d3; if(i3m<0){i3m=0;}
    int i1p = i1+d1; if(i1p>=n1){i1p=n1-1;}
    int i2p = i2+d2; if(i2p>=n2){i2p=n2-1;}
    int i3p = i3+d3; if(i3p>=n3){i3p=n3-1;}
    xmin[0] = i1m;
    xmin[1] = i2m;
    xmin[2] = i3m;
    xmax[0] = i1p;
    xmax[1] = i2p;
    xmax[2] = i3p;
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



}
