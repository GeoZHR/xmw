/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package mef;

import java.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Set weights and constraints for flattening: 
 * set zero weights at fault, 
 * set hard constraints using fault slip vectors.
 * @author Xinming Wu
 * @version 2014.02.09
 */

public class FaultSlipConstraints2 {

  public FaultSlipConstraints2(FaultCurve[] cvs) {
    _cvs = cvs;
  }


  public float[][][] screenPointsX(float[][] wp){
    int n2 = wp.length;
    int n1 = wp[0].length;
    computeUnfaultShifts(_cvs);
    setCells(n1,n2);
    flNormalization();
    ArrayList<float[][]> cl = new ArrayList<float[][]>();
    for (FaultCurve sk:_cvs) {
      for (FaultPoint fp:sk) {
        float fl = fp.getFl();
        float[] h1 = new float[2];
        float[] f1 = new float[2];
        float[] c1 = new float[2];
        float[] h2 = new float[2];
        float[] f2 = new float[2];
        int i1i  = fp.i1;
        int i2i  = fp.i2;
        int i2m1 = bound(fp.i2m,n2);
        int i2p1 = bound(fp.i2p,n2);
        FaultPoint fc = fp.fc;
        int c1i = fc.i1;
        int c2i = fc.i2;
        int c2m1 = bound(fc.i2m,n2);
        int c2p1 = bound(fc.i2p,n2);

        f1[0] = i1i;
        f1[1] = i2m1;
        h1[0] = c1i;  
        h1[1] = c2p1; 
        wp[i2m1][i1i] = 0.0f;
        wp[i2p1][i1i] = 0.0f;

        int d2mi = i2m1-i2i;
        int d2pi = i2p1-i2i;
        int d2mc = c2m1-c2i;
        int d2pc = c2p1-c2i;

        int i2m2 = i2m1;
        int c2p2 = c2p1; 
        if(abs(d2mi)>0){i2m2+=d2mi;}
        else           {i2m2-=d2pi;}
        i2m2=bound(i2m2,n2); 

        if(abs(d2mc)>0){c2p2-=d2mc;}
        else           {c2p2+=d2pc;}
        c2p2=bound(c2p2,n2); 

        h2[0] = c1i;  f2[0] = i1i;
        h2[1] = c2p2; f2[1] = i2m2;

        float s1 = fp.s1;
        float s2 = fp.s2;
        float[] ch = new float[2];
        float[] cf = new float[2];
        float[] ts = new float[]{s1,s2};
        c1[0] = fl; c1[1] = fl;
        cl.add(new float[][]{h1,f1,ts,c1});
        cl.add(new float[][]{h2,h1,ch,c1});
        cl.add(new float[][]{f2,f1,cf,c1});
      }
    }
    int ns = cl.size();
    System.out.println("sets of control points:"+ns);
    float[][][] cs = new float[4][3][ns];
    for (int is=0; is<ns; ++is) {
      float[][] ps = cl.get(is);
      cs[0][0][is] = ps[0][0];
      cs[0][1][is] = ps[0][1];

      cs[1][0][is] = ps[1][0];
      cs[1][1][is] = ps[1][1];

      cs[2][0][is] = ps[2][0];
      cs[2][1][is] = ps[2][1];

      cs[3][0][is] = ps[3][0];
      cs[3][1][is] = ps[3][1];
    }
    return cs;
  }


  public float[][][] screenPoints(float[][] wp){
    int n2 = wp.length;
    int n1 = wp[0].length;
    computeUnfaultShifts(_cvs);
    setCells(n1,n2);
    flNormalization();
    ArrayList<float[][]> cl = new ArrayList<float[][]>();
    for (FaultCurve sk:_cvs) {
      for (FaultPoint fp:sk) {
        float fl = fp.getFl();
        float[] h1 = new float[2];
        float[] f1 = new float[2];
        float[] c1 = new float[2];
        float[] h2 = new float[2];
        float[] f2 = new float[2];
        int i1i  = fp.i1;
        int i2i  = fp.i2;
        int i2m1 = bound(fp.i2m,n2);
        int i2p1 = bound(fp.i2p,n2);

        h1[0] = i1i;  f1[0] = i1i;
        h1[1] = i2p1; f1[1] = i2m1;
        wp[i2m1][i1i] = 0.0f;
        wp[i2p1][i1i] = 0.0f;

        int d2m = i2m1-i2i;
        int d2p = i2p1-i2i;
        int i2m2 = i2m1, i2p2=i2p1; 
        if(abs(d2m)>0){i2p2-=d2m;i2m2+=d2m;}
        else          {i2m2-=d2p;i2p2+=d2p;}
        i2m2=bound(i2m2,n2); 
        i2p2=bound(i2p2,n2);
        h2[0] = i1i;  f2[0] = i1i;
        h2[1] = i2p2; f2[1] = i2m2;

        float t1 = fp.t1;
        float t2 = fp.t2;
        float[] ch = new float[2];
        float[] cf = new float[2];
        float[] ts = new float[]{t1,t2};
        c1[0] = fl; c1[1] = fl;
        cl.add(new float[][]{h1,f1,ts,c1});
        cl.add(new float[][]{h2,h1,ch,c1});
        cl.add(new float[][]{f2,f1,cf,c1});
      }
    }
    int ns = cl.size();
    System.out.println("sets of control points:"+ns);
    float[][][] cs = new float[4][3][ns];
    for (int is=0; is<ns; ++is) {
      float[][] ps = cl.get(is);
      cs[0][0][is] = ps[0][0];
      cs[0][1][is] = ps[0][1];

      cs[1][0][is] = ps[1][0];
      cs[1][1][is] = ps[1][1];

      cs[2][0][is] = ps[2][0];
      cs[2][1][is] = ps[2][1];

      cs[3][0][is] = ps[3][0];
      cs[3][1][is] = ps[3][1];
    }
    return cs;
  }

  private void computeUnfaultShifts(FaultCurve[] curves) {
    int ik = 0;
    for (FaultCurve curve:curves) {
      ik++;
      System.out.println("ik="+ik);
      final FaultPoint[] points = curve.getPoints();
      final int nc = points.length;
      for (int ic=0; ic<nc; ++ic) {
        FaultPoint fci = points[ic];
        FaultPoint fcd = points[ic];
        FaultPoint fcu = points[ic];
        float s1 = fci.s1;
        float s2 = fci.s2;
        float[] yd = new float[2];
        float[] yu = new float[2];
        yd[0] = fci.x1;
        yd[1] = fci.x2;
        yu[0] = fci.x1;
        yu[1] = fci.x2;
        float st = round(s1);
        if (st>0) {
          for (; st>=1.0f; st-=1.0f) {
            fcd = fcd.walkDownDipFrom(yd);
            fcu = fcu.walkUpDipFrom(yu);
          }
        } else {
          for (; st<=-1.0f; st+=1.0f) {
            fcu = fcu.walkDownDipFrom(yu);
            fcd = fcd.walkUpDipFrom(yd);
          }
        }
        float du = s1-fcu.s1;
        float dd = fcd.s1-s1;
        s1 -= (du+dd)*0.5f;
        fci.setUnfaultShifts(new float[]{s1,s2});
      }
    }
  }


  private int bound(int i, int n) {
    if(i<0) {i=0;}
    if(i>=n){i=n-1;}
    return i;
  }

  private void flNormalization() {
    for (FaultCurve curve:_cvs){
      FaultPoint[] fcs = FaultCurve.getPoints(new FaultCurve[]{curve});
      int nc = fcs.length;
      float[] fls = new float[nc];
      for (int ic=0; ic<nc; ++ic) {
        fls[ic] = fcs[ic].getFl();
      }
      div(fls,max(fls),fls);
      for (int ic=0; ic<nc; ++ic) {
        fcs[ic].setFl(fls[ic]);
      }
    }
  }

  private void setCells(int n1, int n2) {
    FaultPoint[][] points = new FaultPoint[n2][n1];
    for (FaultCurve cv:_cvs) {
      for (FaultPoint fp:cv) {
        int i1i = fp.i1;
        int i2i = fp.i2;
        int i2m1 = bound(fp.i2m,n2);
        int i2p1 = bound(fp.i2p,n2);

        int d2m = i2m1-i2i;
        int d2p = i2p1-i2i;
        int i2m2 = i2m1, i2p2=i2p1; 
        if(abs(d2m)>0){i2p2-=d2m;i2m2+=d2m;}
        else          {i2m2-=d2p;i2p2+=d2p;}
        i2m2 = bound(i2m2,n2);
        i2p2 = bound(i2p2,n2);
        FaultPoint cm1 = points[i2m1][i1i];
        FaultPoint cp1 = points[i2p1][i1i];
        FaultPoint cm2 = points[i2m2][i1i];
        FaultPoint cp2 = points[i2p2][i1i];
        FaultPoint[] cls = new FaultPoint[]{cm1,cp1,cm2,cp2};
        int[] i2s = new int[]{i2m1,i2p1,i2m2,i2p2};
        int ik = -1;
        int mc = cls.length;
        for (int im=0; im<mc; ++im) {
          FaultPoint point=cls[im];
          ik ++;
          int i1 = i1i;
          int i2 = i2s[ik];
          if(point==null) {
            points[i2][i1] = fp;
          } else {
            float s1i = abs(fp.t1);
            float s1c = abs(point.t1);
            if(s1i>s1c) {
              point.setFl(0.0f);
              points[i2][i1] = fp;
            }else{fp.setFl(0.0f);}
          }
        }
      }
    }
  }

  private FaultCurve[] _cvs;
}

