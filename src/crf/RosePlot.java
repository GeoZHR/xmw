/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package crf;

import java.awt.*;
import java.util.*;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.ArrayMath.*;

import mef.*;

/**
 * Make a rosette-strike plot for faults/fault cells.
 * @author Xinming Wu, University of Texas at Austin.
 * @version 2016.05.09
 */

public class RosePlot {

  public float[][] removeSignature(float phi, float[][] fps) {
    int np = fps[0].length;
    float[][] fpc = new float[4][np];
    int ik = 0;
    for (int ip=0; ip<np; ++ip) {
      float fpi = fps[3][ip];
      if(abs(fpi-phi)>10f) {
        fpc[0][ik] = fps[0][ip];
        fpc[1][ik] = fps[1][ip];
        fpc[2][ik] = fps[2][ip];
        fpc[3][ik] = fps[3][ip];
        ik++;
      }
    }
    return copy(ik,4,0,0,fpc);
  }

  public void rotate(float phi, float[] fps) {
    int np = fps.length;
    for (int ip=0; ip<np; ++ip) {
      float fpi = fps[ip];
      if(fpi>=0.0f) {
        fpi += phi;
        if (fpi>=360f) fpi-=360f;
        fps[ip] = fpi;
      }
    }
  }

  public void rotateX(float phi, float[] fps) {
    int np = fps.length;
    for (int ip=0; ip<np; ++ip) {
      float fpi = fps[ip];
      if(fpi>=0.0f) {
        fpi = phi-fpi;
        if (fpi<0f) fpi+=360f;
        fps[ip] = fpi;
      }
    }
  }


  public void convert(float[] fps) {
    int np = fps.length;
    for (int ip=0; ip<np; ++ip) {
      float fpi = fps[ip];
      if(fpi>180.0f) {
        fpi -= 180f;
        fps[ip] = fpi;
      }
    }
  }


  public PlotPanel applyForRosePlotsN(
    float alpha, float b1, float e1, int c2, int c3, int n2, int n3, 
    int nbin, float[][] fp, float[][] ob) 
  {
    int np = fp[0].length;
    int d2 = round(n2/c2);
    int d3 = round(n3/c3);
    PlotPanel.AxesPlacement axes = PlotPanel.AxesPlacement.NONE;
    PlotPanel.Orientation orient = PlotPanel.Orientation.X1RIGHT_X2UP;
    PlotPanel pp = new PlotPanel(c3,c2,orient,axes);
    for (int i3=0; i3<c3; ++i3) {
    //for (int i3=c3-1; i3>=0; --i3) {
      int b3 = i3*d3;
      int e3 = (i3+1)*d3;
      if(i3==(c3-1)){e3 = max(e3,n3);}
      System.out.println("b3="+b3);
      System.out.println("e3="+e3);
    for (int i2=0; i2<c2; ++i2) {
      int b2 = i2*d2;
      int e2 = (i2+1)*d2;
      if(i2==(c2-1)){e2 = max(e2,n2);}
      System.out.println("b2="+b2);
      System.out.println("e2="+e2);
      ArrayList<Float> fpa = new ArrayList<Float>();
      for (int ip=0; ip<np; ++ip) {
        int k2 = (int)fp[1][ip];
        int k3 = (int)fp[2][ip];
        if(k2<b2 || k2>e2) {continue;}
        if(k3<b3 || k3>e3) {continue;}
        float obi = ob[k3][k2];
        int k1 = (int)fp[0][ip];
        if(k1>=(b1+obi) && k1<=(e1+obi)) {
        //if(k1>=b1 && k1<=e1) {
          float fpi = fp[3][ip];
          fpa.add(fpi);
        }
      }
      int npk = fpa.size();
      System.out.println("npk="+npk);
      float[] fpk = new float[npk];
      for (int ik=0; ik<npk; ++ik) {
        fpk[ik] = fpa.get(ik);
      }
      fpa.clear();
      roseN(i3,i2,alpha,fpk,nbin,pp);
    }}
    return pp;
  }

  public PlotPanel applyForRosePlotsN(
    float alpha, float[][] b1, float[][] e1, int c2, int c3, int n2, int n3, 
    int nbin, float[][] fp) 
  {
    int np = fp[0].length;
    int d2 = round(n2/c2);
    int d3 = round(n3/c3);
    PlotPanel.AxesPlacement axes = PlotPanel.AxesPlacement.NONE;
    PlotPanel.Orientation orient = PlotPanel.Orientation.X1RIGHT_X2UP;
    PlotPanel pp = new PlotPanel(c3,c2,orient,axes);
    for (int i3=0; i3<c3; ++i3) {
    //for (int i3=c3-1; i3>=0; --i3) {
      int b3 = i3*d3;
      int e3 = (i3+1)*d3;
      if(i3==(c3-1)){e3 = max(e3,n3);}
      System.out.println("b3="+b3);
      System.out.println("e3="+e3);
    for (int i2=0; i2<c2; ++i2) {
      int b2 = i2*d2;
      int e2 = (i2+1)*d2;
      if(i2==(c2-1)){e2 = max(e2,n2);}
      System.out.println("b2="+b2);
      System.out.println("e2="+e2);
      ArrayList<Float> fpa = new ArrayList<Float>();
      for (int ip=0; ip<np; ++ip) {
        int k2 = (int)fp[1][ip];
        int k3 = (int)fp[2][ip];
        if(k2<b2 || k2>e2) {continue;}
        if(k3<b3 || k3>e3) {continue;}
        float b1i = b1[k3][k2];
        float e1i = e1[k3][k2];
        int k1 = (int)fp[0][ip];
        if(k1>b1i && k1<e1i) {
        //if(k1>=b1 && k1<=e1) {
          float fpi = fp[3][ip];
          fpa.add(fpi);
        }
      }
      int npk = fpa.size();
      System.out.println("npk="+npk);
      float[] fpk = new float[npk];
      for (int ik=0; ik<npk; ++ik) {
        fpk[ik] = fpa.get(ik);
      }
      fpa.clear();
      roseN(i3,i2,alpha,fpk,nbin,pp);
    }}
    return pp;
  }



  public PlotPanel applyForRosePlotsX(
    float alpha, float npm, float[][] st, float[][] sb, 
    int c2, int c3, int nbin, float[][] fp) 
  {
    int np = fp[0].length;
    int n3 = st.length;
    int n2 = st[0].length;
    int d2 = round(n2/c2);
    int d3 = round(n3/c3);
    System.out.println("npm="+npm);
    PlotPanel.AxesPlacement axes = PlotPanel.AxesPlacement.NONE;
    PlotPanel.Orientation orient = PlotPanel.Orientation.X1RIGHT_X2UP;
    PlotPanel pp = new PlotPanel(c3,c2,orient,axes);
    for (int i3=0; i3<c3; ++i3) {
    //for (int i3=c3-1; i3>=0; --i3) {
      int b3 = i3*d3;
      int e3 = (i3+1)*d3;
      if(i3==(c3-1)){e3 = max(e3,n3);}
      System.out.println("b3="+b3);
      System.out.println("e3="+e3);
    for (int i2=0; i2<c2; ++i2) {
      int b2 = i2*d2;
      int e2 = (i2+1)*d2;
      if(i2==(c2-1)){e2 = max(e2,n2);}
      System.out.println("b2="+b2);
      System.out.println("e2="+e2);
      ArrayList<Float> fpa = new ArrayList<Float>();
      for (int ip=0; ip<np; ++ip) {
        int k2 = (int)fp[1][ip];
        int k3 = (int)fp[2][ip];
        if(k2<b2 || k2>e2) {continue;}
        if(k3<b3 || k3>e3) {continue;}
        float sti = st[k3][k2];
        float sbi = sb[k3][k2];
        int k1 = (int)fp[0][ip];
        if(k1>=sti && k1<=sbi) {
          float fpi = fp[3][ip];
          fpa.add(fpi);
        }
      }
      int npk = fpa.size();
      float scale = npk/npm;
      System.out.println("npk="+npk);
      float[] fpk = new float[npk];
      for (int ik=0; ik<npk; ++ik) {
        fpk[ik] = fpa.get(ik);
      }
      fpa.clear();
      rose(i3,i2,alpha,scale,fpk,nbin,pp);
    }}
    return pp;
  }
  /*
  public int findMaxSamples(
    float b1, float e1, float d1, int c2, int c3, int n2, int n3, 
    float[][] fp, float[][] ob) {
    int npm = 0;
    for (float x1=b1; x1<e1; x1+=d1) {
      int npk = findMaxSamples(x1,x1+d1,c2,c3,n2,n3,fp,ob);
      if(npk>npm) npm=npk;
    }
    return npm;
  }

  private int findMaxSamples(
    float b1, float e1, int c2, int c3, int n2, int n3, 
    float[][] fp, float[][] ob) 
  {
    int npm = 0;
    int d2 = round(n2/c2);
    int d3 = round(n3/c3);
    int np = fp[0].length;
    for (int i3=0; i3<c3; ++i3) {
      int b3 = i3*d3;
      int e3 = (i3+1)*d3;
      if(i3==(c3-1)){e3 = max(e3,n3);}
    for (int i2=0; i2<c2; ++i2) {
      int npk = 0;
      int b2 = i2*d2;
      int e2 = (i2+1)*d2;
      if(i2==(c2-1)){e2 = max(e2,n2);}
      for (int ip=0; ip<np; ++ip) {
        int k2 = (int)fp[1][ip];
        int k3 = (int)fp[2][ip];
        if(k2<b2 || k2>e2) {continue;}
        if(k3<b3 || k3>e3) {continue;}
        float obi = ob[k3][k2];
        int k1 = (int)fp[0][ip];
        if(k1>=(b1+obi) && k1<=(e1+obi)) {
          npk++;
        }
      }
      if(npk>npm) npm=npk;
      System.out.println("npk="+npk);
    }}
    return npm;
  }
  */

  public int findMaxSamples(
    int c2, int c3, float[][] hu, float[][] hm, float[][] hl, float[][] fp) 
  {
    int np1 = findMaxSamples(c2,c3,hu,hm,fp);
    int np2 = findMaxSamples(c2,c3,hm,hl,fp);
    return max(np1,np2);
  }

  public int findMaxSamples(
    int c2, int c3, float[][] st, float[][] sb, float[][] fp) 
  {
    int npm = 0;
    int n3 = st.length;
    int n2 = st[0].length;
    int d2 = round(n2/c2);
    int d3 = round(n3/c3);
    int np = fp[0].length;
    for (int i3=0; i3<c3; ++i3) {
      int b3 = i3*d3;
      int e3 = (i3+1)*d3;
      if(i3==(c3-1)){e3 = max(e3,n3);}
    for (int i2=0; i2<c2; ++i2) {
      int npk = 0;
      int b2 = i2*d2;
      int e2 = (i2+1)*d2;
      if(i2==(c2-1)){e2 = max(e2,n2);}
      for (int ip=0; ip<np; ++ip) {
        int k2 = (int)fp[1][ip];
        int k3 = (int)fp[2][ip];
        if(k2<b2 || k2>e2) {continue;}
        if(k3<b3 || k3>e3) {continue;}
        float sti = st[k3][k2];
        float sbi = sb[k3][k2];
        int k1 = (int)fp[0][ip];
        if(k1>=sti && k1<=sbi) {
          npk++;
        }
      }
      if(npk>npm) npm=npk;
      System.out.println("npk="+npk);
    }}
    return npm;
  }

  public float[][] faultPoints(float[][] bt, FaultSkin[] sks) {
    ArrayList<Float> fpa = new ArrayList<Float>();
    ArrayList<Float> x1a = new ArrayList<Float>();
    ArrayList<Float> x2a = new ArrayList<Float>();
    ArrayList<Float> x3a = new ArrayList<Float>();
    ArrayList<Float> fva = new ArrayList<Float>();
    for (FaultSkin ski:sks) {
      int na = 0;
      int nb = 0;
      int nc = ski.size();  
      for (FaultCell fci:ski) {
        int i1 = fci.getI1();
        int i2 = fci.getI2();
        int i3 = fci.getI3();
        float bti = bt[i3][i2];
        if (i1<bti) {na++;}
        else        {nb++;}
      }
      float sa = (float)na/(float)nc;
      float sb = (float)nb/(float)nc;
      for (FaultCell fci:ski) {
        int i1 = fci.getI1();
        int i2 = fci.getI2();
        int i3 = fci.getI3();
        fpa.add(fci.getFp());
        x1a.add((float)i1);
        x2a.add((float)i2);
        x3a.add((float)i3);
        float bti = bt[i3][i2];
        if (i1<bti) {fva.add(sa);}
        else        {fva.add(sb);}
      }
    }
    int np = fpa.size();
    float[][] fps = new float[5][np];
    for (int ip=0; ip<np; ++ip) {
      fps[0][ip] = x1a.get(ip);
      fps[1][ip] = x2a.get(ip);
      fps[2][ip] = x3a.get(ip);
      fps[3][ip] = fpa.get(ip);
      fps[4][ip] = fva.get(ip);
    }
    return fps;
  }



  public float[][] faultPoints(float[][] ob, float[][][] fp) {
    int n3 = fp.length;
    int n2 = fp[0].length;
    int n1 = fp[0][0].length;
    ArrayList<Float> fpa = new ArrayList<Float>();
    ArrayList<Float> x1a = new ArrayList<Float>();
    ArrayList<Float> x2a = new ArrayList<Float>();
    ArrayList<Float> x3a = new ArrayList<Float>();
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=round(ob[i3][i2]); i1<n1; ++i1) {
      float fpi = fp[i3][i2][i1];
      if (fpi>=0.0f) {
        fpa.add(fpi);
        x1a.add((float)i1);
        x2a.add((float)i2);
        x3a.add((float)i3);
      }
    }}}
    int np = fpa.size();
    float[][] fps = new float[4][np];
    for (int ip=0; ip<np; ++ip) {
      fps[0][ip] = x1a.get(ip);
      fps[1][ip] = x2a.get(ip);
      fps[2][ip] = x3a.get(ip);
      fps[3][ip] = fpa.get(ip);
    }
    return fps;
  }


  public float[][] faultPoints(float[][][] fp) {
    int n3 = fp.length;
    int n2 = fp[0].length;
    int n1 = fp[0][0].length;
    ArrayList<Float> fpa = new ArrayList<Float>();
    ArrayList<Float> x1a = new ArrayList<Float>();
    ArrayList<Float> x2a = new ArrayList<Float>();
    ArrayList<Float> x3a = new ArrayList<Float>();
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float fpi = fp[i3][i2][i1];
      if (fpi>=0.0f) {
        fpa.add(fpi);
        x1a.add((float)i1);
        x2a.add((float)i2);
        x3a.add((float)i3);
      }
    }}}
    int np = fpa.size();
    float[][] fps = new float[4][np];
    for (int ip=0; ip<np; ++ip) {
      fps[0][ip] = x1a.get(ip);
      fps[1][ip] = x2a.get(ip);
      fps[2][ip] = x3a.get(ip);
      fps[3][ip] = fpa.get(ip);
    }
    return fps;
  }


  public void rose(
    int i3, int i2, float alpha, float scale, float[] phi, int nbin, PlotPanel pp) 
  {
    float[] bins = applyHistogram(nbin,phi);
    mul(bins,scale,bins);
    applyForRose(i3,i2,alpha,bins,pp);
  }

  public void roseN(
    int i3, int i2, float alpha, float[] phi, int nbin, PlotPanel pp) 
  {
    float[] bins = applyHistogram(nbin,phi);
    applyForRoseN(i3,i2,alpha,bins,pp);
  }





  private void applyForRoseN(
    int i3, int i2, float alpha, float[] bins,PlotPanel pp) 
  {
    for (float rdi=0.2f; rdi<=1.0f; rdi+=0.2f) {
      float[][] cps = makeCirclePoints(1000,rdi);
      PointsView pvc = pp.addPoints(i3,i2,cps[0],cps[1]);
      pvc.setLineStyle(PointsView.Line.DASH);
      if(rdi==1.0f) {
        pvc.setLineColor(Color.RED);
        pvc.setLineWidth(3.f);
      } else {
        pvc.setLineColor(Color.BLACK);
        pvc.setLineWidth(2.f);
      }
    }
    addRadials(i3,i2,alpha, pp,Color.BLACK);
    addBinsN(i3,i2,alpha, pp,bins);
  }


  private void applyForRose(
    int i3, int i2, float alpha, float[] bins,PlotPanel pp) 
  {
    for (float rdi=0.2f; rdi<=1.0f; rdi+=0.2f) {
      float[][] cps = makeCirclePoints(1000,rdi);
      PointsView pvc = pp.addPoints(i3,i2,cps[0],cps[1]);
      pvc.setLineStyle(PointsView.Line.DASH);
      if(rdi==1.0f) {
        pvc.setLineColor(Color.RED);
        pvc.setLineWidth(3.f);
      } else {
        pvc.setLineColor(Color.BLACK);
        pvc.setLineWidth(2.f);
      }
    }
    addRadials(i3,i2,alpha,pp,Color.BLACK);
    addBinsF(i3,i2,alpha,pp,bins);
  }



  private float[] applyHistogram(int nbins, float[] phi) {
    int np = phi.length;
    float ps = 1f/(float)np;
    float dp = 180f/(float)nbins;
    float[] bins = new float[nbins];
    for (int ip=0; ip<np; ++ip) {
      int ib = (int)(phi[ip]/dp);
      ib = min(ib,nbins-1);
      bins[ib] += ps;
    }
    return bins;
  }


  private void addRadials(PlotPanel pp, Color color) {
    int np = 8;
    float dp = (float)(2*DBL_PI/np);
    for (int ip=0; ip<np; ++ip) {
      float phi = dp*ip;
      float[][] rps = makeRadialPoints(1f,phi);
      PointsView pv = pp.addPoints(rps[0],rps[1]);
      pv.setLineStyle(PointsView.Line.DASH);
      pv.setLineColor(color);
      pv.setLineWidth(4.f);
    }
  }

  private void addRadials(
    int i3, int i2, float alpha, PlotPanel pp, Color color) 
  {
    int np = 4;
    float dp = 360f/np;
    for (int ip=0; ip<np; ++ip) {
      float phi = dp*ip+alpha;
      if(phi>360f) phi-=360f;
      phi = (float)toRadians(phi);
      float[][] rps = makeRadialPoints(1f,phi);
      PointsView pv = pp.addPoints(i3,i2,rps[0],rps[1]);
      pv.setLineStyle(PointsView.Line.DASH);
      pv.setLineColor(color);
      pv.setLineWidth(2.f);
    }
  }

    /**
   * Converts an angle measured in degrees to radians.
   * @param angdeg an angle, in degrees.
   * @return the angle in radians.
   */
  public static double toRadians(double angdeg) {
    return Math.toRadians(angdeg);
  }

  private void addBins(PlotPanel pp, float[] bins) {
    int nb = bins.length;
    float dp = (float)(2*DBL_PI/nb);
    float rmax = max(bins);
    mul(bins,0.8f/rmax,bins);
    for (int ib=0; ib<nb; ++ib) {
      float phi = ib*dp;
      float rdi = bins[ib];
      Color rgb = getNextColor();
      float[][] rp1 = makeRadialPoints(rdi,phi);
      float[][] rp2 = makeRadialPoints(rdi,phi+dp);
      /*
      float[][] cps = makeCirclePoints(1000,rdi);
      PointsView pvc = pp.addPoints(cps[0],cps[1]);
      pvc.setLineStyle(PointsView.Line.DASH);
      pvc.setLineColor(rgb);
      pvc.setLineWidth(3.f);
      */
      float[] rx3 = new float[]{rp1[0][1],rp2[0][1]};
      float[] ry3 = new float[]{rp1[1][1],rp2[1][1]};
      float[][] rp3 = new float[][]{rx3,ry3};
      PointsView pvr1 = pp.addPoints(rp1[0],rp1[1]);
      PointsView pvr2 = pp.addPoints(rp2[0],rp2[1]);
      PointsView pvr3 = pp.addPoints(rp3[0],rp3[1]);
      pvr1.setLineColor(rgb);
      pvr2.setLineColor(rgb);
      pvr3.setLineColor(rgb);
      pvr1.setLineWidth(4.f);
      pvr2.setLineWidth(4.f);
      pvr3.setLineWidth(4.f);
    }
  }

  private void addBinsB(PlotPanel pp, float[] bins) {
    int nb = bins.length;
    float dp = (float)(2*DBL_PI/nb);
    float rmax = max(bins);
    mul(bins,0.8f/rmax,bins);
    for (int ib=0; ib<nb; ++ib) {
      float phi = ib*dp;
      float rdi = bins[ib];
      Color rgb = getNextColor();
      float[][] rp1 = makeRadialPoints(rdi,phi);
      float[][] rp2 = makeRadialPoints(rdi,phi+dp);
      float dx = rp2[0][1] - rp1[0][1];
      float dy = rp2[1][1] - rp1[1][1];
      float ds = sqrt(dx*dx+dy*dy);
      for (float di=0f; di<=ds; di+=0.01f) {
        float xi = rp1[0][1]+dx*di/ds;
        float yi = rp1[1][1]+dy*di/ds;
        float[] xs = new float[]{0f,xi};
        float[] ys = new float[]{0f,yi};
        PointsView pvr = pp.addPoints(xs,ys);
        pvr.setLineColor(rgb);
        pvr.setLineWidth(4.f);
      }
      /*
      float[][] cps = makeCirclePoints(1000,rdi);
      PointsView pvc = pp.addPoints(cps[0],cps[1]);
      pvc.setLineStyle(PointsView.Line.DASH);
      pvc.setLineColor(rgb);
      pvc.setLineWidth(3.f);
      float[] rx3 = new float[]{rp1[0][1],rp2[0][1]};
      float[] ry3 = new float[]{rp1[1][1],rp2[1][1]};
      float[][] rp3 = new float[][]{rx3,ry3};
      PointsView pvr1 = pp.addPoints(rp1[0],rp1[1]);
      PointsView pvr2 = pp.addPoints(rp2[0],rp2[1]);
      PointsView pvr3 = pp.addPoints(rp3[0],rp3[1]);
      pvr1.setLineColor(rgb);
      pvr2.setLineColor(rgb);
      pvr3.setLineColor(rgb);
      pvr1.setLineWidth(4.f);
      pvr2.setLineWidth(4.f);
      pvr3.setLineWidth(4.f);
      */
    }
  }


  private float addBinsN(int i3, int i2, float alpha, 
    PlotPanel pp, float[] bins) {
    int nb = bins.length;
    float dp = 180f/nb;
    float rmax = max(bins);
    float rc = round(rmax*120f)/100f;
    mul(bins,1.0f/rc,bins);
    for (int ib=0; ib<nb; ++ib) {
      float phi = ib*dp;
      float phm = phi+alpha;
      float php = phi+alpha+dp;
      if (phm>180f) phm-=180f;
      if (php>180f) php-=180f;
      phi = (float)toRadians(phi);
      phm = (float)toRadians(phm);
      php = (float)toRadians(php);
      float rdi = bins[ib];
      float aci = (float)(phi/DBL_PI);
      Color rgb = Color.getHSBColor(aci,1f,1f);
      float[][] rp1 = makeRadialPoints(rdi,phm);
      float[][] rp2 = makeRadialPoints(rdi,php);
      float dx = rp2[0][1] - rp1[0][1];
      float dy = rp2[1][1] - rp1[1][1];
      float ds = sqrt(dx*dx+dy*dy);
      float dd = 0.001f;
      float dk = dd*5f;
      float de = ds-dk;
      dx /= ds;
      dy /= ds;
      float xi,yi;
      float[] xs, ys;
      PointsView pvr;
      float di = dk;
      for ( ; di<de; di+=dd) {
        xi = rp1[0][1]+dx*di;
        yi = rp1[1][1]+dy*di;
        xs = new float[]{-xi,xi};
        ys = new float[]{-yi,yi};
        pvr = pp.addPoints(i3,i2,xs,ys);
        pvr.setLineColor(rgb);
        pvr.setLineWidth(3f);
      }
      // add black bounds
      /*
      xi = rp2[0][1];
      yi = rp2[1][1];
      xs = new float[]{-xi,xi};
      ys = new float[]{-yi,yi};
      pvr = pp.addPoints(i3,i2,xs,ys);
      pvr.setLineColor(Color.BLACK);
      pvr.setLineWidth(2f);

      xi = rp1[0][1];
      yi = rp1[1][1];
      xs = new float[]{-xi,xi};
      ys = new float[]{-yi,yi};
      pvr = pp.addPoints(i3,i2,xs,ys);
      pvr.setLineColor(Color.BLACK);
      pvr.setLineWidth(2f);

      float[][] rt1 = makeRadialPoints(rdi*1.01f,phi);
      float[][] rt2 = makeRadialPoints(rdi*1.01f,phi+dp);
      xs = new float[]{rt1[0][1],rt2[0][1]};
      ys = new float[]{rt1[1][1],rt2[1][1]};
      float[][] rp3 = new float[][]{xs,ys};
      pvr = pp.addPoints(i3,i2,rp3[0],rp3[1]);
      pvr.setLineColor(Color.BLACK);
      pvr.setLineWidth(2.f);
      pvr = pp.addPoints(i3,i2,mul(rp3[0],-1),mul(rp3[1],-1));
      pvr.setLineColor(Color.BLACK);
      pvr.setLineWidth(2.f);
      */
    }
    return rc;
  }


  private float addBinsF(int i3, int i2, float alpha, PlotPanel pp, float[] bins) {
    int nb = bins.length;
    float dp = 180f/nb;
    float rth = 0.3f;
    float rmax = max(bins);
    if (rmax>rth) {
      mul(rth/rmax,bins,bins);
    }
    /*
    for (int ib=0; ib<nb; ++ib) {
      if(bins[ib]>rmax){
        bins[ib]=rmax;
      }
    }
    */
    float rc = rth+0.01f;//round(rmax*120f)/100f;
    mul(bins,1.0f/rc,bins);
    for (int ib=0; ib<nb; ++ib) {
      float phi = ib*dp;
      float phm = phi+alpha;
      float php = phi+alpha+dp;
      if (phm>180f) phm-=180f;
      if (php>180f) php-=180f;
      phi = (float)toRadians(phi);
      phm = (float)toRadians(phm);
      php = (float)toRadians(php);
      float rdi = bins[ib];
      float aci = (float)(phi/DBL_PI);
      Color rgb = Color.getHSBColor(aci,1f,1f);
      float[][] rp1 = makeRadialPoints(rdi,phm);
      float[][] rp2 = makeRadialPoints(rdi,php);
      float dx = rp2[0][1] - rp1[0][1];
      float dy = rp2[1][1] - rp1[1][1];
      float ds = sqrt(dx*dx+dy*dy);
      float dd = 0.001f;
      float dk = dd*5f;
      float de = ds-dk;
      dx /= ds;
      dy /= ds;
      float xi,yi;
      float[] xs, ys;
      PointsView pvr;
      float di = dk;
      for ( ; di<de; di+=dd) {
        xi = rp1[0][1]+dx*di;
        yi = rp1[1][1]+dy*di;
        xs = new float[]{-xi,xi};
        ys = new float[]{-yi,yi};
        pvr = pp.addPoints(i3,i2,xs,ys);
        pvr.setLineColor(rgb);
        pvr.setLineWidth(3f);
      }
    }
    return rc;
  }



  // Makes array of points on the unit circle in a plane.
  private float[][] makeCirclePoints(int nt, float r) {
    double dt = 2.0*DBL_PI/(nt-1);
    float[] x = new float[nt];
    float[] y = new float[nt];
    for (int it=0; it<nt; ++it) {
      float t = (float)(it*dt);
      x[it] = r*sin(t);
      y[it] = r*cos(t);
    }
    return new float[][]{x,y};
  }

  private float[][] makeRadialPoints(float r, float phi) {
    float[] x = new float[2];
    float[] y = new float[2];
    x[0] = 0.0f; x[1] = r*sin(phi); 
    y[0] = 0.0f; y[1] = r*cos(phi);
    return new float[][]{x,y};
  }

  private Color getNextColor() {
    float hue = (float)Math.random();
    float sat = (float)Math.random()*0.6f+0.4f;
    float bri = (float)Math.random()*0.3f+0.7f;
    int c = Color.HSBtoRGB(hue,sat,bri);
    int r = c&0x00ff0000;
    int g = c&0x0000ff00;
    int b = c&0x000000ff;
    return new Color(r|g|b);
  }
  
}
