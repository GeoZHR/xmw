package mef;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.Parallel.*;
import static edu.mines.jtk.util.ArrayMath.*;

import slt.*;

/**
 * Construct smooth and single-valued fault surface using 
 * the screen poisson method.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.02.29
 */

public class FaultReskin {
  
  public float[][][] faultIndicator(FaultSkin skin) {
    FaultCell[] fcs = skin.getCells();
    setCells(fcs);
    System.out.println("n1="+_n1);
    System.out.println("n2="+_n2);
    System.out.println("n3="+_n3);
    System.out.println("fault setting done...");
    float[][][] fl  = new float[_n3][_n2][_n1];
    float[][][] ws  = new float[_n3][_n2][_n1];
    float[][][] u1  = new float[_n3][_n2][_n1];
    float[][][] u2  = new float[_n3][_n2][_n1];
    float[][][] u3  = new float[_n3][_n2][_n1];
    float[][][] g11 = new float[_n3][_n2][_n1];
    float[][][] g12 = new float[_n3][_n2][_n1];
    float[][][] g13 = new float[_n3][_n2][_n1];
    float[][][] g22 = new float[_n3][_n2][_n1];
    float[][][] g23 = new float[_n3][_n2][_n1];
    float[][][] g33 = new float[_n3][_n2][_n1];
    initialTensors(fl,ws,g11,g12,g13,g22,g23,g33);
    System.out.println("tensors done...");
    solveEigenproblems(g11,g12,g13,g22,g23,g33,u1,u2,u3);
    System.out.println("normals done...");
    ScreenPoissonSurfer sps = new ScreenPoissonSurfer();
    sps.setSmoothings(20,20,20);
    mul(ws,u1,u1);
    mul(ws,u2,u2);
    mul(ws,u3,u3);
    float[][][] sf = sps.saltIndicator(fl,u1,u2,u3);
    System.out.println("fault indicator done...");
    return sf;
  }

  private void setCells(FaultCell[] cells) {
    int i1min = Integer.MAX_VALUE;
    int i2min = Integer.MAX_VALUE;
    int i3min = Integer.MAX_VALUE;
    int i1max = -i1min;
    int i2max = -i2min;
    int i3max = -i3min;
    for (FaultCell cell:cells) {
      if (cell.i1<i1min) i1min = cell.i1;
      if (cell.i2<i2min) i2min = cell.i2;
      if (cell.i3<i3min) i3min = cell.i3;
      if (cell.i1>i1max) i1max = cell.i1;
      if (cell.i2>i2max) i2max = cell.i2;
      if (cell.i3>i3max) i3max = cell.i3;
    }
    _j1 = i1min;
    _j2 = i2min;
    _j3 = i3min;
    _n1 = 1+i1max-i1min;
    _n2 = 1+i2max-i2min;
    _n3 = 1+i3max-i3min;
    _cells = new FaultCell[_n3][_n2][_n1];
    for (FaultCell cell:cells)
      set(cell);
  }

  private void set(FaultCell cell) {
    int i1 = cell.i1-_j1;
    int i2 = cell.i2-_j2;
    int i3 = cell.i3-_j3;
    i1 = max(i1,0); i1 = min(i1,_n1-1);
    i2 = max(i2,0); i2 = min(i2,_n2-1);
    i3 = max(i3,0); i3 = min(i3,_n3-1);
    if(_cells[i3][i2][i1]==null){
      _cells[i3][i2][i1] = cell;
    }

  }

  private void setNull(FaultCell cell) {
    int i1 = cell.i1-_j1;
    int i2 = cell.i2-_j2;
    int i3 = cell.i3-_j3;
    i1 = max(i1,0); i1 = min(i1,_n1-1);
    i2 = max(i2,0); i2 = min(i2,_n2-1);
    i3 = max(i3,0); i3 = min(i3,_n3-1);
    _cells[i3][i2][i1]=null;
  }

  private void initialTensors(
    float[][][] fls, float[][][] wss,
    float[][][] g11, float[][][] g12, float[][][] g13,
    float[][][] g22, float[][][] g23, float[][][] g33) 
  {
    int[][][] mk = new int[_n3][_n2][_n1];
    for (int i3=0; i3<_n3; ++i3) {
    for (int i2=0; i2<_n2; ++i2) {
    for (int i1=0; i1<_n1; ++i1) {
      FaultCell cell = _cells[i3][i2][i1];
      if (cell!=null&&mk[i3][i2][i1]!=1) {
        FaultCell cm = cell;
        FaultCell[] cells = findOverlapCells(i1,i2,i3,cell);
        int nc = cells.length;
        if (nc>1) {
        int nbm = nabors(cell);
        for (int ic=0; ic<nc; ++ic) {
          FaultCell fci = cells[ic];
          if(notNearbyCells(cell,fci)) {
          int nbi = nabors(fci);
          if(nbi>nbm) {cm = fci;nbm = nbi;} 
          else {setNull(fci);}
          }
        }}
        float w1 = cm.w1;
        float w2 = cm.w2;
        float w3 = cm.w3;
        float fl = cm.fl;
        int k1 = cm.i1-_j1;
        int k2 = cm.i2-_j2;
        int k3 = cm.i3-_j3;
        mk[k3][k2][k1] = 1;
        fls[k3][k2][k1] = fl; 
        wss[k3][k2][k1] = fl; 
        g11[k3][k2][k1] = w1*w1*fl; 
        g12[k3][k2][k1] = w1*w2*fl; 
        g13[k3][k2][k1] = w1*w3*fl; 
        g22[k3][k2][k1] = w2*w2*fl; 
        g23[k3][k2][k1] = w2*w3*fl; 
        g33[k3][k2][k1] = w3*w3*fl; 
      }
    }}}
    RecursiveGaussianFilterP rgf1 = new RecursiveGaussianFilterP(8.0);
    RecursiveGaussianFilterP rgf2 = new RecursiveGaussianFilterP(16.0);
    float[][][] h = new float[_n3][_n2][_n1];
    float[][][][] gs = {wss,g11,g22,g33,g12,g13,g23};
    for (float[][][] g:gs) {
      rgf1.apply0XX(g,h); copy(g,h);
      rgf2.applyX0X(h,g); copy(h,g);
      rgf2.applyXX0(g,h); copy(h,g);
    }
  }

  private void nearestInterp(final float[][][] f) {
    final int n3 = f.length;
    final int n2 = f[0].length;
    final int n1 = f[0][0].length;
    final Sampling s2 = new Sampling(n2);
    final Sampling s3 = new Sampling(n3);
    loop(n1,new LoopInt() {
    public void compute(int i1) {
      System.out.println("i1="+i1);
      ArrayList<Float> fxl = new ArrayList<Float>();
      ArrayList<Float> x2l = new ArrayList<Float>();
      ArrayList<Float> x3l = new ArrayList<Float>();
      for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        float fi = f[i3][i2][i1];
        if(fi!=0.0f) {
          fxl.add(fi);
          x2l.add((float)i2);
          x3l.add((float)i3);
        }
      }}
      int np = fxl.size();
      float[] fx = new float[np];
      float[] x2 = new float[np];
      float[] x3 = new float[np];
      for (int ip=0; ip<np; ++ip) {
        fx[ip] = fxl.get(ip);
        x2[ip] = x2l.get(ip);
        x3[ip] = x3l.get(ip);
      }
      NearestGridder2 ng = new NearestGridder2(fx,x2,x3);
      float[][] f1i = ng.grid(s2,s3);
      for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        f[i3][i2][i1]=f1i[i3][i2];
      }}
    }});
  }

  private void nearestInterp(int nc, final float[][][] g11) {
    /*
    ArrayList<Float> w12 = new ArrayList<Float>();
    ArrayList<Float> w13 = new ArrayList<Float>();
    ArrayList<Float> w22 = new ArrayList<Float>();
    ArrayList<Float> w23 = new ArrayList<Float>();
    ArrayList<Float> w33 = new ArrayList<Float>();
    */
    int ic = 0;
    float[][] xc = new float[3][nc];
    final float[] w11 = new float[nc];
    for (int i3=0; i3<_n3; ++i3) {
    for (int i2=0; i2<_n2; ++i2) {
    for (int i1=0; i1<_n1; ++i1) {
      if(g11[i3][i2][i1]!=0.0f) {
        xc[0][ic] = i1;
        xc[1][ic] = i2;
        xc[2][ic] = i3;
        w11[ic]=g11[i3][i2][i1];
        ic ++;
      }
    }}}
    final KdTree kt = new KdTree(xc);
    loop(_n3,new LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<_n2; ++i2) {
      for (int i1=0; i1<_n1; ++i1) {
      float[] x = new float[]{i1,i2,i3};
      int kc = kt.findNearest(x);
      g11[i3][i2][i1] = w11[kc];
    }}
    }});
  }

  private boolean notNearbyCells(FaultCell c1, FaultCell c2) {
    int d2 = c1.i2-c2.i2;
    int d3 = c1.i3-c2.i3;
    float ds = d2*d2+d3*d3;
    if(ds<=8f) {return false;}
    else {return true;}
  }

  public FaultCell[] findOverlapCells(int c1, int c2, int c3, FaultCell cell) {
    ArrayList<FaultCell> fcl = new ArrayList<FaultCell>();
    if (abs(cell.w2)>abs(cell.w3)) {
      for (int i2=0; i2<_n2;  ++i2) {
        FaultCell fc = _cells[c3][i2][c1];
        if(fc!=null) fcl.add(fc);
      }
    } else {
      for (int i3=0; i3<_n3;  ++i3) {
        FaultCell fc = _cells[i3][c2][c1];
        if(fc!=null) fcl.add(fc);
      }
    }
    return fcl.toArray(new FaultCell[0]);
  }

  private int nabors(FaultCell cell) {
    final int c1 = cell.i1-_j1;
    final int c2 = cell.i2-_j2;
    final int c3 = cell.i3-_j3;
    final int[] nc = new int[1];
    final int b2 = max(c2-100,0);
    final int b3 = max(c3-100,0);
    final int e2 = min(c2+100,_n2-1);
    final int e3 = min(c3+100,_n3-1);
    final float fp = cell.fp;
    Parallel.loop(b3,e3+1,1,new Parallel.LoopInt() {
    public void compute(int k3) {
      for (int k2=b2; k2<=e2; ++k2) {
      FaultCell fc = _cells[k3][k2][c1];
      if (fc!=null) {
        float del = fc.fp-fp;
        float fp = min(abs(del),abs(del+360.0f),abs(del-360.0f));
        if(fp<=15f) {nc[0] +=1;}
      }
    }
    }});
    return nc[0];
  }


  private int cellsLR(FaultCell cell) {
    int nc = 0;
    FaultCell c = cell;
    for (FaultCell cl=c.cl; cl!=null && cl!=cell; cl=c.cl)
      c = cl;
    FaultCell cLeft = c;
    for (c=c.cr; c!=null && c!=cLeft; c=c.cr) 
      nc ++;
    return nc;
  }


  private void solveEigenproblems(
    final float[][][] g11, final float[][][] g12, final float[][][] g13,
    final float[][][] g22, final float[][][] g23, final float[][][] g33,
    final float[][][] u1, final float[][][] u2, final float[][][] u3) 
  {
    final int n3 = g11.length;
    final int n2 = g11[0].length;
    final int n1 = g11[0][0].length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        double[][] a = new double[3][3];
        double[][] z = new double[3][3];
        double[] e = new double[3];
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            a[0][0] = g11[i3][i2][i1];
            a[0][1] = g12[i3][i2][i1];
            a[0][2] = g13[i3][i2][i1];
            a[1][0] = g12[i3][i2][i1];
            a[1][1] = g22[i3][i2][i1];
            a[1][2] = g23[i3][i2][i1];
            a[2][0] = g13[i3][i2][i1];
            a[2][1] = g23[i3][i2][i1];
            a[2][2] = g33[i3][i2][i1];
            Eigen.solveSymmetric33(a,z,e);
            float u1i = (float)z[0][0];
            float u2i = (float)z[0][1];
            float u3i = (float)z[0][2];
            if (u1i<0.0f) {
              u1i = -u1i;
              u2i = -u2i;
              u3i = -u3i;
            }
            if (u1!=null) u1[i3][i2][i1] = u1i;
            if (u2!=null) u2[i3][i2][i1] = u2i;
            if (u3!=null) u3[i3][i2][i1] = u3i;
          }
        }
      }
    });
  }


  private int _j1,_j2,_j3; // min cell indices
  private int _n1,_n2,_n3; // numbers of cells
  private FaultCell[][][] _cells; // array of cells

}


