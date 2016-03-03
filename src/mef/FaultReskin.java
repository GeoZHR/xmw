package mef;

import util.*;
import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

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
    float[][][] g11 = new float[_n3][_n2][_n1];
    float[][][] g12 = new float[_n3][_n2][_n1];
    float[][][] g13 = new float[_n3][_n2][_n1];
    float[][][] g22 = new float[_n3][_n2][_n1];
    float[][][] g23 = new float[_n3][_n2][_n1];
    float[][][] g33 = new float[_n3][_n2][_n1];
    initialTensors(g11,g12,g13,g22,g23,g33);
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(10.0);
    rgf.apply000(g11,g11);
    rgf.apply000(g12,g12);
    rgf.apply000(g13,g13);
    rgf.apply000(g22,g22);
    rgf.apply000(g23,g23);
    rgf.apply000(g33,g33);
    return null;
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
    float[][][] g11, float[][][] g12, float[][][] g13,
    float[][][] g22, float[][][] g23, float[][][] g33) 
  {
    short[][][] mk = new short[_n3][_n2][_n1];
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
        g11[k3][k2][k1] = w1*w1*fl; 
        g12[k3][k2][k1] = w1*w2*fl; 
        g13[k3][k2][k1] = w1*w3*fl; 
        g22[k3][k2][k1] = w2*w2*fl; 
        g23[k3][k2][k1] = w2*w3*fl; 
        g33[k3][k2][k1] = w3*w3*fl; 
      }
    }}}
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

  private int _j1,_j2,_j3; // min cell indices
  private int _n1,_n2,_n3; // numbers of cells
  private FaultCell[][][] _cells; // array of cells

}


