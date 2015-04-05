package ipfx;

import java.util.*;

public class FaultSkinSub {

  public FaultSkinSub(int f1, int f2, int f3, int l1, int l2, int l3) {
    _f1 = f1;
    _f2 = f2;
    _f3 = f3;
    _l1 = l1;
    _l2 = l2;
    _l3 = l3;
  }

  public void setValueOnFaults(float v, FaultSkin[] skins, float[][][] x) {
    for (FaultSkin skin:skins) {
      for (FaultCell cell:skin) {
        int i1 = cell.i1;
        int i2 = cell.i2;
        int i3 = cell.i3;
        x[i3][i2][i1] = v;
      }
    }
  }

  public void setValueOnFaults(float v, FaultCell[] cells, float[][][] x) {
    for (FaultCell cell:cells) {
      int i1 = cell.i1;
      int i2 = cell.i2;
      int i3 = cell.i3;
      x[i3][i2][i1] = v;
    }
  }


  public FaultSkin[] getSubSkin(FaultSkin[] sks) {
    ArrayList<FaultSkin> sl = new ArrayList<FaultSkin>();
    ArrayList<FaultCell> cl = new ArrayList<FaultCell>();
    for (FaultSkin sk:sks){
      for (FaultCell fc:sk) {
        int i1 = fc.i1;
        int i2 = fc.i2;
        int i3 = fc.i3;
        if(i1<_f1||i1>_l1){continue;}
        if(i2<_f2||i2>_l2){continue;}
        if(i3<_f3||i3>_l3){continue;}
        float fl = fc.fl;
        float fp = fc.fp;
        float ft = fc.ft;
        float x1 = fc.x1-_f1;
        float x2 = fc.x2-_f2;
        float x3 = fc.x3-_f3;
        FaultCell cell = new FaultCell(x1,x2,x3,fl,fp,ft);
        cell.s1 = fc.s1;
        cell.s2 = fc.s2;
        cell.s3 = fc.s3;
        cl.add(cell);
      }
      FaultSkinner fs = new FaultSkinner();
      fs.setGrowLikelihoods(0.01,0.2);
      fs.setMinSkinSize(10);
      FaultSkin[] skins = fs.findSkins(cl.toArray(new FaultCell[0]));
      for (FaultSkin ski:skins){sl.add(ski);}
    }
    return sl.toArray(new FaultSkin[0]);
  }

  public FaultCell[] getSubCell(FaultCell[] cells) {
    ArrayList<FaultCell> cl = new ArrayList<FaultCell>();
    float x2m = 0f;
    float x3m = 0f;
    for (FaultCell fc:cells) {
      int i1 = fc.i1;
      int i2 = fc.i2;
      int i3 = fc.i3;
      if(i1<_f1||i1>_l1){continue;}
      if(i2<_f2||i2>_l2){continue;}
      if(i3<_f3||i3>_l3){continue;}
      float fl = fc.fl;
      float fp = fc.fp;
      float ft = fc.ft;
      float x1 = fc.x1-_f1;
      float x2 = fc.x2-_f2;
      float x3 = fc.x3-_f3;
      if(x2>x2m){x2m=x2;}
      if(x3>x3m){x3m=x3;}
      FaultCell cell = new FaultCell(x1,x2,x3,fl,fp,ft);
      cl.add(cell);
    }
    return cl.toArray(new FaultCell[0]);
  }

  public void setForNewCells(float fl, FaultSkin[] skold, FaultSkin sknew) {
    int n1 = _l1-_f1+1;
    int n2 = _l2-_f2+1;
    int n3 = _l3-_f3+1;
    FaultCell[][][] fcg = new FaultCell[n3][n2][n1];
    for (FaultSkin ski:skold) {
    for (FaultCell cell:ski) {
      int i1i = cell.i1;
      int i2i = cell.i2;
      int i3i = cell.i3;
      int i2m = cell.i2m;
      int i3m = cell.i3m;
      int i2p = cell.i2p;
      int i3p = cell.i3p;
      fcg[i3i][i2i][i1i] = cell;
      fcg[i3m][i2m][i1i] = cell;
      fcg[i3p][i2p][i1i] = cell;
    }}
    for (FaultCell cell:sknew) {
      int i1i = cell.i1;
      int i2i = cell.i2;
      int i3i = cell.i3;
      int i2m = cell.i2m;
      int i3m = cell.i3m;
      int i2p = cell.i2p;
      int i3p = cell.i3p;
      FaultCell fci = fcg[i3i][i2i][i1i];
      FaultCell fcm = fcg[i3m][i2m][i1i];
      FaultCell fcp = fcg[i3p][i2p][i1i];
      if(fci==null&&fcm==null&&fcp==null) {
        cell.setFl(fl);
      }
    }
  }

  private int _f1,_l1;
  private int _f2,_l2;
  private int _f3,_l3;
}
