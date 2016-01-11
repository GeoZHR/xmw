package f3d;

import java.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import ipfx.*;

public class FaultSelect {

  public void setValueOnFaults(float v, FaultSkin[] skins, float[][][] x) {
    for (FaultSkin skin:skins) {
      for (FaultCell cell:skin) {
        int i1 = cell.getI1();
        int i2 = cell.getI2();
        int i3 = cell.getI3();
        x[i3][i2][i1] = v;
      }
    }
  }

  public void setValueOnFaults(float v, FaultCell[] cells, float[][][] x) {
    for (FaultCell cell:cells) {
      int i1 = cell.getI1();
      int i2 = cell.getI2();
      int i3 = cell.getI3();
      x[i3][i2][i1] = v;
    }
  }


  public FaultSkin[] getSubSkin(FaultSkin[] sks) {
    ArrayList<FaultSkin> sl = new ArrayList<FaultSkin>();
    for (FaultSkin sk:sks){
      int i1m = 0;
      for (FaultCell fc:sk) {
        int i1 = fc.getI1();
        if(i1>i1m) {i1m=i1;}
      }
      if(i1m>200) sl.add(sk);
    }
    return sl.toArray(new FaultSkin[0]);
  }

  public void setUncValues(int f1, int l1, float[][][] uc) {
    int n3 = uc.length;
    int n2 = uc[0].length;
    int n1 = uc[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      if (i1<f1||i1>l1) {
        uc[i3][i2][i1] = 0f;
      }
    }}}
  }

}
