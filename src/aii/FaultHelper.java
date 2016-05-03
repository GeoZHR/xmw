package aii;

import java.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import ipfx.*;

public class FaultHelper {

  public void mask(float[][] sf, float[][][] gx) {
    int n3 = gx.length;
    int n2 = gx[0].length;
    float ga = 0.5f*(max(gx)-min(gx));
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      int k1 = round(sf[i3][i2]);
      for (int i1=0; i1<=k1; ++i1)
        gx[i3][i2][i1] = ga;
    }}
  }

  public void mask(float[][] sf1, float[][] sf2, float[][][] gx) {
    int n3 = gx.length;
    int n2 = gx[0].length;
    float ga = 0.5f*(max(gx)-min(gx));
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      int k1 = round(sf1[i3][i2]);
      int k2 = round(sf2[i3][i2]);
      for (int i1=k1; i1<=k2; ++i1)
        gx[i3][i2][i1] = ga;
    }}
  }

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

  public void setValuesOnFaults(float v, FaultSkin[] skins, float[][][] x) {
    for (FaultSkin skin:skins) {
      for (FaultCell cell:skin) {
        int i1 = cell.getI1();
        int m2 = cell.getM2();
        int m3 = cell.getM3();
        int p2 = cell.getP2();
        int p3 = cell.getP3();
        x[m3][m2][i1] = v;
        x[p3][p2][i1] = v;
      }
    }
  }


  public void setValueOnFaultsInt(float v, FaultSkin[] skins, float[][][] x) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    for (FaultSkin skin:skins) {
      for (FaultCell cell:skin) {
        int i1 = cell.getI1();
        int m2 = cell.getM2();
        int m3 = cell.getM3();
        int p2 = cell.getP2();
        int p3 = cell.getP3();
        int i2 = cell.getI2();
        int i3 = cell.getI3();
        if(m2<0) {m2=0;}
        if(p2<0) {p2=0;}
        if(m2>=n2-1) {m2=n2-1;}
        if(p2>=p2-1) {p2=n2-1;}
        if(m3<0) {m3=0;}
        if(p3<0) {p3=0;}
        if(m3>=n3-1) {m3=n3-1;}
        if(p3>=p3-1) {p3=n3-1;}
        int k1 = i1*8;
        for (int k=0; k<8; k++) {
          int p1 = k1+k;
          if (p1<n1) {
            x[m3][m2][k1+k] = v;
            x[p3][p2][k1+k] = v;
            x[i3][i2][k1+k] = v;
          }
        }
      }
    }
  }

  public void getFlOnFaults(FaultSkin[] skins, float[][][] x) {
    int n3 = x.length;
    int n2 = x[0].length;
    for (FaultSkin skin:skins) {
      for (FaultCell cell:skin) {
        int i1 = cell.getI1();
        int m2 = cell.getM2();
        int m3 = cell.getM3();
        int p2 = cell.getP2();
        int p3 = cell.getP3();
        int i2 = cell.getI2();
        int i3 = cell.getI3();
        float fl = cell.getFl();
        if(m2<0) {m2=0;}
        if(p2<0) {p2=0;}
        if(m2>=n2-1) {m2=n2-1;}
        if(p2>=p2-1) {p2=n2-1;}
        if(m3<0) {m3=0;}
        if(p3<0) {p3=0;}
        if(m3>=n3-1) {m3=n3-1;}
        if(p3>=p3-1) {p3=n3-1;}
        x[m3][m2][i1] = fl;
        x[p3][p2][i1] = fl;
        x[i3][i2][i1] = fl;
      }
    }
  }


  public void getFlOnFaultsInt(FaultSkin[] skins, float[][][] x) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    for (FaultSkin skin:skins) {
      for (FaultCell cell:skin) {
        int i1 = cell.getI1();
        int m2 = cell.getM2();
        int m3 = cell.getM3();
        int p2 = cell.getP2();
        int p3 = cell.getP3();
        int i2 = cell.getI2();
        int i3 = cell.getI3();
        float fl = cell.getFl();
        if(m2<0) {m2=0;}
        if(p2<0) {p2=0;}
        if(m2>=n2-1) {m2=n2-1;}
        if(p2>=p2-1) {p2=n2-1;}
        if(m3<0) {m3=0;}
        if(p3<0) {p3=0;}
        if(m3>=n3-1) {m3=n3-1;}
        if(p3>=p3-1) {p3=n3-1;}
        int k1 = i1*8;
        for (int k=0; k<8; k++) {
          int p1 = k1+k;
          if (p1<n1) {
            x[m3][m2][k1+k] = fl;
            x[p3][p2][k1+k] = fl;
            x[i3][i2][k1+k] = fl;
          }
        }
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
