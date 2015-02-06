package ifs;

import java.util.*;

public class FaultSlipTest {

  public FaultSlipTest(FaultSkin sk) {
    _sk = sk;
  }

  public FaultSkin[] getSubSkin() {
    ArrayList<FaultCell> cl = new ArrayList<FaultCell>();
    for (FaultCell fc:_sk) {
      int i1 = fc.i1;
      int i3 = fc.i3;
      if(i1>=68&&i1<88&&i3>=79&&i3<99) {
        fc.ca=null;
        fc.cb=null;
        fc.cl=null;
        fc.cr=null;
        fc.skin=null;
        if(i1==68&&i3==79){fc.fl=0.0f;}
        //if(i1==78&&i3==96){fc.fl=0.0f;}
        cl.add(fc);
      }
    }
    FaultSkinner fs = new FaultSkinner();
    fs.setGrowLikelihoods(0.1,0.2);
    fs.setMinSkinSize(10);
    return fs.findSkins(cl.toArray(new FaultCell[0]));
  }

  public FaultSkin[] applySlipper(FaultSkin[] skins, 
    float[][][] gsx, float[][][] p2, float[][][] p3) {
    FaultSlipper fsl = new FaultSlipper(gsx,p2,p3);
    fsl.setOffset(2.0);
    fsl.setZeroSlope(false);
    fsl.computeDipSlips(skins,0.01,15.0);
    FaultSkinner fsk = new FaultSkinner();
    fsk.setGrowLikelihoods(0.1,0.2);
    fsk.setMinSkinSize(10);
    fsk.setMinMaxThrow(0.01,15.0);
    return fsk.reskin(skins);
  }

  public void checkCellArrays(FaultSkin sk) {
    FaultCell[][] cab = sk.getCellsAB();
    FaultCell[][] clr = sk.getCellsLR();
    int nab = cab.length;
    int nlr = clr.length;
    System.out.println("nab="+nab);
    for (int iab=0; iab<nab; ++iab) {
      System.out.println("iab="+cab[iab].length);
    }
    System.out.println("nlr="+nlr);
    for (int ilr=0; ilr<nlr; ++ilr) {
      System.out.println("ilr="+clr[ilr].length);
    }

  }

  private FaultSkin _sk;
}
