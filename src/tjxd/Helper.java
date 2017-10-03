package tjxd;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class Helper {

  public float[][] fault(int n1, int n2, float[][] p1s, float[][] p2s) {
    float[][] fs = new float[n2][n1];
    int np = p1s.length;
    for (int ip=0; ip<np; ip++) {
      int nk = p1s[ip].length;
      for (int ik=0; ik<nk-1; ik++) {
        int p1i = round(p1s[ip][ik]);
        int p2i = round(p2s[ip][ik]);
        int p1p = round(p1s[ip][ik+1]);
        int p2p = round(p2s[ip][ik+1]);
        float dk = (float)(p1p-p1i)/(float)(p2p-p2i);
        int d = 1; if(p2p<p2i) d=-1;
        int i2 = p2i;
        while(i2!=p2p) {
          int i1 = round(p1i+dk*(i2-p2i));
          fs[i2][i1] = 1;
          i2 +=d;
        }
      }
    }
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(4);
    rgf.apply00(fs,fs);
    fs = sub(fs,min(fs));
    fs = div(fs,max(fs));
    for (int i2=0; i2<n2; i2++) {
    for (int i1=0; i1<n1; i1++) {
      if(fs[i2][i1]>0.85f) fs[i2][i1] = 1.0f;
    }}
    return fs;
  }

}
