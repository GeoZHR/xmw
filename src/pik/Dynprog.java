package pik;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;


/**
 * Local similarity attribute
 * @author Xinming Wu 
 * @version 2016.03.23
 */

public class Dynprog {
 
  public Dynprog(int gate, float an) {
    _gate = gate;
    _an = an;
    _prob = new float[gate*2-1];
  }

  public float[][] applyForWeight(float[][] scan) {
    int n2 = scan.length;
    int n1 = scan[0].length;
    float[][] weight = new float[n1][n2];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      weight[i1][i2]  = exp(-scan[i2][i1]);
    }}
    return weight;
  }

  public float[] applyForPick2(int i0, float[][] scan, float[][] ttime) {
    float[][] wt = applyForWeight(scan);
    int n2 = wt.length;
    int n1 = wt[0].length;
    float[][] wr = new float[n2][n1];
    float[][] tt = new float[n1][n2];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      wr[n2-i2-1][i1] = wt[i2][i1];
    }}
    float[] jt = applyForPick(i0,wr,tt);
    float[] jr = new float[n2];
    for (int i1=0; i1<n1; ++i1) {
    for (int i2=0; i2<n2; ++i2) {
      ttime[i1][n2-i2-1] = tt[i1][i2];
    }}
    for (int i2=0; i2<n2; ++i2) {
      jr[n2-i2-1] = jt[i2];
    }
    return jr;
  }

  public float[] applyForPick(int i0, float[][] weight, float[][] ttime) {
    int n2 = weight.length;
    int n1 = weight[0].length;
    float[] traj = new float[n2];
    float[][] time = new float[n2][n1];
    float[][] what = new float[n2][n1];
    float[] prev = new float[n1];
    float[] next = new float[n1];
    float[] dist = new float[n1];
    for (int i1=0; i1<n1; ++i1) {
      dist[i1] = sqrt(i1*i1+_an*_an); 
    }
	  for (int i1=0; i1<n1; i1++) {
	    float w = 0.5f*(weight[0][i1]+weight[0][i0]);
	    time[0][i1] = abs(i1-i0)*w;
	  }
    
    for (int i1=0; i1<n1; i1++) {
	    float w = 0.5f*(weight[1][i1]+weight[0][i0]);
	    prev[i1] = dist[abs(i1-i0)]*w;
	    what[1][i1] = i0;
	    time[1][i1]=prev[i1];
    }

    for (int i2=2; i2<n2; i2++) {
	    for (int i1=0; i1<n1; i1++) {
	      float w = weight[i2][i1];
	      int ib = max(i1-_gate,-1);
	      int ie = min(i1+_gate,n1);
	      float c = FLT_MAX;
	      int ic = -1;
	      for (int i=ib+1; i<ie; i++) {
		      float w2 = 0.5f*(w+weight[i2-1][i]);
		      float d = dist[abs(i1-i)]*w2+prev[i];
		      int it = i-ib-1;
		      if (d < c) {
		        c =	d;
		        ic = it;
		      }
		      _prob[it]=d;
	      }
        float[] vs = find_minimum(ic,ie-ib-1,ib+1,c,what[i2][i1]);
	      next[i1]= vs[0];
        what[i2][i1] = vs[1];
	    }
	    for (int i1=0; i1<n1; i1++) {
	      prev[i1]=next[i1];
	      time[i2][i1]=prev[i1];
	    }
    }
    pick(traj, next, what);
    for (int i2=0; i2<n2;++i2) {
    for (int i1=0; i1<n1;++i1) {
      ttime[i1][i2]  = time[i2][i1];
    }}
    return traj;
  }


  public float[] applyForPick(int i0, float[][] weight) {
    int n2 = weight.length;
    int n1 = weight[0].length;
    float[] traj = new float[n2];
    float[][] time = new float[n2][n1];
    float[][] what = new float[n2][n1];
    float[] prev = new float[n1];
    float[] next = new float[n1];
    float[] dist = new float[n1];
    for (int i1=0; i1<n1; ++i1) {
      dist[i1] = sqrt(i1*i1+_an*_an); 
    }
	  for (int i1=0; i1<n1; i1++) {
	    float w = 0.5f*(weight[0][i1]+weight[0][i0]);
	    time[0][i1] = abs(i1-i0)*w;
	  }
    
    for (int i1=0; i1<n1; i1++) {
	    float w = 0.5f*(weight[1][i1]+weight[0][i0]);
	    prev[i1] = dist[abs(i1-i0)]*w;
	    what[1][i1] = i0;
	    time[1][i1]=prev[i1];
    }

    for (int i2=2; i2<n2; i2++) {
	    for (int i1=0; i1<n1; i1++) {
	      float w = weight[i2][i1];
	      int ib = max(i1-_gate,-1);
	      int ie = min(i1+_gate,n1);
	      float c = FLT_MAX;
	      int ic = -1;
	      for (int i=ib+1; i<ie; i++) {
		      float w2 = 0.5f*(w+weight[i2-1][i]);
		      float d = dist[abs(i1-i)]*w2+prev[i];
		      int it = i-ib-1;
		      if (d < c) {
		        c =	d;
		        ic = it;
		      }
		      _prob[it]=d;
	      }
        float[] vs = find_minimum(ic,ie-ib-1,ib+1,c,what[i2][i1]);
	      next[i1]= vs[0];
        what[i2][i1] = vs[1];
	    }
	    for (int i1=0; i1<n1; i1++) {
	      prev[i1]=next[i1];
	      time[i2][i1]=prev[i1];
	    }
    }
    pick(traj, next, what);
    return traj;
  }


  public float[][] applyForTime(int i0, float[][] weight) {
    int n2 = weight.length;
    int n1 = weight[0].length;
    float[][] time = new float[n2][n1];
    float[][] what = new float[n2][n1];
    float[] prev = new float[n1];
    float[] next = new float[n1];
    float[] dist = new float[n1];
    for (int i1=0; i1<n1; ++i1) {
      dist[i1] = sqrt(i1*i1+_an*_an); 
    }
	  for (int i1=0; i1<n1; i1++) {
	    float w = 0.5f*(weight[0][i1]+weight[0][i0]);
	    time[0][i1] = abs(i1-i0)*w;
	  }
    
    for (int i1=0; i1<n1; i1++) {
	    float w = 0.5f*(weight[1][i1]+weight[0][i0]);
	    prev[i1] = dist[abs(i1-i0)]*w;
	    what[1][i1] = i0;
	    time[1][i1]=prev[i1];
    }

    for (int i2=2; i2<n2; i2++) {
	    for (int i1=0; i1<n1; i1++) {
	      float w = weight[i2][i1];
	      int ib = max(i1-_gate,-1);
	      int ie = min(i1+_gate,n1);
	      float c = FLT_MAX;
	      int ic = -1;
	      for (int i=ib+1; i<ie; i++) {
		      float w2 = 0.5f*(w+weight[i2-1][i]);
		      float d = dist[abs(i1-i)]*w2+prev[i];
		      int it = i-ib-1;
		      if (d < c) {
		        c =	d;
		        ic = it;
		      }
		      _prob[it]=d;
	      }
        float[] vs = find_minimum(ic,ie-ib-1,ib+1,c,what[i2][i1]);
	      next[i1]= vs[0];
        what[i2][i1] = vs[1];
	    }
	    for (int i1=0; i1<n1; i1++) {
	      prev[i1]=next[i1];
	      time[i2][i1]=prev[i1];
	    }
    }
    float[][] ttime = new float[n1][n2];
    for (int i2=0; i2<n2;++i2) {
    for (int i1=0; i1<n1;++i1) {
      ttime[i1][i2]  = time[i2][i1];
    }}
    return ttime;
  }

  private float[] find_minimum(
    int ic, int nc, int jc, float c, float pick)
  {
    float fm, f0, fp, a, b;
    if (0==ic) {
	    ic++;
	    fm=c;
	    f0=_prob[ic];
	    fp=_prob[ic+1];
    } else if (nc-1==ic) {
	    ic--;
	    fm=_prob[ic-1];
	    f0=_prob[ic];
	    fp=c;
    } else {
	    fm=_prob[ic-1];
	    f0=c;
	    fp=_prob[ic+1];
    }

    ic += jc;
    a = fm+fp-2f*f0;
    if (a <= 0.) { /* no minimum */
	    if (fm < f0 && fm < fp) {
	      pick = ic-1;
	      return new float[]{fm,pick};
	    } 
	    if (fp < f0 && fp < fm) {
	      pick = ic+1;
	      return new float[]{fp,pick};
	    } 
	    pick = ic;
	    return new float[]{f0,pick};
    }

    b = 0.5f*(fm-fp);
    a = b/a;
    if (a > 1.) {
	    pick = ic+1;
	    return new float[]{fp,pick};
    }

    if (a < -1.) {
	    pick = ic-1;
	    return new float[]{fm,pick};
    }

    if (f0 < 0.5*b*a) {
	    pick = ic;
	    return new float[]{f0,pick};
    }

    f0 -= 0.5*b*a;
    pick=ic+a;
    return new float[]{f0,pick};
  }

  private void pick(float[] traj, float[] next, float[][] what)
   {
    float c, d, fc;
    int n1 = next.length;
    int n2 = traj.length;
    c = FLT_MAX;
    fc = 0;
    /* minimum at the bottom */
    for (int i1=0; i1 < n1; i1++) {
	    d = next[i1];
	    if (d < c) {
	      c = d;
	      fc = i1;
	    }
    }
    /* coming up */
    for (int i2=n2-1; i2 >= 0; i2--) {
	    traj[i2]=fc;
	    fc = interpolate(fc,i2,what);
    }
  }

  private float interpolate(float fc, int i2, float[][] what) {
    int n1 = what[0].length;
    int ic = round(fc-0.5f);
    fc -= ic;
    if (n1-1 <= ic) return what[i2][n1-1];
    if (0 > ic) return what[i2][0];
    fc = what[i2][ic]*(1f-fc)+what[i2][ic+1]*fc;
    return fc;
  }


  ///////////////////////////////////////////////////////////////////////////
  // private
  private int _gate;
  private float _an;
  private float[] _prob;
}
