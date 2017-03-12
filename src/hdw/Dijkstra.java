package hdw;

import static edu.mines.jtk.util.ArrayMath.*;
import java.util.*;
public class Dijkstra {

  public Dijkstra(float[][] ws) {
    _gg = new GridGraph(ws);
  }
   
  // Dijkstra's algorithm to find shortest path from s to all other nodes
  public void apply(int s, float[] dist, int[] pred) {
    final boolean [] visited = new boolean [_gg.size()]; // all false initially
    for (int i=0; i<dist.length; i++)
      dist[i] = Integer.MAX_VALUE;
    dist[s] = 0;
    for (int i=0; i<dist.length; i++) {
      System.out.println("i="+i);
      final int next = minVertex(dist, visited);
      visited[next] = true;
      //The shortest path to next is dist[next] and via pred[next].
      final int [] n = _gg.neighbors (next);
      for (int j=0; j<n.length; j++) {
        final int v = n[j];
        final float d = dist[next]+_gg.getWeight(next,v);
        if (dist[v] > d) {
          dist[v] = d;
          pred[v] = next;
        }
      }
    }
  }
  public int coord2Index(int i1, int i2) {
    return _gg.coord2Index(i1,i2);
  }


  private GridGraph _gg = null;

  private static int minVertex (float [] dist, boolean [] v) {
    int y = -1;   // graph not connected, or no unvisited vertices
    float x = Float.MAX_VALUE;
    for (int i=0; i<dist.length; i++) {
      if (!v[i] && dist[i]<x) {y=i; x=dist[i];}
    }
    return y;
  }
  
  public void printPath (int n2, int [] pred, int s, int e) {
    int x = e;
    while (x!=s) {
      int c2 = (int)floor(x/n2);
      int c1 = x-c2*n2;
      System.out.println("c1="+c1+"; c2="+c2);
      x = pred[x];
    }
    int c2 = (int)floor(s/n2);
    int c1 = s-c2*n2;
    System.out.println("c1="+c1+"; c2="+c2);
  }

  public class GridGraph {
    int _n1,_n2;
    int _size;
    float[][] _w;

    // construct a regular grid graph with a weigthing map.
    public GridGraph(float[][] w) {
      _n2 = w.length;
      _n1 = w[0].length;
      _size = _n1*_n2;
      _w = w;
    }

    public int size() {
      return this._size;
    }

    public int[] neighbors(int id) {
      ArrayList<Integer> bs = new ArrayList<Integer>();
      int[] ps = index2Coord(id);
      int p1 = ps[0];
      int p2 = ps[1];
      int a1 = p1-1;
      int b1 = p1+1;
      int l2 = p2-1;
      int r2 = p2+1;
      if(a1>=0) 
        bs.add(coord2Index(a1,p2));
      if(b1<_n1) 
        bs.add(coord2Index(b1,p2));
      if(l2>=0) 
        bs.add(coord2Index(p1,l2));
      if(r2<_n2) 
        bs.add(coord2Index(p1,r2));
      if(a1>=0&&l2>=0)
        bs.add(coord2Index(a1,l2));
      if(a1>=0&&r2<_n2)
        bs.add(coord2Index(a1,r2));
      if(b1<_n1&&l2>=0)
        bs.add(coord2Index(b1,l2));
      if(b1<_n1&&r2<_n2)
        bs.add(coord2Index(b1,r2));
      int ns = bs.size();
      int[] ids = new int[ns];
      for (int is=0; is<ns; ++is)
        ids[is] = bs.get(is);
      return ids;
    }

    public float getWeight(int id1, int id2) {
      int[] a = index2Coord(id1);
      int[] b = index2Coord(id2);
      int a1 = a[0];
      int a2 = a[1];
      int b1 = b[0];
      int b2 = b[1];
      int d1 = b1-a1;
      int d2 = b2-a2;
      float ds = sqrt(d1*d1+d2*d2);
      float wa = _w[a2][a1];
      float wb = _w[b2][b1];
      return 0.5f*(wa+wb)*ds;
    }

    public int coord2Index(int i1, int i2) {
      int id = i2*_n2+i1;
      return id;
    }

    public int[] index2Coord(int id) {
      int i2 = (int)floor(id/_n2);
      int i1 = id-i2*_n2;
      return new int[]{i1,i2};
    }
  }
  
}
