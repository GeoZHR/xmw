package sse;

import java.util.HashSet;
import java.util.Set;

import vec.Point3d;
import vec.Vector3d;

/**
 * Structural element of a {@link ROI3DMesh 3D mesh}, representing a position in 3D space. This
 * class can be overridden to provide additional functionalities, under the condition that
 * subclasses also override the {@link #clone()} method accordingly.
 * 
 * @author Alexandre Dufour
 */
public class Node {

    public final Point3d position = new Point3d();
    public final Vector3d normal = new Vector3d();
    
    public final Set<Integer> neighbors;
    
    public Node(Node v) {
      this(v.position, v.neighbors);
    }
    
    public Node(Point3d position) {
        this(position, 0);
    }
    
    public Node(Point3d position, int nbNeighbors) {
      this.position.set(position);
      this.neighbors = new HashSet<Integer>(nbNeighbors);
    }
    
    public Node(Point3d position, Set<Integer> neighbors){
      this(position, neighbors.size());
      for (Integer i:neighbors)
        this.neighbors.add(i.intValue());
    }
    
    public Node clone(){
      return new Node(this);
    }
    
    public double distanceTo(Node v){
      return position.distance(v.position);
    }
    
    public String toString(){
    return "Vertex {" + position.x + "," + position.y + "," 
      + position.z + "}, " + neighbors.size() + " neighbors";
    }
}
