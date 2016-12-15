package sse;

import java.util.List;

import vec.Point3d;
import vec.Vector3d;



public class Face {

  public int [] nodeIds=new int[3];
  /**
   * Creates a new polygon from the specified vertex indices. 
   * @param vertexIndices
   * the vertex indices, in counter-clockwise order
   */
  public Face(int[] ids) {
    nodeIds[0] = ids[0];
    nodeIds[1] = ids[1];
    nodeIds[2] = ids[2];
  }
    
  public Face clone(){
    return new Face(this.nodeIds);
  }
    
  public void setIndices(int[] ids){
    nodeIds[0] = ids[0];
    nodeIds[1] = ids[1];
    nodeIds[2] = ids[2];
  }

  /**
   * Calculate the area of this face, using the specified vertex list to fetch the vertex
   * positions
   * 
   * @param vertices
   *            the vertex list where to fetch the positions from
   * @return the surface of this polygon
   */
  public double getArea(List<Node> nodes) {
    Vector3d v12 = new Vector3d();
    Vector3d v13 = new Vector3d();
    Vector3d cross = new Vector3d();
    Vector3d v1 = new Vector3d(nodes.get(nodeIds[0]).position);
    Point3d v2 = nodes.get(nodeIds[1]).position;
    Point3d v3 = nodes.get(nodeIds[2]).position;
    v12.sub(v2, v1);
    v13.sub(v3, v1);
    cross.cross(v12, v13);
    return cross.length() * 0.5f;
  }

  public void replace(int oldIndex, int newIndex) {
    int position = indexOf(oldIndex);
    if (position >= 0) nodeIds[position] = newIndex;
  }


  public boolean isEdgeOrdered(int v1, int v2) throws IllegalArgumentException {
    int i1 = indexOf(v1);
    if (i1 == -1) throw new IllegalArgumentException("Vertex index " + i1 + " does not belong to this face");
    int i2 = indexOf(v2);
    if (i2 == -1) throw new IllegalArgumentException("Vertex index " + i2 + " does not belong to this face");
    return (i1 + 1)%3 == i2;
  }

  public boolean contains(int id) {
    for(int i : nodeIds)
      if (i==id) return true;
    return false;
  }


  public int indexOf(int nodeId){
    for (int i = 0; i <3; i++)
      if (nodeIds[i]==nodeId) return i;
    return -1;
  }

}
