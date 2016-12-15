package sse;

import vec.Point3d;
import vec.Tuple3d;
import vec.Vector3d;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;


public class TriMesh3D {
    
  public List<Face> faces = new ArrayList<Face>();
  public ArrayList<Node> nodes = new ArrayList<Node>();
  public Point3d center = new Point3d();
  public boolean maskUpToDate = false;
  public boolean maskUpdating = false;
  public Tuple3d pixelSize = new Point3d(1.0, 1.0, 1.0);
    
  public void setNodes(Collection<Node> newNodes) {
    synchronized (nodes){
      nodes.clear();
      nodes.ensureCapacity(newNodes.size());
      nodes.addAll(newNodes);
    }
  }
    
  public void setFaces(Collection<Face> newFaces) {
    synchronized (this) {
      faces.clear();
      if (faces instanceof ArrayList) {
        ((ArrayList<Face>) faces).ensureCapacity(newFaces.size());
      }
      faces.addAll(newFaces);
    }
  }
    
  public void addFace(Face face) {
    synchronized (this) {
      faces.add(face);
    }
  }
    
  public void replaceCell(Face oldFace, Face newFace) {
    synchronized (this) {
      faces.set(faces.indexOf(oldFace), newFace);
    }
  }
    
  /**
   * @param index
   *            the index of the cell to retrieve
   * @return the cell stored at the specified index
   */
  public Face getFace(int index) {
    return faces.get(index);
  }
    
  public int getNumberOfCells() {
    return faces.size();
  }
    
  public List<Face> getFaces() {
    List<Face> ft = new ArrayList<Face>(faces.size());
    ft.addAll(faces);
    return ft;
  }

  public float[] getTriGroup() {
    int nf = faces.size();
    float[] xyz = new float[nf*3*3];
    int k = 0;
    for (Face face:faces) {
    for (int id:face.nodeIds) {
      Node node = nodes.get(id);
      xyz[k++] = (float)node.position.x;
      xyz[k++] = (float)node.position.y;
      xyz[k++] = (float)node.position.z;
    }}
    return xyz;
  }
    
  /**
   * Adds the specified vertex to this mesh (if no vertex already exists at this vertex location)
   * 
   * @param v
   *            the vertex to add
   * @return the index of the new vertex, or that of the old vertex at that location (if any)
   */
  public int addNode(Node node) {
    return addNode(node, true, true);
  }
    
  public int addNode(Node node, boolean tidy, boolean mergeDuplicates) {
    if (!tidy && !mergeDuplicates) {
      nodes.add(node);
      return nodes.size() - 1;
    }
    Integer index, nullIndex = -1;
    for (index=0; index<nodes.size(); index++) {
      Node existingVertex = nodes.get(index);
      if (existingVertex == null) {
        nullIndex = index;
      } else if (mergeDuplicates && existingVertex.position.epsilonEquals(node.position, 0.00001)) {
          return index;
      }
    }
    //if there is a free spot in the list, use it
    if (tidy && nullIndex >= 0) {
      index = nullIndex;
      nodes.set(index, node);
    } else {
      nodes.add(node);
    }
    return index;
  }
    
  public Node createNode(Point3d position) {
    return new Node(position);
  }
    
  /**
   * @param index
   * @return the vertex at the specified index in the buffer
   */
  public Node getNode(int index) {
    return nodes.get(index);
  }
    
  /**
   * @return the size of the vertex buffer
   * @param excludeNullElements
   * set to <code>true</code> if <code>null</code> elements should not be counted
   * (yielding the actual number of vertices, which may be smaller than the buffersize)
   */
  public int getNumberOfNodes(boolean excludeNullElements) {
    if (!excludeNullElements) return nodes.size();
    int count = 0;
    for (Node node : nodes)
      if (node != null) count++;
    return count;
  }
    
  /**
   * @return a safe copy of the list of vertices (changes to this list will not affect the mesh)
   */
  public List<Node> getNodes() {
    List<Node> nt = new ArrayList<Node>(nodes.size());
    nt.addAll(nodes);
    return nt;
  }
    
  public void setVertex(int index, Node node) {
    nodes.set(index, node);
  }
    
  /**
   * Checks mesh integrity and will print out to the console any missing vertex
   * @return <code>true</code> if the check succeeded.
   */
  public boolean checkMeshIntegrity() {
    boolean success = true;
    for (Face face:faces) {
    for (int i:face.nodeIds) {
      if (nodes.get(i)==null){
        success = false;
        System.err.println("missing vertex : " + i);
      }
    }}
    return success;
  }
    
  /**
   * Optimizes the vertex buffer by shifting all elements such that there are no <code>null</code>
   * elements intermingled between vertices.
   * 
   * @param tighten
   * <ul>
   * <li>if <code>true</code>: the vertex buffer is resized after optimization to
   * ensure it does not contain empty values anymore</li>
   * <li>if <code>false</code>: the size of the vertex buffer remains unchanged,
   * therefore <code>null</code> elements may be found at the end of the buffer after
   * optimization.</li>
   * </ul>
   */
  public void optimizeNodeBuffer() {
    int nds = nodes.size();
    int nRealNds = getNumberOfNodes(true);
    if (nRealNds==nds) return;
    ArrayList<Node> newNodes = new ArrayList<Node>(nRealNds);
    int blankSpaces = 0;
    for (int i=0; i<nds; i++) {
      Node node = nodes.get(i);
      // count the number of consecutive blank spaces
      if (node==null) {
        blankSpaces++;
        continue;
      }
      newNodes.add(node);
      if (blankSpaces > 0) {
      // first available space = i - blankSpaces
      // redirect all cells to this space
        for (Face face : faces)
          face.replace(i,i-blankSpaces);
      }
    }
    setNodes(newNodes);
  }
  
  /**
   * @return The major axis of this contour, i.e. an unnormalized vector linking the two most
   *         distant contour points (the orientation of this vector is arbitrary)
   */
  public Vector3d getMajorAxis() {
    Vector3d axis = new Vector3d();
    int nbPoints = nodes.size();
    // TODO this is not optimal, geometric moments should be used
    double maxDistSq = 0;
    Vector3d vec = new Vector3d();
    for (int i = 0; i < nbPoints; i++) {
      Node v1 = nodes.get(i);
      if (v1 == null) continue;
      for (int j=i+1; j<nbPoints; j++) {
        Node v2 = nodes.get(j);
        if (v2 == null) continue;
        vec.sub(v1.position, v2.position);
        double dSq = vec.lengthSquared();
        if (dSq > maxDistSq) {
          maxDistSq = dSq;
          axis.set(vec);
        }
      }
    }
    return axis;
  }

  /**
   * @param pixelSize
   *  the size of a pixel in metric space. This information is necessary to perform
   *  conversions from metric to image space.
   */
  public void setPixelSize(Tuple3d pixelSize) {
    this.pixelSize.set(pixelSize);
  }
  /**
   * @return the pixel size associated to this mesh
   */
  public Tuple3d getPixelSize() {
    return new Point3d(pixelSize);
  }

  public Point3d getCenter(boolean convertToImageSpace) {
    double dx = pixelSize.x;
    double dy = pixelSize.y;
    double dz = pixelSize.z;
    return convertToImageSpace ? new Point3d(center.x/dx, center.y/dy, center.z/dz) : new Point3d(center);
  }

  /**
   * Calculates Returns a 3D cuboid representing the contour's bounding box. The bounding box is
   * defined as the smallest cuboid that entirely contains the contour.
   * 
   * @return a {@link icy.type.rectangle.Rectangle3D.Double} object containing the bounding box of
   *         this mesh, expressed in voxel (image) space
   */
  public Rectangle3D computeBounds3D() {
    if (pixelSize.x == 0 || pixelSize.y == 0 || pixelSize.z == 0) {
      throw new RuntimeException("Invalid pixel size: (" + pixelSize.x + "," + pixelSize.y + "," + pixelSize.z + ')');
    }
    double minX = Double.MAX_VALUE, minY = Double.MAX_VALUE, minZ = Double.MAX_VALUE;
    double maxX = 0, maxY = 0, maxZ = 0;
    for (Node v : getNodes()) {
      if (v == null) continue;
      if (v.position.x < minX) minX = v.position.x;
      if (v.position.x > maxX) maxX = v.position.x;
      if (v.position.y < minY) minY = v.position.y;
      if (v.position.y > maxY) maxY = v.position.y;
      if (v.position.z < minZ) minZ = v.position.z;
      if (v.position.z > maxZ) maxZ = v.position.z;
    }
    // switch from real to voxel space
    minX /= pixelSize.x;
    minY /= pixelSize.y;
    minZ /= pixelSize.z;
    maxX /= pixelSize.x;
    maxY /= pixelSize.y;
    maxZ /= pixelSize.z;
    return new Rectangle3D.Double(minX, minY, minZ, maxX-minX+1, maxY-minY+1, maxZ-minZ+1);
  }
    
}
