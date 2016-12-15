/****************************************************************************
Copyright 2006, Colorado School of Mines and others.
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****************************************************************************/
package sse;

import vec.Point3d;
import vec.Vector3d;
import quickhull.QuickHull3D;

import java.util.List;
import java.util.Stack;
import java.util.Iterator;
import java.util.ArrayList;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;


/** 
 * A 3-D triangulated manifold oriented surface, possibly with boundary.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.12.11
 */
public class ActiveMesh extends TriMesh3D {

  public ActiveMesh(float sampling, float c1, float c2, float c3, float r) {
    float dt = 1f;
    float dp = 1f;
    float ft = 0.f;
    float fp = 0.f;
    int nt = round((180f-ft)/dt);
    int np = round((360f-fp)/dt);
    Sampling st = new Sampling(nt,dt,ft);
    Sampling sp = new Sampling(np,dp,fp);
    double[] xyz = new double[np*nt*3];
    int k = 0;
    for (int it=0; it<nt; ++it) {
      float ti = (float)toRadians(st.getValue(it));
    for (int ip=0; ip<np; ++ip) {
      float pi = (float)toRadians(sp.getValue(ip));
      xyz[k++] = c3+r*sin(ti)*cos(pi);
      xyz[k++] = c2+r*sin(ti)*sin(pi);
      xyz[k++] = c1+r*cos(ti);
    }}
    try{ 
      QuickHull3D qh3 = new QuickHull3D(xyz);
      qh3.triangulate();
      for (Point3d pt : qh3.getVertices()){
        pt.x *= pixelSize.x;
        pt.y *= pixelSize.y;
        pt.z *= pixelSize.z;
        addNode(createNode(pt), false, false);
      }
      for (int[] face : qh3.getFaces()) {
        addFace(face);
      }
    } catch (Exception e) {
      nodes.clear();
      faces.clear();
    }
    System.out.println("Qhull triangulation done!");
    /*
    try{
    // finish with a proper sampling
      reSampleToAverageDistance(sampling, 0.4);
    } catch (MeshTopologyException e){
      System.err.println("Warning: couldn't initialize contour. Here is the stack trace: ");
      e.printStackTrace();
    }
    */
  }

  public ActiveMesh(double sampling, double[] pts_dbl) {
    try{ 
      QuickHull3D qh3 = new QuickHull3D(pts_dbl);
      qh3.triangulate();
      for (Point3d pt : qh3.getVertices()){
        pt.x *= pixelSize.x;
        pt.y *= pixelSize.y;
        pt.z *= pixelSize.z;
        addNode(createNode(pt), false, false);
      }
      int k = 0;
      for (int[] face : qh3.getFaces()) {
        addFace(face);
        k++;
      }
      System.out.println("k="+k);
    } catch (Exception e) {
      nodes.clear();
      faces.clear();
    }
    try{
    // finish with a proper sampling
      reSampleToAverageDistance(sampling, 0.4);
    } catch (MeshTopologyException e){
      System.err.println("Warning: couldn't initialize contour. Here is the stack trace: ");
      e.printStackTrace();
    }
  }


  public void updateU(int niter) {
    for (int iter=0; iter<niter; ++iter) {
      updateNormals();
      for (Node node:nodes) {
        node.position.x += node.normal.x;
        node.position.y += node.normal.y;
        node.position.z += node.normal.z;
      }
    }
  }

  public void updateT(int niter, float[][][] v1, float[][][] v2, float[][][] v3) {
    int n3 = v1.length;
    int n2 = v1[0].length;
    int n1 = v1[0][0].length;
    for (int iter=0; iter<niter; ++iter) {
      updateNormals();
      for (Node node:nodes) {
        double xi = node.position.x;
        double yi = node.position.y;
        double zi = node.position.z;
        //double dc = abs(xi-296)+abs(yi-282)+abs(zi-182);
        //if(dc==0.0) {continue;}
        double xs = xi;
        double ys = yi;
        double zs = zi;
        double np = 1.;
        for (int id:node.neighbors) {
          xs += nodes.get(id).position.x;
          ys += nodes.get(id).position.y;
          zs += nodes.get(id).position.z;
          np += 1.0;
        }
        double xm = xs/np-xi; 
        double ym = ys/np-yi; 
        double zm = zs/np-zi; 
        double ux = node.normal.x;
        double uy = node.normal.y;
        double uz = node.normal.z;
        double mn = xm*ux+ym*uy+zm*uz;
        double nx = mn*ux;
        double ny = mn*uy;
        double nz = mn*uz;
        double tx = xm-nx;
        double ty = ym-ny;
        double tz = zm-nz;
        int p1 = (int)round(zi);
        int p2 = (int)round(yi);
        int p3 = (int)round(xi);
        p1 = min(p1,n1-1);
        p2 = min(p2,n2-1);
        p3 = min(p3,n3-1);
        p1 = max(p1,0);
        p2 = max(p2,0);
        p3 = max(p3,0);

        double vx = v3[p3][p2][p1];
        double vy = v2[p3][p2][p1];
        double vz = v1[p3][p2][p1];
        double uv = vx*ux+vy*uy+vz*uz;
        double xu  = xi+tx+nx+vx*10;
        double yu  = yi+ty+ny+vy*10;
        double zu  = zi+tz+nz+vz*10;
        zu = min(zu,n1-1);
        node.position.x = xu;
        node.position.y = yu;
        node.position.z = zu;
      }
    }

  }



  /**
   * Creates a new cell for this mesh (without adding it), using the specified vertex indices.
   * <br/>
   * <br/>
   * By default, this method will return a {@link Polygon3D} object, but it can be overridden to
   * provide a custom implementation.
   * 
   * @param vertexIndices
   *            the indices of the vertices forming the mesh, in counter-clockwise order
   * @return the newly created face
   */
  public Face createFace(int[] nodeIds){
    return new Face(nodeIds);
  }

  public Face createFace(int id1, int id2, int id3){
    return new Face(new int[]{id1,id2,id3});
  }



  /**
   * Adds a new face to this mesh based on the specified vertex indices. Note that these indices
   * are not checked for existence in the vertex buffer.
   * 
   * @param vertexIndices
   */
  public void addFace(int[] nodeIds){
    addFace(createFace(nodeIds));
  }

  public void addFace(int id1, int id2, int id3){
    addFace(createFace(new int[]{id1,id2,id3}));
  }
    
  /**
   * Adds a new face to this mesh. Note that the vertex indices of the specified face are not
   * checked for existence in the vertex buffer.<br/>
   * For optimal compatibility, it is recommended that the specified face is obtained via the
   * {@link #createCell(int...) createCompatibleCell} method
   * 
   * @param face
   *            the new face to add
   * @return the newly added face
   * @throws MeshTopologyException
   */
  public void addFace(Face face) {
    super.addFace(face);
    // let vertices know about their new neighbors
    for (int i=0; i<3; i++) {
      int vi = face.nodeIds[i];
      int vj = face.nodeIds[(i+1)%3];
      nodes.get(vi).neighbors.add(vj);
      nodes.get(vj).neighbors.add(vi);
    }
  }

  public void reSampleToAverageDistance(double resolution, double tolerance) 
    throws MeshTopologyException {
    double minAllowedLength = resolution * (1.0 - tolerance);
    double maxAllowedLength = resolution * (1.0 + tolerance);
    // if there are 2 faces only in the mesh, it should be destroyed
    if (faces.size() < 10) {
      System.out.println("vanished (too small)");
      throw new MeshTopologyException(ActiveMesh.this, new ActiveMesh[0]);
    }
    int cpt = -1;
    boolean noChange = false;
    while (noChange == false) {
      noChange = true;
      cpt++;
      // we are looking for 2 faces f1 = a-b-c1 and f2 = b-a-c2
      // such that they share an edge a-b that is either
      // - shorter than the low-threshold (resolution * min)
      // or
      // - longer than the high-threshold (resolution * max)
      boolean split = false, merge = false;
      int e1 = 0, e2 = 0, f1v3 = -1, f2v3 = -1;
      Face f1 = null, f2 = null;
      faceSearch:
      for (int i = 0; i<faces.size(); i++) {
        f1 = getFace(i);
        // Check all edges of this face
        int[] f1v123 = { f1.nodeIds[0], f1.nodeIds[1], f1.nodeIds[2] };
        int[] f1v231 = { f1.nodeIds[1], f1.nodeIds[2], f1.nodeIds[0] };
        int[] f1v312 = { f1.nodeIds[2], f1.nodeIds[0], f1.nodeIds[1] };
        double minEdgeLength = maxAllowedLength;
        int merge_e1 = -1, merge_e2 = -1, merge_f1v3 = -1;
        double maxEdgeLength = minAllowedLength;
        int split_e1 = -1, split_e2 = -1, split_f1v3 = -1;
        // find the extreme edge sizes
        for (int v = 0; v < 3; v++) {
          double edgeLength = nodes.get(f1v123[v]).position.distance(nodes.get(f1v231[v]).position);
          if (edgeLength < minEdgeLength) {
            minEdgeLength = edgeLength;
            merge_e1 = f1v123[v];
            merge_e2 = f1v231[v];
            merge_f1v3 = f1v312[v];
          }
          if (edgeLength > maxEdgeLength) {
            maxEdgeLength = edgeLength;
            split_e1 = f1v123[v];
            split_e2 = f1v231[v];
            split_f1v3 = f1v312[v];
          }
        }
        // favor merging over splitting
        if (merge = minEdgeLength < minAllowedLength) {
          e1 = merge_e1;
          e2 = merge_e2;
          f1v3 = merge_f1v3;
        } else if (split = maxEdgeLength > maxAllowedLength){
          e1 = split_e1;
          e2 = split_e2;
          f1v3 = split_f1v3;
        }
        if (!split && !merge) continue faceSearch;
        // a change will occur
        noChange = false;
        // => we need the second associated face for that edge
        // start from the current face => O(N)
        for (int j=i + 1; j<faces.size(); j++) {
          f2 = getFace(j);
          // check if f2 contains [v1-v2]
          if (e1 == f2.nodeIds[0] && e2 == f2.nodeIds[2]) {
            f2v3 = f2.nodeIds[1];
            break;
          } else if (e1 == f2.nodeIds[1] && e2 == f2.nodeIds[0]) {
            f2v3 = f2.nodeIds[2];
            break;
          } else if (e1 == f2.nodeIds[2] && e2 == f2.nodeIds[1]) {
            f2v3 = f2.nodeIds[0];
            break;
          }
        }
        break faceSearch;
      }
      // if everything is fine, return happily!
      if (noChange) break;
      // CASE 0: THE MESH IS INCONSISTENT //
      if (f2v3 == -1) {
      // should never happen:
      // if f2v3 does not exist, then [v1-v2] only belongs to a single face (f1)
      // => this means the mesh is inconsistent (not closed)
        System.err.println("[MESH RESAMPLING ERROR] Problematic edge: " + e1 + "-" + e2 + ":");
        System.err.print(" Node " + f1.nodeIds[0] + " has neighbors: ");
        for (Integer nn : getNode(f1.nodeIds[0]).neighbors)
             System.err.print(nn.intValue() + "  ");
        System.err.println();
        System.err.print(" Node " + f1.nodeIds[1] + " has neighbors: ");
        for (Integer nn : getNode(f1.nodeIds[1]).neighbors)
             System.err.print(nn.intValue() + "  ");
        System.err.println();
        System.err.print(" Node " + f1.nodeIds[2] + " has neighbors: ");
        for (Integer nn : getNode(f1.nodeIds[2]).neighbors)
             System.err.print(nn.intValue() + "  ");
        System.err.println();
        System.err.println("The mesh will be removed from further computations");
        throw new MeshTopologyException(ActiveMesh.this, null);
      } else
      // we're about to change the mesh structure
      // => lock everything to prevent nasty bugs
        synchronized (this) {
        // CASE 1: MERGE //
        if (merge) {
          // Check first if the edge to merge is at the base of a tetrahedron
          // if so, delete the whole tetrahedron
          if (nodes.get(f1v3).neighbors.size()==3) {
            // deleteTetrahedron(f1v3, f1v1, f1v2);
            deleteTetrahedron(f1v3, e1, e2);
          } else if (nodes.get(f2v3).neighbors.size() == 3) {
            // deleteTetrahedron(f2v3, f1v2, f1v1);
            deleteTetrahedron(f2v3, e2, e1);
          } else {
            Node v1 = nodes.get(e1);
            Node v2 = nodes.get(e2);
            // if v1 and v2 have a 3rd neighbor n in addition to f1v3 and f2v3,
            // then the mesh has reached a tubular structure.
            // => split the mesh using the virtual face [v1,v2,n]
            for (int n : v1.neighbors)
              if (v2.neighbors.contains(n) && n != f1v3 && n != f2v3) {
                splitContourAtVertices(e1, e2, n);
                // don't go further
                return;
              }
            // Now, we can confidently merge v1 and v2:
            // 1) Remove f1, f2 (they used to [v1-v2])
            // 2) Move v1 to the center of [v1-v2]
            // 3) Remove v2 from its neighborhood
            // 4) Faces pointing to v2 should now point to v1
            // 5) Delete v2
             // 1) remove the 2 old faces sharing [v1-v2]
            faces.remove(f1);
            faces.remove(f2);
            // 2a) Move v1 to the middle of v1-v2
            v1.position.interpolate(v2.position, 0.5);
            // 3) Remove v2 from its neighborhood...
            for (int n : v2.neighbors) {
              Node vn = nodes.get(n);
              vn.neighbors.remove(e2);
              // additionally, add v2's neighbors to v1
              // (except the existing ones: e1, f1v3, f2v3)
              if (n != e1 && n != f1v3 && n != f2v3) {
                v1.neighbors.add(n);
                vn.neighbors.add(e1);
              }
            }
            // 4) All faces pointing to v2 should point to v1
            for (Face f : faces) {
              f.replace(e2, e1);
            }
            // 5) delete the old vertex and notify its neighbors
            nodes.set(e2, null);
          }
        } else if (split) {
        // CASE 2: INVERT or SPLIT //
         // 1) remove the old faces
          faces.remove(f1);
          faces.remove(f2);
          // 2) the vertices won't be neighbors anymore
          nodes.get(e1).neighbors.remove(e2);
          nodes.get(e2).neighbors.remove(e1);
          // 3) invert or split?
          Node v1 = nodes.get(f1v3);
          Node v2 = nodes.get(f2v3);
          if (v1.distanceTo(v2) < maxAllowedLength && v1.neighbors.contains(f2v3)) {
            // INVERT
            // 3) create the two new faces
            addFace(f1v3, e1, f2v3);
            addFace(f2v3, e2, f1v3);
          } else {
            // SPLIT
            // 3) create a vertex in the middle of the edge
            int c = addVertexBetween(e1, e2);
            // 4) create 4 new faces around the new vertex
            addFace(e1, c, f1v3);
            addFace(f1v3, c, e2);
            addFace(e1, f2v3, c);
            addFace(c, f2v3, e2);
          }
        }
      }
      // prevent infinite loop
      if (cpt > getNumberOfNodes(false) * 2) noChange = true;
    }
  }

  /**
   * Deletes a tetrahedron from the mesh, and fill the hole with a new face
   * 
   * @param topVertex
   *            the vertex at the top of the tetrahedron
   * @param v1
   *            one of the three vertices at the base of the tetrahedron
   * @param v2
   *            another of the vertices at the base of the tetrahedron
   */
  private void deleteTetrahedron(int topNode, int v1, int v2) {
    // find the third node at the base of the tetrahedron
    int v3 = -1;
    for (int n : nodes.get(topNode).neighbors) {
      if (n != v1 && n != v2) v3 = n;
      // take the opportunity to update the neighborhood
      nodes.get(n).neighbors.remove(topNode);
    }
    for (int i=0; i<faces.size(); i++) {
      Face face = faces.get(i);
      if (face.contains(topNode)) {
        faces.remove(i--);
      }
    }
    // fill the hole
    addFace(v1, v2, v3);
    // delete the top vertex for good
    nodes.set(topNode, null);
  }



  /**
   * Splits the current contour using the 'cutting' face defined by the given vertices. <br>
   * 
   * <pre>
   * How this works:
   *  - separate all vertices on each side of the cutting face (without considering 
   *  the vertices of the cutting face), 
   *  - separate all faces touching at least one vertex of each group 
   *  (will include faces touching the cutting face),
   *  - create a contour with each group of vertices and faces,
   *  - add the cutting face and its vertices to each created contour
   * </pre>
   * 
   * @param v1
   * @param v2
   * @param v3
   * @throws MeshTopologyException
   */
  private void splitContourAtVertices(int v1, int v2, int v3) throws MeshTopologyException {
    ActiveMesh[] children = new ActiveMesh[2];
    ArrayList<Node> visitedNodes = new ArrayList<Node>(getNumberOfNodes(false));
    // mark the vertices from the cutting triangle as visited
    // => that should separate the mesh into 2 (open) meshes
    visitedNodes.add(nodes.get(v1));
    visitedNodes.add(nodes.get(v2));
    visitedNodes.add(nodes.get(v3));
    // the mesh split into two components, extract them via connected component analysis
    for (int child = 0; child < 2; child++) {
      // use the first non-visited vertex as seed
      Node seed = null;
      for (Node v : nodes)
        if (v != null && !visitedNodes.contains(v)){
          seed = v;
          break;
        }
      if (seed == null) {
        System.err.println("While splitting mesh at vertices (" + v1 + "," + v2 + "," + v3 + "):");
        System.err.println("Couldn't create child mesh #" + (child + 1) + ": no more seeds");
        throw new MeshTopologyException(ActiveMesh.this, null);
      }
      ArrayList<Face> newFaces = new ArrayList<Face>();
      ArrayList<Node> newNodes = new ArrayList<Node>();
      // create a null empty list that will be used to clone the source vertices
      for (int i=0; i<nodes.size(); i++)
        newNodes.add(null);
      extractNodes(seed, visitedNodes, nodes, newNodes);
      extractFaces(newNodes, faces, newFaces);
      // Fill the hole using the cutting face
      // First, add the 3 vertices from the cut
      for (int v:new int[]{v1,v2,v3 }) {
        Node newV = nodes.get(v).clone();
        newNodes.set(v, newV);
        // the neighbors have also been cloned
        // remove those on the wrong side of the cut
        // => they should point to null in the child's vertex list
        Iterator<Integer> it = newV.neighbors.iterator();
        while (it.hasNext()) {
          int n = it.next();
          if (n != v1 && n != v2 && n != v3 && newNodes.get(n) == null) {
            // newV should not point to n anymore
            it.remove();
          }
        }
      }
      // Add the new face used for the cut
      // => find any edge (e.g. v1-v2) to check its ordering
      for (Face f : newFaces)
        if (f.contains(v1) && f.contains(v2)) {
          // if the edge v1-v2 is counter-clockwise in f,
          // the new face must be clockwise and vice-versa
          newFaces.add(f.isEdgeOrdered(v1,v2) ? createFace(v1,v3,v2):createFace(v1,v2,v3));
          break;
        }
      try {
        ActiveMesh newMesh = getClass().newInstance();
        newMesh.setPixelSize(pixelSize);
        newMesh.setNodes(newNodes);
        newMesh.setFaces(newFaces);
        // keep meshes with more than 10 vertices
        if (newMesh.getNumberOfCells() > 10) {
          children[child] = newMesh;
         }
      } catch (Exception e) {
        throw new RuntimeException(e);
      }
    }
    if (children[0] == null) {
      if (children[1] == null) throw new MeshTopologyException(ActiveMesh.this, null);
    } else {
      if (children[1] != null) throw new MeshTopologyException(ActiveMesh.this, children);
    }
  }



  private int addVertexBetween(int v1, int v2) {
    // Create the middle vertex using v1
    Point3d p = new Point3d(nodes.get(v1).position);
    // move it half way towards v2
    p.interpolate(nodes.get(v2).position, 0.5);
    return addNode(createNode(p));
  }


  private static void extractNodes( Node seed, List<Node> visitedNodes, 
    List<Node> oldNodes, List<Node> newNodes) {
    Stack<Node> seeds = new Stack<Node>();
    seeds.add(seed);
    while (!seeds.isEmpty()) {
      Node currentNode = seeds.pop();
      // don't process a visited vertex
      if (visitedNodes.contains(currentNode)) continue;
      // mark the vertex as visited
      visitedNodes.add(currentNode);
      // extract the vertex into the new list (at the same position!!)
      newNodes.set(oldNodes.indexOf(currentNode), currentNode);
      // add the neighbors to the list of seeds
      for (int n:currentNode.neighbors) {
        try {
          seeds.push(oldNodes.get(n));
        } catch (IndexOutOfBoundsException e) {
          n = n + 1;
        }
      }
    }
  }

  private static void extractFaces(List<Node> nodesList, List<Face> oldFaces, 
    List<Face> newFaces) {
    for (Face face : oldFaces) {
    for (int i:face.nodeIds) {
      if (nodesList.get(i) != null) {
        newFaces.add(face);
        break;
      }
    }}
  }

  public void updateNormals() {
    Vector3d v31 = new Vector3d();
    Vector3d v12 = new Vector3d();
    Vector3d v23 = new Vector3d();
    for (Face f : faces) {
      // Accumulate face normals in each vertex
      // All face vertices are in the same plane
      // => use any 3 consecutive vertices
      Node v1 = nodes.get(f.nodeIds[0]);
      Node v2 = nodes.get(f.nodeIds[1]);
      Node v3 = nodes.get(f.nodeIds[2]);
      v31.sub(v1.position, v3.position);
      v12.sub(v2.position, v1.position);
      v23.sub(v3.position, v2.position);
      // normal at v1 = [v1 v2] ^ [v1 v3] = [v3 v1] ^ [v1 v2]
      v1.normal.x += v31.y * v12.z - v31.z * v12.y;
      v1.normal.y += v31.z * v12.x - v31.x * v12.z;
      v1.normal.z += v31.x * v12.y - v31.y * v12.x;
      // normal at v2 = [v2 v3] ^ [v2 v1] = [v1 v2] ^ [v2 v3]
      v2.normal.x += v12.y * v23.z - v12.z * v23.y;
      v2.normal.y += v12.z * v23.x - v12.x * v23.z;
      v2.normal.z += v12.x * v23.y - v12.y * v23.x;
      // normal at v3 = [v3 v1] ^ [v3 v2] = [v2 v3] ^ [v3 v1]
      v3.normal.x += v23.y * v31.z - v23.z * v31.y;
      v3.normal.y += v23.z * v31.x - v23.x * v31.z;
      v3.normal.z += v23.x * v31.y - v23.y * v31.x;
    }
    // Normalize the accumulated normals
    for (Node v : nodes)
      if (v != null) v.normal.normalize();
  }

}
