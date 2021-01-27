#include <float.h>
#include <assert.h>
#include "meshEdit.h"
#include "mutablePriorityQueue.h"
#include "error_dialog.h"

namespace CMU462 {

VertexIter HalfedgeMesh::splitEdge(EdgeIter e0) {
  // TODO: (meshEdit)
  // This method should split the given edge and return an iterator to the
  // newly inserted vertex. The halfedge of this vertex should point along
  // the edge that was split, rather than the new edges.

  // HALFEDGES
    if (e0->isBoundary()) return e0->halfedge()->vertex();
    if (e0->halfedge()->face()->degree() != 3 || e0->halfedge()->twin()->face()->degree() != 3) return e0->halfedge()->vertex();
    HalfedgeIter h = e0->halfedge(), t = h->twin();
    FaceIter f0 = h->face(), f1 = t->face();
    VertexIter v0 = h->vertex(), v1 = t->vertex(), v2 = h->next()->next()->vertex(),v3 = t->next()->next()->vertex();
    VertexIter v = newVertex();
    EdgeIter e1 = newEdge(), e2 = newEdge(),e3 = newEdge();
    FaceIter f2 = newFace(), f3 = newFace();
    v->position = e0->centroid();
    HalfedgeIter h1 = newHalfedge(), t1 = newHalfedge(),h2 = newHalfedge(), t2 = newHalfedge(),h3 = newHalfedge(), t3 = newHalfedge();
    h1->setNeighbors(h->next()->next(), t1, v, e1, f0);
    t1->setNeighbors(h2, h1, v2, e1, f2);
    h2->setNeighbors(h->next(), t2, v, e2, f2);
    t2->setNeighbors(h3, h2, v1, e2, f3);
    h3->setNeighbors(t->next()->next(), t3, v, e3, f3);
    t3->setNeighbors(t, h3, v3, e3, f1);
    v->halfedge() = h2;
    e1->halfedge() = h1;
    e2->halfedge() = h2;
    e3->halfedge() = h3;
    f0->halfedge() = h;
    f1->halfedge() = t;
    f2->halfedge() = t1;
    f3->halfedge() = h3;
    h2->next()->next() = t1;
    h3->next()->next() = t2;
    t->next()->next() = t3;

    h->next() = h1;

    t->vertex() = v;
    return v;
}

VertexIter HalfedgeMesh::collapseEdge(EdgeIter e) {
  // TODO: (meshEdit)
  // This method should collapse the given edge and return an iterator to
  // the new vertex created by the collapse.
   /*if (e->isBoundary()) { return e->halfedge()->vertex(); }

    HalfedgeIter h = e->halfedge();
    HalfedgeIter old_h = h; // h might change after we delete halfedge from triagnles
    HalfedgeIter h0 = h->twin();
    FaceIter f0 = h->face();
    FaceIter f1 = h0->face();
    //check if edge is in triangle
    bool right = false;
    bool left = false;
    if (f0->degree() == 3) { right = true; }
    if (f1->degree() == 3) { left = true; }

    //set the center vertex
    VertexIter v = h->vertex();
    VertexIter v0 = h0->vertex();
    v->position = e->centroid();
    v->halfedge() = h0->next()->twin()->next();
    v0->halfedge() = h->next()->twin()->next();

    if (right) {
        if (h->next()->isBoundary() || h->next()->next()->isBoundary()) {
            return e->halfedge()->vertex();
        }
        //get related halfedges
        HalfedgeIter r1 = h->next();
        HalfedgeIter r2 = r1->next();
        HalfedgeIter r3 = r2->twin();
        HalfedgeIter r4 = r1->twin();
        //get related edge
        EdgeIter er1 = r2->edge();
        EdgeIter er2 = r1->edge();
        //get related vertex
        VertexIter v1 = r2->vertex();

        //reassign value
        v1->halfedge() = r3->next();
        r3->twin() = r4;
        r4->twin() = r3;
        r4->edge() = er1;
        er1->halfedge() = r3;

        h = r4->next();

        deleteHalfedge(r1);
        deleteHalfedge(r2);
        deleteEdge(er2);
    }
    else {
        //get related halfedges
        HalfedgeIter r1 = h->next();
        HalfedgeIter r2 = h;
        while (r2->next() != h) {
            r2 = r2->next();
        }
        //reassign value
        r2->next() = r1;
        r1->vertex() = v;
        r1->face()->halfedge() = r1;

        h = r1->twin()->next();
    }

    //update the bottom vertex to the returned vertex
    //HalfedgeIter temp1 = h->next()->twin()->next();
    cout << "num of edges; " << h->vertex()->degree() << endl;
    HalfedgeIter temp1 = h;
    while (temp1 != h0) {
        temp1->vertex() = v;
        temp1 = temp1->twin()->next();
    }
    cout << "collapse 1" << endl;
    if (left) {
        if (h0->next()->isBoundary() || h0->next()->next()->isBoundary()) {
            return e->halfedge()->vertex();
        }
        //get related halfedges
        HalfedgeIter l1 = h0->next();
        HalfedgeIter l2 = l1->next();
        HalfedgeIter l3 = l2->twin();
        HalfedgeIter l4 = l1->twin();
        //get related edge
        EdgeIter el1 = l1->edge();
        EdgeIter el2 = l2->edge();
        //get related vertex
        VertexIter v2 = l2->vertex();

        //reassign value
        v2->halfedge() = l3->next();
        l3->twin() = l4;
        l4->twin() = l3;
        l3->edge() = el1;
        el1->halfedge() = l4;
        //l3->vertex() = v;
        deleteHalfedge(l1);
        deleteHalfedge(l2);
        deleteEdge(el2);
    }
    else {
        //get related halfedges
        HalfedgeIter l1 = h0->next();
        HalfedgeIter l2 = h0;
        while (l2->next() != h0) {
            l2 = l2->next();
        }
        //reassign value
        l2->next() = l1;
        l1->face()->halfedge() = l1;
    }
    //deleteHalfedge(h);
    deleteHalfedge(old_h);
    deleteHalfedge(h0);
    deleteEdge(e);
    deleteVertex(v0);
    if (right) { deleteFace(f0); }
    if (left) { deleteFace(f1); }
    cout << "collapse complete" << endl;
    return v;*/
    if (e->isBoundary()) { return e->halfedge()->vertex(); }
    HalfedgeIter h1 = e->halfedge(), h2 = h1->twin(), oh = e->halfedge();
    VertexIter v1 = h1->vertex(), v2 = h2->vertex();
    FaceIter f1 = h1->face(), f2 = h2->face() ;
    bool b1 = f1->degree() == 3, b2 = f2->degree() == 3;
    v1->position = e->centroid();
    v1->halfedge() = h1->next()->twin()->next();
    v2->halfedge() = h2->next()->twin()->next();

    //update left side face
    if (b1) {
        HalfedgeIter r1 = h1->next(), r2 = r1->twin(), r3 = r1->next(), r4 = r3->twin();
        
        EdgeIter e1 = r1->edge(), e2 = r3->edge();
        FaceIter f3 = r2->face();
        r3->next() = r2->next();
        r3->face() = f3;
        h1 = r2->next();

        f3->halfedge() = r3;

        deleteHalfedge(r1);
        deleteHalfedge(r2);
        deleteEdge(e1);
    }
    else {
        HalfedgeIter r1 = h1->next(), r2 = h1;
        while (r2->next() != h1) {
            r2 = r2->next();
        }
        r1->vertex() = v1;
        r2->next() = r1;
        r1->face()->halfedge() = r1;
        
        h1 = r1->twin()->next();
    }

    //update upper part
    while (h1 ->twin()-> next()!= h2) {
        h1->vertex() = v1;
        h1 = h1->twin()->next();
    }
    h1->vertex() = v1;
    h1 = h1->twin();

    //update right part
    if (b2) {
        HalfedgeIter r1 = h2->next(), r2 = r1->twin(), r3 = r2->next() ;
        EdgeIter e1 = r1->edge();
        h1->next() = r3;
        h1->face() = r3->face();
        
        r3->face()->halfedge() = h1;
        deleteHalfedge(r1);
        deleteHalfedge(r2);
        deleteEdge(e1);
    }
    else {
        h1->next() = h2->next();
        h1->face()->halfedge() = h2->next();
    }

    //delete 
    deleteHalfedge(oh);
    deleteHalfedge(h2);
    deleteEdge(e);
    deleteVertex(v2);
    if (b1) deleteFace(f1);
    if (b2) deleteFace(f2);
    return v1;
    
}

VertexIter HalfedgeMesh::collapseFace(FaceIter f) {
  // TODO: (meshEdit)
  // This method should collapse the given face and return an iterator to
  // the new vertex created by the collapse.
    showError("collapseface() not implemented.");
    return VertexIter();
}

FaceIter HalfedgeMesh::eraseVertex(VertexIter v) {
  // TODO: (meshEdit)
  // This method should replace the given vertex and all its neighboring
  // edges and faces with a single face, returning the new face.
    HalfedgeIter h = v->halfedge(), t = h;
    int d = v->degree();
    vector<HalfedgeIter> hfs;
    vector<EdgeIter> egs;
    vector<FaceIter> fs;
    vector<HalfedgeIter> ori_h;
    FaceIter f = newFace();
    for (int i = 0; i < d; i++) {
        hfs.push_back(h);
        hfs.push_back(h->twin());
        egs.push_back(h->edge());
        fs.push_back(h->face());
        do {
            ori_h.push_back(t->next());
            t = t->next();
        } while (t->next()->next() != h);
        h = t->next()->twin();
        t = h;
    }
    
    f->halfedge() = ori_h[0];
    for (int i = 0; i < ori_h.size(); i++) {
        ori_h[i]->next() = ori_h[(i + 1) % ori_h.size()];
        ori_h[i]->face() = f;
        ori_h[i]->vertex()->halfedge() = ori_h[i];
            
    }
    for(int i = 0; i<hfs.size(); i++) deleteHalfedge(hfs[i]);
    for (int i = 0; i < egs.size(); i++) deleteEdge(egs[i]);
    for (int i = 0; i < fs.size(); i++) deleteFace(fs[i]);
    deleteVertex(v);
    
    return f;

    


    
}

FaceIter HalfedgeMesh::eraseEdge(EdgeIter e) {
  // TODO: (meshEdit)
  // This method should erase the given edge and return an iterator to the
  // merged face.
    vector<HalfedgeIter> hf1;
    vector<HalfedgeIter> hf2;
    HalfedgeIter h0 = e->halfedge(), tmp = h0;
    FaceIter f1 = h0->face();
    do {
        hf1.push_back(tmp);
        tmp = tmp->next();
    } while (tmp != h0);
    HalfedgeIter h1 = h0->twin();
    tmp = h1;
    FaceIter f2 = h1->face();

    do {
        hf2.push_back(tmp);
        tmp = tmp->next();
    } while (tmp != h1);
    size_t s1 = hf1.size() - 1, s2 = hf2.size() - 1;
    hf1[s1]->next() = hf2[1];
    hf2[s2]->next() = hf1[1];

    for (HalfedgeIter& h : hf2)  h->face() = f1;
    deleteFace(f2);
    deleteHalfedge(hf1[0]);
    deleteHalfedge(hf2[0]);
    deleteEdge(e);
    return f1;
}

EdgeIter HalfedgeMesh::flipEdge(EdgeIter e0) {
  // TODO: (meshEdit)
  // This method should flip the given edge and return an iterator to the
  // flipped edge.
    if (e0->isBoundary()) return e0;

    HalfedgeIter h0 = e0->halfedge();
    HalfedgeIter h1 = h0->next();
    HalfedgeIter h2 = h1->next();
    HalfedgeIter h3 = h0->twin();
    HalfedgeIter h4 = h3->next();
    HalfedgeIter h5 = h4->next();
    HalfedgeIter h6 = h1->twin();
    HalfedgeIter h7 = h2->twin();
    HalfedgeIter h8 = h4->twin();
    HalfedgeIter h9 = h5->twin();

    // VERTICES
    VertexIter v0 = h0->vertex();
    VertexIter v1 = h3->vertex();
    VertexIter v2 = h8->vertex();
    VertexIter v3 = h2->vertex();
    // ...you fill in the rest!...

    // EDGES
    EdgeIter e1 = h5->edge();
    EdgeIter e2 = h4->edge();
    EdgeIter e3 = h2->edge();
    EdgeIter e4 = h1->edge();

    // ...you fill in the rest!...

    // FACES
    FaceIter f0 = h0->face();
    FaceIter f1 = h3->face();

    h0->next() = h1;
    h0->twin() = h3;
    h0->vertex() = v2;
    h0->edge() = e0;
    h0->face() = f0;

    h1->next() = h2;
    h1->twin() = h7;
    h1->vertex() = v3;
    h1->edge() = e3;
    h1->face() = f0;

    h2->next() = h0;
    h2->twin() = h8;
    h2->vertex() = v0;
    h2->edge() = e2;
    h2->face() = f0;

    h3->next() = h4;
    h3->twin() = h0;
    h3->vertex() = v3;
    h3->edge() = e0;
    h3->face() = f1;

    h4->next() = h5;
    h4->twin() = h9;
    h4->vertex() = v2;
    h4->edge() = e1;
    h4->face() = f1;

    h5->next() = h3;
    h5->twin() = h6;
    h5->vertex() = v1;
    h5->edge() = e4;
    h5->face() = f1;

    h6->next() = h6->next();
    h6->twin() = h5;
    h6->vertex() = v3;
    h6->edge() = e4;
    h6->face() = h6->face();

    h7->next() = h7->next();
    h7->twin() = h1;
    h7->vertex() = v0;
    h7->edge() = e3;
    h7->face() = h7->face();

    h8->next() = h8->next();
    h8->twin() = h2;
    h8->vertex() = v2;
    h8->edge() = e2;
    h8->face() = h8->face();

    h9->next() = h9->next();
    h9->twin() = h4;
    h9->vertex() = v1;
    h9->edge() = e1;
    h9->face() = h9->face();


    // VERTICES
    v0->halfedge() = h2;
    v1->halfedge() = h5;
    v2->halfedge() = h4;
    v3->halfedge() = h1;
    // ...you fill in the rest!...

    // EDGES
    e0->halfedge() = h0;
    e1->halfedge() = h4;
    e2->halfedge() = h2;
    e3->halfedge() = h1;
    e4->halfedge() = h5;


    // ...you fill in the rest!...

    // FACES
    f0->halfedge() = h0;
    f1->halfedge() = h3;
    
    return e0;
}

void HalfedgeMesh::subdivideQuad(bool useCatmullClark) {
  // Unlike the local mesh operations (like bevel or edge flip), we will perform
  // subdivision by splitting *all* faces into quads "simultaneously."  Rather
  // than operating directly on the halfedge data structure (which as you've
  // seen
  // is quite difficult to maintain!) we are going to do something a bit nicer:
  //
  //    1. Create a raw list of vertex positions and faces (rather than a full-
  //       blown halfedge mesh).
  //
  //    2. Build a new halfedge mesh from these lists, replacing the old one.
  //
  // Sometimes rebuilding a data structure from scratch is simpler (and even
  // more
  // efficient) than incrementally modifying the existing one.  These steps are
  // detailed below.

  // TODO Step I: Compute the vertex positions for the subdivided mesh.  Here
  // we're
  // going to do something a little bit strange: since we will have one vertex
  // in
  // the subdivided mesh for each vertex, edge, and face in the original mesh,
  // we
  // can nicely store the new vertex *positions* as attributes on vertices,
  // edges,
  // and faces of the original mesh.  These positions can then be conveniently
  // copied into the new, subdivided mesh.
  // [See subroutines for actual "TODO"s]
  if (useCatmullClark) {
    computeCatmullClarkPositions();
  } else {
    computeLinearSubdivisionPositions();
  }

  // TODO Step II: Assign a unique index (starting at 0) to each vertex, edge,
  // and
  // face in the original mesh.  These indices will be the indices of the
  // vertices
  // in the new (subdivided mesh).  They do not have to be assigned in any
  // particular
  // order, so long as no index is shared by more than one mesh element, and the
  // total number of indices is equal to V+E+F, i.e., the total number of
  // vertices
  // plus edges plus faces in the original mesh.  Basically we just need a
  // one-to-one
  // mapping between original mesh elements and subdivided mesh vertices.
  // [See subroutine for actual "TODO"s]
  assignSubdivisionIndices();

  // TODO Step III: Build a list of quads in the new (subdivided) mesh, as
  // tuples of
  // the element indices defined above.  In other words, each new quad should be
  // of
  // the form (i,j,k,l), where i,j,k and l are four of the indices stored on our
  // original mesh elements.  Note that it is essential to get the orientation
  // right
  // here: (i,j,k,l) is not the same as (l,k,j,i).  Indices of new faces should
  // circulate in the same direction as old faces (think about the right-hand
  // rule).
  // [See subroutines for actual "TODO"s]
  vector<vector<Index> > subDFaces;
  vector<Vector3D> subDVertices;
  buildSubdivisionFaceList(subDFaces);
  buildSubdivisionVertexList(subDVertices);

  // TODO Step IV: Pass the list of vertices and quads to a routine that clears
  // the
  // internal data for this halfedge mesh, and builds new halfedge data from
  // scratch,
  // using the two lists.
  rebuild(subDFaces, subDVertices);
}

/**
 * Compute new vertex positions for a mesh that splits each polygon
 * into quads (by inserting a vertex at the face midpoint and each
 * of the edge midpoints).  The new vertex positions will be stored
 * in the members Vertex::newPosition, Edge::newPosition, and
 * Face::newPosition.  The values of the positions are based on
 * simple linear interpolation, e.g., the edge midpoints and face
 * centroids.
 */
void HalfedgeMesh::computeLinearSubdivisionPositions() {
  // TODO For each vertex, assign Vertex::newPosition to
  // its original position, Vertex::position.
    for (VertexIter v = verticesBegin(); v != verticesEnd(); ++v) {
        v->newPosition = v->position;
    }
  // TODO For each edge, assign the midpoint of the two original
  // positions to Edge::newPosition.
    for (EdgeIter e = edgesBegin(); e != edgesEnd(); ++e) {
        e->newPosition = e->centroid();
    }
  // TODO For each face, assign the centroid (i.e., arithmetic mean)
  // of the original vertex positions to Face::newPosition.  Note
  // that in general, NOT all faces will be triangles!
    for (FaceIter f = facesBegin(); f != facesEnd(); ++f) {
        f->newPosition = f->centroid();
    }
}

/**
 * Compute new vertex positions for a mesh that splits each polygon
 * into quads (by inserting a vertex at the face midpoint and each
 * of the edge midpoints).  The new vertex positions will be stored
 * in the members Vertex::newPosition, Edge::newPosition, and
 * Face::newPosition.  The values of the positions are based on
 * the Catmull-Clark rules for subdivision.
 */
void HalfedgeMesh::computeCatmullClarkPositions() {
  // TODO The implementation for this routine should be
  // a lot like HalfedgeMesh::computeLinearSubdivisionPositions(),
  // except that the calculation of the positions themsevles is
  // slightly more involved, using the Catmull-Clark subdivision
  // rules. (These rules are outlined in the Developer Manual.)

  // TODO face
    for (FaceIter f = facesBegin(); f != facesEnd(); f++) {
        f->newPosition = f->centroid();
    }
  // TODO edges
    for (EdgeIter e = edgesBegin(); e != edgesEnd(); e++) {
        HalfedgeIter h = e->halfedge(), t = h->twin();
        VertexIter v1 = h->vertex(), v2 = t->vertex();
        FaceIter f1 = h->face(), f2 = t->face();
        e->newPosition = (v1->position + v2->position + f1->centroid() + f2->centroid()) / 4;
    }
  // TODO vertices

    for (VertexIter v = verticesBegin(); v != verticesEnd(); v++) {
        HalfedgeIter h = v->halfedge();
        Vector3D Q = Vector3D(), R = Vector3D();
        int n = v->degree();
        for (int i = 0; i < n; i++) {
            Q += h->face()->centroid();
            R += h->edge()->centroid();
            h = h->twin()->next();
        }
        Q /= n;
        R /= n;
        v->newPosition = (Q + 2 * R + (n - 3) * v->position) / n;
    }
}

/**
 * Assign a unique integer index to each vertex, edge, and face in
 * the mesh, starting at 0 and incrementing by 1 for each element.
 * These indices will be used as the vertex indices for a mesh
 * subdivided using Catmull-Clark (or linear) subdivision.
 */
void HalfedgeMesh::assignSubdivisionIndices() {
  // TODO Start a counter at zero; if you like, you can use the
  // "Index" type (defined in halfedgeMesh.h)
    int idx = 0;
  // TODO Iterate over vertices, assigning values to Vertex::index
    for (VertexIter v = verticesBegin(); v != verticesEnd(); ++v) {
        v->index = idx;
        idx++;
    }
  // TODO Iterate over edges, assigning values to Edge::index
    for (EdgeIter e = edgesBegin(); e != edgesEnd(); ++e) {
        e->index = idx;
        idx++;
    }
  // TODO Iterate over faces, assigning values to Face::index
    for (FaceIter f = facesBegin(); f != facesEnd(); ++f) {
        f->index = idx;
        idx++;
    }
}

/**
 * Build a flat list containing all the vertex positions for a
 * Catmull-Clark (or linear) subdivison of this mesh.  The order of
 * vertex positions in this list must be identical to the order
 * of indices assigned to Vertex::newPosition, Edge::newPosition,
 * and Face::newPosition.
 */
void HalfedgeMesh::buildSubdivisionVertexList(vector<Vector3D>& subDVertices) {
  // TODO Resize the vertex list so that it can hold all the vertices.
    
    
  // TODO Iterate over vertices, assigning Vertex::newPosition to the
  // appropriate location in the new vertex list.
    
  // TODO Iterate over edges, assigning Edge::newPosition to the appropriate
  // location in the new vertex list.
    
    
  // TODO Iterate over faces, assigning Face::newPosition to the appropriate
  // location in the new vertex list.
    for (VertexIter v = verticesBegin(); v != verticesEnd(); ++v) {
        subDVertices.push_back(v->newPosition);
    }

    for (EdgeIter e = edgesBegin(); e != edgesEnd(); ++e) {
        subDVertices.push_back(e->newPosition);
    }

    for (FaceIter f = facesBegin(); f != facesEnd(); ++f) {
        subDVertices.push_back(f->newPosition);
    }
  
}

/**
 * Build a flat list containing all the quads in a Catmull-Clark
 * (or linear) subdivision of this mesh.  Each quad is specified
 * by a vector of four indices (i,j,k,l), which come from the
 * members Vertex::index, Edge::index, and Face::index.  Note that
 * the ordering of these indices is important because it determines
 * the orientation of the new quads; it is also important to avoid
 * "bowties."  For instance, (l,k,j,i) has the opposite orientation
 * of (i,j,k,l), and if (i,j,k,l) is a proper quad, then (i,k,j,l)
 * will look like a bowtie.
 */
void HalfedgeMesh::buildSubdivisionFaceList(vector<vector<Index> >& subDFaces) {
  // TODO This routine is perhaps the most tricky step in the construction of
  // a subdivision mesh (second, perhaps, to computing the actual Catmull-Clark
  // vertex positions).  Basically what you want to do is iterate over faces,
  // then  for each face, append N quads to the list (where N is the
  // degree of the face).  For this routine, it may be more convenient to simply
  // append quads to the end of the list (rather than allocating it ahead of
  // time), though YMMV.  You can of course iterate around a face by starting
  // with its first halfedge and following the "next" pointer until you get
  // back to the beginning.  The tricky part is making sure you grab the right
  // indices in the right order---remember that there are indices on vertices,
  // edges, AND faces of the original mesh.  All of these should get used.  Also
  // remember that you must have FOUR indices per face, since you are making a
  // QUAD mesh!

  // TODO iterate over faces
  // TODO loop around face
  // TODO build lists of four indices for each sub-quad
  // TODO append each list of four indices to face list
    vector<Index> l(4); 
    HalfedgeIter h1, h2, oh;
    for (FaceIter f = facesBegin(); f != facesEnd(); f++) {
        h1 = f->halfedge(), h2 = h1->next(), oh = h1;
        do{
            EdgeIter e1 = h1->edge(), e2 = h2->edge();
            VertexIter v = h2->vertex();
            l[0] = f->index;
            l[1] = e1->index;
            l[2] = v->index;
            l[3] = e2->index;
            subDFaces.push_back(l);
            h1 = h1->next();
            h2 = h2->next();   
        } while (h1 != oh);
    }

}

FaceIter HalfedgeMesh::bevelVertex(VertexIter v) {
  // TODO This method should replace the vertex v with a face, corresponding to
  // a bevel operation. It should return the new face.  NOTE: This method is
  // responsible for updating the *connectivity* of the mesh only---it does not
  // need to update the vertex positions.  These positions will be updated in
  // HalfedgeMesh::bevelVertexComputeNewPositions (which you also have to
  // implement!)
    int d = v->degree();
    vector<HalfedgeIter> ori_h;
    vector<VertexIter> ori_v;
    vector<FaceIter> ori_f;
    HalfedgeIter h = v->halfedge(),oh= h;
    do {
        HalfedgeIter h1 = h->twin();
        ori_h.push_back(h);
        ori_v.push_back(h1->vertex());
        ori_f.push_back(h->face());
        h = h1->next();
    } while (h != oh);
    vector<HalfedgeIter> hfs;
    vector<HalfedgeIter> twins;
    vector<EdgeIter> egs;
    vector<VertexIter> vs;
    FaceIter f = newFace();
    for(int i = 0; i<d ; i++){
        vs.push_back(newVertex());
        hfs.push_back(newHalfedge());
        twins.push_back(newHalfedge());
        egs.push_back(newEdge());
    }
    for (int i = 0; i < d; i++) {
        hfs[i]->twin() = twins[i];
        hfs[i]->vertex() = vs[i];
        hfs[i]->next() = hfs[(i + d - 1) % d];
        hfs[i]->face() = f;
        hfs[i]->edge() = egs[i];
        twins[i]->twin() = hfs[i];
        twins[i]->vertex() = vs[(i + d - 1) % d];
        twins[i]->next() = ori_h[i];
        twins[i]->face() = ori_f[i]; 
        twins[i]->edge() = egs[i];
        ori_h[i]->twin()->next() = twins[(i + 1) % d];
        ori_h[i]->vertex() = vs[i];
        egs[i]->halfedge() = hfs[i];
        vs[i]->halfedge() = ori_h[i];
        vs[i]->position = v->position;
    }
    f->halfedge() = hfs[0];
    deleteVertex(v);
    return f;
}

FaceIter HalfedgeMesh::bevelEdge(EdgeIter e) {
  // TODO This method should replace the edge e with a face, corresponding to a
  // bevel operation. It should return the new face.  NOTE: This method is
  // responsible for updating the *connectivity* of the mesh only---it does not
  // need to update the vertex positions.  These positions will be updated in
  // HalfedgeMesh::bevelEdgeComputeNewPositions (which you also have to
  // implement!)

  showError("bevelEdge() not implemented.");
  return facesBegin();
}

FaceIter HalfedgeMesh::bevelFace(FaceIter f) {
  // TODO This method should replace the face f with an additional, inset face
  // (and ring of faces around it), corresponding to a bevel operation. It
  // should return the new face.  NOTE: This method is responsible for updating
  // the *connectivity* of the mesh only---it does not need to update the vertex
  // positions.  These positions will be updated in
  // HalfedgeMesh::bevelFaceComputeNewPositions (which you also have to
  // implement!)

  showError("bevelFace() not implemented.");
  return facesBegin();
}


void HalfedgeMesh::bevelFaceComputeNewPositions(
    vector<Vector3D>& originalVertexPositions,
    vector<HalfedgeIter>& newHalfedges, double normalShift,
    double tangentialInset) {
  // TODO Compute new vertex positions for the vertices of the beveled face.
  //
  // These vertices can be accessed via newHalfedges[i]->vertex()->position for
  // i = 1, ..., newHalfedges.size()-1.
  //
  // The basic strategy here is to loop over the list of outgoing halfedges,
  // and use the preceding and next vertex position from the original mesh
  // (in the originalVertexPositions array) to compute an offset vertex
  // position.
  //
  // Note that there is a 1-to-1 correspondence between halfedges in
  // newHalfedges and vertex positions
  // in orig.  So, you can write loops of the form
  //
  //for( int i = 0; i < newHalfedges.size(); hs++ )
  // {
  //    Vector3D pi = originalVertexPositions[i]; // get the original vertex
  //    position correponding to vertex i
  // }
  //
    FaceIter f = newHalfedges[0]->face();
    Vector3D normal = f->normal();
    for (int i = 0; i < newHalfedges.size(); i++) {
        Vector3D pi = originalVertexPositions[i];
        Vector3D p = newHalfedges[i]->twin()->vertex()->position;
        newHalfedges[i]->vertex()->position = (pi * (1 - tangentialInset) + p * tangentialInset)+normalShift*normal;
    }

}

void HalfedgeMesh::bevelVertexComputeNewPositions(
    Vector3D originalVertexPosition, vector<HalfedgeIter>& newHalfedges,
    double tangentialInset) {
  // TODO Compute new vertex positions for the vertices of the beveled vertex.
  //
  // These vertices can be accessed via newHalfedges[i]->vertex()->position for
  // i = 1, ..., hs.size()-1.
  //
  // The basic strategy here is to loop over the list of outgoing halfedges,
  // and use the preceding and next vertex position from the original mesh
  // (in the orig array) to compute an offset vertex position.
    for (int i = 0; i < newHalfedges.size(); i++) {
        Vector3D pi = originalVertexPosition;
        Vector3D p = newHalfedges[i]->twin()->vertex()->position;
        newHalfedges[i]->vertex()->position += (p - pi) * tangentialInset;
    }
}

void HalfedgeMesh::bevelEdgeComputeNewPositions(
    vector<Vector3D>& originalVertexPositions,
    vector<HalfedgeIter>& newHalfedges, double tangentialInset) {
  // TODO Compute new vertex positions for the vertices of the beveled edge.
  //
  // These vertices can be accessed via newHalfedges[i]->vertex()->position for
  // i = 1, ..., newHalfedges.size()-1.
  //
  // The basic strategy here is to loop over the list of outgoing halfedges,
  // and use the preceding and next vertex position from the original mesh
  // (in the orig array) to compute an offset vertex position.
  //
  // Note that there is a 1-to-1 correspondence between halfedges in
  // newHalfedges and vertex positions
  // in orig.  So, you can write loops of the form
  //
  // for( int i = 0; i < newHalfedges.size(); i++ )
  // {
  //    Vector3D pi = originalVertexPositions[i]; // get the original vertex
  //    position correponding to vertex i
  // }
  //

}

void HalfedgeMesh::splitPolygons(vector<FaceIter>& fcs) {
  for (auto f : fcs) splitPolygon(f);
}

void HalfedgeMesh::splitPolygon(FaceIter f) {
  // TODO: (meshedit) 
  // Triangulate a polygonal face
    HalfedgeIter h = f->halfedge(), next = h->next(), oh = h;
    VertexIter v = h->vertex();
    int d = f->degree();
    
    for (int i = 0; i < d-3; i++) {
        HalfedgeIter n = newHalfedge();
        HalfedgeIter t = newHalfedge();
        FaceIter fn = newFace();
        EdgeIter e = newEdge();
        n->twin() = t;
        n->next() = h;
        n->vertex() = next->twin()->vertex();
        n->edge() = e;
        n->face() = fn;

        t->twin() = n;
        t->next() = next->next();
        t->vertex() = v;
        t->edge() = e;
        t->face() = next->next()->face();

        fn->halfedge() = n;
        e->halfedge() = n;
        h->face() = fn;
        next->face() = fn;
        next->next() = n;
        h = t;
        next = h->next(); 
    }
    next->next()->next() = h;
    f->halfedge() = h;
}

EdgeRecord::EdgeRecord(EdgeIter& _edge) : edge(_edge) {
  // TODO: (meshEdit)
  // Compute the combined quadric from the edge endpoints.
  // -> Build the 3x3 linear system whose solution minimizes the quadric error
  //    associated with these two endpoints.
  // -> Use this system to solve for the optimal position, and store it in
  //    EdgeRecord::optimalPoint.
  // -> Also store the cost associated with collapsing this edg in
  //    EdgeRecord::Cost.
}

void MeshResampler::upsample(HalfedgeMesh& mesh)
// This routine should increase the number of triangles in the mesh using Loop
// subdivision.
{
  // TODO: (meshEdit)
  // Compute new positions for all the vertices in the input mesh, using
  // the Loop subdivision rule, and store them in Vertex::newPosition.
  // -> At this point, we also want to mark each vertex as being a vertex of the
  //    original mesh.
  // -> Next, compute the updated vertex positions associated with edges, and
  //    store it in Edge::newPosition.
  // -> Next, we're going to split every edge in the mesh, in any order.  For
  //    future reference, we're also going to store some information about which
  //    subdivided edges come from splitting an edge in the original mesh, and
  //    which edges are new, by setting the flat Edge::isNew. Note that in this
  //    loop, we only want to iterate over edges of the original mesh.
  //    Otherwise, we'll end up splitting edges that we just split (and the
  //    loop will never end!)
  // -> Now flip any new edge that connects an old and new vertex.
  // -> Finally, copy the new vertex positions into final Vertex::position.

  // Each vertex and edge of the original surface can be associated with a
  // vertex in the new (subdivided) surface.
  // Therefore, our strategy for computing the subdivided vertex locations is to
  // *first* compute the new positions
  // using the connectity of the original (coarse) mesh; navigating this mesh
  // will be much easier than navigating
  // the new subdivided (fine) mesh, which has more elements to traverse.  We
  // will then assign vertex positions in
  // the new mesh based on the values we computed for the original mesh.

  // Compute updated positions for all the vertices in the original mesh, using
  // the Loop subdivision rule.

  // Next, compute the updated vertex positions associated with edges.

  // Next, we're going to split every edge in the mesh, in any order.  For
  // future
  // reference, we're also going to store some information about which
  // subdivided
  // edges come from splitting an edge in the original mesh, and which edges are
  // new.
  // In this loop, we only want to iterate over edges of the original
  // mesh---otherwise,
  // we'll end up splitting edges that we just split (and the loop will never
  // end!)

  // Finally, flip any new edge that connects an old and new vertex.

  // Copy the updated vertex positions to the subdivided mesh.
  showError("upsample() not implemented.");
}

void MeshResampler::downsample(HalfedgeMesh& mesh) {
  // TODO: (meshEdit)
  // Compute initial quadrics for each face by simply writing the plane equation
  // for the face in homogeneous coordinates. These quadrics should be stored
  // in Face::quadric
  // -> Compute an initial quadric for each vertex as the sum of the quadrics
  //    associated with the incident faces, storing it in Vertex::quadric
  // -> Build a priority queue of edges according to their quadric error cost,
  //    i.e., by building an EdgeRecord for each edge and sticking it in the
  //    queue.
  // -> Until we reach the target edge budget, collapse the best edge. Remember
  //    to remove from the queue any edge that touches the collapsing edge
  //    BEFORE it gets collapsed, and add back into the queue any edge touching
  //    the collapsed vertex AFTER it's been collapsed. Also remember to assign
  //    a quadric to the collapsed vertex, and to pop the collapsed edge off the
  //    top of the queue.
  showError("downsample() not implemented.");
}

void MeshResampler::resample(HalfedgeMesh& mesh) {
  // TODO: (meshEdit)
  // Compute the mean edge length.
  // Repeat the four main steps for 5 or 6 iterations
  // -> Split edges much longer than the target length (being careful about
  //    how the loop is written!)
  // -> Collapse edges much shorter than the target length.  Here we need to
  //    be EXTRA careful about advancing the loop, because many edges may have
  //    been destroyed by a collapse (which ones?)
  // -> Now flip each edge if it improves vertex degree
  // -> Finally, apply some tangential smoothing to the vertex positions
  showError("resample() not implemented.");
}

}  // namespace CMU462
