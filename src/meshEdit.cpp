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
    HalfedgeIter h = e0->halfedge(), t = h->twin(), r1 = h->next(), r2 = r1->next(),r3= t->next(), r4 = r3->next();
    FaceIter f0 = h->face(), f1 = t->face();
    VertexIter v0 = h->vertex(), v1 = t->vertex(), v2 = r2->vertex(),v3 = r4->vertex();
    VertexIter v = newVertex();
    EdgeIter e1 = newEdge(), e2 = newEdge(),e3 = newEdge();
    v->isNew = true;
    e1->isNew = true;
    e2->isNew = true;
    e3->isNew = true;
    FaceIter f2 = newFace(), f3 = newFace();
    v->position = e0->centroid();
    HalfedgeIter h1 = newHalfedge(), t1 = newHalfedge(),h2 = newHalfedge(), t2 = newHalfedge(),h3 = newHalfedge(), t3 = newHalfedge();
    h1->setNeighbors(r2, t1, v, e1, f0);
    t1->setNeighbors(h2, h1, v2, e1, f2);
    h2->setNeighbors(r1, t2, v, e2, f2);
    t2->setNeighbors(h3, h2, v1, e2, f3);
    h3->setNeighbors(r4, t3, v, e3, f3);
    t3->setNeighbors(t, h3, v3, e3, f1);
    v->halfedge() = t;
    e1->halfedge() = h1;
    e2->halfedge() = h2;
    e3->halfedge() = h3;
    f0->halfedge() = h;
    f1->halfedge() = t;
    f2->halfedge() = t1;
    f3->halfedge() = h3;
    r1->next() = t1;
    r1->face() = f2;
    r4->next() = t2;
    r4->face() = f3;
    r3->next() = t3;
    h->next() = h1;
   
    t->vertex() = v;
    return v;
}

VertexIter HalfedgeMesh::collapseEdge(EdgeIter e) {
  // TODO: (meshEdit)
  // This method should collapse the given edge and return an iterator to
  // the new vertex created by the collapse.
    
    if (e->isBoundary()) { return e->halfedge()->vertex(); }
    
    HalfedgeIter h1 = e->halfedge(), h2 = h1->twin(), tmp = h1;
    VertexIter v1 = h1->vertex(), v2 = h2->vertex();
    FaceIter f1 = h1->face(), f2 = h2->face() ;
    bool b1 = f1->degree() == 3, b2 = f2->degree() == 3;
    
    HalfedgeIter h = h1->next()->twin()->next();
    //update left side face
    if (b1) {
        
        HalfedgeIter r1 = tmp->next(), r2 = r1->twin(), r3 = r1->next(), r4 = r3->twin();
       
        HalfedgeIter r5 = r2->next();
        while (r5->next() != r2) {
            r5 = r5->next();
        }
        
        
        EdgeIter e1 = r1->edge(), e2 = r3->edge();
        FaceIter f3 = r2->face();
       
        r3->next() = r2->next();
        r5->next() = r3;
        r3->face() = f3;
        tmp = r3->next();
        
        f3->halfedge() = r3;
        
        deleteHalfedge(r1);
        deleteHalfedge(r2);
        deleteEdge(e1);
        
        

    }
    else {
        HalfedgeIter r1 = tmp->next(), r2 = tmp;
        while (r2->next() != tmp) {
            r2 = r2->next();
        }
        r1->vertex() = v1;
        r2->next() = r1;
        r1->face()->halfedge() = r1;
        
        tmp = r1->twin()->next();
    }
    
    
    //update upper part
    while (tmp ->twin()-> next()!= h2) {
        tmp->vertex() = v1;
        tmp = tmp->twin()->next();
    }
    //update right part
    if (b2) {
        HalfedgeIter r1 = h2->next(), r2 = r1->twin(), r3 = tmp, r4 = r3->twin();
        HalfedgeIter r5 = r3->next();
        FaceIter f4 = r5->face();
        while (r5->next() != r3) {
            r5 = r5->next();
        }
        if (r5->next() == r3) cout << "sdsdsfsdfsdfs" << endl;
        EdgeIter e2 = r3->edge();
        r1->next() = r3->next();
        r5->next() = r1;
        r1->face() = f4;
        f4->halfedge() = r1;
        checkConsistency();
        deleteHalfedge(r3);
        deleteHalfedge(r4);
        deleteEdge(e2);
        
    }
    else {
        tmp->vertex() = v1;
        tmp = tmp->twin();
        tmp->next() = h2->next();
        tmp->face()->halfedge() = h2->next();
    }
    v1->halfedge() = h;
    v1->position = e->centroid();
    if(b1) deleteFace(f1);
    if(b2) deleteFace(f2);
    
    //delete 
    deleteHalfedge(h1);
    deleteHalfedge(h2);
    deleteEdge(e);
    deleteVertex(v2);
    checkConsistency();
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
    HalfedgeIter h0 = e->halfedge(), h1 = h0->twin(), h0_prev = h0->next(), h1_prev = h1->next(), h2 = h0->next(), h3 = h1->next();
    VertexIter v0 = h0->vertex(), v1 = h1->vertex();
    while (h0_prev->next() != h0) 
        h0_prev = h0_prev->next();
    while (h1_prev->next() != h1) 
        h1_prev = h1_prev->next();
    FaceIter f1 = h0->face(), f2 = h1->face();
    h1_prev->next() = h2;
    h0_prev->next() = h3;
    HalfedgeIter tmp = h3;
    while (tmp->next() != h3) {
        tmp->face() = f1;
        tmp = tmp->next();
    }
    v0->halfedge() = h3;
    v1->halfedge() = h2;
    f1->halfedge() = h2;
    //tmp->face() = f1;
    deleteFace(f2);
    deleteHalfedge(h0);
    deleteHalfedge(h1);
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
    VertexIter v2 = h5->vertex();
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

    h0->next() = h2;
    h0->vertex() = v2;
    h3->next() = h5;
    h3->vertex() = v3;
    h2->next() = h4;
    h4->next() = h0;
    h1->next() = h3;
    h1->face() = f1;
    h4->face() = f0;
    h5->next() = h1;
    v0->halfedge() = h4;
    v1->halfedge() = h1;
    v2->halfedge() = h5;
    v3->halfedge() = h2;

    f1->halfedge() = h3;
    f0->halfedge() = h0;
            
    return e0;
    /*if (e0->isBoundary()) {
        return e0;
    }

    HalfedgeIter h0 = e0->halfedge();
    HalfedgeIter h1 = h0->next();
    HalfedgeIter h2 = h1->next();
    HalfedgeIter h3 = h0->twin();
    HalfedgeIter h4 = h3->next();
    HalfedgeIter h5 = h4->next();
    HalfedgeIter h6 = h1->twin();
    //HalfedgeIter h7 = h2->twin();
    HalfedgeIter h8 = h4->twin();
    //HalfedgeIter h9 = h5->twin();
    HalfedgeIter h10 = h3;
    while (h10->next() != h3) {
        h10 = h10->next();
    }
    HalfedgeIter h11 = h0;
    while (h11->next() != h0) {
        h11 = h11->next();
    }

    VertexIter v0 = h0->vertex();
    VertexIter v1 = h3->vertex();
    VertexIter v2 = h6->vertex();
    VertexIter v3 = h8->vertex();



    FaceIter f0 = h0->face();
    FaceIter f1 = h3->face();


    v0->halfedge() = h4;
    v1->halfedge() = h1;
    v2->halfedge() = h2;
    v3->halfedge() = h5;
    f0->halfedge() = h2;
    f1->halfedge() = h5;

    h0->next() = h2;
    h0->twin() = h3;
    h0->vertex() = v3;
    h0->face() = f0;
    h1->next() = h3;
    //h1->twin() = h6;
    //h1->vertex() = v1;
    h1->face() = f1;
    if (h2->next() == h0) {
        h2->next() = h4;
    }
    //h2->twin() = h7;
    //h2->vertex() = v2;
    //h2->face() = f0;
    h3->next() = h5;
    h3->twin() = h0;
    h3->vertex() = v2;
    h3->face() = f1;
    h4->next() = h0;
    //h4->twin() = h8;
    //h4->vertex() = v0;
    h4->face() = f0;
    h10->next() = h1;
    //h10->twin()
    //h10->vertex()
    //h10->face()
    h11->next() = h4;

    // TODO: (meshEdit)
    // This method should flip the given edge and return an iterator to the
    // flipped edge.


    return e0;*/
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
    HalfedgeIter h = _edge->halfedge(), tmp = h;
    VertexIter i = h->vertex(), j = h->twin()->vertex();
    Matrix4x4 kij = i->quadric + j->quadric;
    Matrix3x3 A = Matrix3x3();
    A.zero();
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            A(i, j) = kij(i, j);
        }
    }
    Vector3D b = Vector3D(-kij(0, 3), -kij(1, 3), -kij(2, 3));
    optimalPoint = A.inv() * b;
    Vector4D x = Vector4D(optimalPoint, 1);
    Vector4D kx = kij * x;
    score = dot(x,kx);
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

    for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
        v->isNew = false;
        v->halfedge()->edge()->isNew = false;
        HalfedgeIter h = v->halfedge();
        int n = v->degree();
        double u = n==3? 3/16 :3 / (8 * n);

        v->newPosition = (1. - n * u) * v->position + (u * n) * v->neighborhoodCentroid();
    }

  
  // Next, compute the updated vertex positions associated with edges.
    for(EdgeIter e = mesh.edgesBegin(); e!= mesh.edgesEnd(); e++){
        e->isNew = false;
        HalfedgeIter h = e->halfedge(), t = h->twin();
        if (h->face()->degree() != 3) continue;
        VertexIter v1 = h->vertex(), v2 = t->vertex(), v3 = h->next()->next()->vertex(), v4 = t->next()->next()->vertex();
        e->newPosition = (3./8.) * (v1->position + v2->position)  + (1./8.)*(v3->position + v4->position);
    }
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
    int n = mesh.nEdges();
    EdgeIter e = mesh.edgesBegin();
    for (int i = 0; i < n; i++) {

        // get the next edge NOW!
        EdgeIter nextEdge = e;
        nextEdge++;

        // now, even if splitting the edge deletes it...
        if(!e->isNew){
            VertexIter v_temp = mesh.splitEdge(e);
            v_temp->newPosition = e->newPosition;
            v_temp->isNew = true;
            HalfedgeIter he_temp = v_temp->halfedge()->twin()->next();
            while (he_temp != v_temp->halfedge()) {
                he_temp->edge()->isNew = true;
                he_temp = he_temp->twin()->next();
            }

            v_temp->halfedge()->edge()->isNew = false;
            v_temp->halfedge()->twin()->next()->twin()->next()->edge()->isNew = false;
        }

        // ...we still have a valid reference to the next edge.
        e = nextEdge;
    }
  // Finally, flip any new edge that connects an old and new vertex.
    n = mesh.nEdges();
    e = mesh.edgesBegin();
    for (int i = 0; i < n; i++) {
        HalfedgeIter h = e->halfedge();
        // get the next edge NOW!
        EdgeIter nextEdge = e;
        nextEdge++;

        // now, even if splitting the edge deletes it...
        if (e->isNew) {
            if ((h->vertex()->isNew && !h->twin()->vertex()->isNew) ||
                (!h->vertex()->isNew && h->twin()->vertex()->isNew)){
                mesh.flipEdge(e);
            }
        }

        // ...we still have a valid reference to the next edge.
        e = nextEdge;
    }
        
     for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
       if (!v->isNew) v->position = v->newPosition;
     }

      // Copy the updated vertex positions to the subdivided mesh.
    
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
    MutablePriorityQueue<EdgeRecord> q;
    int nedges = mesh.nEdges();
    for (FaceIter f = mesh.facesBegin(); f != mesh.facesEnd(); f++) {
        HalfedgeIter h = f->halfedge();
        VertexIter v = h->vertex();
        Vector4D u = Vector4D(v->position, 1);
        Vector4D d = Vector4D(f->normal(), -dot(v->position, f->normal()));
        f->quadric = outer(u, d);
    }
    for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
        HalfedgeIter h = v->halfedge();
        for (int i = 0; i < v->degree(); i++) {
            v->quadric += h->face()->quadric;
            h = h->twin()->next();
        }
    }
    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
        e->record = EdgeRecord(e);
        q.insert(e->record);
    }
    while (mesh.nEdges() > nedges * 0.25) {
        EdgeRecord r = q.top();
        q.pop();
        EdgeIter e = r.edge;

        HalfedgeIter h1 = e->halfedge(), h2 = h1->twin();
        VertexIter i = h1->vertex();
        VertexIter j = h2->vertex();
        Matrix4x4 _new = i->quadric + j->quadric;
        for (int a = 0; a < i->degree(); a++) {
            q.remove(h1->edge());
            h1 = h1->twin()->next();
        }
        for (int a = 0; a < j->degree(); a++) {
            q.remove(h2->edge());
            h2 = h2->twin()->next();
        }
        Vector3D p = e->centroid();
        VertexIter v = mesh.collapseEdge(e);
        v->position = p;
        HalfedgeIter h = v->halfedge();
        for (int a = 0; a < v->degree(); a++) {
            EdgeIter e = h->edge();
            e->record = EdgeRecord(e);
            q.insert(e->record);
            h = h->twin()->next();
        }
    }
    
  
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
