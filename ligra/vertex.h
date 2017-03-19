#ifndef VERTEX_H
#define VERTEX_H
#include "vertexSubset.h"
#include "trace.hpp"
using namespace std;

namespace decode_uncompressed {
  template <class V, class F>
  inline void decodeInNghBreakEarly(V* v, long i, bool* vertexSubset, F &f, bool* next, bool parallel = 0) {
    uintE d = v->getInDegree();
    if (!parallel || d < 1000) {
      for (uintE j=0; j<d; j++) {
        uintE ngh = v->getInNeighbor(j);
#ifndef WEIGHTED
        if (vertexSubset[ngh] && f.update(ngh,i))
#else
        if (vertexSubset[ngh] && f.update(ngh,i,v->getInWeight(j)))
#endif
          next[i] = 1;
        if(!f.cond(i)) break;
      }
    } else {
      parallel_for(uintE j=0; j<d; j++) {
        uintE ngh = v->getInNeighbor(j);
#ifndef WEIGHTED
        if (vertexSubset[ngh] && f.updateAtomic(ngh,i))
#else
        if (vertexSubset[ngh] && f.updateAtomic(ngh,i,v->getInWeight(j)))
#endif
          next[i] = 1;
      }
    }
  }

  template <class V, class F>
  inline void decodeOutNgh(V* v, long i, bool* vertexSubset, F &f, bool* next) {
    uintE d = v->getOutDegree();
    if(d < 1000) {
      for(uintE j=0; j<d; j++) {
        uintE ngh = v->getOutNeighbor(j);
#ifndef WEIGHTED
        if (f.cond(ngh) && f.updateAtomic(i,ngh))
#else 
        if (f.cond(ngh) && f.updateAtomic(i,ngh,v->getOutWeight(j))) 
#endif
          next[ngh] = 1;
      }
    } else {
      parallel_for(uintE j=0; j<d; j++) {
        uintE ngh = v->getOutNeighbor(j);
#ifndef WEIGHTED
        if (f.cond(ngh) && f.updateAtomic(i,ngh)) 
#else
          if (f.cond(ngh) && f.updateAtomic(i,ngh,v->getOutWeight(j)))
#endif
        next[ngh] = 1;
      }
    }
  }

  template <class V, class F>
  inline void decodeOutNghSparse(V* v, long i, uintT o, F &f, uintE* outEdges) {
    uintE d = v->getOutDegree();
    if(d < 1000) {
      for (uintE j=0; j < d; j++) {
        uintE ngh = v->getOutNeighbor(j);
#ifndef WEIGHTED
        if(f.cond(ngh) && f.updateAtomic(i,ngh)) 
#else
        if(f.cond(ngh) && f.updateAtomic(i,ngh,v->getOutWeight(j)))
#endif
          outEdges[o+j] = ngh;
        else outEdges[o+j] = UINT_E_MAX;
      }
    } else {
      parallel_for (uintE j=0; j < d; j++) {
        uintE ngh = v->getOutNeighbor(j);
#ifndef WEIGHTED
        if(f.cond(ngh) && f.updateAtomic(i,ngh)) 
#else
        if(f.cond(ngh) && f.updateAtomic(i,ngh,v->getOutWeight(j)))
#endif
          outEdges[o+j] = ngh;
        else outEdges[o+j] = UINT_E_MAX;
      }
    }
  }
}

struct symmetricVertex {
#if defined(DO_TRACE)
  uintT id;
#endif
#ifndef WEIGHTED
  uintE* neighbors;
#else 
  intE* neighbors;
#endif
  uintT degree;
  void del() {free(neighbors); }
#if defined(DO_TRACE)
#ifndef WEIGHTED
symmetricVertex(uintT id, uintE* n, uintT d)
#else 
symmetricVertex(uintT id, intE* n, uintT d)
#endif
: id(id), neighbors(n), degree(d) {}
#else
#ifndef WEIGHTED
symmetricVertex(uintE* n, uintT d)
#else
symmetricVertex(intE* n, uintT d)
#endif
: neighbors(n), degree(d) {}
#endif
#ifndef WEIGHTED
  uintE* getInNeighbors () { return neighbors; }
  uintE* getOutNeighbors () { return neighbors; }
  uintE getInNeighbor(uintT j) { TRACE_EDGE_READ(id, neighbors[j], &neighbors[j]); return neighbors[j]; }
  uintE getOutNeighbor(uintT j) { TRACE_EDGE_READ(id, neighbors[j], &neighbors[j]); return neighbors[j]; }
  void setInNeighbors(uintE* _i) { neighbors = _i; }
  void setOutNeighbors(uintE* _i) { neighbors = _i; }
#else
  //weights are stored in the entry after the neighbor ID
  //so size of neighbor list is twice the degree
  intE* getInNeighbors () { return neighbors; }
  intE* getOutNeighbors () { return neighbors; }
  intE getInNeighbor(intT j) { TRACE_EDGE_READ(id, neighbors[2*j], &neighbors[2*j]); return neighbors[2*j]; }
  intE getOutNeighbor(intT j) { TRACE_EDGE_READ(id, neighbors[2*j], &neighbors[2*j]); return neighbors[2*j]; }
  intE getInWeight(intT j) { TRACE_EDGE_READ(id, neighbors[2*j], &neighbors[2*j+1]); return neighbors[2*j+1]; }
  intE getOutWeight(intT j) { TRACE_EDGE_READ(id, neighbors[2*j], &neighbors[2*j+1]); return neighbors[2*j+1]; }
  void setInNeighbors(intE* _i) { neighbors = _i; }
  void setOutNeighbors(intE* _i) { neighbors = _i; }
#endif

  uintT getInDegree() { TRACE_VERTEX_READ(id, &this->degree); return degree; }
  uintT getOutDegree() { TRACE_VERTEX_READ(id, &this->degree); return degree; }
  void setInDegree(uintT _d) { TRACE_VERTEX_WRITE(id, &this->degree); degree = _d; }
  void setOutDegree(uintT _d) { TRACE_VERTEX_WRITE(id, &this->degree); degree = _d; }
  void flipEdges() {}

  template <class F>
  inline void decodeInNghBreakEarly(long i, bool* vertexSubset, F &f, bool* next, bool parallel = 0) {
    decode_uncompressed::decodeInNghBreakEarly<symmetricVertex, F>(this, i, vertexSubset, f, next, parallel);
  }

  template <class F>
  inline void decodeOutNgh(long i, bool* vertexSubset, F &f, bool* next) {
     decode_uncompressed::decodeOutNgh<symmetricVertex, F>(this, i, vertexSubset, f, next);
  }

  template <class F>
  inline void decodeOutNghSparse(long i, uintT o, F &f, uintE* outEdges) {
    decode_uncompressed::decodeOutNghSparse<symmetricVertex, F>(this, i, o, f, outEdges);
  }
};

struct asymmetricVertex {
#if defined(DO_TRACE)
  uintT id;
#endif
#ifndef WEIGHTED
  uintE* inNeighbors, *outNeighbors;
#else
  intE* inNeighbors, *outNeighbors;
#endif
  uintT outDegree;
  uintT inDegree;
  void del() {free(inNeighbors); free(outNeighbors);}
#if defined(DO_TRACE)
#ifndef WEIGHTED
asymmetricVertex(uintT id, uintE* iN, uintE* oN, uintT in, uintT od)
#else
asymmetricVertex(uintT id, intE* iN, intE* oN, uintT in, uintT od)
#endif
: id(id), inNeighbors(iN), outNeighbors(oN), inDegree(in), outDegree(od) {}
#else
#ifndef WEIGHTED
asymmetricVertex(uintE* iN, uintE* oN, uintT in, uintT od)
#else
asymmetricVertex(intE* iN, intE* oN, uintT in, uintT od)
#endif
: inNeighbors(iN), outNeighbors(oN), inDegree(in), outDegree(od) {}
#endif
#ifndef WEIGHTED
  uintE* getInNeighbors () { return inNeighbors; }
  uintE* getOutNeighbors () { return outNeighbors; }
  uintE getInNeighbor(uintT j) {  TRACE_EDGE_READ(id, inNeighbors[j], &inNeighbors[j]); return inNeighbors[j]; }
  uintE getOutNeighbor(uintT j) { TRACE_EDGE_READ(id, outNeighbors[j], &outNeighbors[j]); return outNeighbors[j]; }
  void setInNeighbors(uintE* _i) { inNeighbors = _i; }
  void setOutNeighbors(uintE* _i) { outNeighbors = _i; }
#else 
  intE* getInNeighbors () { return inNeighbors; }
  intE* getOutNeighbors () { return outNeighbors; }
  intE getInNeighbor(uintT j) {  TRACE_EDGE_READ(id, inNeighbors[2*j], &inNeighbors[2*j]);   return inNeighbors[2*j]; }
  intE getOutNeighbor(uintT j) { TRACE_EDGE_READ(id, outNeighbors[2*j], &outNeighbors[2*j]); return outNeighbors[2*j]; }
  intE getInWeight(uintT j) {  TRACE_EDGE_READ(id, inNeighbors[2*j], &inNeighbors[2*j+1]);   return inNeighbors[2*j+1]; }
  intE getOutWeight(uintT j) { TRACE_EDGE_READ(id, outNeighbors[2*j], &outNeighbors[2*j+1]); return outNeighbors[2*j+1]; }
  void setInNeighbors(intE* _i) { inNeighbors = _i; }
  void setOutNeighbors(intE* _i) { outNeighbors = _i; }
#endif

  uintT getInDegree() { TRACE_VERTEX_READ(id, &this->inDegree); return inDegree; }
  uintT getOutDegree() { TRACE_VERTEX_READ(id, &this->outDegree); return outDegree; }
  void setInDegree(uintT _d) {  TRACE_VERTEX_WRITE(id, &this->inDegree); inDegree = _d; }
  void setOutDegree(uintT _d) { TRACE_VERTEX_WRITE(id, &this->outDegree); outDegree = _d; }
  void flipEdges() { TRACE_VERTEX_RW(id, this); swap(inNeighbors,outNeighbors); swap(inDegree,outDegree); }

  template <class F>
  inline void decodeInNghBreakEarly(long i, bool* vertexSubset, F &f, bool* next, bool parallel = 0) {
    decode_uncompressed::decodeInNghBreakEarly<asymmetricVertex, F>(this, i, vertexSubset, f, next, parallel);
  }

  template <class F>
  inline void decodeOutNgh(long i, bool* vertexSubset, F &f, bool* next) {
    decode_uncompressed::decodeOutNgh<asymmetricVertex, F>(this, i, vertexSubset, f, next);
  }

  template <class F>
  inline void decodeOutNghSparse(long i, uintT o, F &f, uintE* outEdges) {
    decode_uncompressed::decodeOutNghSparse<asymmetricVertex, F>(this, i, o, f, outEdges);
  }
};

#endif
