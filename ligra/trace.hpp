#ifndef TRACE_HPP
#define TRACE_HPP

#define DO_TRACE

#if defined(DO_TRACE)

#define TRACE_OPEN(filename) trace_open(filename)
#define TRACE_FINISH() trace_finish()
#define TRACE_LOAD_GRAPH_VERTEX_COUNT(n) trace_prop("vertex_count", n)
#define TRACE_LOAD_GRAPH_EDGE_COUNT(m) trace_prop("edge_count", m)

#include <iostream>
#include <fstream>


static std::ofstream trace_out;

static inline void trace_open(const char* filename) {
	trace_out.open(filename, std::ofstream::out);
}

static inline void trace_finish() {
	trace_out.close();
}

template<typename T>
static inline void trace_prop(const char* name, const T& value) {
	trace_out << name << "=" << value << std::endl;
}

#else

#define TRACE_OPEN(filename)
#define TRACE_FINISH()
#define TRACE_LOAD_GRAPH_VERTEX_COUNT(n)
#define TRACE_LOAD_GRAPH_EDGE_COUNT(m)

#endif




#endif // TRACE_HPP
