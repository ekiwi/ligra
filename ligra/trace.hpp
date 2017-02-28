#ifndef TRACE_HPP
#define TRACE_HPP

#define DO_TRACE

#if defined(DO_TRACE)

// essential macros
#define TRACE_OPEN(filename) trace_open(filename)
#define TRACE_FINISH() trace_finish()
#define TRACE_PROP(name, value) trace_prop((name), (value))
#define TRACE_VERTEX_READ(id)  trace_access(DataType::Vertex,   AccessType::Read, (id))
#define TRACE_VERTEX_RW(id)    trace_access(DataType::Vertex,   AccessType::ReadWrite, (id))
#define TRACE_VERTEX_WRITE(id) trace_access(DataType::Vertex,   AccessType::Write, (id))
#define TRACE_EDGE_READ(id)    trace_access(DataType::Edge,     AccessType::Read, (id))
#define TRACE_EDGE_WRITE(id)   trace_access(DataType::Edge,     AccessType::Write, (id))
#define TRACE_EDGE_RW(id)      trace_access(DataType::Edge,     AccessType::ReadWrite, (id))
#define TRACE_PROP_READ(id)    trace_access(DataType::Property, AccessType::Read, (id))
#define TRACE_PROP_WRITE(id)   trace_access(DataType::Property, AccessType::Write, (id))
#define TRACE_PROP_RW(id)      trace_access(DataType::Property, AccessType::ReadWrite, (id))

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

enum class DataType   : uint8_t { Edge, Vertex, Property };
static inline const char* to_string(DataType tt) {
	switch(tt) {
	case DataType::Edge: return "E";
	case DataType::Vertex: return "V";
	case DataType::Property: return "P";
	default: return "?";
}}

enum class AccessType : uint8_t { Read, ReadWrite, Write };
static inline const char* to_string(AccessType tt) {
	switch(tt) {
	case AccessType::Read: return "r";
	case AccessType::ReadWrite: return "rw";
	case AccessType::Write: return "w";
	default: return "?";
}}
using IdType = uint64_t;
// TODO: trace address range
static inline void trace_access(DataType data, AccessType access, IdType id) {
	trace_out << to_string(data) << "," << to_string(access) << "," << id << std::endl;
}

#else

#define TRACE_OPEN(filename)
#define TRACE_FINISH()
#define TRACE_PROP(name, value)
#define TRACE_VERTEX_READ(id)
#define TRACE_VERTEX_WRITE(id)
#define TRACE_EDGE_READ(id)
#define TRACE_EDGE_WRITE(id)

#endif


// short hands / constants
#define TRACE_LOAD_GRAPH_VERTEX_COUNT(n) TRACE_PROP("vertex_count", n)
#define TRACE_LOAD_GRAPH_EDGE_COUNT(m) TRACE_PROP("edge_count", m)


#endif // TRACE_HPP