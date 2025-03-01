#include "utility/min_degree_pq.h"

// data structure for dynamically finding the min degree vertex
// graph of n vertices and k non-terminals
// not nice but assume that isolated vertices have been handled/removed
MinDegreePQ::MinDegreePQ(int n, int k)
        : minDegree(1), degreeThreshold(n / degreeThresholdRatio), vertices(k), buckets(degreeThreshold) {}

void MinDegreePQ::add(int vertexId, int degree){
    vertices[vertexId] = new Vertex(vertexId, std::min(degree, degreeThreshold-1));
    insert(vertices[vertexId]);
}

void MinDegreePQ::eliminate(int vertexId){
    remove(vertices[vertexId]);
}

void MinDegreePQ::update(int vertexId, int newDegree) {
    Vertex* vertex = vertices[vertexId];
    int degree = std::min(newDegree, degreeThreshold-1);
    if(degree == vertex->degree) return;

    remove(vertex);
    vertex->degree = degree;
    insert(vertex);

    // Update minDegree
    if (newDegree < minDegree) {
        minDegree = newDegree;
    }
}

// only guaranteed to return the correct value if the last operation was pop
int MinDegreePQ::get_min_deg() {
    return minDegree;
}

int MinDegreePQ::pop(){
    while (!buckets[minDegree]) {
        ++minDegree;
    }
    Vertex* vertex = buckets[minDegree];

    remove(vertex);
    return vertex->id;
}

void MinDegreePQ::remove(Vertex* vertex) {
    if (vertex->prev) {
        vertex->prev->next = vertex->next;
    } else {
        buckets[vertex->degree] = vertex->next;
    }
    if (vertex->next) {
        vertex->next->prev = vertex->prev;
    }
    //vertex->prev = vertex->next = nullptr;
}

void MinDegreePQ::insert(Vertex* vertex) {
    auto degree = vertex->degree;
    vertex->prev = nullptr;
    vertex->next = buckets[degree];
    if (buckets[degree]) {
        buckets[degree]->prev = vertex;
    }
    buckets[degree] = vertex;
}
