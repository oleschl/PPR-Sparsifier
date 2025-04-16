#ifndef UNTITLED2_MIN_DEGREE_PQ_H
#define UNTITLED2_MIN_DEGREE_PQ_H

#include <vector>

/*
https://student.cs.uwaterloo.ca/~cs466/notes/Notes_5.1.1_DynamicMinDegree.pdf
*/

class MinDegreePQ {
private:

    struct Vertex{
        int id;
        int degree;
        Vertex* prev;
        Vertex* next;

        Vertex() : id(-1), degree(0), prev(nullptr), next(nullptr) {}
        Vertex(int id, int degree) : id(id), degree(degree), prev(nullptr), next(nullptr) {}
    };

    int minDegree;
    // TODO find appropriate value
    int degreeThresholdRatio = 1;
    int degreeThreshold;
    std::vector<Vertex*> vertices;
    std::vector<Vertex*> buckets;

    void remove(Vertex* vertex);
    void insert(Vertex* vertex);

public:
    MinDegreePQ(int n, int k);
    void add(int vertexId, int degree);
    void eliminate(int vertexId);
    void update(int vertexId, int newDegree);
    void increase(int vertexId);
    void decrease(int vertexId);
    int get_min_deg();
    int pop();
};

#endif //UNTITLED2_MIN_DEGREE_PQ_H
