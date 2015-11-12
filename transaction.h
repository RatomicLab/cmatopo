#ifndef __CMA_TRANSACTION_H
#define __CMA_TRANSACTION_H

#include <vector>

#include <types.h>

namespace cma {

class Topology;

class TopologyTransaction
{
public:
    TopologyTransaction(Topology& topology) : _topology(topology) {};
    virtual ~TopologyTransaction() {};

    virtual void commit() {};
    virtual void rollback() = 0;

protected:
    Topology& _topology;
};

class EdgeTransaction : public TopologyTransaction
{
public:
    EdgeTransaction(Topology& topology, edge* e);
    ~EdgeTransaction();

    void commit();
    void rollback();

private:
    edge* _edge;
};

class NodeTransaction : public TopologyTransaction
{
public:
    NodeTransaction(Topology& topology, node* n);
    ~NodeTransaction();

    void commit();
    void rollback();

private:
    node* _node;
};

class FaceTransaction : public TopologyTransaction
{
public:
    FaceTransaction(Topology& topology, face* f);
    ~FaceTransaction();

    void commit();
    void rollback();

private:
    face* _face;
};

class AddFaceIndexTransaction : public TopologyTransaction
{
public:
    AddFaceIndexTransaction(
        Topology& topology,
        std::vector<edgeid_set_ptr>* index,
        int faceId,
        int edgeId);
    ~AddFaceIndexTransaction();

    void rollback();

private:
    std::vector<edgeid_set_ptr>* _index;
    int _faceId;
    int _edgeId;
};

class RemoveFaceIndexTransaction : public TopologyTransaction
{
public:
    RemoveFaceIndexTransaction(
        Topology& topology,
        std::vector<edgeid_set_ptr>* index,
        int faceId,
        int edgeId);
    ~RemoveFaceIndexTransaction();

    void rollback();

private:
    std::vector<edgeid_set_ptr>* _index;
    int _faceId;
    int _edgeId;
};

} // namespace cma

#endif // __CMA_TRANSACTION_H
