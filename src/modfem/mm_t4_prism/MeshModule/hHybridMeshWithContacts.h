#ifndef H_HYBRID_MESH_WITH_CONTACTS_H
#define H_HYBRID_MESH_WITH_CONTACTS_H

#include "hHybridMesh.h"

class hHybridMeshWithContacts : public hHybridMesh
{
public:
    int n_Groups();
    int getGroups(int array[]);

    bool normalizationProcessor();
protected:
    void    collectGroupIds(std::vector<int>& groupsIds);
    std::vector<int>    groups_ids;
};

#endif //#ifndef H_HYBRID_MESH_WITH_CONTACTS_H
