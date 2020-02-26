#include "hHybridMeshWithContacts.h"

void    hHybridMeshWithContacts::collectGroupIds(std::vector<int>& groupsIds)
{
    std::vector<int>    group_count;
    group_count.clear();

    for(hHybridMesh::ElemPool::Iterator<hHybridMesh::ElemPool::allObj> it(& elements_);
        !it.done(); ++it)
    {
        int group_id = it->flags(GROUP_ID);
        if(group_id >= group_count.size()) {
            group_count.resize(group_id+1);
            group_count[group_id]=1;
        }
        else {
            ++group_count[group_id];
        }
    }

    for(int i=0; i < group_count.size(); ++i) {
        if(group_count[i] > 0) {
            groupsIds.push_back(i);
        }
    }
}

bool hHybridMeshWithContacts::normalizationProcessor()
{
    collectGroupIds(this->groups_ids);
    return hHybridMesh::normalizationProcessor();
}

int hHybridMeshWithContacts::n_Groups()
{
    return this->groups_ids.size();
}

int hHybridMeshWithContacts::getGroups(int array[])
{
    int rc=-1;
    if(array != NULL) {
        memcpy(array,this->groups_ids.data(),groups_ids.size()*sizeof(int));
        rc = groups_ids.size();
    }
    return rc;
}
