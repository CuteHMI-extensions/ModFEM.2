#ifndef MESHFILEWRITER_H
#define MESHFILEWRITER_H

#include "IMeshWriter.h"

#include <string>
#include <fstream>

namespace MeshWrite
{

class MeshFileWriter : public MeshWrite::IMeshWriter
{
    public:
        MeshFileWriter();
        MeshFileWriter(const std::string & file_name);
        virtual ~MeshFileWriter();

        void    Free();

        bool    Init();
        virtual bool    Init(const std::string & file_name);

		virtual bool	doWrite(const hHybridMesh * mesh);
    protected:
        std::ofstream   _file;
        std::string     _file_name;
    private:
};

};
#endif // MESHFILEWRITER_H
