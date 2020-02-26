#include "MeshFileWriter.h"

using MeshWrite::MeshFileWriter;

MeshFileWriter::MeshFileWriter(const std::string & file_name) : _file_name(file_name)
{
}

MeshFileWriter::MeshFileWriter()
{
    //ctor
}

MeshFileWriter::~MeshFileWriter()
{
    Free();
}

void    MeshFileWriter::Free()
{
    if(_file.is_open())
    {
        _file.close();
    }
}

bool    MeshFileWriter::Init()
{
    return  Init(_file_name);
}

bool    MeshFileWriter::Init(const std::string & file_name)
{
    if(!file_name.empty())
    {
        _file_name = file_name;
        _file.open(_file_name.c_str(), std::ios::out|std::ios::trunc);
    }
    return _file.is_open();
}

bool	MeshFileWriter::doWrite(const hHybridMesh * mesh)
{
	_file.close();
	return false;
}