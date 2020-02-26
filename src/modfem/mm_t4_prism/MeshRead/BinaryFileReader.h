#ifndef BINARYFILEREADER_H_INCLUDED
#define BINARYFILEREADER_H_INCLUDED

#include "IMeshReader.h"

#include <fstream>

namespace MeshRead
{

  class BinaryFileReader : public           IMeshReader
{
public:
  
  BinaryFileReader();
  BinaryFileReader(const std::string & fileName);
	  
  virtual ~BinaryFileReader();

    bool Init();
    bool Init(const std::string &);
    void Free();

    bool doRead(hHybridMesh * mesh);
	std::string Name() const;
    //interface implementation
    int GetCoordinatesDimension() const;
    int GetVerticesCount();
    bool GetNextVertex(double *);
    void GetElementCount(int*);
    bool GetNextElement(Tind*, Tind*,Tind*,Tind&,Tind&,int8_t&,int8_t&);
    bool GetBoundaryConditions(double **,int&); 
	
	bool needsSetup() const ;
	
    template <typename T>
    friend BinaryFileReader& operator>>(BinaryFileReader& reader,T & value);

    template <typename T>
    friend BinaryFileReader& operator>>(BinaryFileReader& reader,T * value);

private:
    std::string     _file_name;
    std::ifstream   _file;
};

template <typename T>
BinaryFileReader& operator>> (BinaryFileReader& reader, T & value)
{
    reader._file.read(reinterpret_cast<char *>(&value), sizeof(value));
    return reader;
}

template <typename T>
BinaryFileReader& operator>> (BinaryFileReader& reader, T * value)
{
    reader._file.read(reinterpret_cast<char *>(value), sizeof(T));
    return reader;
}


};

#endif // BINARYFILEREADER_H_INCLUDED
