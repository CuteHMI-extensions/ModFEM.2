#ifndef _IN_FILE_IMPORTER_H_
#define  _IN_FILE_IMPORTER_H_

#include <string>

#include "MeshFileImporter.h"

namespace MeshRead
{
  
  class InFileImporter : public MeshFileImporter
  {
  public:
	InFileImporter();
	InFileImporter(const std::string & fn);
	
	virtual ~InFileImporter();

	void Free();
	
	bool Init();
	bool Init(const std::string& fn);
	bool doRead(hHybridMesh * m);
	
  protected:
	static const std::string nodesBegin;
	static const std::string elemsBegin;
	static const std::string nodesEnd;
	static const std::string cmBlockBegin;
	static const char eol = '\n';

	struct EBlockLine
	{
	  static const int maxNodes=8;
	  
	  
	EBlockLine()
	  : nNodes(0),good(false)
	  {};

	  
	  int nNodes,N,
		field[9];
	uTind	nodes[maxNodes];
	  bool good;
	};

	friend std::istream& operator>>(std::istream & ifs,InFileImporter::EBlockLine & b);
	
  };  
}; // namespace MeshRead
  
#endif //  _IN_FILE_IMPORTER_H_
