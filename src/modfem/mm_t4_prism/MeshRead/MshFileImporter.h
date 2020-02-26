#ifndef _MSH_FILE_IMPORTER_H_
#define  _MSH_FILE_IMPORTER_H_

#include <string>

#include "MeshFileImporter.h"

namespace MeshRead
{
  
  class MshFileImporter : public MeshFileImporter
  {
  public:
	MshFileImporter();
	MshFileImporter(const std::string & fn);
	
	virtual ~MshFileImporter();

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

	friend std::istream& operator>>(std::istream & ifs,EBlockLine & b);
	
  };  
}; // namespace MeshRead
  
#endif //  _MSH_FILE_IMPORTER_H_
