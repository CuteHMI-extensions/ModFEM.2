
#include "AMVisualExporter.h"

AMVisualExporter::AMVisualExporter(const std::string & file_name)
{
	_file_name = file_name;
}

AMVisualExporter::~AMVisualExporter(void)
{
	Free();
}

bool	AMVisualExporter::Init()
{
	_file.open(_file_name.c_str(), std::ios::app);
	return _file.is_open();
}

void	AMVisualExporter::Free()
{
	_file.close();
}

void	AMVisualExporter::WriteVerticesCount(const int noVert)
{
	_file << noVert << std::endl;
}

void	AMVisualExporter::WriteVertex(double coords[])
{
	// 3D coordinates
	_file << coords[0] << std::endl
		  << coords[1] << std::endl
		  << coords[2] << std::endl;
}

void	AMVisualExporter::WriteElementCount(const int noElems)
{
	_file << noElems << std::endl;
}

void	AMVisualExporter::WriteElement(int vertices[], int neighbours[])
{
	_file << vertices[0] << std::endl
		  << vertices[1] << std::endl
		  << vertices[2] << std::endl
		  << vertices[2] << std::endl;
	_file << neighbours[0] << std::endl
		  << neighbours[1] << std::endl
		  << neighbours[2] << std::endl
		  << neighbours[3] << std::endl;
}
