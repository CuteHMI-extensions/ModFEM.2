
#include "AdinaDatImporter.h"
#include <string.h>

namespace MeshRead
{

AdinaDatImporter::AdinaDatImporter(const std::string & file_name)
: MeshFileImporter(file_name)
{
}

AdinaDatImporter::~AdinaDatImporter(void)
{
	MeshFileImporter::Free();
}

bool	AdinaDatImporter::Init(const std::string & file_name)
{
	if(MeshFileImporter::Init(file_name))
	{
		// First run over file.
		char	token[256];
		do
		{
			_file.getline(token,255,'\n');
		} while (strcmp(token,"C***  NODAL POINT DATA") != 0);

		// we find nodal data, now check it


		// we find element data, now ckeck it
	}
	return _file.is_open();
}

int     AdinaDatImporter::GetCoordinatesDimension() const
{
    // TODO: add Adina Dat file reading
    return 0;
}

int     AdinaDatImporter::GetVerticesCount()
{
    return 0;
}

bool    AdinaDatImporter::GetNextVertex(double coords[])
{
    return false;
}

int     AdinaDatImporter::GetElementCount()
{
    return 0;
}

bool    AdinaDatImporter::GetNextElement(int vertices[], int neighbours[])
{
    return false;
}


bool    AdinaDatImporter::GetBoundaryConditions(double ** bc, int & param)
{
    return false;
}

}// namespace
