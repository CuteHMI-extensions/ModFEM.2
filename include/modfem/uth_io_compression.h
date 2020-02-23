#ifndef _UTH_IO_COMPRESSION_H_
#define _UTH_IO_COMPRESSION_H_

/**
 \defgroup UTM_IO_COMPRESSION IO Compression Utilities
 \ingroup UTM


The IO Compression Utilities provides automatic ZIP compression of all dumped mesh and field files.
General assumptions:
    -   compression and decompression is applied for all I/O operations executed with usage of \link UTM utilities \endlink
    -   decompression seeks for given filename with or without ".zip" ending
        (e.g. if in problem definition file a "mesh.dat" filename is given, a "mesh.dat.zip" file is accepted as correct data file and decompressed to "mesh.dat")
    -   if decompressed file is found with matching name, it will be selected for read
    -   decompression temporarily creates files in working directory
    -   dumped files are compressed right after file write is complete
    -   original (not-compressed) files are auto-deleted at program exit
    -   if program fails to exit correctly both compressed and uncompressed file can be found in work directory
    -   compressed files are about 15-30% of the size of not compressed files
    -   compression is recommended for files > 1MB, and strongly recommended for files > 10MB


The IO Compression Utilities are enabled by default.
  @{
 */
#ifdef __cplusplus
extern "C"
{
#endif



/** \brief This functions tries to find compressed file and decompress it, in given directory.
 *
 * \param Work_dir const char* Working directory - in this path given \link Filename \endlink will be searched.
 * \param Filename const char* filename with or without ".zip" ending to decompress.
 * \param Decompressed_filename OUT: filename to proceed with, after decompression.
 * \return int
 *
 */
int utr_io_decompress_file(const char* Work_dir, const char* Filename, char* Decompressed_filename);

/** \brief This functions tires to find given file and compress it, saving compressed file under the same name, but with ".zip" postfix
 *  (e.g. "field0.dmp" will become "field0.dmp.zip").
 *
 * \param Work_dir const char* Working directory - in this path given \link Filename \endlink will be searched.
 * \param Filename const char* file to compress
 * \return int
 *
 */
int utr_io_compress_file(const char* Work_dir, const char* Filename);


#ifdef __cplusplus
}
#endif

 /** @} */ // end of group
#endif //_UTH_IO_COMPRESSION_H_
