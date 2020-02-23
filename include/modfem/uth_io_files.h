#ifndef UTH_IO_FILES_H_
#define UTH_IO_FILES_H_
#ifdef __cplusplus
extern "C"
{
#endif

/**
 \defgroup UTM_IO Input-Ouput Utilities
 \ingroup UTM
 @{
 */

/**---------------------------------------------------------
utr_io_initialize_mesh - to initialize mesh of a specified type
---------------------------------------------------------*/
int utr_io_initialize_mesh( ///> returns mesh_id */
  FILE *Interactive_output, ///> file or stdout to write messages */
  const char* Work_dir, ///> path to working directory
  char Mesh_type, ///> letter symbol denoting mesh_type (j, p, t or h) */
  const char* Mesh_file ///> mesh file name - conforming to naming convention */
  );

/**--------------------------------------------------------
  utr_io_read_mesh - to read mesh with given filename(s)
  Reads all mesh files from current working directory
  that matches with passed regular expression Mesh_files_regex.
  All files HAVE TO be the same type (ex. nas,dat,in etc.).
  NOTE: single filename ex. "mesh.dat" is also a valid regular expression.
//--------------------------------------------------------*/
int utr_io_read_mesh(
   /// returns: > 0 - number of files correctly read; <=0 - error
  const int Mesh_id,  ///  IN: id of mesh to read to
  const char * Working_dir,  /// IN: directory with mesh file(s)
  const char * Mesh_files_regex,  /// IN: regular expression pattern
  ///  NOTE: using regex 'POSIX grep' standard
  const char Mesh_type  /// IN: type of mesh
              );




int utr_io_write_img_to_pbm( /// returns 0 if all ok.
        const char * Work_dir,   /// Directory to write in.
        const char * Filename,  ///  without extension
        const char * Comment, ///  written into img file header, can be NULL if no comment.
        const int Width,   /// >0
        const int Height,  /// >0
        const int Max_color_component_value,
        const int Magic_number,  /// {4,5,6}
        const unsigned char* Img_data,  /// pointer to the array with image data
        FILE* Interactive_output
        );






int utr_io_write_mesh_info_to_PAM(int Mesh_id
        , const char *Work_dir, const char *Filename, const char *Comment, FILE *Interactive_output);
/**
 utr_io_write_img_to_pam - to write data into PAM image file.
 http://en.wikipedia.org/wiki/Netpbm
 Allowed combinations of parameters:
 -----------------------------------------------------------
 TUPLTYPE           |MAXVAL |DEPTH  |comment
 -----------------------------------------------------------
 BLACKANDWHITE        1       1       special case of GRAYSCALE
 GRAYSCALE            2…65535	1       2 bytes per pixel for MAXVAL > 255
 RGB                  1…65535	3       6 bytes per pixel for MAXVAL > 255
 BLACKANDWHITE_ALPHA	1       2       2 bytes per pixel
 GRAYSCALE_ALPHA      2…65535	2       4 bytes per pixel for MAXVAL > 255
 RGB_ALPHA            1…65535	4       8 bytes per pixel for MAXVAL > 255s of parameters:
 -----------------------------------------------------------*/
int utr_io_write_img_to_pam(  /// returns 0 if all ok.
        const char * Work_dir,   /// Directory to write in.
        const char * Filename,   /// without extension
        const int Width,   /// >0
        const int Height,  /// >0
        const int Depth,   /// <1:4>
        const int Maxval,  /// <1:65535>,
        const char* TUPLTYPE,
        const char* Img_data,  /// pointer to the array with image data
        FILE* Interactive_output
        );
/**
 utr_io_write_img_to_pnm - to write data into PNM image file.
 The Portable Bit/Grey/PixMap formats PBM, PGM, PPM.
 They are collectively referred to as PNM (Portable aNy Map).
 http://en.wikipedia.org/wiki/Netpbm_format
 Allowed combinations:
 -----------------------------------------------------------
 Type                 Magic number	Extension	Colors
 -----------------------------------------------------------
 Portable BitMap		P4	binary      .pbm        0–1 (black & white)
 Portable GrayMap		P5	binary      .pgm        0–255 (gray scale)
 Portable PixMap		P6	binary      .ppm        0–255 (RGB)
 -----------------------------------------------------------*/
int utr_io_write_img_to_pnm( /// returns 0 if all ok.
        const char * Work_dir,   /// Directory to write in.
        const char * Filename,   /// without extension
        const char * Comment,  /// written into img file header, can be NULL if no comment.
        const int Width,   /// >0
        const int Height,  /// >0
        const int Max_color_component_value,  /// see Colors above
        const int Magic_number,  /// {4,5,6}
        const char* Img_data,  /// pointer to the array with image data
        FILE* Interactive_output
        );

/** @}*/

#ifdef __cplusplus
}
#endif

/** @} */ // end of group


#endif //UTH_IO_FILES_H_
