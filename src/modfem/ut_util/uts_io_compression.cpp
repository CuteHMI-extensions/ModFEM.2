#include <set>
#include <string>
#include <boost/filesystem.hpp>

#include <modfem/uth_intf.h>
#include <modfem/uth_system.h>
#include <modfem/uth_io_compression.h>
#define MINIZ_NO_ZLIB_COMPATIBLE_NAMES //< to avoid naming clashes (miniz supports subset of zlib).
#include "miniz.h"
#include <modfem/uth_log.h>

std::set<std::string> utv_io_tmp_files;

void utr_io_remove_file(const std::string & filename)
{
#ifndef PARALLEL
	{

		boost::filesystem::path file(filename);
		if (boost::filesystem::exists(file)) {
			boost::filesystem::remove(file);
		}
	}
#endif
}

void utr_io_remove_tmp_files()
{
	if (! utv_io_tmp_files.empty()) {
		std::set<std::string>::iterator it = utv_io_tmp_files.begin();

		for (; it != utv_io_tmp_files.end(); ++it) {
			if (!it->empty()) {
				// at this point utv_log_out is already freed.
				// mf_log_info("Removing temporary file %s",it->c_str());
				utr_io_remove_file(*it);
			}
		}

		utv_io_tmp_files.clear();
	}
}

#ifdef __cplusplus
extern "C"
{
#endif

/**
/** \defgroup UTM_IO_COMPRESSION IO Compression Utilities
/** \ingroup UTM
/** @{
/** */

///
int utr_io_decompress_file(const char * Work_dir, const char * Filename, char * Decompressed_filename)
{

	double time_begin = time_clock();


	mf_check_mem(Filename);
	mf_check_mem(Decompressed_filename);

	if (utv_io_tmp_files.empty()) {
		atexit(utr_io_remove_tmp_files);
	}

	std::string filepath, zippath, filename(Filename);

	if (Work_dir != NULL) {
		filepath.append(Work_dir);
	}


	// Remove ".zip" ending if given file is a zip.
	if (filename.substr(filename.size() - 4, 4) == ".zip") {
		filename = filename.substr(0, filename.size() - 4);
	}

	strcpy(Decompressed_filename, filename.c_str());

	// for "." work dir
	if (filepath.size() == 1) {
		filepath.clear();
	}
	else if (!filepath.empty()) {
		filepath.append("/");
	}

	filepath.append(filename);
	zippath = filepath + ".zip";

	filename = boost::filesystem::path(filepath).filename().string();

	mz_zip_archive  zip_archive;
	memset(&zip_archive, 0, sizeof(zip_archive));

	mz_bool success = mz_zip_reader_init_file(& zip_archive, zippath.c_str(), MZ_ZIP_FLAG_DO_NOT_SORT_CENTRAL_DIRECTORY);

	mf_log_info("Decompressing file %s from %s", filename.c_str(), zippath.c_str());

	if ( success == MZ_TRUE ) {

		success =  mz_zip_reader_extract_file_to_file(& zip_archive, filename.c_str(), filepath.c_str(), 0 );

		if ( success ==  MZ_TRUE ) {
			utv_io_tmp_files.insert(filepath);
			mf_log_info("Temporary file %s will be deleted at exit.", filepath.c_str());
		}
		else {
			if (boost::filesystem::exists(filepath)) {
				mf_log_warn("Unable to extract file %s from archive %s, but this file already exist.", filename.c_str(), zippath.c_str());
			}
			else {
				mf_log_err("Unable to extract file %s from archive %s", filename.c_str(), zippath.c_str());
			}
		}

		mz_zip_reader_end(& zip_archive);

	}
	else if (boost::filesystem::exists(filepath)) {
		mf_log_warn("Unable to init decompression of archive  %s, but file %s already exist.", zippath.c_str(), filepath.c_str() );
	}
	else {
		mf_log_err("Unable to init decompression of archive %s", zippath.c_str());
	}

	utv_time.decompress_files += time_clock() - time_begin;


	return success;
}

int utr_io_compress_file(const char * Work_dir, const char * Filename)
{

	double time_begin = time_clock();

	mf_check_mem(Filename);

	if (utv_io_tmp_files.empty()) {
		atexit(utr_io_remove_tmp_files);
	}

	std::string filepath, zippath, filename(Filename);

	if (Work_dir != NULL) {
		filepath.append(Work_dir);
	}

	// Remove ".zip" ending if given file is a zip.
	if (filename.substr(filename.size() - 4, 4) == ".zip") {
		filename = filename.substr(0, filename.size() - 4);
	}

	// for "." work dir
	if (filepath.size() == 1) {
		filepath.clear();
	}
	else if (!filepath.empty()) {
		filepath.append("/");
	}

	filepath.append(filename);
	zippath = filepath + ".zip";

	mz_bool success = MZ_FALSE;


	filepath = boost::filesystem::path(filepath).generic_string();  //make_preferred().string();
	zippath = boost::filesystem::path(zippath).generic_string();
	filename = boost::filesystem::path(filename).filename().generic_string();

	if (boost::filesystem::exists(filepath)) {

		mf_log_info("Compressing file %s into archive %s", filepath.c_str(), zippath.c_str());

		mz_zip_archive  zip_archive;
		memset(&zip_archive, 0, sizeof(zip_archive));

		success = mz_zip_writer_init_file(& zip_archive, zippath.c_str(), 32);

		if (success == MZ_FALSE ) {
			mf_log_err("Unable to init compression or file %s", zippath.c_str());
		}
		else {
			success = mz_zip_writer_add_file(& zip_archive, filename.c_str(), filepath.c_str(), NULL, 0, MZ_BEST_COMPRESSION);

			if (success == MZ_FALSE ) {
				mf_log_err("Unable to compress file %s into archive %s", filepath.c_str(), zippath.c_str());
			}
			else {
				// Here are 2 options:
				// 1. Delete original files at exit - safer, but both files stay until end.
//                utv_io_tmp_files.insert(filepath.c_str());
//                mf_log_info("Temporary file %s will be deleted at exit.", filepath.c_str());

				// 2. Delete original file at once - cleaner, faster and more elegant, but if zip is corrupted
				// or something, there will be no backup.
				utr_io_remove_file(filepath.c_str());
				mf_log_info("Temporary file %s deleted.", filepath.c_str());
			}
		}

		mz_zip_writer_finalize_archive(& zip_archive);
		mz_zip_writer_end(& zip_archive );
	}
	else {
		mf_log_err("Unable to locate file %s for compression into archive %s", filepath.c_str(), zippath.c_str() );
	}

	utv_time.compress_files += time_clock() - time_begin;

	return success;
}

/**@}*/

#ifdef __cplusplus
}
#endif
