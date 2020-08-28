import qbs

import cutehmi

Project {
	name: "ModFEM.2"

	cutehmi.CppExtension {
		name: "ModFEM.Interfaces.2"

		vendor: "ModFEM"

		domain: "modfem.agh.edu.pl"

		friendlyName: "ModFEM interfaces"

		description: "ModFEM interfaces."

		files: [
			"include/modfem/aph_intf.h",
			"include/modfem/apph_intf.h",
			"include/modfem/ddh_intf.h",
			"include/modfem/gitversion.h",
			"include/modfem/info.h",
			"include/modfem/lin_alg_intf.h",
			"include/modfem/mf_version.h",
			"include/modfem/mmh_intf.h",
			"include/modfem/mmph_intf.h",
			"include/modfem/mod_fem_viewer.h",
			"include/modfem/pch_intf.h",
			"include/modfem/pdh_control_intf.h",
			"include/modfem/pdh_intf.h",
			"include/modfem/sih_intf.h",
			"include/modfem/svnversion.h",
			"include/modfem/tmh_intf.h",
			"include/modfem/uth_bc.h",
			"include/modfem/uth_err.h",
			"include/modfem/uth_hash.h",
			"include/modfem/uth_intf.h",
			"include/modfem/uth_io_compression.h",
			"include/modfem/uth_io_files.h",
			"include/modfem/uth_io_results.h",
			"include/modfem/uth_log.h",
			"include/modfem/uth_mat.h",
			"include/modfem/uth_mesh.h",
			"include/modfem/uth_system.h",
		]
	}

	cutehmi.CppExtension {
		name: "ModFEM.2"

		vendor: "ModFEM"

		domain: "modfem.agh.edu.pl"

		friendlyName: "ModFEM modules"

		description: "ModFEM modules."

		Export {

			Depends { name: "modfem.config" }

			Depends { name: "modfem.libs.libconfig" }

			Depends { name: "modfem.libs.libm" }

			Depends { name: "modfem.libs.libpthread" }

			Depends { name: "modfem.libs.mkl" }

			Depends { name: "modfem.libs.zlib" }

			Depends { name: "modfem.libs.boost.system" }

			Depends { name: "modfem.libs.boost.filesystem" }

			Depends { name: "modfem.libs.boost.regex" }

			Depends { name: "modfem.libs.iomp5"; condition: modfem.config.acceleration === "openmp" && modfem.config.openMPLibrary === "iomp5" }
			Depends { name: "modfem.libs.gomp"; condition: modfem.config.acceleration === "openmp" && modfem.config.openMPLibrary === "gomp" }

			Depends { name: "ModFEM.Interfaces.2" }

			cpp.defines: {
				var result = []

				if (modfem.config.renumbering)
					result.push("RENUMBERING")

				if (modfem.config.internalRenumbering)
					result.push("INTERNAL_RENUMBERING")

				if (modfem.config.mixedApproximationMatrixStorage === "dofByDof")
					result.push("MIXED_DOF_BY_DOF")
				else if (modfem.config.mixedApproximationMatrixStorage === "fieldByField")
					result.push("MIXED_FIELD_BY_FIELD")

				if (modfem.config.acceleration !== "none")
					result.push("MULTITHREADED")
				else if (modfem.config.acceleration !== "openmp") {
					if (modfem.config.mkbDirectSolver === "viennacl")
						result.push("VIENNACL_WITH_OPENMP")
				}

				if (!modfem.config.mpi)
					result.push("NPARALLEL")

				return result
			}

			cpp.cFlags: {
				var result =[]

				if ((modfem.config.acceleration === "openmp") && (modfem.config.openMPLibrary != "iomp5"))
					result.push("-fopenmp")

				return result
			}

			cpp.cppFlags: {
				var result =[]

				if ((modfem.config.acceleration === "openmp") && (modfem.config.openMPLibrary != "iomp5"))
					result.push("-fopenmp")

				return result
			}

			Group {
				name: "ut_system"

				files: [
					"include/modfem/ut_system/utd_threads.h",
					"src/modfem/ut_system/utd_threads.c",
				]

				Group {
					name: "ut_system - Windows"

					condition: qbs.targetOS.contains("windows")

					files: [
						"src/modfem/ut_system/WIN/uts_time.cpp",
					]
				}

				Group {
					name: "ut_system - Linux"

					condition: qbs.targetOS.contains("linux")

					files: [
						"src/modfem/ut_system/UNIX/uts_time.c",
					]
				}
			}

			Group {
				name: "ut_util"

				files: [
					"src/modfem/ut_util/base64.cpp",
					"src/modfem/ut_util/base64.h",
					"src/modfem/ut_util/uts_adapt.c",
					"src/modfem/ut_util/uts_coloring.cpp",
					"src/modfem/ut_util/uts_util.c",
					"src/modfem/ut_util/uts_ls_intf.c",
					"src/modfem/ut_util/uts_write_paraview.cpp",
				]

				Group {
					name: "ut_util - OpenCL dependent files"

					condition: modfem.config.acceleration === "opencl"

					files: [
						"src/modfem/ut_util/uts_accel_intf.c"
					]
				}

				Group {
					name: "ut_bc"

					files: [
						"src/modfem/ut_util/uts_bc.cpp"
					]
				}

				Group {
					name: "ut_hash"

					files: [
						"src/modfem/ut_util/uts_lookup3.c",
					]
				}

				Group {
					name: "ut_io"

					files: [
						"src/modfem/ut_util/miniz.c",
						"src/modfem/ut_util/miniz.h",
						"src/modfem/ut_util/uts_io_compression.cpp",
						"src/modfem/ut_util/uts_io_intf.cpp",
					]
				}

				Group {
					name: "ut_io_result"

					files: [
						"src/modfem/ut_util/uts_io_results.cpp",
					]
				}
				Group {
					name: "ut_log"

					files: [
						"src/modfem/ut_util/uts_log.cpp",
					]
				}

				Group {
					name: "ut_mat"

					files: [
						"src/modfem/ut_util/fastmathparser/exprtk.hpp",
						"src/modfem/ut_util/uts_mat.cpp",
					]
				}

				Group {
					name: "ut_mesh"

					files: [
						"src/modfem/ut_util/uts_mesh.cpp"
					]
				}
			}

			Group {
				name: "mm_prism"

				condition: modfem.config.mm === "prism"

				files: [
					"include/modfem/mm_prism/mmh_prism.h",
					"src/modfem/mm_prism/mms_prism_datstr.c",
					"src/modfem/mm_prism/mms_prism_input_gradmesh.c",
					"src/modfem/mm_prism/mms_prism_intf.c",
					"src/modfem/mm_prism/mms_prism_io_dump.c",
					"src/modfem/mm_prism/mms_prism_ref.c",
					"src/modfem/mm_prism/mms_prism_util.c",
					"src/modfem/mm_prism/mms_util.c",
				]
			}

			Group {
				name: "mm_t4_prism"

				condition: modfem.config.mm === "mm_t4"

				files: [
					"src/modfem/mm_t4_prism/Common.h",
					"src/modfem/mm_t4_prism/Field.hpp",
					"src/modfem/mm_t4_prism/mms_prism_t4_intf.cpp",
					"src/modfem/mm_t4_prism/ProfTimer.hpp",
//					"src/modfem/mm_t4_prism/HybridMeshUnitTest.cpp",
					"src/modfem/mm_t4_prism/MeshModule/ArrayPool.hpp",
					"src/modfem/mm_t4_prism/MeshModule/Edge.cpp",
					"src/modfem/mm_t4_prism/MeshModule/Edge.h",
					"src/modfem/mm_t4_prism/MeshModule/ElemPrism.cpp",
					"src/modfem/mm_t4_prism/MeshModule/ElemPrism.h",
					"src/modfem/mm_t4_prism/MeshModule/ElemT4.cpp",
					"src/modfem/mm_t4_prism/MeshModule/ElemT4.h",
					"src/modfem/mm_t4_prism/MeshModule/ElemT4tables.hpp",
					"src/modfem/mm_t4_prism/MeshModule/EntityAttributes.hpp",
					"src/modfem/mm_t4_prism/MeshModule/Face3.cpp",
					"src/modfem/mm_t4_prism/MeshModule/Face3.h",
					"src/modfem/mm_t4_prism/MeshModule/Face4.cpp",
					"src/modfem/mm_t4_prism/MeshModule/Face4.h",
					"src/modfem/mm_t4_prism/MeshModule/FixedSizeAllocator.hpp",
					"src/modfem/mm_t4_prism/MeshModule/MeshEntity.hpp",
					"src/modfem/mm_t4_prism/MeshModule/StaticPool.hpp",
					"src/modfem/mm_t4_prism/MeshModule/SuperFastHash.cpp",
					"src/modfem/mm_t4_prism/MeshModule/SuperFastHash.h",
					"src/modfem/mm_t4_prism/MeshModule/Vertex.cpp",
					"src/modfem/mm_t4_prism/MeshModule/Vertex.h",
					"src/modfem/mm_t4_prism/MeshModule/VtsSqId.hpp",
					"src/modfem/mm_t4_prism/MeshModule/hHybridMesh.cpp",
					"src/modfem/mm_t4_prism/MeshModule/hHybridMesh.h",
					"src/modfem/mm_t4_prism/MeshModule/hHybridMeshWithContacts.cpp",
					"src/modfem/mm_t4_prism/MeshModule/hHybridMeshWithContacts.h",
					"src/modfem/mm_t4_prism/MeshModule/hObj.cpp",
					"src/modfem/mm_t4_prism/MeshModule/hObj.h",
					"src/modfem/mm_t4_prism/MeshModule/hParent.cpp",
					"src/modfem/mm_t4_prism/MeshModule/hParent.hpp",
					"src/modfem/mm_t4_prism/MeshModule/mesh_configuration.cpp",
					"src/modfem/mm_t4_prism/MeshModule/mesh_configuration.h",
					"src/modfem/mm_t4_prism/MeshModule/mmh_vec3.h",
					"src/modfem/mm_t4_prism/MeshModule/mmr_vec3.cpp",
					"src/modfem/mm_t4_prism/MeshRead/AdinaDatImporter.cpp",
					"src/modfem/mm_t4_prism/MeshRead/AdinaDatImporter.h",
					"src/modfem/mm_t4_prism/MeshRead/BinaryFileReader.cpp",
					"src/modfem/mm_t4_prism/MeshRead/BinaryFileReader.h",
					"src/modfem/mm_t4_prism/MeshRead/DmpFileImporter.cpp",
					"src/modfem/mm_t4_prism/MeshRead/DmpFileImporter.h",
//					"src/modfem/mm_t4_prism/MeshRead/GrdMeshBuilder3D.cpp",
//					"src/modfem/mm_t4_prism/MeshRead/GrdMeshBuilder3D.h",
//					"src/modfem/mm_t4_prism/MeshRead/IMeshBuilder.h",
					"src/modfem/mm_t4_prism/MeshRead/IMeshReader.cpp",
					"src/modfem/mm_t4_prism/MeshRead/IMeshReader.h",
					"src/modfem/mm_t4_prism/MeshRead/InFileImporter.cpp",
					"src/modfem/mm_t4_prism/MeshRead/InFileImporter.h",
					"src/modfem/mm_t4_prism/MeshRead/KazFileImporter.cpp",
					"src/modfem/mm_t4_prism/MeshRead/KazFileImporter.h",
					"src/modfem/mm_t4_prism/MeshRead/MeshFileImporter.cpp",
					"src/modfem/mm_t4_prism/MeshRead/MeshFileImporter.h",
//					"src/modfem/mm_t4_prism/MeshRead/MshFileImporter.cpp",
//					"src/modfem/mm_t4_prism/MeshRead/MshFileImporter.h",
					"src/modfem/mm_t4_prism/MeshRead/NasFileImporter.cpp",
					"src/modfem/mm_t4_prism/MeshRead/NasFileImporter.h",
					"src/modfem/mm_t4_prism/MeshRead/SimplePointReader.cpp",
					"src/modfem/mm_t4_prism/MeshRead/SimplePointReader.h",
					"src/modfem/mm_t4_prism/MeshWrite/AMVisualExporter.cpp",
					"src/modfem/mm_t4_prism/MeshWrite/AMVisualExporter.h",
					"src/modfem/mm_t4_prism/MeshWrite/BinaryFileWriter.cpp",
					"src/modfem/mm_t4_prism/MeshWrite/BinaryFileWriter.h",
					"src/modfem/mm_t4_prism/MeshWrite/DmpFileExporter.cpp",
					"src/modfem/mm_t4_prism/MeshWrite/DmpFileExporter.h",
//					"src/modfem/mm_t4_prism/MeshWrite/HDF5Exporter.cpp",
//					"src/modfem/mm_t4_prism/MeshWrite/HDF5Exporter.h",
					"src/modfem/mm_t4_prism/MeshWrite/IMeshWriter.h",
					"src/modfem/mm_t4_prism/MeshWrite/KazFileExporter.cpp",
					"src/modfem/mm_t4_prism/MeshWrite/KazFileExporter.h",
					"src/modfem/mm_t4_prism/MeshWrite/MeshFileWriter.cpp",
					"src/modfem/mm_t4_prism/MeshWrite/MeshFileWriter.h",
					"src/modfem/mm_t4_prism/MeshWrite/MileninExporter.cpp",
					"src/modfem/mm_t4_prism/MeshWrite/MileninExporter.h",
					"src/modfem/mm_t4_prism/GeometryModule/GeometryModule.hpp",
				]
			}

			Group {
				name: "mm_remesh"

				files: [ "include/modfem/mm_remesh/delanouy3d.h",
					"include/modfem/mm_remesh/figury.h",
					"include/modfem/mm_remesh/losowanie.h",
					"include/modfem/mm_remesh/plansza.h",
					"include/modfem/mm_remesh/rozrost.h",
					"include/modfem/mm_remesh/texty.h",
					"include/modfem/mm_remesh/ticooMesh3D.h",
					"include/modfem/mm_remesh/wektor.h",
					"include/modfem/mm_remesh/ziarno.h",
					"src/modfem/mm_remesh/delanouy3d.cpp",
					"src/modfem/mm_remesh/figury.cpp",
					"src/modfem/mm_remesh/losowanie.cpp",
					"src/modfem/mm_remesh/plansza.cpp",
					"src/modfem/mm_remesh/rozrost.cpp",
					"src/modfem/mm_remesh/teksty.cpp",
					"src/modfem/mm_remesh/ticooMesh3D.cpp",
					"src/modfem/mm_remesh/wektor.cpp",
					"src/modfem/mm_remesh/ziarno.cpp",
				]

				Group {
					name: "mm_remesh - prism_t4"

					condition: modfem.config.mm === "remesh"

					files: [
						"src/modfem/mm_remesh/mms_prism_t4_remesh_intf.cpp",
					]
				}
			}

			Group {
				name: "ap_std_lin"

				condition: modfem.config.ap === "std_lin"

				files: [
					"include/modfem/ap_std_lin/aph_std_lin.h",
					"src/modfem/ap_std_lin/aps_gauss_util_std.c",
					"src/modfem/ap_std_lin/aps_std_lin_intf.c",
					"src/modfem/ap_std_lin/aps_std_util.c",
				]
			}

			Group {
				name: "tm_opencl"

				condition: modfem.config.acceleration === "opencl"

				files: [
					"include/modfem/tm_opencl/tmh_ocl.h",
					"include/modfem/tm_opencl/tmh_ocl_num_int.h",
					"src/modfem/tm_opencl/tms_ocl_intf.c",
					"src/modfem/tm_opencl/tms_ocl_num_int.c",
				]
			}

			Group {
				name: "si_mkb"

				files: [
					"include/modfem/si_mkb/sih_mkb.h",
					"include/modfem/si_mkb/sih_mkb_fem_intf.h",
					"src/modfem/si_mkb/sis_mkb_fem_intf.c",
					"src/modfem/si_mkb/sis_mkb_intf.c",
				]
			}

			Group {
				name: "si_lapack"

				files: [
					"include/modfem/si_lapack/sih_lapack.h",
					"src/modfem/si_lapack/sis_lapack_intf.c",
				]
			}

			Group {
				name: "ls_mkb"

				files: [
					"include/modfem/ls_mkb/lah_intf.h",
					"include/modfem/ls_mkb/lsh_mkb_intf.h",
					"src/modfem/ls_mkb/las_intf.c",
					"src/modfem/ls_mkb/lss_mkb_fem_intf.c",
					"src/modfem/ls_mkb/lss_mkb_intf.c",
				]

				Group {
					name: "ls_mkb - lad_amg_ext"

					files: [
						"src/modfem/ls_mkb/amg_ext/amg_ext.c",
						"src/modfem/ls_mkb/amg_ext/amg_ext.h",
						"src/modfem/ls_mkb/amg_ext/lah_petsc_interface.h",
					]
				}

				Group {
					name: "ls_mkb - lad_bcrs"

					files: [
						"src/modfem/ls_mkb/lad_bcrs/lah_bcrs.h",
						"src/modfem/ls_mkb/lad_bcrs/las_bcrs_intf.c",
					]
				}

				Group {
					name: "ls_mkb - lad_crs"

					files: [
						"src/modfem/ls_mkb/lad_crs/lah_crs.h",
						"src/modfem/ls_mkb/lad_crs/las_crs_intf.c",
					]
				}

				Group {
					name: "ls_mkb - lad_crs_generic"

					files: [
						"src/modfem/ls_mkb/lad_crs_generic/lah_crs_generic.h",
						"src/modfem/ls_mkb/lad_crs_generic/las_crs_generic_intf.c",
					]
				}

				Group {
					name: "ls_mkb - lad_block"

					files: [
						"src/modfem/ls_mkb/lad_block/lah_block.h",
						"src/modfem/ls_mkb/lad_block/las_block_intf.c",
						"src/modfem/ls_mkb/lad_block/las_block_util.c",
					]
				}

				Group {
					name: "ls_mkb - lsd_mkb_core"

					files: [
						"src/modfem/ls_mkb/lsd_mkb_core/lsh_mkb_core.h",
						"src/modfem/ls_mkb/lsd_mkb_core/lsh_mkb_core_fem_intf.h",
						"src/modfem/ls_mkb/lsd_mkb_core/lss_mkb_core.c",
					]
				}

				Group {
					name: "ls_mkb - lsd_mkb_mumps"

					condition: false

					files: [
						"src/modfem/ls_mkb/lsd_mkb_mumps/lsh_mkb_mumps.h",
						"src/modfem/ls_mkb/lsd_mkb_mumps/lss_mkb_mumps.c",
					]
				}

				Group {
					name: "ls_mkb - lsd_mkb_pardiso"

					condition: modfem.config.mkbDirectSolver === "pardiso"

					files: [
						"src/modfem/ls_mkb/lsd_mkb_pardiso/lsh_mkb_pardiso.h",
						"src/modfem/ls_mkb/lsd_mkb_pardiso/lss_mkb_pardiso.c",
					]
				}

				Group {
					name: "ls_mkb - lsd_mkb_superlu"

					condition: false

					files: [
						"src/modfem/ls_mkb/lsd_mkb_superlu/README",
						"src/modfem/ls_mkb/lsd_mkb_superlu/lsh_mkb_superlu.h",
						"src/modfem/ls_mkb/lsd_mkb_superlu/lss_mkb_superlu.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/colamd.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/colamd.h",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dcolumn_bmod.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dcolumn_dfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dcomplex.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dcopy_to_ucol.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ddiagonal.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dgscon.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dgsequ.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dgsisx.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dgsitrf.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dgsrfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dgssv.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dgssvx.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dgstrf.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dgstrs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dlacon.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dlacon2.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dlamch.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dlangs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dlaqgs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dldperm.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dmemory.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dmyblas2.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dpanel_bmod.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dpanel_dfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dpivotL.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dpivotgrowth.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dpruneL.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dreadhb.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dreadrb.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dreadtriple.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dsnode_bmod.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dsnode_dfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dsp_blas2.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dsp_blas3.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dutil.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/dzsum1.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/get_perm_c.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/heap_relax_snode.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/html_mainpage.h",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/icmax1.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_ccolumn_dfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_ccopy_to_ucol.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_cdrop_row.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_cpanel_dfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_cpivotL.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_csnode_dfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_dcolumn_dfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_dcopy_to_ucol.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_ddrop_row.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_dpanel_dfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_dpivotL.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_dsnode_dfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_heap_relax_snode.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_relax_snode.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_scolumn_dfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_scopy_to_ucol.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_sdrop_row.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_spanel_dfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_spivotL.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_ssnode_dfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_zcolumn_dfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_zcopy_to_ucol.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_zdrop_row.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_zpanel_dfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_zpivotL.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ilu_zsnode_dfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/input_error.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/izmax1.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/lsame.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/mark_relax.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/mc64ad.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/memory.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/mmd.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/qselect.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/relax_snode.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/scolumn_bmod.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/scolumn_dfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/scomplex.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/scopy_to_ucol.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/scsum1.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/sdiagonal.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/sgscon.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/sgsequ.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/sgsisx.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/sgsitrf.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/sgsrfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/sgssv.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/sgssvx.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/sgstrf.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/sgstrs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/slacon.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/slacon2.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/slamch.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/slangs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/slaqgs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/sldperm.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/slu_Cnames.h",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/slu_cdefs.h",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/slu_dcomplex.h",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/slu_ddefs.h",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/slu_scomplex.h",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/slu_sdefs.h",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/slu_util.h",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/slu_zdefs.h",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/smemory.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/smyblas2.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/sp_coletree.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/sp_ienv.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/sp_preorder.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/spanel_bmod.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/spanel_dfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/spivotL.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/spivotgrowth.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/spruneL.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/sreadhb.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/sreadrb.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/sreadtriple.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ssnode_bmod.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ssnode_dfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ssp_blas2.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/ssp_blas3.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/superlu_enum_consts.h",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/superlu_timer.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/supermatrix.h",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/sutil.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/util.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_seq/xerbla.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/await.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/cholnzcnt.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/colamd.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/colamd.h",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/dclock.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/dgscon.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/dgsequ.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/dgsrfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/dgstrs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/dlacon.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/dlamch.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/dlangs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/dlaqgs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/dmatgen.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/dmyblas2.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/dpivotgrowth.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/dreadhb.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/dreadmt.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/dreadrb.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/dsp_blas2.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/dsp_blas3.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/get_perm_c.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/heap_relax_snode.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/lsame.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/mmd.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pdgssv.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pdgssvx.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pdgstrf.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pdgstrf_bmod1D.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pdgstrf_bmod1D_mv2.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pdgstrf_bmod2D.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pdgstrf_bmod2D_mv2.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pdgstrf_column_bmod.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pdgstrf_column_dfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pdgstrf_copy_to_ucol.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pdgstrf_factor_snode.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pdgstrf_init.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pdgstrf_panel_bmod.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pdgstrf_panel_dfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pdgstrf_pivotL.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pdgstrf_pruneL.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pdgstrf_snode_bmod.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pdgstrf_snode_dfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pdgstrf_thread.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pdgstrf_thread_finalize.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pdgstrf_thread_init.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pdmemory.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pdsp_defs.h",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pdutil.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pmemory.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pxgstrf_finalize.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pxgstrf_mark_busy_descends.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pxgstrf_pruneL.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pxgstrf_relax_snode.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pxgstrf_scheduler.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pxgstrf_super_bnd_dfs.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pxgstrf_synch.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/pxgstrf_synch.h",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/qrnzcnt.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/slu_mt_Cnames.h",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/slu_mt_cdefs.h",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/slu_mt_ddefs.h",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/slu_mt_machines.h",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/slu_mt_util.h",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/sp_coletree.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/sp_colorder.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/sp_ienv.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/superlu_timer.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/supermatrix.h",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/util.c",
						"src/modfem/ls_mkb/lsd_mkb_superlu/superlu_threads/xerbla.c",
					]
				}

				Group {
					name: "ls_mkb - lsd_mkb_viennacl"

					condition: false

					files: [
						"src/modfem/ls_mkb/lsd_mkb_viennacl/lsh_mkb_viennacl.h",
						"src/modfem/ls_mkb/lsd_mkb_viennacl/lss_mkb_viennacl.cpp",
					]
				}

				Group {
					name: "ls_mkb - lsd_ns_supg_ext"

					condition: false

					files: [
						"src/modfem/ls_mkb/lsd_ns_supg_ext/lsh_ns_supg_ext_intf.h",
						"src/modfem/ls_mkb/lsd_ns_supg_ext/lss_ns_supg_ext_intf.c",
					]
				}
			}
		}

	}
}
