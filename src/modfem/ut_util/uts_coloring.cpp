#ifndef _uts_coloring_
#define _uts_coloring_

#include <vector>
#include <algorithm>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/shared_array_property_map.hpp>
#include <boost/graph/smallest_last_ordering.hpp>
#include <boost/graph/sequential_vertex_coloring.hpp>
#include <modfem/mmh_intf.h>

#include <modfem/uth_log.h>
#include <modfem/uth_system.h>


int utr_generate_mesh_coloring( // return number of colors
		const int Mesh_id, // IN: mesh id to color
		int * Elem_colors) // OUT: array of size number of active elems containing colors (numbers)
{

	mf_check_mem(Elem_colors);

	// Constructing and filling the graph.
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
	const int nElems = mmr_get_nr_elem(Mesh_id);
	Graph g(nElems);

	int nel = 0;
	while ( (nel = mmr_get_next_act_elem(Mesh_id, nel)) != 0) {
		int neigs[MMC_MAXELFAC + 1] = {0};
		mmr_el_eq_neig(Mesh_id, nel, neigs, NULL);
		for (int i = 1; i < neigs[0]; ++i) {
			add_edge(nel, neigs[i], g);
		}
	}

	std::fill(Elem_colors, Elem_colors + nElems, 0);
	boost::smallest_last_vertex_ordering(g, Elem_colors);

	std::vector<int> color_counts;
	int cur_color = -1;
	do {
		++cur_color;
		color_counts.push_back( std::count(Elem_colors, Elem_colors + nElems, cur_color) );
	}
	while (color_counts[cur_color] > 0);

	mf_check_mem(Elem_colors);

	return cur_color;
}

#ifdef __cplusplus
extern "C"
{
#endif

struct element {
	int id;
	element * next;
};


int utr_first_fit_coloring_with_front( // return number of colors
		const int nElems, // IN: number of elements to color
		const std::vector<int> * elem2elem,  // IN: array elements neighborhood
		//int ** elem2elem, // IN: array elements neighborhood
		int * Elem_colors) // OUT: array of size number of active elems containing colors (numbers)
{
	const int MAX_COLOR = 500;
	int i, color_id, flag;
	int current_elem = 0, nr_color = 0, colored_elems = 0, nr_colored;
	element * head = NULL;
	element * tail = NULL;
	bool * avail_elem = (bool *)malloc(nElems * sizeof(bool));

	for (i = 0; i < nElems; ++i) {
		avail_elem[i] = true;
		Elem_colors[i] = -1;
	}

	avail_elem[current_elem] = false;
	Elem_colors[current_elem] = 0;
	nr_colored = 1;
	nr_color++;
	bool avail_color[MAX_COLOR];
	for (i = 0; i < MAX_COLOR; ++i) {
		avail_color[i] = true;
	}

	for (i = 1; i <= elem2elem[current_elem][0]; ++i) { //add neighbors to check list
		if (avail_elem[elem2elem[current_elem][i]] == true) {
			element * nowy = new element;
			nowy->id = elem2elem[current_elem][i];
			nowy->next = NULL;
			avail_elem[elem2elem[current_elem][i]] = false;
			if (head == NULL) {
				head = nowy;
				tail = nowy;
			}
			else {
				tail->next = nowy;
				tail = nowy;
			}
		}
	}



	while (nr_colored < nElems) {

		if (head != NULL) { //is neighbor
			current_elem = head->id;
			element * nowy = head;
			head = head->next;
			if (head == NULL)tail = NULL;
			delete nowy;
			avail_elem[current_elem] = false;
			for (i = 1; i <= elem2elem[current_elem][0]; ++i) { //checking neighbors colors
				color_id = Elem_colors[elem2elem[current_elem][i]];
				if (color_id != -1) {
					avail_color[color_id] = false;
				} else if (avail_elem[elem2elem[current_elem][i]] == true) { //add neighbors to check list
					element * nowy = new element;
					nowy->id = elem2elem[current_elem][i];
					nowy->next = NULL;
					avail_elem[elem2elem[current_elem][i]] = false;
					if (head == NULL) {
						head = nowy;
						tail = nowy;
					}
					else {
						tail->next = nowy;
						tail = nowy;
					}
				}
			}
			flag = 0;
			for (color_id = 0; color_id < nr_color; ++color_id) { //set color with minimal id
				if (avail_color[color_id] == true) {
					Elem_colors[current_elem] = color_id;
					nr_colored++;
					flag = 1;
					break;
				}
			}
			if (flag == 0) {                      //set color with new id
				Elem_colors[current_elem] = nr_color;
				nr_color++;
				nr_colored++;
				/* realloc if too much colors
				if(nr_color>MAX_COLOR){
				  MAX_COLOR=nr_color+10;
				  avail_color = (bool*) realloc (avail_color, MAX_COLOR * sizeof(bool));
				}
				*/
			}
			for (color_id = 0; color_id < nr_color; ++color_id)avail_color[color_id] = true;
		}
		else { //isnt neighbor
			for (current_elem = 0; current_elem < nElems && Elem_colors[current_elem] != -1; ++current_elem);
			if (Elem_colors[current_elem] == -1) {
				avail_elem[current_elem] = false;
				for (i = 1; i <= elem2elem[current_elem][0]; ++i) { //checking neighbors colors
					color_id = Elem_colors[elem2elem[current_elem][i]];
					if (color_id != -1) {
						avail_color[color_id] = false;
					} else if (avail_elem[elem2elem[current_elem][i]] == true) { //add neighbors to check list
						element * nowy = new element;
						nowy->id = elem2elem[current_elem][i];
						nowy->next = NULL;
						avail_elem[elem2elem[current_elem][i]] = false;
						if (head == NULL) {
							head = nowy;
							tail = nowy;
						}
						else {
							tail->next = nowy;
							tail = nowy;
						}
					}
				}
				flag = 0;
				for (color_id = 0; color_id < nr_color; ++color_id) { //set color with minimal id
					if (avail_color[color_id] == true) {
						Elem_colors[current_elem] = color_id;
						flag = 1;
						nr_colored++;
						break;
					}
				}
				if (flag == 0) {                    //set color with new id
					Elem_colors[current_elem] = nr_color;
					nr_color++;
					nr_colored++;
					/* realloc if too much colors
					if(nr_color>MAX_COLOR){
					  MAX_COLOR=nr_color+10;
					  avail_color = (bool*) realloc (avail_color, MAX_COLOR * sizeof(bool));
					}
					*/
				}
				for (color_id = 0; color_id < nr_color; ++color_id)avail_color[color_id] = true;
			}
		}
	}//while

	free(avail_elem);

	return nr_color;
}

int utr_first_fit_coloring_without_front( // return number of colors
		const int nElems, // IN: number of elements to color
		const std::vector<int> * elem2elem, // IN: array elements neighborhood
		//int ** elem2elem, // IN: array elements neighborhood
		int * Elem_colors) // OUT: array of size number of active elems containing colors (numbers)
{
	const int MAX_COLOR = 500;
	int i, color_id, flag;
	int current_elem = 0, nr_color = 0, nr_colored;
	bool avail_color[MAX_COLOR];


	for (current_elem = 0; current_elem < nElems; ++current_elem) {
		Elem_colors[current_elem] = -1;
	}

	Elem_colors[0] = 0;
	nr_colored = 1;
	nr_color++;


	for (current_elem = 0; current_elem < nElems; ++current_elem) {
		for (color_id = 0; color_id < MAX_COLOR; ++color_id) {
			avail_color[color_id] = true;
		}
		for (i = 1; i <= elem2elem[current_elem][0]; ++i) { //get neighbors color
			color_id = Elem_colors[elem2elem[current_elem][i]];
			if (color_id != -1) {
				avail_color[color_id] = false;
			}
		}
		for (color_id = 0, flag = 0; color_id < nr_color; ++color_id) { //set color with minimal id
			if (avail_color[color_id] == true && flag == 0) {
				Elem_colors[current_elem] = color_id;
				nr_colored++;
				flag = 1;
				break;
			}
		}
		if (flag == 0) { //set color with new id
			Elem_colors[current_elem] = nr_color;
			nr_color++;
			nr_colored++;
			/* realloc if too much colors
			if(nr_color>MAX_COLOR){
			  MAX_COLOR=nr_color+10;
			  avail_color = (bool*) realloc (avail_color, MAX_COLOR * sizeof(bool));
			}
			*/
		}
	}

	return nr_color;
}


int utr_boost_sequential_coloring( // return number of colors
		const int nElems, // IN: number of elements to color
		int ** elem2elem, // IN: array elements neighborhood
		int * Elem_colors) // OUT: array of size number of active elems containing colors (numbers)
{

	int i, j;
	int nr_color = 0;
	mf_check_mem(Elem_colors);

	// Constructing and filling the graph.
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
	//const int nElems = mmr_get_nr_elem(Mesh_id);
	Graph g(nElems);


	for (i = 0; i < nElems; i++) {
		for (j = 1; j <= elem2elem[i][0]; j++)
			add_edge(i, elem2elem[i][j], g);
	}


	nr_color = boost::sequential_vertex_coloring(g, Elem_colors);

	return nr_color;
}


int utr_generate_int_ent_coloring( // return number of colors
		const int nElems, // IN: number of elements to color
		const std::vector<int> * elem2elem, // IN: array elements neighborhood
		//int ** elem2elem, // IN: array elements neighborhood
		int * Elem_colors) // OUT: array of size number of active elems containing colors (numbers)
{
	//int nr_colors=utr_boost_sequential_coloring(nElems,elem2elem,Elem_colors);
	//int nr_colors=utr_first_fit_coloring_without_front(nElems,elem2elem,Elem_colors);
	int nr_colors = utr_first_fit_coloring_with_front(nElems, elem2elem, Elem_colors);


	/*kbw
	printf("nr_colors=%d\n",nr_colors);
	for(i=0;i<nElems;i++){
	  printf("%d ",Elem_colors[i]);
	}
	/*kew*/


	return nr_colors;
}

int utr_color_int_ent_for_assembly_with_std_vector(
		int Problem_id,           /* in: Problem id */
		int Level_id,             /* in: Level id */
		int Nr_elems,             /* in: number of elements */
		int Nr_faces,             /* in: number of faces */
		int * L_int_ent_type,      /* in: integration entities type,
							   out: coloured  */
		int * L_int_ent_id,        /* in: integration entities id,
							   out: coloured  */
		int Nr_dof_ent,           /* number of dof entities dof_ents */
		int * Nr_int_ent_loc,     /* number of int_ents for each dof_ent */
		int ** L_int_ent_index,   /* list of indices in int_ent table for each dof_ent */
		int * Nr_colors_elems,    /* out: number of colours for elements  */
		int ** L_color_index_elems, /* out: crs table with index of L_int_ent.. where colour starts and ends  */
		int * Nr_colors_faces,    /* out: number of colours for faces  */
		int ** L_color_index_faces /* out: crs table with index of L_int_ent.. where colour starts and ends  */
)
{

	int i, j, k, l, flag;

	int nr_int_ent = Nr_elems + Nr_faces;

	std::vector<int> * int_ent_neig = new std::vector<int>[nr_int_ent + 1];

	for (i = 0; i < nr_int_ent; i++) {
		int_ent_neig[i].reserve(40);
		int_ent_neig[i].push_back(0);
	}

#ifdef TIME_TEST
	double time_create_graph = 0.0;
	double time_coloring = 0.0;
	double time_tmp;
	time_create_graph = time_clock();
#endif

	for (i = 0; i < Nr_dof_ent; i++) { // for each dof entity

		for (j = 0; j < Nr_int_ent_loc[i]; j++) { // for each integration entity containing dof entity

			int int_ent_current = L_int_ent_index[i][j]; // index of integration ent

			//if(L_int_ent_type[int_ent_current]==PDC_ELEMENT){
			if (int_ent_current < Nr_elems) { //element


				for (k = 0; k < Nr_int_ent_loc[i]; k++) {
					//for(k=j+1;k<Nr_int_ent_loc[i];k++){

					int neig_index = L_int_ent_index[i][k];

					//if(k!=j && L_int_ent_type[neig_index]==PDC_ELEMENT){
					//if(L_int_ent_type[neig_index]==PDC_ELEMENT){
					if (neig_index < Nr_elems) { //element
						(int_ent_neig[int_ent_current]).push_back(neig_index);
					}

				}

			}
			else { // if face


				for (k = 0; k < Nr_int_ent_loc[i]; k++) {
					//for(k=j+1;k<Nr_int_ent_loc[i];k++){


					int neig_index = L_int_ent_index[i][k];

					//if(k!=j && L_int_ent_type[neig_index]==PDC_FACE){
					//if(L_int_ent_type[neig_index]==PDC_FACE){
					if (neig_index >= Nr_elems) { // face
						(int_ent_neig[int_ent_current]).push_back(neig_index - Nr_elems);
					}
				}

			}

		}

	}


	#pragma omp parallel for
	for (i = 0; i < nr_int_ent; i++) {
		//unique sort - we dont need
		//std::sort( (int_ent_neig[i]).begin()+1, (int_ent_neig[i]).end() );
		//(int_ent_neig[i]).erase( std::unique( (int_ent_neig[i]).begin()+1, (int_ent_neig[i]).end() ),(int_ent_neig[i]).end());
		int_ent_neig[i][0] = int_ent_neig[i].size() - 1;
	}


	/*kbw
	for(i=0;i<(Nr_faces+Nr_elems);i++){
	  printf("%d=",i);
	  for(j=0;j<=int_ent_neig[i][0];j++){
		printf("%d ",int_ent_neig[i][j]);

	  }printf("\n");
	}
	/*kew*/

#ifdef TIME_TEST
	time_tmp = time_clock();
	time_create_graph = time_tmp - time_create_graph;
	time_coloring = time_clock();
#endif

	int * elem_color = (int *)malloc(Nr_elems * sizeof(int));
	int * face_color = (int *)malloc(Nr_faces * sizeof(int));

	*Nr_colors_elems = utr_generate_int_ent_coloring(Nr_elems, int_ent_neig, elem_color);

	printf("\nNr_colors_for_elems= %d\n", *Nr_colors_elems);

	*Nr_colors_faces = utr_generate_int_ent_coloring(Nr_faces, &int_ent_neig[Nr_elems], face_color);

	printf("Nr_colors_for_faces= %d\n", *Nr_colors_faces);

#ifdef TIME_TEST
	time_tmp = time_clock();
	time_coloring = time_tmp - time_coloring;
	printf("\nGraph creation for coloring: %lf\nColoring time: %lf\n", time_create_graph, time_coloring);
#endif

	//memory free
	for (i = 0; i < nr_int_ent; i++)(int_ent_neig[i]).clear();
	delete [] int_ent_neig;


	*L_color_index_elems = (int *)malloc((*Nr_colors_elems + 1) * sizeof(int));
	*L_color_index_faces = (int *)malloc((*Nr_colors_faces + 1) * sizeof(int));

	int * int_ent_list_id = (int *)malloc(nr_int_ent * sizeof(int));
	int * int_ent_list_type = (int *)malloc(nr_int_ent * sizeof(int));
	for (i = 0; i < nr_int_ent; i++) {
		int_ent_list_id[i] = L_int_ent_id[i];
		int_ent_list_type[i] = L_int_ent_type[i];
	}


	//reorder L_int_ent... list of elements using colors
	(*L_color_index_elems)[0] = 0;
	for (i = 0, k = 0; i < (*Nr_colors_elems); i++) {
		for (j = 0; j < Nr_elems; j++) {
			if (elem_color[j] == i) {
				L_int_ent_id[k] = int_ent_list_id[j];
				L_int_ent_type[k] = int_ent_list_type[j];

				//int_ent_list_id[k]=L_int_ent_id[j];
				//int_ent_list_type[k]=L_int_ent_type[j];

				k++;
			}
		}
		(*L_color_index_elems)[i + 1] = k;
	}

	//reorder L_int_ent... list of faces using colors
	(*L_color_index_faces)[0] = Nr_elems;
	for (i = 0, k = Nr_elems; i < (*Nr_colors_faces); i++) {
		for (j = 0; j < Nr_faces; j++) {
			if (face_color[j] == i) {
				L_int_ent_id[k] = int_ent_list_id[j + Nr_elems];
				L_int_ent_type[k] = int_ent_list_type[j + Nr_elems];

				//int_ent_list_id[k]=L_int_ent_id[j+Nr_elems];
				//int_ent_list_type[k]=L_int_ent_type[j+Nr_elems];

				k++;
			}
		}
		(*L_color_index_faces)[i + 1] = k;
	}

	free(elem_color);
	free(face_color);
	free(int_ent_list_id);
	free(int_ent_list_type);


	return (0);
}

int utr_first_fit_coloring_with_front_old( // return number of colors
		const int nElems, // IN: number of elements to color
		int ** elem2elem, // IN: array elements neighborhood
		int * Elem_colors) // OUT: array of size number of active elems containing colors (numbers)
{
	int MAX_COLOR = 500;
	int i, color_id, flag;
	int current_elem = 0, nr_color = 0, colored_elems = 0, nr_colored;
	element * head = NULL;
	element * tail = NULL;
	bool * avail_elem = (bool *)malloc(nElems * sizeof(bool));

	for (i = 0; i < nElems; ++i) {
		avail_elem[i] = true;
		Elem_colors[i] = -1;
	}

	avail_elem[current_elem] = false;
	Elem_colors[current_elem] = 0;
	nr_colored = 1;
	nr_color++;
	bool * avail_color = (bool *)malloc(MAX_COLOR * sizeof(bool));
	for (i = 0; i < MAX_COLOR; ++i) {
		avail_color[i] = true;
	}

	while (nr_colored < nElems) {
		for (i = 1; i <= elem2elem[current_elem][0]; ++i) { //add neighbors to check list
			if (avail_elem[elem2elem[current_elem][i]] == true) {
				element * nowy = new element;
				nowy->id = elem2elem[current_elem][i];
				nowy->next = NULL;
				avail_elem[elem2elem[current_elem][i]] = false;
				if (head == NULL) {
					head = nowy;
					tail = nowy;
				}
				else {
					tail->next = nowy;
					tail = nowy;
				}
			}
		}
		if (head != NULL) { //is neighbor
			current_elem = head->id;
			element * nowy = head;
			head = head->next;
			if (head == NULL)tail = NULL;
			delete nowy;
			avail_elem[current_elem] = false;
			for (i = 1; i <= elem2elem[current_elem][0]; ++i) { //checking neighbors colors
				color_id = Elem_colors[elem2elem[current_elem][i]];
				if (color_id != -1) {
					avail_color[color_id] = false;
				}
			}
			flag = 0;
			for (color_id = 0; color_id < nr_color; ++color_id) { //set color with minimal id
				if (avail_color[color_id] == true && flag == 0) {
					Elem_colors[current_elem] = color_id;
					nr_colored++;
					flag = 1;
					break;
				}
				avail_color[color_id] == true;
			}
			if (flag == 0) { //set color with new id
				Elem_colors[current_elem] = nr_color;
				nr_color++;
				nr_colored++;
				if (nr_color > MAX_COLOR) {
					avail_color = (bool *) realloc (avail_color, nr_color * sizeof(bool));
					for (color_id = 0; color_id < nr_color; ++color_id)avail_color[color_id] = true;
				}
			} else {
				for (color_id = 0; color_id < nr_color; ++color_id)avail_color[color_id] = true;
			};
		}
		else { //isnt neighbor
			for (current_elem = 0; current_elem < nElems && Elem_colors[current_elem] != -1; ++current_elem);
			if (Elem_colors[current_elem] == -1) {
				avail_elem[current_elem] = false;
				for (i = 1; i <= elem2elem[current_elem][0]; ++i) {
					color_id = Elem_colors[elem2elem[current_elem][i]];
					if (color_id != -1) {
						avail_color[color_id] = false;
					}
				}
				flag = 0;
				for (color_id = 0; color_id < nr_color; ++color_id) {
					if (avail_color[color_id] == true && flag == 0) {
						Elem_colors[current_elem] = color_id;
						flag = 1;
						nr_colored++;
						break;
					}
					avail_color[color_id] == true;
				}
				if (flag == 0) {
					Elem_colors[current_elem] = nr_color;
					nr_color++;
					nr_colored++;
					if (nr_color > MAX_COLOR) {
						avail_color = (bool *) realloc (avail_color, nr_color * sizeof(bool));
						MAX_COLOR = nr_color;
						for (color_id = 0; color_id < nr_color; ++color_id)avail_color[color_id] = true;
					}
				} else {
					for (color_id = 0; color_id < nr_color; ++color_id)avail_color[color_id] = true;
				}
			}
		}
	}//while

	free(avail_elem);
	free(avail_color);

	return nr_color;
}


int utr_generate_int_ent_coloring_old( // return number of colors
		const int nElems, // IN: number of elements to color
		int ** elem2elem, // IN: array elements neighborhood
		int * Elem_colors) // OUT: array of size number of active elems containing colors (numbers)
{
	//int nr_colors=utr_boost_sequential_coloring(nElems,elem2elem,Elem_colors);
	//int nr_colors=utr_first_fit_coloring_without_front(nElems,elem2elem,Elem_colors);
	int nr_colors = utr_first_fit_coloring_with_front_old(nElems, elem2elem, Elem_colors);


	/*kbw
	printf("nr_colors=%d\n",nr_colors);
	for(i=0;i<nElems;i++){
	  printf("%d ",Elem_colors[i]);
	}
	/*kew*/


	return nr_colors;
}




#ifdef __cplusplus
}
#endif

#endif // uts_coloring_
