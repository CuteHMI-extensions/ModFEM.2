/*! \mainpage ModFEM Index Page
 *
 * \section intro_sec Introduction
 *
 *The name of the framework - ModFEM stands for two things. The first is an abbre-
viation from the Modular Finite Element Method. This highlights the main modular
principle used in ModFEM. The second is an abbreviation from the Modifiable Finite
Element Method. This indicates main functional advantage of the framework.
ModFEM falls into the category of open-source scientific codes and libraries. It is
a computational framework for generating a parallel adaptive finite element modelling
applications. By its’ modular architecture the framework allows for introducing new
mesh algorithms, approximation spaces, adaptation algorithms or weak formulations.
 The Moors’ Law shows clearly, that the computing power grows, and will grow
further. For large scale applications the dominant role falls to parallel distributed
memory environments. We address this requirement by introducing parallel over-
lay modules in ModFEM. The modules are capable of parallelize sequential scientific
codes without the need of a separate process of developing a distributed application.
The multi-scale and multi-physics modelling is also supported in ModFEM ar-
chitecture. Modules are responsible for modelling a single physics phenomena. The
super-modules combine two or more problem modules, thus allowing the development
of a strongly coupled solution for complicated multi-physics phenomena. Same prin-
ciple can be used for multi-scale modelling. Using a super-module, one can benefit,
when simulating a complicated problem, by directly passing values between differ-
ent scale levels. The wrapper modules also support interfacing with standalone FEM
codes and applications using other methods than FEM (e.g. cellular automata).
 *
 * \section mod_sec Modular structure
 *
 *The ModFEM framework provides modular structure, with all interfaces defined and
several example modules already implemented. Modules are independent of each oth-
er. Each module can be easily exchanged with another module with the same func-
tionality, as long as it properly implements the module interface.
The whole program is split into several, currently less then ten, modules (Fig. 1).
Modules interact only through strictly defined interfaces. The interfaces have been
created also to maintain the similarity to the classic FEM codes formulation:
1. Defining a problem and choosing FEM weak formulation.
2. Selecting appropriate approximation function space
3. Selecting conforming geometric discretization with mesh type and properties.
4. Adjusting solver and its preconditioner.
5. Thus you have well-defined application for well-defined problem.
All modules are implemented for cross platform usage and all mainstream compil-
ers support. Process of application building is provided by cross platform compilation
and linking framework. Whole FEM-application logic have been hidden behind the
seven orthogonal interfaces.
 *
 *
 */
