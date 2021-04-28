# HHD_and_Trait_Performance

2021, Alexander Strang
All rights reserved.

Contents:

1. Building the Curl: contains a series of example codes that show how to build the curl given an arbitrary network topology. 
   These are meant as examples and are not optimized. Make_Fundamental_Basis builds the curl associated with a fundamental cycle basis,
   while Make_Weakly_Fundamental_Basis builds the curl associated with a weakly fundamental basis. Both are based on search procedures.
   The weakly fundamental basis search typically produces more natural basis with smaller loops, fewer reused edges, and fewer nested loops.
2. Performing the HHD: contains three tools for performing some of the key steps in the HHD. Get_Topology builds a variety of data structures 
   representing the topology of the network given a set of observed events. These structures could also be built by hand, and used as input to
   the other codes. Get_Operators builds the gradient, and an overdetermined version of the curl (all triangles) if the network is complete. Otherwise
   the curl should either be built by hand (if the graph structure is known), or using one of the search procedures. Perform_HHD performs the HHD 
   given an edge flow and the operators.
3. Replicating Figures and Sampling Networks: contains the code used to produce figure 7 (scatter plot representing the distribution of random
   three competitor networks on the transitive intransitive plane), figure 10 (heat maps showing the distribution of randomly sampled networks
   using the proposed tunable null models), and figure SM3 in the supplement (comparison of the correlation rho between the fair fight and 
   press your advantage models for uniform, exponential, and Gaussian distributed traits). The code to generate figure 10 used the tight_subplot.mat
   function available on the MathWorks file exchange. The license for tight_subplot.mat is provided in the folder, and all credit for the command
   goes to Pekka Kumpulainen. 

Excluding tight_subplot.mat all codes were written and tested by Alexander Strang. 
