# An Implementation of the CHAMP Algorithm
Working with Dartmouth Applied Mathematics Professor Peter Mucha as a URAD Scholar, 
I implemented the CHAMP algorithm that was outlined in his paper 
[*Post-Processing Partitions to Identify Domains of Modularity Optimization*](https://www.mdpi.com/1999-4893/10/3/93).

For those unfamiliar with the paper and the evolution of community detection in R: Community detection is the process of grouping a network into communities often on the basis of similarity or connectivity. The ongoing and most common practice for community detection analysis in R is centered around calling ```cluster_leiden``` on a network. Often, the function is called using the default value for the ```resolution_parameter``` ($\gamma = 1$). This default value, most of the time, is incorrect. At the simplest level, the goal of the CHAMP algorithm is to vary the ```resolution_parameter``` in order to find its optimal value(s).


### Example 1: Weighted Karate Club Network
Calling ```champ = run_champ(karate, weighted = TRUE, name = "Karate Club")```, you get the following summary

<img width="697" alt="Screen Shot 2022-06-16 at 7 23 44 PM" src="https://user-images.githubusercontent.com/54069119/174203768-912a896c-9979-42cf-b99e-63a5ef9f35d0.png">

This highlights the performance of distinct partitions, where performance is measured by stability to variations in the values of gamma (```gamma_range```). Based on this metric, we find that in the Weighted Karate Club Network, the 2 and 4 cluster solutions perform the best.

The partition of the 2 cluster solution can be represented by the following two figures (```communities(champ$partitions[[1]])``` and ```figures\partition1```)

<img width="763" alt="Screen Shot 2022-06-17 at 1 52 00 PM" src="https://user-images.githubusercontent.com/54069119/174392590-ddcfbe77-d52d-449f-bfa6-d77ffe9d721e.png">
<img src="figures/partition1.png"/>

The partition of the 4 cluster solution can be represented by the following two figures (```communities(champ$partitions[[4]])``` and ```figures\partition1```)

<img width="754" alt="Screen Shot 2022-06-17 at 1 52 46 PM" src="https://user-images.githubusercontent.com/54069119/174392655-1751a5db-cc04-4a95-9ebb-8fb8b42b0b87.png">
<img src="figures/partition4.png"/>

In the ```figures``` folder, we also find three other figures. They display three different visualizations of the CHAMP algorithm (modularity lines for the top partitions).

<img src="figures/figure1.png"/>
<img src="figures/figure2.png"/>
<img src="figures/figure3.png"/>

These figures provide us with a visual of which partitions are on the upper envelope as modularity varies. Those partitions that are on the upper envelope for the longest (partition 1 and 4) may be the most stable solutions.
