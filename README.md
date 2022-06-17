# An Implementation of the CHAMP Algorithm
Working with Dartmouth Applied Mathematics Professor Peter Mucha as a URAD Scholar, 
I implemented the CHAMP algorithm that was outlined in his paper 
[*Post-Processing Partitions to Identify Domains of Modularity Optimization*](https://www.mdpi.com/1999-4893/10/3/93).
The goal of the CHAMP Algorithm is to make community detection analysis for networks more straightforward.

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
