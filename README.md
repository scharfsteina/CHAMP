# An Implementation of the CHAMP Algorithm
(in progress)
Working with Dartmouth Applied Mathematics Professor Peter Mucha as a URAD Scholar, 
I implemented the CHAMP algorithm that was outlined in his paper 
[*Post-Processing Partitions to Identify Domains of Modularity Optimization*](https://www.mdpi.com/1999-4893/10/3/93).
The goal of the CHAMP Algorithm is to make community detection analysis for networks more straightforward.

### Example 1: Weighted Karate Club Network
Calling ```run_champ(karate, weighted = TRUE, name = "Karate Club")```, you get the following summary

<img width="690" alt="Screen Shot 2022-06-11 at 1 55 42 PM" src="https://user-images.githubusercontent.com/54069119/173202940-2fb4da61-41c0-4d09-a1e2-988323116812.png">

This highlights the performance of distinct partitions, where performance is measured by stability to variations in the values of gamma (```gamma_range```). Based on this metric, we find that in the Weighted Karate Club Network, the __ and __ cluster solutions perform the best.
The partitions of the __ cluster solution are as follows ()
insert image
The partitions of the __ cluster solution are as follows ()
image

In the ```figures``` folder, we get five different figures. Figures 1,2,and 3 display three different visualizations of the CHAMP algorithm (modularity lines for the top partitions). The other two figures are the top two partitions determined by the stability metric.

![figure1](https://user-images.githubusercontent.com/54069119/173202965-fb22e3fa-8bf7-42f4-840f-0949badc0ff8.png)
![figure2](https://user-images.githubusercontent.com/54069119/173202966-76981446-f3d7-436e-98ec-9c8ca5c71b1b.png)
![figure3](https://user-images.githubusercontent.com/54069119/173202967-3dab71ac-c285-4a70-bebd-510fcf15a407.png)

add figures here!
