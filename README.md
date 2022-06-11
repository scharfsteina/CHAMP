# An Implementation of the CHAMP Algorithm
Working with Dartmouth Applied Mathematics Professor Peter Mucha as a URAD Scholar, 
I implemented the CHAMP algorithm that was outlined in his paper 
[*Post-Processing Partitions to Identify Domains of Modularity Optimization*](https://www.mdpi.com/1999-4893/10/3/93).
The goal of the CHAMP Algorithm is to make community detection analysis for networks more straightforward.

### Example 1: Karate Club Network
Calling ```run_champ(karate, weighted = TRUE, name = "Karate Club")```, you get the following summary

<img width="690" alt="Screen Shot 2022-06-11 at 1 55 42 PM" src="https://user-images.githubusercontent.com/54069119/173202940-2fb4da61-41c0-4d09-a1e2-988323116812.png">

and figures saved in the figures folder:

![figure1](https://user-images.githubusercontent.com/54069119/173202965-fb22e3fa-8bf7-42f4-840f-0949badc0ff8.png)
![figure2](https://user-images.githubusercontent.com/54069119/173202966-76981446-f3d7-436e-98ec-9c8ca5c71b1b.png)
![figure3](https://user-images.githubusercontent.com/54069119/173202967-3dab71ac-c285-4a70-bebd-510fcf15a407.png)
