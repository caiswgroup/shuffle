We do fine-grained comprehensive assessments for the tricks of the kissat-MAB-sat based parallel solvers in our GitHub repository.
The SC20 and SC21 instances are divided into several categories according to the benchmark documents of SAT Competitions.
The experiment settings are the same as in the paper.
In the attachments named "sat2020-detailed-results.txt" and "sat2021-detailed-results.txt", we report the comparison of kissat-MAB-sat with "+d", "+LA/RS", "+LA/RS+s", "+RS+d", "+RS+d+s" for each category.
From the summary table in the paper, we can learn that the crucial variable sampling method named LA shows better performance on satisfiable instances than diversification-based methods and SOTA solvers, while LA is weak on unsatisfiable instances. To gain a competitive solver with overall ability, we combine LA with clause sharing and diversification. For better adaption with clause sharing, we propose RS as the generalized version of LA, which further pushes the SOTA to the edge when combined with clause sharing.
From the comparison of kissat+LA/RS with pakis(kissat-MAB), i.e., kissat+d, we can learn that the LA/RS method shows better performance than diversification on some fine-grained categories, such as antibandwidth problem, school timetable problem, circuit multiplier, Stedman triples problem, hypertree decomposition, and so on.
For most of the subcategories, "+s" has a greater influence on "RS" than "LA", which confirms the adaptivity of the generalization. 

In the attachments called "sat2020-detailed-speedups.txt" and "sat2021-detailed-speedups.txt", we calculate the speedups with different CPU cores ranging from 4 to 64 and report the average speedups for each small class.
From the tables, we notice that the RS+s method shows better performance than LA+s when the core number is large for most categories.
Our methods show great potential to accelerate the solving processes in categories like circuit multiplier, abstract argumentation, minimal superpermutation, rbsat, fermat, coloring, and discrete logarithm.

