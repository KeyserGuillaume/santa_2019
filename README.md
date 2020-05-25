# santa_2019
My work on the [2019 Kaggle Santa Problem](https://www.kaggle.com/c/santa-workshop-tour-2019/overview/description) 

Best score 69085.55 - late submission so no rank to reflect it. I worked on this for several months after the end of the competition on my free time between my master's courses.

## Heuristic

LocalSearch.h/.cpp contain the implementation of a localsearch (restart was coded as an afterthought, there's nothing to hightlight it in the code. I only found disappointing solutions with this method though (score around 73000).

## Flow method

Flow.h/.cpp uses flow theory to compute lower bounds on the problem (max flow min cost). The lower bounds are loose on the day costs, and those bounds are pretty bad overall. They still yielded solutions down to 71400 with some basic branching algorithms, additional constraints and code optimization.

## Lagrangian Relaxation

day_sub_problem.h / .cpp use dynamic programing to compute an excellent lower bound, a few dozen of euros short of the optimal solution. However I struggled to make use of it because the computations are relatively slow (around 1sec). The principle is very simple: in a simple (non-linear) modelization of the problem there should be a constraint that states that each family has only exactly one assigned day. Dualizing this constraint enables one to solve the problem with dynamic programing. We just need to compute a knapsack problem for every day and every possible number of people assigned to the day. By solving these problems with dynamic programing, it's only one problem per day. Then more dynamic programing takes the values we just computed and returns the lower bound. 

The lower bounds depends on the dual values we set when we dualized the constraint. I obviously did not invent anything, this technique is named Lagrangian Relaxation. The highest lower bound is found with a subgradient descent, as can often be done in Lagrangian Relaxation.

My best solutions were found with this method, thanks to a greedy algorithm which is guided by the lower bound.

## There are files in which I tried to implement Dantzig-Wolfe decomposition, with the same solving process as the Lagrangian Relaxation. Unfortunately I never got it to work. I found it rather more difficult to implement. I hoped I could find an optimal solution with it, but I had never worked with this algorithm before and would have needed more time to master it.
