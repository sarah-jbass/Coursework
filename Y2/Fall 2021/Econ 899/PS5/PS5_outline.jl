using Random, Interpolation, Optim

#=
Value function iteration
    When normally we would have determined budget and use the grid search to look at all k' and pick the best one. Now what we do is a search over the continuous approximation using the linear interpolation. Then use optimization routine to numerically solve for the k'.

First forecast k'
Conditional on productivity, that gives l today
Given k, l, can calculate r, w
Given k', get_index for interpolation
