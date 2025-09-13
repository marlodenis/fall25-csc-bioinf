# Week 1

## Insertion Sort

### pseudocode
```
insertion_sort(A, n)
	for j <- 2 to n
		do key <- A[j]
			i <- j - 1
			while i > 0 and A[i] > key
				do 
					A[i + 1] <- A[i]
					i <- i - 1
			A[i + i] <- key
```

### Validation of a loop (数学归纳法)

How to check if a loop is right?

0. loop invariant: a statement, which is true before the loop and remains true all the way even when the loop is ended
	* `A[1 ... j - 1]` is already sorted
1. initialization: initializaiton fits loop invariant
	* `A[1]` is sorted
2. maintenance: loop invariant remains true after each change
	* insert `A[j]` into `A[1 ... j - 1]` at the right spot 
3. termination: after the loop, loop invariant is still true
	* `A[1 ... n]` is sorted

## Analysis

1. worst-case
	* to avoid worst-case for insertion sort: randomly shuffle `A` first
2. average-case
	* depend on the distribution of the data
3. best-case (not cencerned)

## Big Idea of Algorithm Analysis

e.g. `T(n) = an^2 + bn + c` => `O(n^2)`, the order of growth
1. look at the growth of `T(n)` as `n -> infinity` (leading term: `an^2`)
2. ignore machine-dependent constant (`a`)

## merge sort

### analysis

`T(n) = D(n) + a*T(n/b) + C(n)` => `O(nlgn)`
* 计算：递归树

### idea: divide and conquer

1. divide a sub-problem `n/b`
2. conquer each sub-problem `T(n/b)` by solving them recursively
3. combine solution

### pseudocode

```
merge_sort(A, p, r)
	if p < r
		then q <- floor((p + r) / 2)
			merge_sort(A, p, q)
			merge_sort(A, q + 1, r)
			merge(A, p, q, r)
```

```
merge(A, p, q, r)
	n1 <- q - p + 1
	n2 <- r - q
	create arrays L[1 .. n1 + 1] and R[1 .. n2 + 1]
	for i <- 1 to n1
		do L[i] <- A[p + i - 1]
	for j <- 1 to n2
		do R[j] <- A[q + j]
	L[n1 + 1] <- infinity
	R[n2 + 1] <- infinity
	i <- 1
	j <- 1
	for k <- p to r
		do if L[i] <= R[j]
			then A[k] <- L[i]
				i <- i + 1
			else A[k] <- R[j]
				j <- j + 1
```

### validation

loop invariant: `A[p ... k - 1]` is already sorted

### application

* counting inversions: count inversions (i < j but ai > aj, this is one inversion) when merging

## notations 渐近记号

### big O: upper bound 渐近上界

1. f(n) = O(g(n)) 
	* means ther are constant c > 0, exist n0 > 0, such that `0 <= f(n) <= c*g(n)` for all n >= n0
	* e.g. 2(n^2) = O(n^3)
	* e.g. 2(n^2) = O(n^2)
	* e.g. 2(n^2) = O(n) is wrong!
2. O(g(n)) = { f(n), where `0 <= f(n) <= c*g(n)` for all n > n0, ther is constant c, exist n0 > 0 }
	* e.g. f(n) = n^3 + O(n^2) means that there is a function h(n) belongs to O(n^2) such that f(n) = n^3 + h(n)
	* e.g. n^2 + O(n) = O(n^2) means for any f(n) belonging to O(n), there is an h(n) belonging to O(n^2) such that n^2 + f(n) = h(n)

### big omega: lower bound 渐近下界

1. omega(g(n)) = { f(n) there is const c, exist n0 > 0, such that `0 <= c*g(n) <= f(n)` for all n > n0 }
	* e.g. sqrt(n) = omega(lgn)
	* e.g. 2n^2 = omega(n)
	* e.g. 2n^2 = omega(n^2)
	* e.g. 2n^2 = omega(n^3) is wrong!

### big theta 渐近紧确界

1. theta(g(n)) = O(g(n)) and omega(g(n))
2. theta(g(n)) = { f(n) there is const c1, c2 exist n0 > 0, such that `c1*g(n) <= f(n) <= c2*g(n)` for all n > n0 }

## find a peak (anyone is ok)

1. 1-d: O(lgn)
2. 2-d: O(nlgn)

## Summary

1. What is algorithm
2. Why is it important
3. correctness
	* loop invariant: inition, maintanence, termination
	* e.g. insertion sort, merge sort
4. time complexity
	* drop low order, ignore constant
	* worst case, average case, best case(not concerned)
5. notations
	* big O: upper bound (<= big O)
	* big theta: (~~)
	* big omega: lower bound (>= big omega)
	* small o: upper bound, not asymptotically tight
