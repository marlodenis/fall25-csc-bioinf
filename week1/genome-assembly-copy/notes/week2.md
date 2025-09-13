# Week 2

## establish bounds for recurrence

3 methods

### substitution method

1. guess the form of the solution
2. verify by induction
3. solve for constants

### recursion tree

1. draw a recursion tree
2. sum every row
3. get the height
4. total sum

### master theorem

* only solve specific recursive structures

## recursion

### idea

* divide the problem int subproblem
* conquer each subproblem recursively
* combine the solution

### examples

1. merge sort
	* T(n)=2T(n/2)+O(n) 
	* -> T(n)=theta(nlgn)
2. binary search
	* T(n)=T(n/2)+theta(1) 
	* -> T(n)=theta(lgn)
3. powering a number: `x*x* .. *x`
	* T(n)=T(n/2)+theta(1)
	* -> T(n)=theta(lgn)
4. Fibonacci numbers: F(n)=F(n-1)+F(n-2), if n>=2
	* naive recursion: 
		- T(n)=omega(thai^n)
		- thai=(1+sqrt(5))/2
	* bottom up
		- O(n)
	* formula
		- F(n)=thai^n/sqrt(5)
		- theta(lgn) <- powering thai n times
		- but the answer is a float...
	* powering the matrix
		- [[F(n+1), F(n)],[F(n), F(n-1)]]=[[1,1],[1,0]]^n
		- theta(lgn) <- powering matrix n times
5. matrix multiplication
	* square matrix multiplication: O(n^3)
	* simple recursion: 8 times, O(n^3)
	* strassen recursion: 7 times, O(n^(log(2)7))~O(n^2.81)
