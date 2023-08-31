al form' see,

   * Donald E. Knuth, The Art of Computer Programming: Fundamental
     Algorithms (Vol 1, 3rd Ed, 1997), Addison-Wesley, ISBN 0201896850.
     Section 1.3.3, An Unusual Correspondence, p.178–179.


File: gsl-ref.info,  Node: Combinations,  Next: Multisets,  Prev: Permutations,  Up: Top

10 Combinations
***************

This chapter describes functions for creating and manipulating
combinations.  A combination c is represented by an array of k integers
in the range 0 to n - 1, where each value c_i occurs at most once.  The
combination c corresponds to indices of k elements chosen from an n
element vector.  Combinations are useful for iterating over all
k-element subsets of a set.

The functions described in this chapter are defined in the header file
‘gsl_combination.h’.

* Menu:

* The Combination struct::
* Combination allocation::
* Accessing combination elements::
* Combination properties::
* Combination functions::
* Reading and writing combinations::
* Examples: Examples<5>.
* References and Further Reading: References and Further Reading<6>.


File: gsl-ref.info,  Node: The Combination struct,  Next: Combination allocation,  Up: Combinations

10.1 The Combination struct
===========================

 -- Type: gsl_combination

     A combination is defined by a structure containing three
     components, the values of n and k, and a pointer to the combination
     array.  The elements of the combination array are all of type
     ‘size_t’, and are stored in increasing order.  The *note
     gsl_combination: 41f. structure looks like this:

          typedef struct
          {
            size_t n;
            size_t k;
            size_t *data;
          } gsl_combination;


File: gsl-ref.info,  Node: Combination allocation,  Next: Accessing combination elements,  Prev: The Combination struct,  Up: Combinations

10.2 Combination allocation
===========================

 -- Function: gsl_combination * gsl_combination_alloc (size_t n,
          size_t k)

     This function allocates memory for a new combination with
     parameters ‘n’, ‘k’.  The combination is not initialized and its
     elements are undefined.  Use the function *note
     gsl_combination_calloc(): 422. if you want to create a combination
     which is initialized to the lexicographically first combination.  A
     null pointer is returned if insufficient memory is available to
     create the combination.

 -- Function: gsl_combination * gsl_combination_calloc (size_t n,
          size_t k)

     This function allocates memory for a new combination with
     parameters ‘n’, ‘k’ and initializes it to the lexicographically
     first combination.  A null pointer is returned if insufficient
     memory is available to create the combination.

 -- Function: void gsl_combination_init_first (gsl_co