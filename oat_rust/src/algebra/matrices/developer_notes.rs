//! Notes for developers
//! 
//! # Design principles for matrix oracles
//!
//!  * Principle: define >=6 distinct oracle traits for ascending/descending, major/minor, ordered/unordered views
//!    * Justification
//!      * each variant serves a distinct purpose.  Ascending and descending views both make appearances in
//!        the computatino of persistent homology, for example, and unordered views are especially convenient
//!        for some applications because they have fewer "safety" requirements to ensure.
//!      * an individual project typically only requires a subset of the oracle traits to function; by keeping the
//!        traits separate we avoid forcing users to implement methods that they will not use.
//!    * Why it is preferred to alternative approaches
//!      * Rather than creating a single matrix struct that implements several different oracle traits,
//!        one could construct multiple different structs that all implement the same oracle trait (one struct
//!        for ascending major views, one for ascending minor views, etc.).  We prefer the current design over
//!        this alternative because (i) the alternative seems prone to forcing the user to define a large number
//!        of structs for the same underlying matrix, (ii) traits are designed to capture "shared behavior."  
//!        although several variants of the oracle traits exist, we expect each variant to be used extensively
//!        across domains.  Defining the behavior of each variant in a trait gives a simple, effective way to 
//!        formalize the shared behavior of each variant, and export this behavior safely across domains.  
//!        Contrast this with the alternative approach: suppose you are writing a function that takes an 
//!        oracle as an argument, which must implement `ViewMajorAscend`.  If you wanted to ensure that 
//!        this argument met the desired criteria without the trait, you would need to carefully define 
//!        the desired behavior in the docstrings for the function.  The existence of the trait relieves 
//!        this burden from documentation, since it gives a definition for behavior in explicit detail.
//!   * Principle: use different names for the method associated with each trait.
//!     * Justification
//!       * we often wish for a single struct to implement multiple oracle traits.  If the trait methods
//!         were all named identically, then we would have to disambiguate them using "fully qualified syntax", 
//!         which is highly verbose.
//!   * Principle: use associated types for ColIndex, RowIndex, Coefficient, and Ascending/Descending Major/Minor views.
//!     * Justification
//!       * helps the compiler with type inference; in the past, when we did not use associated types, the compiler would sometimes find types ambiguous and request explicit type declarations, which were very long and very non-human-readable
//!       * reduces the number of type parameters required for many matrix-related structs, since the type parameters for any `impl MatrixOracle for T ...` implementation must usually come either from the type parameters of the struct or from associated types
//!         * this makes a difference in practice when we switched to associated types:
//!           * many of our structs dropped 3-5 type parameters
//!           * we were able to drop many "PhantomData" fields
//!           * in the case of U-match factorization, it solved two problems:
//!             * an issue with type inference (as mentioned above)
//!             * a problem with the COMB structures; without associated types we were forced to include `ViewMinorDescend` types as
//!               generic type parameters in our COMBs.  However, we wanted the flexibility to run U-match decomposition on matrices that
//!               did not implement ViewCol traits (only ViewRow ones).  In the end, it seemed to us that we would have to create
//!               two different structs for each COMB -- one with a type parameter for minor views and one without.  This was unsatisfactory
//!               as it lead to huge redundancy of code, and all the problems with support/debugging that accompanies redundancy.  Switching to
//!               associated types addressed this problem.
//!       * in effect, provides "optional" type parameters
//!         * cf the discussion of type parameters for minor views of COMBs, above
//!       * there are advantages to associated types, in terms or human-readability (if you see Matrix::ViewMajorAscend, you know it's the type of the ascending major view of the type Matrix)
//!        
//! 
//!     * Added bonus
//!       * one can implement oracle traits on references automatically when the trait has no type parameters (c.f. [this stack overflow exchange](https://stackoverflow.com/questions/72368327/conflicting-implementation-error-that-only-appears-when-type-parameters-are-intr/72368471#72368471))
