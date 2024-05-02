//! Notes on the design of Umatch factorization, for developers.
//! 
//! # Row-major factorization.
//! 
//! ## Using `(usize, SnzVal)` for the entries of `R_{\rho \rho}^{-1}`
//! 
//! ### Arguments in favor
//! 
//! - minor views of this matrix will naturally be integer-indexed, anyway, since the matrix is stored in vec-of-vec format; so storing entries with major keys will not remove all the instances that require type conversion
//! - integer indices can easily be converted to *either* major or minor indices (converting from major or minor keys to another format is harder)
//!   - counter-point: in practice, it's not clear that we ever need to convert to anything other than major keys
//! 
//! ### Arguments against
//! 
//! - forces us to use various type conversions to convert entries into a format compatible with major keys
//!   - counter-point: not sure that we ever need to convert integer ordinals into anything other than major keys
//! 
//! ## Use of order comparators for the entries of minor descending views of the mapping array
//! 
//! ### Use 1: obtaining major views
//! 
//! The `UmatchRowMajor` struct in this module has a field named `order_comparator_viewminordescendentry`; this 
//! piece of data is necessary for several column look-up operations (it's what we use to simplify linear combinations of
//! column vectors), but it's also important for certain row look-up operations, because many of formulae for 
//! rows given by the "Inner Identities" of Umatch factorization involve linear combinations of the rows of
//! `R_{\rho \rho}^{-1}`.  To simplify these linear combinations you typically need an order comparator for the entries.
//! At present, the entries of `R_{\rho \rho}^{-1}` have type `(usize, SnzVal)` where the `usize` corresponds to the
//! ordinal of a matched major key; no special order comparator is necessary for this because the major
//! ordinals encode the order on matched major keys via the usual total order on integers; however, in future
//! we may use a different data structure to store `R_{\rho \rho}^{-1}`, which uses major keys as incides (this is
//! motivated by the case where we do not store apparent pairs in memory); in that case it seems more likely that
//! we will need to store an order comparator for key-val pairs where the key has type KeyMaj.
//! 
// //! In the current implementation of Umatch on oat_rust you can technically side-step this requirement because
// //! we store a bijection between matched row indices and their ordinals (in fact, entries of R_{\rho \rho})
// //! are stored with ordinals for indices, rather than major keys).  However in future iterations we may no 
// //! longer store such a bijection (e.g., if we eventually stop storing apparent pairs).  In that case we may
// //! have little choice but to hang on to an order comparator for key-val pairs with major keys.
// //! 
// //! 
// //! may seem unnecessary for the purpose of running row-lookup operations on the Umatch object, 
// //! because extraneous to the Umatch factorization, because (i) 
// //! 
// //! 
// //! since we only need to compute the matched block of
// //! R^{-1}, for some codomain COMB R.  However, in order to reconstruct rows and columns of the matrices involved in the
// //! factorization, we have to *solve systems of linear equations* that involve the matched block of R^{-1}.  Those 
// //! systems are hard to solve without leveraging the fact that R^{-1} is upper triangular, so we rely heavily on the
// //! fact that it is.  However, it is upper triangular *with respect to the order on **major** indices*.  Our 
// //! triangular solvers require access to information about the order of indices in order to function correctly, so it might seem 
// //! unavoidable that a Umatch struct would have to contain order information about major keys in order to function correctly.
// //! This isn't quite true, as we can store R^{-1}, as the factorization algorithm makes it easy to store the matched block of R^{-1}
// //! as a square matrix indexed by *integers*, where the order on integers agrees with the order on major keys.  However, in
// //! future versions of the algorithm, it may no longer be easy to store an integer-indexed matrix (for example, this comes up in contexts where we do not wish to store
// //! apparent pairs in memory).  So in future it may truly be unavoidable to store order information about major keys.
//! 
//! ### Use 2: obtaining minor views
//! 
//! Storing order information is unavoidable if one wishes to recover minor views of matrices in the factorization.
//! 
//! ### Consequenes for type parameters
//! 
//! Give the notes above, and given that our strongest motivating use case (computing cycle representatives for persistent homology)
//! involves minor views of the factorization, it seems to make sense to store order information about major keys in the 
//! most commonly used umatch struct.  Currently, that means adding the field `order_comparator_viewminordescendentry` to
//! `UmatchRowMajor`.  Consequently, we need to store a type parameter for `order_comparator_viewminordescendentry` in `UmatchRowMajor`.
//! 
//! ### Solutions for potential future problems.
//! 
//! This may cause problems in the future for users who don't have anatural order on major keys.  To accomodate such users, it may make
//! sense to
//!
//! - create a (highly unsafe) order comparator that returns false for every pair; this order comparator requires no type parameters
//! - use this comparator to generate a `UmatchRowMajor` struct, and wrap the resulting struct as a private field in a new wrapper struct
//! - implement only row look-ups on the wrapper struct
//! - if implemented, this method should be tested heavily in order to ensure that no order information about major keys is inadvertently used to access rows, resulting in errors
//! - note that this method also depends on storing the matched block of R^{-1} as an integer-indexed matrix
//!   - it might be possible, alternatively,to store this block as a minor-key-indexed matrix, but it doesn't seem clear that this is feasible
//! 
//! 
//! ## Infeasibility of naive "antitranspose" method for computing minor views
//! 
//! There are sections of the oat_rust library where methods for minor views can be implemented by a simple trick: 
//! applying a method for major views to the (anti)transpose of a matrix.  One might wish to apply such a trick here
//! in order to get the minor views of a Umatch decomposition.  However this seems difficult, because an important
//! part of storing a (compressed) Umatch decomposition is storing the matched block of R^{-1}, where R is the codomain COMB.  If we 
//! antitransposed the underlying matrix, the corresponding storage scheme would involve storing the matched block of C,
//! where C is the domain COMB.  **We cannot simply antitranspose the matched block of R.**