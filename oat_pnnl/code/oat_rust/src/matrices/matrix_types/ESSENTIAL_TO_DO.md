


* update the oracle traits with lifetimes as learned in the crate "iter_forward_and_back_traits"
* update the reduction algorithms to work with oracles using these updated traits (recall the issue was getting iterators that run backward)
* write mini-tutorial (just a few paragraphs plus example) explaining how to deal with lifetime parameters in this context -- eg in cases where you want to get an iterator now and save it for later, but later you'll also want to change the internal state of the matrix oracle (is that even possible?)