//! Tools that use the Gurobi solver
//! 
//! Gurobi is a cutting-edge optimization engine with excellent performance.  Using Gurobi takes some extra work (getting a license, downloading, etc.) but if performance is critical, then it may be worth the work.
//!
//! # Help installing the Gurobi optimizer
//! 
//! This module uses the Gurobi optimization engine.  To use it, you'll need to install Gurobi.
//! 
//! #### Get a free license
//! 
//! Gurobi provides free licenses for published research and education.
//! 
//! #### Update your system path
//! 
//! Once you've installed Gurobi, you may need to add some variables to your system path.  This isn't as hard as it sounds!
//! 
//! - See [grb](https://docs.rs/grb/latest/grb/) for additional hints
//! - First you'll need to find where some Gurobi files are located.  The Gurobi website will give some tips, e.g. for most Mac users the files will be at `/Library/gurobi1001/macos_universal2/bin/gurobi_cl`.  BUT BE AWARE that Macs have many folders called `Library`.  To get to the right one, open Finder and from the menu bar select `Go > Computer > Library`, not `Go > Library`.
//! - Second you'll need to add some information to your system path.  On mac, go to your user folder and find `.profile`, `.bash_profile`, `.zshrc`, or a variant thereof
//!   - this file is sometimes hidden; to reveal hidden files select "COMMAND + SHIFT + ."
//!   - if you don't have one you can create one; just open an editor like sublime, atom, or vs code, and save a file with name `.profile`
//!   - open the file, and add the following lines, replacing `user_name` with your user name
//!     ```text
//!     export GRB_LICENSE.md_FILE="/Users/user_name/gurobi.lic" # the license file
//!     export GUROBI_HOME="/Library/gurobi1001/macos_universal2" # /Library/gurobi1001/macos_universal2/bin/gurobi_cl
//!     ```
//! - Save your changes
//! - Finally, you may need to refresh your terminal in order for your changes to take effect.  You can do this be executing either `source path/to/.profile` or `. path/to/.profile` in a terminal
//! 
//! #### Test your installation
//! 
//!   
//! To test the gurobi installation, see [this link](https://www.gurobi.com/documentation/9.5/quickstart_mac/testing_your_license.html)

pub mod l1_norm;