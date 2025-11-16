use std::env;
use std::path::PathBuf;

fn main() {
    // Tell cargo to look for the library in the build directory
    let manifest_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
    let lib_dir = PathBuf::from(&manifest_dir)
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("build");

    println!("cargo:rustc-link-search=native={}", lib_dir.display());
    println!("cargo:rustc-link-lib=dylib=yasat");

    // Also check LD_LIBRARY_PATH
    if let Ok(ld_path) = env::var("LD_LIBRARY_PATH") {
        for path in ld_path.split(':') {
            println!("cargo:rustc-link-search=native={}", path);
        }
    }

    // Check system library paths
    println!("cargo:rustc-link-search=native=/usr/local/lib");
    println!("cargo:rustc-link-search=native=/usr/lib");

    // Rerun build script if library changes
    println!("cargo:rerun-if-changed=../../build/libyasat.so");
    println!("cargo:rerun-if-changed=../../src/c_api/yasat.h");
}
